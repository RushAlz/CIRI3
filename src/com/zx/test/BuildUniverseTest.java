package com.zx.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;

import com.zx.findcircrna.CircRNAUniverseIO;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;

public class BuildUniverseTest {

    /**
     * Builds the circRNA universe from all SCAN1 outputs and writes a universe file.
     *
     * samplesMetaTsv: two-column TSV (samFile path, scan1_meta path) per sample.
     *   The samFile column is accepted for backward compatibility but is not
     *   opened — only the scan1_meta file is read to obtain bsjPrefix,
     *   fileSplitNum, and readLen.
     * faFile:         FASTA reference genome
     * outputPrefix:   prefix for output files; writes {outputPrefix}.universe
     */
    public void build(String samplesMetaTsv, String faFile, String outputPrefix) throws IOException {
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :BUILD_UNIVERSE start");

        ArrayList<String> bsjPrefixList = new ArrayList<String>();
        HashMap<String, Integer> fileSplitNumMap = new HashMap<String, Integer>();
        int globalReadLen = 0;

        BufferedReader tsvBr = new BufferedReader(new FileReader(new File(samplesMetaTsv)));
        String tsvLine = tsvBr.readLine();
        while (tsvLine != null) {
            if (tsvLine.startsWith("#") || tsvLine.equals("")) {
                tsvLine = tsvBr.readLine();
                continue;
            }
            String[] tsvArr = tsvLine.split("\t");
            // Column 0 (samFile) is kept for backward compatibility but not used.
            String metaFilePath = tsvArr[1].trim();

            String bsjPrefix = null;
            int splitNum = 0;

            BufferedReader metaBr = new BufferedReader(new FileReader(new File(metaFilePath)));
            String metaLine = metaBr.readLine();
            while (metaLine != null) {
                if (metaLine.startsWith("readLen=")) {
                    int readLen = Integer.parseInt(metaLine.split("=")[1].trim());
                    if (readLen > globalReadLen) globalReadLen = readLen;
                } else if (metaLine.startsWith("fileSplitNum=")) {
                    splitNum = Integer.parseInt(metaLine.split("=")[1].trim());
                } else if (metaLine.startsWith("bsjPrefix=")) {
                    bsjPrefix = metaLine.split("=", 2)[1].trim();
                }
                metaLine = metaBr.readLine();
            }
            metaBr.close();

            if (bsjPrefix == null || splitNum == 0) {
                System.out.println("WARNING: incomplete scan1_meta at " + metaFilePath + " — skipping sample");
                tsvLine = tsvBr.readLine();
                continue;
            }
            bsjPrefixList.add(bsjPrefix);
            fileSplitNumMap.put(bsjPrefix, splitNum);
            tsvLine = tsvBr.readLine();
        }
        tsvBr.close();
        System.out.println(df.format(System.currentTimeMillis()) + " :Loaded " + bsjPrefixList.size() + " samples");

        buildFromBsjPrefixes(bsjPrefixList, fileSplitNumMap, globalReadLen, faFile, outputPrefix, df);
    }

    /**
     * Builds the circRNA universe from a direct BSJ-list TSV, bypassing scan1_meta
     * files. This is the preferred input format for cloud environments where the
     * paths stored in scan1_meta are VM-local and not accessible across steps.
     *
     * bsjListTsv: three-column TSV per sample:
     *   bsjPrefix<TAB>fileSplitNum<TAB>readLen
     *   - bsjPrefix:    path prefix where the BSJ files (BSJ1, BSJ2, …) are located
     *                   on the current machine (e.g. after staging from object storage)
     *   - fileSplitNum: number of BSJ files written by SCAN1 for this sample
     *   - readLen:      maximum read length reported by SCAN1 for this sample
     * faFile:       FASTA reference genome
     * outputPrefix: prefix for output files; writes {outputPrefix}.universe
     */
    public void buildFromBsjList(String bsjListTsv, String faFile, String outputPrefix) throws IOException {
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :BUILD_UNIVERSE (direct BSJ list) start");

        ArrayList<String> bsjPrefixList = new ArrayList<String>();
        HashMap<String, Integer> fileSplitNumMap = new HashMap<String, Integer>();
        int globalReadLen = 0;

        BufferedReader tsvBr = new BufferedReader(new FileReader(new File(bsjListTsv)));
        String tsvLine = tsvBr.readLine();
        while (tsvLine != null) {
            if (tsvLine.startsWith("#") || tsvLine.equals("")) {
                tsvLine = tsvBr.readLine();
                continue;
            }
            String[] arr = tsvLine.split("\t");
            String bsjPrefix = arr[0].trim();
            int splitNum = Integer.parseInt(arr[1].trim());
            int readLen = Integer.parseInt(arr[2].trim());
            bsjPrefixList.add(bsjPrefix);
            fileSplitNumMap.put(bsjPrefix, splitNum);
            if (readLen > globalReadLen) globalReadLen = readLen;
            tsvLine = tsvBr.readLine();
        }
        tsvBr.close();
        System.out.println(df.format(System.currentTimeMillis()) + " :Loaded " + bsjPrefixList.size() + " samples");

        buildFromBsjPrefixes(bsjPrefixList, fileSplitNumMap, globalReadLen, faFile, outputPrefix, df);
    }

    /**
     * Core universe-building logic shared by build() and buildFromBsjList().
     * Reads all BSJ files, aggregates circRNA candidates across samples, builds
     * spatial-index structures, and writes the universe file.
     */
    private void buildFromBsjPrefixes(ArrayList<String> bsjPrefixList,
            HashMap<String, Integer> fileSplitNumMap, int globalReadLen,
            String faFile, String outputPrefix, SimpleDateFormat df) throws IOException {

        // Load FA file
        ReadFaFile RF = new ReadFaFile();
        RF.readFa(faFile);
        HashMap<String, Integer> chrLenMap = RF.getChrLenMap();
        RF = null;
        System.out.println(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files");

        int seqLen = globalReadLen - 12;
        System.out.println(df.format(System.currentTimeMillis()) + " :seqLen=" + seqLen + " (readLen=" + globalReadLen + ")");

        // Read all BSJ files from all samples to build chrCircSiteMap
        HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
        HashSet<String> circSiteSet = new HashSet<String>();

        for (String bsjPrefix : bsjPrefixList) {
            int splitNum = fileSplitNumMap.get(bsjPrefix);
            for (int j = 1; j <= splitNum; j++) {
                BufferedReader BSJbr = new BufferedReader(new FileReader(new File(bsjPrefix + "BSJ" + j)), 262144);
                String line = BSJbr.readLine();
                while (line != null) {
                    String[] BSJArr = line.split("\t", 5);
                    circSiteSet = chrCircSiteMap.computeIfAbsent(BSJArr[3], k -> new HashSet<>());
                    circSiteSet.add(BSJArr[4]);
                    line = BSJbr.readLine();
                }
                BSJbr.close();
            }
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :BSJ sites loaded from all samples");

        // Build circFSJMap and site-index structures
        HashMap<String, Integer> circFSJMap = new HashMap<String, Integer>();
        HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1 =
                new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        HashMap<Integer, ArrayList<SiteSort>> SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
        ArrayList<SiteSort> siteList1 = new ArrayList<SiteSort>();
        HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2 =
                new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        HashMap<Integer, ArrayList<SiteSort>> SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
        ArrayList<SiteSort> siteList2 = new ArrayList<SiteSort>();

        // LinkedHashMap preserves insertion order so universe file order is stable.
        LinkedHashMap<String, String> universeDataMap = new LinkedHashMap<String, String>();

        for (String chrKey : chrCircSiteMap.keySet()) {
            circSiteSet = chrCircSiteMap.get(chrKey);
            byte[] siteArray1 = new byte[(chrLenMap.get(chrKey) / seqLen) + 1];
            byte[] siteArray2 = new byte[(chrLenMap.get(chrKey) / seqLen) + 1];
            SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
            SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
            for (String siteInfor : circSiteSet) {
                String[] arr = siteInfor.split("\t");
                String temCircRNA = chrKey + "\t" + arr[0] + "\t" + arr[1];
                circFSJMap.put(temCircRNA, 0);
                universeDataMap.put(temCircRNA, siteInfor);
                int site1 = Integer.parseInt(arr[0]) / seqLen;
                int site2 = Integer.parseInt(arr[1]) / seqLen;
                siteArray1[site1] = 1;
                siteArray2[site2] = 1;
                if (!SiteMap1.containsKey(site1)) {
                    siteList1 = new ArrayList<SiteSort>();
                    siteList1.add(new SiteSort(Integer.parseInt(arr[0]), arr));
                    SiteMap1.put(site1, siteList1);
                } else {
                    siteList1 = SiteMap1.get(site1);
                    siteList1.add(new SiteSort(Integer.parseInt(arr[0]), arr));
                    SiteMap1.put(site1, siteList1);
                }
                if (!SiteMap2.containsKey(site2)) {
                    siteList2 = new ArrayList<SiteSort>();
                    siteList2.add(new SiteSort(Integer.parseInt(arr[1]), arr));
                    SiteMap2.put(site2, siteList2);
                } else {
                    siteList2 = SiteMap2.get(site2);
                    siteList2.add(new SiteSort(Integer.parseInt(arr[1]), arr));
                    SiteMap2.put(site2, siteList2);
                }
            }
            chrSiteMap1.put(chrKey, SiteMap1);
            chrSiteMap2.put(chrKey, SiteMap2);
        }
        chrCircSiteMap = null;
        circSiteSet = null;

        // Sort site lists
        for (String chrKey : chrSiteMap1.keySet()) {
            SiteMap1 = chrSiteMap1.get(chrKey);
            for (Integer site : SiteMap1.keySet()) {
                siteList1 = SiteMap1.get(site);
                Collections.sort(siteList1);
                SiteMap1.put(site, siteList1);
            }
            chrSiteMap1.put(chrKey, SiteMap1);
        }
        for (String chrKey : chrSiteMap2.keySet()) {
            SiteMap2 = chrSiteMap2.get(chrKey);
            for (Integer site : SiteMap2.keySet()) {
                siteList2 = SiteMap2.get(site);
                Collections.sort(siteList2);
                SiteMap2.put(site, siteList2);
            }
            chrSiteMap2.put(chrKey, SiteMap2);
        }

        System.out.println(df.format(System.currentTimeMillis()) + " :Universe contains " + circFSJMap.size() + " circRNA candidates");

        String universePath = outputPrefix + ".universe";
        CircRNAUniverseIO.writeUniverse(universePath, seqLen, universeDataMap);
        System.out.println(df.format(System.currentTimeMillis()) + " :Universe written to " + universePath);
    }
}
