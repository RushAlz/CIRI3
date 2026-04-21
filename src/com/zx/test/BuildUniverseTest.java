package com.zx.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import com.zx.findcircrna.CircRNAUniverseIO;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;

public class BuildUniverseTest {

    /**
     * Builds the circRNA universe from all SCAN1 outputs and writes a universe file.
     *
     * samplesMetaTsv: three-column TSV per sample:
     *   samFile  scan1_meta_path  bsjFile1,bsjFile2,...,bsjFileN
     * faFile:         FASTA reference genome
     * outputPrefix:   prefix for output files; writes {outputPrefix}.universe
     *
     * The comma-separated BSJ file list in column 3 is the authoritative set of BSJ
     * files produced by SCAN1 for that sample.  Using an explicit list prevents stale
     * files left over from earlier runs with a different thread count from leaking in.
     */
    public void build(String samplesMetaTsv, String faFile, String outputPrefix) throws IOException {
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :BUILD_UNIVERSE start");

        // Step 1: Read samples TSV and all scan1_meta files
        ArrayList<String> filePathList = new ArrayList<String>();
        HashMap<String, ArrayList<String>> bsjFilesMap = new HashMap<String, ArrayList<String>>();
        int globalReadLen = 0;

        BufferedReader tsvBr = new BufferedReader(new FileReader(new File(samplesMetaTsv)));
        String tsvLine = tsvBr.readLine();
        while (tsvLine != null) {
            if (tsvLine.startsWith("#") || tsvLine.equals("")) {
                tsvLine = tsvBr.readLine();
                continue;
            }
            String[] tsvArr = tsvLine.split("\t");
            String samFilePath = tsvArr[0].trim();
            String metaFilePath = tsvArr[1].trim();
            String bsjFilesCol  = tsvArr[2].trim();
            filePathList.add(samFilePath);

            // Parse explicit BSJ file list from column 3
            ArrayList<String> bsjPaths = new ArrayList<String>(Arrays.asList(bsjFilesCol.split(",")));
            bsjFilesMap.put(samFilePath, bsjPaths);

            // Read scan1_meta for globalReadLen
            BufferedReader metaBr = new BufferedReader(new FileReader(new File(metaFilePath)));
            String metaLine = metaBr.readLine();
            while (metaLine != null) {
                if (metaLine.startsWith("readLen=")) {
                    int readLen = Integer.parseInt(metaLine.split("=")[1].trim());
                    if (readLen > globalReadLen) {
                        globalReadLen = readLen;
                    }
                }
                metaLine = metaBr.readLine();
            }
            metaBr.close();
            tsvLine = tsvBr.readLine();
        }
        tsvBr.close();

        System.out.println(df.format(System.currentTimeMillis()) + " :Loaded " + filePathList.size() + " samples");

        // Step 2: Load FA file
        ReadFaFile RF = new ReadFaFile();
        RF.readFa(faFile);
        HashMap<String, Integer> chrLenMap = RF.getChrLenMap();
        RF = null;
        System.out.println(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files");

        // Step 3: Determine global seqLen (max readLen across samples, minus 12)
        int seqLen = globalReadLen - 12;
        System.out.println(df.format(System.currentTimeMillis()) + " :seqLen=" + seqLen + " (readLen=" + globalReadLen + ")");

        // Step 4: Read all BSJ files from all samples to build chrCircSiteMap
        // (identical to MutFileTest lines 306-327 / MutTest lines 256-273)
        // BSJ file paths come from the explicit list in column 3 of the input TSV,
        // not reconstructed from samFile+fileSplitNum, to prevent stale-file leakage.
        HashMap<String, HashSet<String>> chrCircSiteMap = new HashMap<String, HashSet<String>>();
        HashSet<String> circSiteSet = new HashSet<String>();

        for (int i = 0; i < filePathList.size(); i++) {
            String samFilePath = filePathList.get(i);
            for (String bsjPath : bsjFilesMap.get(samFilePath)) {
                BufferedReader BSJbr = new BufferedReader(new FileReader(new File(bsjPath)));
                String line = BSJbr.readLine();
                while (line != null) {
                    String[] BSJArr = line.split("\t", 5);
                    if (!chrCircSiteMap.containsKey(BSJArr[3])) {
                        circSiteSet = new HashSet<String>();
                        circSiteSet.add(BSJArr[4]);
                        chrCircSiteMap.put(BSJArr[3], circSiteSet);
                    } else {
                        circSiteSet = chrCircSiteMap.get(BSJArr[3]);
                        circSiteSet.add(BSJArr[4]);
                        chrCircSiteMap.put(BSJArr[3], circSiteSet);
                    }
                    line = BSJbr.readLine();
                }
                BSJbr.close();
            }
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :BSJ sites loaded from all samples");

        // Step 5: Build circFSJMap and site-index structures
        // (identical to MutTest lines 282-358)
        HashMap<String, Integer> circFSJMap = new HashMap<String, Integer>();
        HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap1 =
                new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        HashMap<Integer, ArrayList<SiteSort>> SiteMap1 = new HashMap<Integer, ArrayList<SiteSort>>();
        ArrayList<SiteSort> siteList1 = new ArrayList<SiteSort>();
        HashMap<String, HashMap<Integer, ArrayList<SiteSort>>> chrSiteMap2 =
                new HashMap<String, HashMap<Integer, ArrayList<SiteSort>>>();
        HashMap<Integer, ArrayList<SiteSort>> SiteMap2 = new HashMap<Integer, ArrayList<SiteSort>>();
        ArrayList<SiteSort> siteList2 = new ArrayList<SiteSort>();

        // universeDataMap: circFSJMap-key -> siteInfor (for writing the universe file)
        HashMap<String, String> universeDataMap = new HashMap<String, String>();

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

        // Sort site lists (identical to MutTest lines 340-357)
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

        // Step 6: Write universe file
        String universePath = outputPrefix + ".universe";
        CircRNAUniverseIO.writeUniverse(universePath, seqLen, universeDataMap);
        System.out.println(df.format(System.currentTimeMillis()) + " :Universe written to " + universePath);
    }
}
