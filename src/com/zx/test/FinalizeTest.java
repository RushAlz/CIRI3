package com.zx.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.HashMap;

import com.zx.findcircrna.CircRNAUniverseIO;
import com.zx.findcircrna.GetAnnotationInformation;
import com.zx.findcircrna.ReadFaFile;
import com.zx.findcircrna.SiteSort;
import com.zx.findcircrna.Summary;
import com.zx.hg38.Annotation;
import com.zx.hg38.AnnotationIntron;

public class FinalizeTest {
    private int strigency, relExp;
    private boolean intronLable;

    public FinalizeTest(int minMapqUni, int maxCircle, int minCircle, int linear_range_size_min,
            boolean intronLable, int strigency, int relExp, String mitochondrion,
            boolean mlable, boolean spLable) {
        this.intronLable = intronLable;
        this.strigency = strigency;
        this.relExp = relExp;
    }

    /**
     * Merges FSJ counts from all samples, runs Summary + Annotation, and writes matrix outputs.
     *
     * finalizeInputTsv: five-column TSV per sample:
     *   samFile  fsjCountsFile  fileSplitNum  sampleName  bsjPrefix
     *   The samFile column (column 1) is accepted for backward compatibility but
     *   is not used — only bsjPrefix (column 5) is consulted when reading BSJ files.
     * faFile:         FASTA reference genome
     * annotationFile: GTF/GFF3 or "F"
     * outputPrefix:   prefix for output files
     */
    public void finalize(String finalizeInputTsv, String faFile, String annotationFile,
            String outputPrefix) throws IOException {
        finalize(finalizeInputTsv, faFile, annotationFile, outputPrefix, null);
    }

    /**
     * Variant that takes the BUILD_UNIVERSE output as the authoritative circRNA set.
     * When universeFile is non-null, circRowMap is seeded from it rather than from
     * the first sample's .fsj_counts, so FINALIZE is tolerant of partial re-runs.
     */
    public void finalize(String finalizeInputTsv, String faFile, String annotationFile,
            String outputPrefix, String universeFile) throws IOException {
        finalize(finalizeInputTsv, faFile, annotationFile, outputPrefix, universeFile, null);
    }

    /**
     * Variant with a frozen circRNA filter loaded from an existing BSJ_Matrix file.
     * When freezeMatrixFile is non-null, Summary re-evaluation is skipped entirely:
     * the accepted circRNA set is taken verbatim from the first column of that matrix.
     * This is the correct mode when adding a new sample to an existing cohort.
     */
    public void finalize(String finalizeInputTsv, String faFile, String annotationFile,
            String outputPrefix, String universeFile, String freezeMatrixFile) throws IOException {
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :FINALIZE start");

        // Parse five-column TSV: samFile  fsjCountsFile  fileSplitNum  sampleName  [bsjPrefix]
        // samFile (col 1) is vestigial — present only for backward compatibility.
        ArrayList<String> fsjCountsFileList = new ArrayList<String>();
        ArrayList<String> sampleNameList = new ArrayList<String>();
        ArrayList<String> bsjPrefixList = new ArrayList<String>();
        HashMap<String, Integer> bsjPrefixSplitNumMap = new HashMap<String, Integer>();

        BufferedReader tsvBr = new BufferedReader(new FileReader(new File(finalizeInputTsv)));
        String tsvLine = tsvBr.readLine();
        while (tsvLine != null) {
            if (tsvLine.startsWith("#") || tsvLine.equals("")) {
                tsvLine = tsvBr.readLine();
                continue;
            }
            String[] arr = tsvLine.split("\t");
            String samFilePath = arr[0].trim();   // not used, kept for compat
            String fsjCountsPath = arr[1].trim();
            int splitNum = Integer.parseInt(arr[2].trim());
            String sampleName = arr[3].trim();
            String bsjPrefix = (arr.length >= 5 && !arr[4].trim().isEmpty()) ? arr[4].trim() : samFilePath;
            fsjCountsFileList.add(fsjCountsPath);
            sampleNameList.add(sampleName);
            bsjPrefixList.add(bsjPrefix);
            bsjPrefixSplitNumMap.put(bsjPrefix, splitNum);
            tsvLine = tsvBr.readLine();
        }
        tsvBr.close();

        runFinalize(fsjCountsFileList, sampleNameList, bsjPrefixList, bsjPrefixSplitNumMap,
                faFile, annotationFile, outputPrefix, universeFile, freezeMatrixFile, df);
    }

    /**
     * Merges FSJ counts and builds output matrices using a direct BSJ-list TSV,
     * bypassing the samFile column that is present in the legacy format. This is
     * the preferred input format for cloud environments where file paths are not
     * persistent across pipeline steps.
     *
     * finalizeInputTsv: four-column TSV per sample:
     *   bsjPrefix<TAB>fileSplitNum<TAB>fsjCountsFile<TAB>sampleName
     *   - bsjPrefix:    path prefix where the BSJ files are located on the current machine
     *   - fileSplitNum: number of BSJ files for this sample (from SCAN1 metadata)
     *   - fsjCountsFile: path to the .fsj_counts file produced by SCAN2
     *   - sampleName:   column header for this sample in the output matrices
     * faFile:         FASTA reference genome
     * annotationFile: GTF/GFF3 or "F"
     * outputPrefix:   prefix for output files
     */
    public void finalizeFromBsjList(String finalizeInputTsv, String faFile, String annotationFile,
            String outputPrefix) throws IOException {
        finalizeFromBsjList(finalizeInputTsv, faFile, annotationFile, outputPrefix, null);
    }

    public void finalizeFromBsjList(String finalizeInputTsv, String faFile, String annotationFile,
            String outputPrefix, String universeFile) throws IOException {
        finalizeFromBsjList(finalizeInputTsv, faFile, annotationFile, outputPrefix, universeFile, null);
    }

    public void finalizeFromBsjList(String finalizeInputTsv, String faFile, String annotationFile,
            String outputPrefix, String universeFile, String freezeMatrixFile) throws IOException {
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :FINALIZE (direct BSJ list) start");

        ArrayList<String> fsjCountsFileList = new ArrayList<String>();
        ArrayList<String> sampleNameList = new ArrayList<String>();
        ArrayList<String> bsjPrefixList = new ArrayList<String>();
        HashMap<String, Integer> bsjPrefixSplitNumMap = new HashMap<String, Integer>();

        BufferedReader tsvBr = new BufferedReader(new FileReader(new File(finalizeInputTsv)));
        String tsvLine = tsvBr.readLine();
        while (tsvLine != null) {
            if (tsvLine.startsWith("#") || tsvLine.equals("")) {
                tsvLine = tsvBr.readLine();
                continue;
            }
            String[] arr = tsvLine.split("\t");
            String bsjPrefix = arr[0].trim();
            int splitNum = Integer.parseInt(arr[1].trim());
            String fsjCountsPath = arr[2].trim();
            String sampleName = arr[3].trim();
            bsjPrefixList.add(bsjPrefix);
            bsjPrefixSplitNumMap.put(bsjPrefix, splitNum);
            fsjCountsFileList.add(fsjCountsPath);
            sampleNameList.add(sampleName);
            tsvLine = tsvBr.readLine();
        }
        tsvBr.close();
        System.out.println(df.format(System.currentTimeMillis()) + " :Loaded " + sampleNameList.size() + " samples");

        runFinalize(fsjCountsFileList, sampleNameList, bsjPrefixList, bsjPrefixSplitNumMap,
                faFile, annotationFile, outputPrefix, universeFile, freezeMatrixFile, df);
    }

    /**
     * Core FINALIZE logic shared by finalize() and finalizeFromBsjList().
     */
    private void runFinalize(
            ArrayList<String> fsjCountsFileList,
            ArrayList<String> sampleNameList,
            ArrayList<String> bsjPrefixList,
            HashMap<String, Integer> bsjPrefixSplitNumMap,
            String faFile, String annotationFile, String outputPrefix,
            String universeFile, String freezeMatrixFile,
            SimpleDateFormat df) throws IOException {

        int sampleCount = sampleNameList.size();
        System.out.println(df.format(System.currentTimeMillis()) + " :Loaded " + sampleCount + " samples");

        // Load FA file
        ReadFaFile RF = new ReadFaFile();
        RF.readFa(faFile);
        HashMap<String, String> chrTCGAMap = RF.getChrTCGAMap();
        RF = null;
        System.out.println(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files");

        // Load annotation (optional)
        HashMap<String, ArrayList<SiteSort>> geneExonMap = new HashMap<String, ArrayList<SiteSort>>();
        HashMap<String, ArrayList<Integer[]>> exonListMap = new HashMap<String, ArrayList<Integer[]>>();
        HashMap<String, String> chrExonStartMap = new HashMap<String, String>();
        HashMap<String, String> chrExonEndMap = new HashMap<String, String>();
        HashMap<String, ArrayList<String>> chrExonStartTranscriptMap = new HashMap<String, ArrayList<String>>();
        HashMap<String, ArrayList<String>> chrExonEndTranscriptMap = new HashMap<String, ArrayList<String>>();

        if (!annotationFile.equals("F")) {
            GetAnnotationInformation GAI = new GetAnnotationInformation();
            GAI.hand(annotationFile, intronLable);
            if (intronLable) {
                chrExonStartTranscriptMap = GAI.getChrExonStartTranscriptMap();
                chrExonEndTranscriptMap = GAI.getChrExonEndTranscriptMap();
                geneExonMap = GAI.getGeneExonMap();
                exonListMap = GAI.getExonListMap();
            } else {
                chrExonStartMap = GAI.getChrExonStartMap();
                chrExonEndMap = GAI.getChrExonEndMap();
                geneExonMap = GAI.getGeneExonMap();
                exonListMap = GAI.getExonListMap();
            }
            GAI = null;
            if (geneExonMap.size() == 0) {
                System.out.println("please input formatted annotation file");
                return;
            }
            System.out.println(df.format(System.currentTimeMillis()) + " :Successfully imported comment files");
        }

        // Seed the universe: prefer BUILD_UNIVERSE output when supplied, else
        // fall back to the first sample's .fsj_counts for backward compatibility.
        HashMap<String, Integer> circFSJMerged = new HashMap<String, Integer>();
        if (universeFile != null && !universeFile.isEmpty()) {
            int[] seqLenOut = new int[1];
            HashMap<String, String[]> universeDataMap = CircRNAUniverseIO.readUniverse(universeFile, seqLenOut);
            for (String circKey : universeDataMap.keySet()) {
                circFSJMerged.put(circKey, 0);
            }
            System.out.println(df.format(System.currentTimeMillis()) + " :Universe seeded from " + universeFile + ": " + circFSJMerged.size() + " circRNAs");
        } else {
            BufferedReader firstBr = new BufferedReader(new FileReader(new File(fsjCountsFileList.get(0))));
            String line = firstBr.readLine();
            while (line != null) {
                String[] parts = line.split("\t");
                String circKey = parts[0] + "\t" + parts[1] + "\t" + parts[2];
                circFSJMerged.put(circKey, 0);
                line = firstBr.readLine();
            }
            firstBr.close();
        }

        // Build circRowMap (universe index) and per-sample FSJ matrix
        int circNum = circFSJMerged.size();
        int[][] FSJmatrix = new int[circNum][sampleCount];
        HashMap<String, Integer> circRowMap = new HashMap<String, Integer>();
        int rowIdx = 0;
        for (String circKey : circFSJMerged.keySet()) {
            circRowMap.put(circKey, rowIdx++);
        }

        // Read all .fsj_counts: populate per-sample FSJ matrix
        for (int j = 0; j < sampleCount; j++) {
            BufferedReader fsjBr = new BufferedReader(new FileReader(new File(fsjCountsFileList.get(j))));
            String line = fsjBr.readLine();
            while (line != null) {
                String[] parts = line.split("\t");
                String circKey = parts[0] + "\t" + parts[1] + "\t" + parts[2];
                int count = Integer.parseInt(parts[3]);
                if (circRowMap.containsKey(circKey)) {
                    FSJmatrix[circRowMap.get(circKey)][j] = count;
                    circFSJMerged.put(circKey, circFSJMerged.get(circKey) + count);
                }
                line = fsjBr.readLine();
            }
            fsjBr.close();
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :FSJ counts merged: " + circFSJMerged.size() + " circRNAs");

        // Build BSJ matrix from BSJ files
        int[][] BSJmatrix = new int[circNum][sampleCount];
        for (int i = 0; i < sampleCount; i++) {
            HashMap<String, Integer> circMap = new HashMap<String, Integer>();
            String bsjPrefix = bsjPrefixList.get(i);
            int splitNum = bsjPrefixSplitNumMap.get(bsjPrefix);
            long totalLines = 0, tagOneLines = 0, unknownCirc = 0;
            for (int j = 1; j <= splitNum; j++) {
                long fileLines = 0, fileTagOne = 0;
                BufferedReader BSJBr = new BufferedReader(new FileReader(new File(bsjPrefix + "BSJ" + j)), 262144);
                String line = BSJBr.readLine();
                while (line != null) {
                    fileLines++;
                    String[] circLineArr = line.split("\t", 7);
                    if (circLineArr[2].equals("1")) {
                        fileTagOne++;
                        String chrStartEnd = circLineArr[3] + "\t" + circLineArr[4] + "\t" + circLineArr[5];
                        circMap.merge(chrStartEnd, 1, Integer::sum);
                    }
                    line = BSJBr.readLine();
                }
                BSJBr.close();
                totalLines += fileLines;
                tagOneLines += fileTagOne;
                System.out.println(df.format(System.currentTimeMillis()) + " :DIAG sample[" + i + "] BSJ" + j + ": total=" + fileLines + " tag1=" + fileTagOne);
            }
            for (String circKey : circMap.keySet()) {
                if (circRowMap.containsKey(circKey)) {
                    BSJmatrix[circRowMap.get(circKey)][i] = circMap.get(circKey);
                } else {
                    unknownCirc++;
                }
            }
            System.out.println(df.format(System.currentTimeMillis()) + " :DIAG sample[" + i + "] name=" + sampleNameList.get(i) + " bsj_total=" + totalLines + " tag1=" + tagOneLines + " circRNAs_in_universe=" + circMap.size() + " circRNAs_not_in_universe=" + unknownCirc);
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :BSJ matrix built");

        // Determine accepted circRNA set
        ArrayList<String> SummaryCircList = new ArrayList<String>();
        HashMap<String, String> circTrueIdMap;
        if (freezeMatrixFile != null && !freezeMatrixFile.isEmpty()) {
            circTrueIdMap = new HashMap<String, String>();
            BufferedReader freezeBr = new BufferedReader(new FileReader(new File(freezeMatrixFile)));
            freezeBr.readLine(); // skip header
            String freezeLine = freezeBr.readLine();
            while (freezeLine != null) {
                String circId = freezeLine.split("\t")[0];
                circTrueIdMap.put(circId, "");
                freezeLine = freezeBr.readLine();
            }
            freezeBr.close();
            System.out.println(df.format(System.currentTimeMillis()) + " :Frozen circRNA set from " + freezeMatrixFile + ": " + circTrueIdMap.size() + " circRNAs (Summary skipped)");
        } else {
            Summary summary = new Summary(strigency, chrTCGAMap);
            SummaryCircList = summary.summary(bsjPrefixList, bsjPrefixSplitNumMap, circFSJMerged, "");
            circTrueIdMap = summary.getCircMap();
            summary = null;
            System.out.println(df.format(System.currentTimeMillis()) + " :Summary completed");
        }
        chrTCGAMap = null;

        // Write BSJ_Matrix and FSJ_Matrix
        String outPutBSJCountFile = outputPrefix + ".BSJ_Matrix";
        String outPutFSJCountFile = outputPrefix + ".FSJ_Matrix";
        BufferedWriter BSJCount = new BufferedWriter(new FileWriter(new File(outPutBSJCountFile)));
        BufferedWriter FSJCount = new BufferedWriter(new FileWriter(new File(outPutFSJCountFile)));
        BSJCount.write("circRNA_ID");
        FSJCount.write("circRNA_ID");
        for (String sampleName : sampleNameList) {
            BSJCount.write("\t" + sampleName);
            FSJCount.write("\t" + sampleName);
        }
        BSJCount.write("\n");
        FSJCount.write("\n");
        for (String circKey : circRowMap.keySet()) {
            String circId = circKey.replaceFirst("\t", ":").replace("\t", "|");
            if (circTrueIdMap.containsKey(circId)) {
                BSJCount.write(circId);
                FSJCount.write(circId);
                for (int j = 0; j < sampleNameList.size(); j++) {
                    BSJCount.write("\t" + BSJmatrix[circRowMap.get(circKey)][j]);
                    FSJCount.write("\t" + FSJmatrix[circRowMap.get(circKey)][j]);
                }
                BSJCount.write("\n");
                FSJCount.write("\n");
            }
        }
        BSJCount.close();
        FSJCount.close();
        System.out.println(df.format(System.currentTimeMillis()) + " :Matrix files written");

        // Annotation
        if (!annotationFile.equals("F") && geneExonMap.size() > 0) {
            if (intronLable) {
                AnnotationIntron annotation = new AnnotationIntron();
                annotation.annotation(SummaryCircList, geneExonMap, chrExonStartTranscriptMap,
                        chrExonEndTranscriptMap, exonListMap, outputPrefix, relExp);
            } else {
                Annotation annotation = new Annotation();
                annotation.annotation(SummaryCircList, geneExonMap, chrExonStartMap,
                        chrExonEndMap, exonListMap, outputPrefix, relExp);
            }
        }
        System.out.println(df.format(System.currentTimeMillis()) + " :FINALIZE completed");
    }
}
