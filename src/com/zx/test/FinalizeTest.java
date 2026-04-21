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
     * finalizeInputTsv: four-column TSV per sample:
     *   samFile  fsjCountsFile  bsjFile1,bsjFile2,...  sampleName
     *
     * Column 3 is a comma-separated list of the exact BSJ file paths produced by SCAN1
     * (and appended to by SCAN2) for that sample.  Using an explicit list prevents stale
     * files from a previous run leaking into the BSJ matrix or Summary filtering.
     *
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
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :FINALIZE start");

        // Step 1: Read finalize input TSV
        // TSV format: samFile  fsjCountsFile  bsjFile1,bsjFile2,...  sampleName
        ArrayList<String> filePathList = new ArrayList<String>();
        HashMap<String, ArrayList<String>> bsjFilesMap = new HashMap<String, ArrayList<String>>();
        ArrayList<String> fsjCountsFileList = new ArrayList<String>();
        ArrayList<String> sampleNameList = new ArrayList<String>();

        BufferedReader tsvBr = new BufferedReader(new FileReader(new File(finalizeInputTsv)));
        String tsvLine = tsvBr.readLine();
        while (tsvLine != null) {
            if (tsvLine.startsWith("#") || tsvLine.equals("")) {
                tsvLine = tsvBr.readLine();
                continue;
            }
            String[] arr = tsvLine.split("\t");
            String samFilePath = arr[0].trim();
            String fsjCountsPath = arr[1].trim();
            ArrayList<String> bsjPaths = new ArrayList<String>(Arrays.asList(arr[2].trim().split(",")));
            String sampleName = arr[3].trim();
            filePathList.add(samFilePath);
            bsjFilesMap.put(samFilePath, bsjPaths);
            fsjCountsFileList.add(fsjCountsPath);
            sampleNameList.add(sampleName);
            tsvLine = tsvBr.readLine();
        }
        tsvBr.close();
        System.out.println(df.format(System.currentTimeMillis()) + " :Loaded " + filePathList.size() + " samples");

        // Step 2: Load FA file
        ReadFaFile RF = new ReadFaFile();
        RF.readFa(faFile);
        HashMap<String, String> chrTCGAMap = RF.getChrTCGAMap();
        RF = null;
        System.out.println(df.format(System.currentTimeMillis()) + " :Successful import of reference genome files");

        // Step 3: Load annotation (optional)
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

        // Step 4: Read all .fsj_counts files and build merged circFSJMap
        // Seed the universe: prefer the BUILD_UNIVERSE output when supplied, else
        // fall back to the first sample's .fsj_counts for backwards compatibility.
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
        int sampleCount = filePathList.size();
        int[][] FSJmatrix = new int[circNum][sampleCount];
        HashMap<String, Integer> circRowMap = new HashMap<String, Integer>();
        int rowIdx = 0;
        for (String circKey : circFSJMerged.keySet()) {
            circRowMap.put(circKey, rowIdx++);
        }

        // Read all .fsj_counts: populate per-sample FSJ matrix and merged map
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

        // Step 5: Build BSJ matrix from BSJ files (SCAN2 appends additional BSJ reads to the same files).
        // BSJ file paths come from the explicit list in column 3 of the input TSV,
        // not reconstructed from samFile+fileSplitNum, to prevent stale-file leakage.
        int[][] BSJmatrix = new int[circNum][sampleCount];
        for (int i = 0; i < sampleCount; i++) {
            HashMap<String, Integer> circMap = new HashMap<String, Integer>();
            String samFilePath = filePathList.get(i);
            ArrayList<String> bsjPaths = bsjFilesMap.get(samFilePath);
            long totalLines = 0, tagOneLines = 0, unknownCirc = 0;
            int slotIdx = 0;
            for (String bsjPath : bsjPaths) {
                slotIdx++;
                long fileLines = 0, fileTagOne = 0;
                BufferedReader BSJBr = new BufferedReader(new FileReader(new File(bsjPath)));
                String line = BSJBr.readLine();
                while (line != null) {
                    fileLines++;
                    String[] circLineArr = line.split("\t", 7);
                    if (circLineArr[2].equals("1")) {
                        fileTagOne++;
                        String chrStartEnd = circLineArr[3] + "\t" + circLineArr[4] + "\t" + circLineArr[5];
                        if (!circMap.containsKey(chrStartEnd)) {
                            circMap.put(chrStartEnd, 1);
                        } else {
                            circMap.put(chrStartEnd, circMap.get(chrStartEnd) + 1);
                        }
                    }
                    line = BSJBr.readLine();
                }
                BSJBr.close();
                totalLines += fileLines;
                tagOneLines += fileTagOne;
                System.out.println(df.format(System.currentTimeMillis()) + " :DIAG sample[" + i + "] BSJ" + slotIdx + ": total=" + fileLines + " tag1=" + fileTagOne);
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

        // Step 6: Run Summary (using explicit BSJ file lists, not fileSplitNum reconstruction)
        Summary summary = new Summary(strigency, chrTCGAMap);
        ArrayList<String> SummaryCircList = summary.summaryFromBsjFiles(filePathList, bsjFilesMap, circFSJMerged, "");
        HashMap<String, String> circTrueIdMap = summary.getCircMap();
        summary = null;
        chrTCGAMap = null;
        System.out.println(df.format(System.currentTimeMillis()) + " :Summary completed");

        // Write BSJ_Matrix and FSJ_Matrix (identical to MutFileTest lines 509-541)
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

        // Step 7: Annotation
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
