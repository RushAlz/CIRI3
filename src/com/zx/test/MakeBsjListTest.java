package com.zx.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;

public class MakeBsjListTest {

    /**
     * Converts legacy or cloud-staged SCAN1 outputs into the three-column BSJ
     * list TSV required by BUILD_UNIVERSE -IB.
     *
     * Input TSV (auto-detected by column count):
     *
     *   1 col:  scan1_meta_path
     *     Use the bsjPrefix, fileSplitNum, and readLen stored in the meta file
     *     as-is.  Suitable when the pipeline runs on a single machine and paths
     *     have not changed since SCAN1 ran.
     *
     *   2 cols: scan1_meta_path <TAB> staged_bsj_prefix
     *     Override the stale bsjPrefix from the meta file with the path where
     *     BSJ files are accessible on the current machine (e.g. after staging
     *     from object storage in a cloud workflow).  fileSplitNum and readLen
     *     are still read from the meta file.
     *
     *   3 cols: sam_file_path <TAB> staged_bsj_prefix <TAB> read_len
     *     For truly legacy outputs produced by the old joint pipeline
     *     (MutFileTest / MutSTARTest) where no .scan1_meta file was written.
     *     fileSplitNum is auto-detected by probing {staged_bsj_prefix}BSJ1,
     *     BSJ2, … until a file is not found.  read_len must be provided
     *     explicitly because it cannot be recovered from the BSJ files.
     *
     * Output: {outputPrefix}.bsj_list.tsv
     *   Three-column TSV ready for BUILD_UNIVERSE -IB:
     *     bsjPrefix <TAB> fileSplitNum <TAB> readLen
     */
    public void make(String inputTsv, String outputPrefix) throws IOException {
        SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        System.out.println(df.format(System.currentTimeMillis()) + " :MAKE_BSJ_LIST start");

        String outPath = outputPrefix + ".bsj_list.tsv";
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outPath)));

        int rowsOut = 0;
        int rowsSkipped = 0;

        BufferedReader br = new BufferedReader(new FileReader(new File(inputTsv)));
        String line = br.readLine();
        while (line != null) {
            if (line.startsWith("#") || line.trim().isEmpty()) {
                line = br.readLine();
                continue;
            }

            String[] cols = line.split("\t");
            int ncols = cols.length;

            String bsjPrefix;
            int fileSplitNum;
            int readLen;

            if (ncols == 1) {
                // 1-col: read everything from scan1_meta
                String metaPath = cols[0].trim();
                int[] meta = readMeta(metaPath);
                if (meta == null) {
                    System.out.println("WARNING: could not parse scan1_meta: " + metaPath + " — skipping");
                    rowsSkipped++;
                    line = br.readLine();
                    continue;
                }
                bsjPrefix = readMetaBsjPrefix(metaPath);
                fileSplitNum = meta[0];
                readLen = meta[1];

            } else if (ncols == 2) {
                // 2-col: bsjPrefix override, rest from scan1_meta
                String metaPath = cols[0].trim();
                bsjPrefix = cols[1].trim();
                int[] meta = readMeta(metaPath);
                if (meta == null) {
                    System.out.println("WARNING: could not parse scan1_meta: " + metaPath + " — skipping");
                    rowsSkipped++;
                    line = br.readLine();
                    continue;
                }
                fileSplitNum = meta[0];
                readLen = meta[1];

            } else if (ncols >= 3) {
                // 3-col: legacy (no scan1_meta); fileSplitNum auto-detected
                // col[0] = sam_file_path (not opened, kept for user reference)
                bsjPrefix = cols[1].trim();
                readLen = Integer.parseInt(cols[2].trim());
                fileSplitNum = autoDetectSplitNum(bsjPrefix);
                if (fileSplitNum == 0) {
                    System.out.println("WARNING: no BSJ files found at " + bsjPrefix + "BSJ1 — skipping row");
                    rowsSkipped++;
                    line = br.readLine();
                    continue;
                }

            } else {
                System.out.println("WARNING: unexpected column count (" + ncols + ") in row: " + line + " — skipping");
                rowsSkipped++;
                line = br.readLine();
                continue;
            }

            bw.write(bsjPrefix + "\t" + fileSplitNum + "\t" + readLen + "\n");
            rowsOut++;
            line = br.readLine();
        }
        br.close();
        bw.close();

        System.out.println(df.format(System.currentTimeMillis()) + " :MAKE_BSJ_LIST done — " + rowsOut + " rows written to " + outPath
                + (rowsSkipped > 0 ? " (" + rowsSkipped + " rows skipped)" : ""));
    }

    /** Returns [fileSplitNum, readLen] from a .scan1_meta file, or null on parse failure. */
    private int[] readMeta(String metaPath) throws IOException {
        int fileSplitNum = 0;
        int readLen = 0;
        BufferedReader mbr = new BufferedReader(new FileReader(new File(metaPath)));
        String ml = mbr.readLine();
        while (ml != null) {
            if (ml.startsWith("fileSplitNum=")) {
                fileSplitNum = Integer.parseInt(ml.split("=", 2)[1].trim());
            } else if (ml.startsWith("readLen=")) {
                readLen = Integer.parseInt(ml.split("=", 2)[1].trim());
            }
            ml = mbr.readLine();
        }
        mbr.close();
        if (fileSplitNum == 0 || readLen == 0) return null;
        return new int[]{fileSplitNum, readLen};
    }

    /** Reads the bsjPrefix field from a .scan1_meta file. */
    private String readMetaBsjPrefix(String metaPath) throws IOException {
        BufferedReader mbr = new BufferedReader(new FileReader(new File(metaPath)));
        String ml = mbr.readLine();
        while (ml != null) {
            if (ml.startsWith("bsjPrefix=")) {
                mbr.close();
                return ml.split("=", 2)[1].trim();
            }
            ml = mbr.readLine();
        }
        mbr.close();
        return null;
    }

    /** Probes {prefix}BSJ1, BSJ2, … and returns the count of contiguous existing files. */
    private int autoDetectSplitNum(String bsjPrefix) {
        int n = 0;
        while (new File(bsjPrefix + "BSJ" + (n + 1)).exists()) {
            n++;
        }
        return n;
    }
}
