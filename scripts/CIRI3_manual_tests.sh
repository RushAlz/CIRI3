# ============================================================================
# CIRI3 Test Data Analysis (Shell Script)
# ============================================================================
cd /pastel/Github_scripts/CIRI3 # Git clone
THREADS=64
REF_FASTA=test_decoupled_manual/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
GTF_FILE=test_decoupled_manual/ref/gencode.v32.primary_assembly.annotation.gtf
GENOMEDIR=${PWD}
FASTQDIR=${PWD}/test_decoupled_manual/fastq
BWADIR=${PWD}/test_decoupled_manual/bwa
STARDIR=${PWD}/test_decoupled_manual/star

# git clone https://github.com/gyjames/CIRI3.git
CIRI3_JAR_PATH="/pastel/tools/circRNA_tools/CIRI3/CIRI3_Java_18.0.1.jar"

# Decoupled processing
DECOUPLED_JAR="/pastel/Github_scripts/CIRI3/CIRI3_decoupled.jar" # scripts/build_jar.sh

# List of sample prefixes
SAMPLES=(
  "Div_100_S91"
  "Div_101_S92"
  "PARDOS_1_S1"
  "PARDOS_2_S2"
)

### SKIP if already processed
# # BWA only
# # conda create -n bwa_env bwa
# conda activate bwa_env
# 
# # Index reference
# bwa index -a bwtsw ${REF_FASTA}
# 
# for SAMPLE in "${SAMPLES[@]}"; do
#   R1="${FASTQDIR}/${SAMPLE}_R1_001.fastq.gz"
#   R2="${FASTQDIR}/${SAMPLE}_R2_001.fastq.gz"
#   bwa mem -t ${THREADS} -T 19 ${REF_FASTA} ${R1} ${R2} > ${BWADIR}/${SAMPLE}.sam
# done
# 
# # git clone https://github.com/gyjames/CIRI3.git
# # cd CIRI3
# # conda env create -n CIRI3 -f ./environment.yaml
# conda activate CIRI3
# 
# # STAR
# # Index reference
# mkdir -p ${STARDIR}/STAR
# STAR --runThreadN ${THREADS} \
# --runMode genomeGenerate \
# --genomeDir ${STARDIR}/STAR \
# --genomeFastaFiles \
# ${REF_FASTA} \
# --sjdbGTFfile ${GTF_FILE}
# 
# for SAMPLE in "${SAMPLES[@]}"; do
#   R1="${FASTQDIR}/${SAMPLE}_R1_001.fastq.gz"
#   R2="${FASTQDIR}/${SAMPLE}_R2_001.fastq.gz"
#   
#   STAR --runThreadN ${THREADS} \
#   --genomeDir ${STARDIR}/STAR \
#   --outSAMtype SAM \
#   --readFilesIn ${R1} ${R2} \
#   --readFilesCommand zcat \
#   --outFileNamePrefix ${STARDIR}/STAR_output_${SAMPLE}/ \
#   --outReadsUnmapped Fastx \
#   --outSJfilterOverhangMin 15 12 12 12 \
#   --alignSJoverhangMin 15 \
#   --alignSJDBoverhangMin 15 \
#   --outFilterMultimapNmax 20 \
#   --outFilterScoreMin 1 \
#   --outFilterMatchNmin 1 \
#   --outFilterMismatchNmax 2 \
#   --chimSegmentMin 15 \
#   --chimScoreMin 15 \
#   --chimJunctionOverhangMin 15
# done
# 
# conda deactivate
# conda activate bwa_env
# 
# for SAMPLE in "${SAMPLES[@]}"; do
#   bwa mem -t ${THREADS} \
#   -T 19 \
#   ${REF_FASTA} \
#   ${STARDIR}/STAR_output_${SAMPLE}/Unmapped.out.mate1 \
#   ${STARDIR}/STAR_output_${SAMPLE}/Unmapped.out.mate2 > ${STARDIR}/STAR_output_${SAMPLE}/bwa.sam
# done
# 
# conda deactivate

conda activate CIRI3

# Processing CIRI3 by sample (single-bam)
mkdir -p test_decoupled_manual/single_bam
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${STARDIR}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${STARDIR}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${STARDIR}/STAR_output_${SAMPLE}/bwa.sam
  
  $CONDA_PREFIX/bin/java -jar ${CIRI3_JAR_PATH} \
    -I "${ChimericOutJunction},${AlignedOutSam},${BWASam}" \
    -O test_decoupled_manual/single_bam/${SAMPLE}.CIRI3results.txt \
    -F ${REF_FASTA} \
    -A ${GTF_FILE} \
    -Ma 1 \
    -W 0 \
    -T ${THREADS}
done

# Joint processing (multi-bam)
mkdir -p test_decoupled_manual/multi_bam
cat /dev/null > test_decoupled_manual/multi_bam/my_samples.tsv
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${STARDIR}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${STARDIR}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${STARDIR}/STAR_output_${SAMPLE}/bwa.sam
  echo "${ChimericOutJunction},${AlignedOutSam},${BWASam}" >> test_decoupled_manual/multi_bam/my_samples.tsv
done

$CONDA_PREFIX/bin/java -jar ${CIRI3_JAR_PATH} \
  -I test_decoupled_manual/multi_bam/my_samples.tsv \
  -O test_decoupled_manual/multi_bam/CIRI3_joint_results.txt \
  -F ${REF_FASTA} \
  -W 1 \
  -Ma 1 \
  -T ${THREADS}

##### Decoupled processing
# SCAN1
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${STARDIR}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${STARDIR}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${STARDIR}/STAR_output_${SAMPLE}/bwa.sam
  
  mkdir -p test_decoupled_manual/${SAMPLE}.SCAN1
  $CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} SCAN1 \
  -I "${ChimericOutJunction},${AlignedOutSam},${BWASam},${SAMPLE}" \
  -O test_decoupled_manual/${SAMPLE}.SCAN1/ \
  -F ${REF_FASTA} \
  -A ${GTF_FILE} \
  -Ma 1 \
  -W 0 \
  -T ${THREADS}
done

cat /dev/null > samples_scan1.tsv
for SAMPLE in "${SAMPLES[@]}"; do
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  echo "${BWASam}\t${SAMPLE}.SCAN1.scan1_meta" >> samples_scan1.tsv
done

# BUILD_UNIVERSE
$CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} BUILD_UNIVERSE \
-I  samples_scan1.tsv \
-F  ${REF_FASTA} \
-O  decoupled_test

# SCAN2
for SAMPLE in "${SAMPLES[@]}"; do
  ChimericOutJunction=${PWD}/STAR_output_${SAMPLE}/Chimeric.out.junction
  AlignedOutSam=${PWD}/STAR_output_${SAMPLE}/Aligned.out.sam
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  
  $CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} SCAN2 \
  -I "${ChimericOutJunction},${AlignedOutSam},${BWASam}" \
  -CU decoupled_test.universe \
  -O ${SAMPLE}.SCAN2 \
  -F ${REF_FASTA} \
  -Ma 1 \
  -T ${THREADS}
done

# FINALIZE
cat /dev/null > finalize_samples.tsv
# /abs/path/sample1.sam   /abs/path/sample1.fsj_counts   8   sample1
# /abs/path/sample2.sam   /abs/path/sample2.fsj_counts   8   sample2

for SAMPLE in "${SAMPLES[@]}"; do
  BWASam=${PWD}/STAR_output_${SAMPLE}/bwa.sam
  FSJ_counts=${PWD}/${SAMPLE}.SCAN2.fsj_counts
  echo "${BWASam}\t${FSJ_counts}\t8\t${SAMPLE}" >> finalize_samples.tsv
done


$CONDA_PREFIX/bin/java -jar ${DECOUPLED_JAR} FINALIZE \
-I  samples_scan1.tsv \
-F  ${REF_FASTA} \
-O  decoupled_test \
-A ${GTF_FILE}
