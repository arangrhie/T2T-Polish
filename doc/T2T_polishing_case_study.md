## Telomere-2-Telomere polishing workflow

----

### CHM13 chr20 case-study
Polishing of CHM13 chr20 with T2T polishing pipeline

##### Step 1: Download and prepare input data (Requirement: Docker and Samtools)
```bash
BASE="${HOME}/t2t-polishing-case-study"

# Set up input data
INPUT_DIR="${BASE}/input/data"
REF="chm13.draft_v0.9.chr20.fasta"
HIFI_BAM="chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam"
ONT_BAM="chm13.draft_v0.9.ont_guppy_3.6.0.wm_1.11.chr20.bam"
ILMN_BAM="chm13.draft_v0.9.pcrfree.chr20.bam"

# Set the number of CPUs to use
THREADS="64"

# Set up output directory
OUTPUT_DIR="${BASE}/output"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.chr20.fasta
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.chr20.fasta.fai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.ont_guppy_3.6.0.wm_1.11.chr20.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.ont_guppy_3.6.0.wm_1.11.chr20.bam.bai
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.pcrfree.chr20.bam
wget -P ${INPUT_DIR} https://storage.googleapis.com/pepper-deepvariant-public/T2T_case_study/chm13.draft_v0.9.pcrfree.chr20.bam.bai

samtools merge -@"${THREADS}" "${INPUT_DIR}"/chm13.draft_v0.9.illumina_hifi_hybrid.chr20.bam "${INPUT_DIR}"/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam "${INPUT_DIR}"/chm13.draft_v0.9.pcrfree.chr20.bam
samtools index -@"${THREADS}" "${INPUT_DIR}"/chm13.draft_v0.9.illumina_hifi_hybrid.chr20.bam
```

##### Step 3: Find small errors
```bash
# DeepVariant Hybrid Illumina + PacBio
mkdir -p "${OUTPUT_DIR}"/intermediate_results_dir

BIN_VERSION="1.2.0"
OUTPUT_VCF="CHM13_DeepVariant_hybrid.chr20.vcf.gz"

sudo docker pull google/deepvariant:"${BIN_VERSION}"

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
google/deepvariant:"${BIN_VERSION}" \
/opt/deepvariant/bin/run_deepvariant \
--model_type "HYBRID_PACBIO_ILLUMINA" \
--ref "${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta \
--reads "${INPUT_DIR}"/chm13.draft_v0.9.pcrfree.chr20.bam \
--output_vcf "${OUTPUT_DIR}"/"${OUTPUT_VCF}" \
--num_shards "${THREADS}" \
--regions chr20 \
--intermediate_results_dir "${OUTPUT_DIR}"/intermediate_results_dir

# PEPPER-DeepVariant ONT
sudo docker pull kishwars/pepper_deepvariant:r0.4

OUTPUT_PREFIX="CHM13_PEPPER_DeepVariant.chr20"
OUTPUT_VCF_ONT="CHM13_PEPPER_DeepVariant.chr20.vcf.gz"

# Run PEPPER-Margin-DeepVariant
sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}"/chm13.draft_v0.9.ont_guppy_3.6.0.wm_1.11.chr20.bam \
-f "${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta \
-o "${OUTPUT_DIR}"/pepper_deepvariant_output \
-p "${OUTPUT_PREFIX}" \
-t "${THREADS}" \
--ont

cp "${OUTPUT_DIR}"/pepper_deepvariant_output/"${OUTPUT_VCF_ONT}" "${OUTPUT_DIR}"/
tabix -p vcf "${OUTPUT_DIR}"/"${OUTPUT_VCF_ONT}"

# Variant filtration and merge
SMALL_VARIANTS_HYBRID=CHM13_DeepVariant_hybrid.chr20.filtered.vcf.gz
SMALL_VARIANTS_ONT=CHM13_PEPPER_DeepVariant.chr20.filtered.vcf.gz

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools view -f "PASS" -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=30' -Oz "${OUTPUT_DIR}"/"${OUTPUT_VCF}" > "${OUTPUT_DIR}"/"${SMALL_VARIANTS_HYBRID}"

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools view -f "PASS" -v indels -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=25' -Oz "${OUTPUT_DIR}"/"${OUTPUT_VCF_ONT}" > "${OUTPUT_DIR}"/"${SMALL_VARIANTS_ONT}"

# Merge the small variants
HAPPY_OUTPUT_DIR="${OUTPUT_DIR}"/happy_output
HAPPY_OUTPUT_PREFIX="HAPPY_OUTPUT"
HAPPY_OUTPUT_VCF="${HAPPY_OUTPUT_DIR}"/"${HAPPY_OUTPUT_PREFIX}".vcf.gz
mkdir -p "${HAPPY_OUTPUT_DIR}"

sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
"${OUTPUT_DIR}"/"${SMALL_VARIANTS_HYBRID}" \
"${OUTPUT_DIR}"/"${SMALL_VARIANTS_ONT}" \
-r "${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta \
-o "${HAPPY_OUTPUT_DIR}"/"${HAPPY_OUTPUT_PREFIX}" \
--pass-only \
--engine=vcfeval \
--threads="${THREADS}"

MERGE_VCF_OUTPUT="${OUTPUT_DIR}"/CHM13_MERGED_SMALL_VARIANTS.vcf.gz
sudo docker pull kishwars/t2t_polishing:0.1

sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/t2t_polishing:0.1 \
python3 vcf_merge_t2t.py \
-v1 "${OUTPUT_DIR}"/"${SMALL_VARIANTS_HYBRID}" \
-v2 "${OUTPUT_DIR}"/"${SMALL_VARIANTS_ONT}" \
-hv "${HAPPY_OUTPUT_VCF}" \
-o "${MERGE_VCF_OUTPUT}"

tabix -p vcf "${MERGE_VCF_OUTPUT}"

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools stats "${MERGE_VCF_OUTPUT}"
# this will show we have found 79 records that are potential errors in the assembly
```
