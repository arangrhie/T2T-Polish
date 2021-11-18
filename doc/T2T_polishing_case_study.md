## Telomere-2-Telomere polishing workflow

----

### CHM13 chr20 case-study
This case-study is an example of T2T polishing scheme on `CHM13 chr20`. You can copy-paste the commands from here to generate the final polished assembly. This case-study focuses on finding SV and SNV-like errors in the assembly.

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
##### Step 2: Find SV-like errors
```bash
# Run parliament2 (Illumina)
PARLIAMENT_OUTPUT_PREFIX="PARLIAMENT_SV_OUTPUT"
PARLIAMENT_OUTPUT="${OUTPUT_DIR}"/"PARLIAMENT_SV_OUTPUT.combined.genotyped.vcf"
sudo docker pull dnanexus/parliament2:latest

sudo docker run \
-v "${INPUT_DIR}":/home/dnanexus/in \
-v "${OUTPUT_DIR}":/home/dnanexus/out \
dnanexus/parliament2:latest \
--bam chm13.draft_v0.9.pcrfree.chr20.bam \
-r chm13.draft_v0.9.chr20.fasta \
--prefix "${PARLIAMENT_OUTPUT_PREFIX}" \
--bai chm13.draft_v0.9.pcrfree.chr20.bam.bai \
--fai chm13.draft_v0.9.chr20.fasta.fai \
--filter_short_contigs --breakdancer --breakseq --manta --cnvnator \
--lumpy --delly_deletion --genotype --svviz_only_validated_candidates

# Run sniffles (HiFi, ONT)
sudo docker pull quay.io/biocontainers/sniffles:1.0.12--h8b12597_1

# Sniffles on HiFi reads
HIFI_SVS_SNIFFLES="chm13_chr20_sniffles_hifi_svs.vcf"

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
quay.io/biocontainers/sniffles:1.0.12--h8b12597_1 \
sniffles \
-m "${INPUT_DIR}"/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam \
-v "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES}" \
-d 500 -n -1 -s 3

# Sniffles on ONT reads
ONT_SVS_SNIFFLES="chm13_chr20_sniffles_ont_svs.vcf"

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
quay.io/biocontainers/sniffles:1.0.12--h8b12597_1 \
sniffles \
-m "${INPUT_DIR}"/chm13.draft_v0.9.ont_guppy_3.6.0.wm_1.11.chr20.bam \
-v "${OUTPUT_DIR}"/"${ONT_SVS_SNIFFLES}" \
-d 500 -n -1 -s 3

ONT_SVS_SNIFFLES_FILTERED="chm13_chr20_sniffles_ont_svs.filtered.vcf"
HIFI_SVS_SNIFFLES_FILTERED="chm13_chr20_sniffles_hifi_svs.filtered.vcf"
# filter the SV calls
sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/t2t_polishing:0.1 \
python3 filter.py "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES}" > "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED}"

sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/t2t_polishing:0.1 \
python3 filter.py "${OUTPUT_DIR}"/"${ONT_SVS_SNIFFLES}" > "${OUTPUT_DIR}"/"${ONT_SVS_SNIFFLES_FILTERED}"

# refine SVs with IRIS (HiFi)
IRIS_OUTPUT_DIR="IRIS_OUT"
HIFI_SVS_SNIFFLES_FILTERED_IRIS="chm13_chr20_sniffles_hifi_svs.filtered.iris.vcf"

sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
quay.io/biocontainers/irissv:1.0.4--hdfd78af_2 \
iris --hifi --keep_long_variants --keep_files \
genome_in="${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta \
reads_in="${INPUT_DIR}"/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam \
vcf_in="${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED}" \
vcf_out="${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED_IRIS}" \
out_dir="${OUTPUT_DIR}"/"${IRIS_OUT}"

# As we found no SVs after filtering this file will be empty
echo "${PARLIAMENT_OUTPUT}" > "${OUTPUT_DIR}"/SV_filelist.txt
echo "${OUTPUT_DIR}"/"${ONT_SVS_SNIFFLES_FILTERED}" >> "${OUTPUT_DIR}"/SV_filelist.txt
echo "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED_IRIS}" >> "${OUTPUT_DIR}"/SV_filelist.txt

sudo docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
quay.io/biocontainers/jasminesv:1.1.4--hdfd78af_0 \
jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 --output_genotypes \
file_list="${OUTPUT_DIR}"/SV_filelist.txt out_file="${OUTPUT_DIR}"/"SV_like_errors.vcf"
```
After merging we have the finale VCF: `"${OUTPUT_DIR}"/SV_like_errors.vcf` which contains all the errors. We then manually inspected all of these variants and confirmed all of them are either FPs or true hets. So in this case we did not find any SV-like errors.

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
bcftools view -f "PASS" -V indels -e 'FORMAT/VAF<=0.5 | FORMAT/GQ<=25' -Oz "${OUTPUT_DIR}"/"${OUTPUT_VCF_ONT}" > "${OUTPUT_DIR}"/"${SMALL_VARIANTS_ONT}"

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

##### Step 4: Run merfin and generate a polished assembly
```bash
cd "${HOME}"
git clone https://github.com/arangrhie/merfin.git
cd merfin/src
make -j 12
cd "${HOME}"

wget -P ${INPUT_DIR} https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/MERFIN_2021/chm13/evaluation/chm13.k21.gt1.meryl.tar.gz
wget -P ${INPUT_DIR} https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/MERFIN_2021/chm13/evaluation/lookup_table.k21.txt

tar -xvf "${INPUT_DIR}"/"chm13.k21.gt1.meryl.tar.gz" -C "${INPUT_DIR}"/

MERGE_VCF_OUTPUT="${OUTPUT_DIR}"/CHM13_MERGED_SMALL_VARIANTS.vcf.gz
gunzip "${MERGE_VCF_OUTPUT}"
MERGE_VCF_OUTPUT="${OUTPUT_DIR}"/CHM13_MERGED_SMALL_VARIANTS.vcf

READMER_DB="${INPUT_DIR}"/"chm13.k21.gt1.meryl"
LOOKUP_TABLE="${INPUT_DIR}"/"lookup_table.k21.txt"
MERFIN_OUTPUT="${OUTPUT_DIR}"/"CHM13_MERGED_SMALL_VARIANTS.merfin"

"${HOME}"/merfin/build/bin/merfin -polish \
-vcf "${MERGE_VCF_OUTPUT}" \
-threads "${THREADS}" \
-sequence "${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta \
-readmers "${READMER_DB}" \
-prob "${LOOKUP_TABLE}" \
-peak 106.7 \
-output "${MERFIN_OUTPUT}"

MERFIN_OUTPUT="${OUTPUT_DIR}"/"CHM13_MERGED_SMALL_VARIANTS.merfin.polish.vcf"

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools view -Oz "${MERFIN_OUTPUT}" > "${MERFIN_OUTPUT}.gz"

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools index "${MERFIN_OUTPUT}.gz"

POLISHED_FASTA="${OUTPUT_DIR}"/chm13.draft_v0.9.chr20.polished.fasta

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
bcftools consensus \
-f "${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta \
-H 1 "${MERFIN_OUTPUT}.gz" > "${POLISHED_FASTA}"
# Applied 46 variants
```
After this we have `"${OUTPUT_DIR}"/chm13.draft_v0.9.chr20.polished.fasta` which is the polished assembly we have generated.

##### Step 5: Evaluating the assembly quality before and after t2t_polishing
```bash
MERQURY_ASSEMENT_DIR="${OUTPUT_DIR}/merqury_assessement"
mkdir "${MERQURY_ASSEMENT_DIR}"

RAW_ASM_MERQURY_ASSEMENT_DIR="${OUTPUT_DIR}/merqury_assessement/raw_assembly_qv"
mkdir "${RAW_ASM_MERQURY_ASSEMENT_DIR}"

wget -P ${INPUT_DIR} https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/qc/hybrid.k21.meryl.tar.gz
tar -xvf "${INPUT_DIR}"/"hybrid.k21.meryl.tar.gz" -C "${RAW_ASM_MERQURY_ASSEMENT_DIR}"/

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${RAW_ASM_MERQURY_ASSEMENT_DIR}":/data \
juklucas/hpp_merqury:latest \
merqury.sh hybrid.meryl "${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta chm13_chr20_raw_qv

cat "${RAW_ASM_MERQURY_ASSEMENT_DIR}"/chm13_chr20_raw_qv.qv
# The output will show this:
# chm13.draft_v0.9.chr20	170	66210241	69.127	1.22266e-07
# which means the quality of the raw assembly is 69.12 with 170 missing k-mers from the hybrid k-mer db

POLISHED_ASM_MERQURY_ASSEMENT_DIR="${OUTPUT_DIR}/merqury_assessement/polished_assembly_qv"
mkdir "${POLISHED_ASM_MERQURY_ASSEMENT_DIR}"
tar -xvf "${INPUT_DIR}"/"hybrid.k21.meryl.tar.gz" -C "${POLISHED_ASM_MERQURY_ASSEMENT_DIR}"/

sudo docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${POLISHED_ASM_MERQURY_ASSEMENT_DIR}":/data \
juklucas/hpp_merqury:latest \
merqury.sh hybrid.meryl "${OUTPUT_DIR}"/chm13.draft_v0.9.chr20.polished.fasta chm13_chr20_polished_qv

cat "${POLISHED_ASM_MERQURY_ASSEMENT_DIR}"/chm13_chr20_polished_qv.qv
# The output will show this:
# chm13.draft_v0.9.chr20.polished	144	66210231	69.8478	1.03566e-07
# which means the quality of the polished assembly is 69.84 and with 144 missing k-mers
```

In this case-study we improved the raw assembly quality of CHM13 v0.9 from `QV69.12` to `QV69.84`.
