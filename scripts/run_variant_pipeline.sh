#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# NGS Cancer Biomarkers - Tumour/Normal Variant Pipeline
#
# Input:
#   - HCC1395 normal BAM
#   - HCC1395 tumour BAM
#   - GRCh37 chromosome 17 reference
#
# Output:
#   - Raw variants
#   - Tumour candidate variants
#   - Candidate VCF
# ============================================================


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

NORMAL_BAM="data/raw/HCC1395_normal.bam"
TUMOUR_BAM="data/raw/HCC1395_tumour.bam"

REFERENCE="data/reference/GRCh37_chr17.fa"

VARIANT_DIR="results/variants"

RAW_VCF="${VARIANT_DIR}/raw_variants.vcf.gz"
RENAMED_VCF="${VARIANT_DIR}/raw_variants_renamed.vcf.gz"

VARIANT_TABLE="${VARIANT_DIR}/variants_table.tsv"
CANDIDATE_TABLE="${VARIANT_DIR}/tumour_candidates.tsv"

REGIONS_FILE="${VARIANT_DIR}/candidate_regions.tsv"
CANDIDATE_VCF="${VARIANT_DIR}/tumour_candidates.vcf.gz"


echo "=============================================="
echo " NGS Cancer Biomarkers Variant Pipeline"
echo "=============================================="
echo


# ------------------------------------------------------------
# 2. Check required software
# ------------------------------------------------------------

echo "[1/8] Checking software..."

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools is not installed."
    exit 1
}

command -v bcftools >/dev/null 2>&1 || {
    echo "ERROR: bcftools is not installed."
    exit 1
}

echo "samtools: OK"
echo "bcftools: OK"
echo


# ------------------------------------------------------------
# 3. Check BAM files
# ------------------------------------------------------------

echo "[2/8] Checking BAM files..."

samtools quickcheck -v "$NORMAL_BAM"
samtools quickcheck -v "$TUMOUR_BAM"

echo "BAM files passed quickcheck."
echo


# ------------------------------------------------------------
# 4. Create BAM and reference indexes if necessary
# ------------------------------------------------------------

echo "[3/8] Creating indexes..."

if [ ! -f "${NORMAL_BAM}.bai" ]; then
    samtools index "$NORMAL_BAM"
fi

if [ ! -f "${TUMOUR_BAM}.bai" ]; then
    samtools index "$TUMOUR_BAM"
fi

if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx "$REFERENCE"
fi

echo "Indexes ready."
echo


# ------------------------------------------------------------
# 5. Variant calling
# ------------------------------------------------------------

echo "[4/8] Calling variants..."

bcftools mpileup \
    -Ou \
    -f "$REFERENCE" \
    -r 17:7000000-8000000 \
    "$NORMAL_BAM" \
    "$TUMOUR_BAM" \
| bcftools call \
    -mv \
    -Oz \
    -o "$RAW_VCF"

bcftools index -f "$RAW_VCF"

echo "Raw variants created."
echo


# ------------------------------------------------------------
# 6. Rename samples
# ------------------------------------------------------------

echo "[5/8] Renaming samples..."

printf "NORMAL\nTUMOUR\n" > data/sample_names.txt

bcftools reheader \
    -s data/sample_names.txt \
    -o "$RENAMED_VCF" \
    "$RAW_VCF"

bcftools index -f "$RENAMED_VCF"

echo "Samples:"
bcftools query -l "$RENAMED_VCF"
echo


# ------------------------------------------------------------
# 7. Create variant table
# ------------------------------------------------------------

echo "[6/8] Creating variant table..."

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%AD]\n' \
    "$RENAMED_VCF" \
    > "$VARIANT_TABLE"

sed -i \
'1i CHROM\tPOS\tREF\tALT\tQUAL\tNORMAL_GT\tNORMAL_AD\tTUMOUR_GT\tTUMOUR_AD' \
"$VARIANT_TABLE"

echo "Variant table created."
echo


# ------------------------------------------------------------
# 8. Filter tumour-specific candidate variants
#
# Criteria:
#   Normal genotype = 0/0
#   Tumour contains alternative allele
#   QUAL >= 30
#   Tumour depth >= 10
#   Tumour ALT reads >= 3
#   Tumour VAF >= 10%
#   Normal VAF <= 5%
# ------------------------------------------------------------

echo "[7/8] Filtering candidate variants..."

awk -F'\t' '
BEGIN {OFS="\t"}

NR==1 {
    print $0, \
          "NORMAL_DP", \
          "NORMAL_ALT_READS", \
          "NORMAL_VAF", \
          "TUMOUR_DP", \
          "TUMOUR_ALT_READS", \
          "TUMOUR_VAF"
    next
}

{
    split($7, n, ",")

    normal_ref = n[1]
    normal_alt = 0

    for (i=2; i<=length(n); i++)
        normal_alt += n[i]

    normal_dp = normal_ref + normal_alt

    normal_vaf = \
        (normal_dp > 0) ? normal_alt / normal_dp : 0


    split($9, t, ",")

    tumour_ref = t[1]
    tumour_alt = 0

    for (i=2; i<=length(t); i++)
        tumour_alt += t[i]

    tumour_dp = tumour_ref + tumour_alt

    tumour_vaf = \
        (tumour_dp > 0) ? tumour_alt / tumour_dp : 0


   keep = ($6 == "0/0" && $8 ~ /[1-9]/ && $5 >= 30 && tumour_dp >= 10 && tumour_alt >= 3 && tumour_vaf >= 0.10 && normal_vaf <= 0.05)

if (keep) {
    print $0, normal_dp, normal_alt, normal_vaf, tumour_dp, tumour_alt, tumour_vaf
}
}
' "$VARIANT_TABLE" > "$CANDIDATE_TABLE"


# ------------------------------------------------------------
# 9. Create candidate VCF
# ------------------------------------------------------------

echo "[8/8] Creating candidate VCF..."

awk -F'\t' \
    'NR>1 {print $1 "\t" $2 "\t" $2}' \
    "$CANDIDATE_TABLE" \
    > "$REGIONS_FILE"

bcftools view \
    -R "$REGIONS_FILE" \
    "$RENAMED_VCF" \
    -Oz \
    -o "$CANDIDATE_VCF"

bcftools index -f "$CANDIDATE_VCF"


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

RAW_COUNT=$(bcftools view -H "$RAW_VCF" | wc -l)

CANDIDATE_COUNT=$(
    bcftools view -H "$CANDIDATE_VCF" | wc -l
)

echo
echo "=============================================="
echo " Pipeline completed"
echo "=============================================="
echo
echo "Raw variants:       $RAW_COUNT"
echo "Candidate variants: $CANDIDATE_COUNT"
echo
echo "Candidate VCF:"
echo "$CANDIDATE_VCF"
echo
echo "Candidate table:"
echo "$CANDIDATE_TABLE"
echo
