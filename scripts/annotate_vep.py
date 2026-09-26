import csv
import json
import urllib.request
from pathlib import Path


# --------------------------------------------------
# Files
# --------------------------------------------------

input_file = Path("results/variants/tumour_candidates.tsv")
output_file = Path("results/annotation/annotated_candidates.tsv")


# --------------------------------------------------
# Read candidate variants
# --------------------------------------------------

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)


# --------------------------------------------------
# Convert variants to VEP input format
#
# Example:
# 17 7578406 . C T . . .
# --------------------------------------------------

variants = []

for row in rows:
    variant = (
        f"{row['CHROM']} "
        f"{row['POS']} "
        f". "
        f"{row['REF']} "
        f"{row['ALT']} "
        f". . ."
    )

    variants.append(variant)


# --------------------------------------------------
# Send all variants to Ensembl VEP GRCh37
# --------------------------------------------------

url = (
    "https://grch37.rest.ensembl.org/"
    "vep/homo_sapiens/region"
    "?hgvs=1&canonical=1&pick=1"
)

payload = json.dumps({
    "variants": variants
}).encode("utf-8")

request = urllib.request.Request(
    url,
    data=payload,
    headers={
        "Content-Type": "application/json",
        "Accept": "application/json"
    },
    method="POST"
)

with urllib.request.urlopen(request, timeout=60) as response:
    annotations = json.load(response)


# --------------------------------------------------
# Match VEP results with our variants
# --------------------------------------------------

annotation_by_input = {}

for annotation in annotations:
    key = " ".join(annotation["input"].split())
    annotation_by_input[key] = annotation


# --------------------------------------------------
# Annotation columns
# --------------------------------------------------

annotation_columns = [
    "GENE",
    "GENE_ID",
    "TRANSCRIPT",
    "CONSEQUENCE",
    "HGVSC",
    "HGVSP",
    "AMINO_ACIDS",
    "PROTEIN_POSITION"
]


# --------------------------------------------------
# Create final annotated table
# --------------------------------------------------

with output_file.open("w", newline="") as f:

    fieldnames = list(rows[0].keys()) + annotation_columns

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    final_rows = []

    for row, variant in zip(rows, variants):

        key = " ".join(variant.split())
        annotation = annotation_by_input.get(key, {})

        transcripts = annotation.get(
            "transcript_consequences",
            []
        )

        transcript = transcripts[0] if transcripts else {}

        row["GENE"] = transcript.get("gene_symbol", "")
        row["GENE_ID"] = transcript.get("gene_id", "")
        row["TRANSCRIPT"] = transcript.get("transcript_id", "")

        row["CONSEQUENCE"] = ",".join(
            transcript.get(
                "consequence_terms",
                [annotation.get("most_severe_consequence", "")]
            )
        )

        row["HGVSC"] = transcript.get("hgvsc", "")
        row["HGVSP"] = transcript.get("hgvsp", "")
        row["AMINO_ACIDS"] = transcript.get("amino_acids", "")
        row["PROTEIN_POSITION"] = transcript.get("protein_start", "")

        writer.writerow(row)

        final_rows.append(row)


# --------------------------------------------------
# Simple summary
# --------------------------------------------------

print()
print("Annotation completed")
print("--------------------")

for row in final_rows:

    position = (
        f"{row['CHROM']}:{row['POS']} "
        f"{row['REF']}>{row['ALT']}"
    )

    gene = row["GENE"] or "-"
    consequence = row["CONSEQUENCE"] or "-"
    protein = row["HGVSP"] or "-"

    print(
        f"{position:25} "
        f"{gene:12} "
        f"{consequence:25} "
        f"{protein}"
    )

print()
print("Output:")
print(output_file)
