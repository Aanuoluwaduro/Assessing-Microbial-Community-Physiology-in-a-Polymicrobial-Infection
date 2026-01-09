#!/bin/bash
#this script uses the mpa_vOct22_CHOCOPhlAnSGB_202403 database to run metaphlan4
#it uses a statq of 0.01 to ensure the inclusion of low coverage species in our models
#this decision was made after using a series of values to detect the presence of the species while reducing the possibility of false positives
#also, we used the 0.1 default for the human samples


#Prevent MetaPhlAn from downloading databases
export METAPHLAN_DOWNLOAD=NO

#Output directory
OUTDIR="all_samples_metaphlan"
mkdir -p "${OUTDIR}"

echo "Starting MetaPhlAn profiling at: $(date)"

for sample in $(cat all_samples.txt); do
    echo "Processing ${sample} at $(date)..."

    metaphlan "${sample}.fastq.gz" \
        --input_type fastq \
        --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
        --bowtie2db /home/ibbelab/databases/metaphlan_databases \
        --read_min_len 22 \
        -t rel_ab_w_read_stats \
        --nproc 8 \
        --min_mapq_val 5 \
        --min_cu_len 2000 \
        --ignore_eukaryotes \
        --ignore_archaea \
        --stat_q 0.01 \
        --force \
        -o "${OUTDIR}/${sample}_profile.tsv" \
        &> "${OUTDIR}/${sample}_log.txt"

    echo "✅ Finished ${sample}"
done

echo "✅ All samples processed. Finished at: $(date)"
