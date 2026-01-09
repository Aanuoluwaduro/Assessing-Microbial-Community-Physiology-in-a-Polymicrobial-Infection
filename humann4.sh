#!/bin/bash
#####this script uses a reduced prescreen-threshold to include the species with low coverage in our Metaphlan data
#####the metaphlan data was extensively analyzed to see that this threshold allowed the inclusion of the key species in our in vitro samples files
##### This script also shows the collapsing of our raw UniRef ID data to the COG category

echo "Using bowtie2 from: $(which bowtie2)"
echo "Using humann from: $(which humann)"
echo "Start time: $(date)"

mkdir all_samples_humann

for sample in $(cat all_samples.txt); do
    echo "Running HUMAnN on $sample at $(date)"
    
    humann \
        --input ${sample}.fastq.gz \
        --input-format fastq.gz \
        --output all_samples_humann/${sample} \
        --threads 8 \
        --taxonomic-profile all_samples_metaphlan/${sample}_profile.tsv \
        --prescreen-threshold 0.00012 \
        --average-read-length 90 \
        --count-normalization Counts \
        --nucleotide-query-coverage-threshold 80.0 \
        --translated-query-coverage-threshold 80.0 \
        --evalue 1e-3 \
        --nucleotide-database /home/ibbelab/databases/humann4_dbs/chocophlan \
        --protein-database /home/ibbelab/databases/humann4_dbs/uniref \
        --utility-database /home/ibbelab/databases/humann4_dbs/utility_mapping \
        --pathways-database /home/ibbelab/databases/humann4_dbs/utility_mapping/metacyc_reactions_level4ec_only.uniref.bz2,/home/ibbelab/databases/humann4_dbs/utility_mapping/metacyc_pathways_structured_filtered_v24_subreactions
done

echo "All samples finished at: $(date)"

##regrouping to COG
for i in *_genefamilies.tsv; do humann_regroup_table --input ${i} --output /home/ibbelab/aanuoluwa/model_comm/RNAseq/all_samples_humann/gene_families_only/human_cog/cog_${i} --custom /home/ibbelab/databases/humann_dbs/utility_mapping/map_eggnog_uniref90.txt.gz; done
