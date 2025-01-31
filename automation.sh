#!/bin/zsh

cd /Users/mattkahler/Desktop/cowpea_test 

# Filtering Ensembl data set for strictly gene regions of chromosomes 

awk '$3 == "gene" && $1 ~ /^[1-7]/{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $9}' Vigna_unguiculata.ASM411807v1.60.chr.gff3 > filtered_assembly.gff3

# Removing the header from cowpea_test.tsv and Compara

tail -n +2 cowpea_test.tsv > new.tsv

tail -n +2 Compara_cowpea_homologs_ara.tsv > new2.tsv

# Input files

data1="new.tsv"
data2="filtered_assembly.gff3"
data3="new2.tsv" 

output_dir="output_data"
mkdir -p "$output_dir"

while IFS=$'\t' read -r line; do
 
    chr1=$(echo "$line" | awk '{print $1}')
    start1=$(echo "$line" | awk '{print $2}')
    end1=$(echo "$line" | awk '{print $3}')
 
   chromosome_dir="${output_dir}/chr${chr1}"
   mkdir -p "$chromosome_dir"

   output_file="${chromosome_dir}/${chr1}_${start1}_${end1}.tsv"

   matched_genes=$(awk -v chr="$chr1" -v start="$start1" -v end="$end1" \
   '$1 == chr && $3 >= start && $4 <= end {print $5}' "$data2")

   if [[ -n "$matched_genes" ]]; then
    for gene_id_mod in $matched_genes; do
      core=$(echo "$gene_id_mod" | awk -F '\t' '{if ($1 ~/ID=gene:/){match($1, /Vigun[^;]+/); print substr($1, RSTART, RLENGTH)}}')

    matching_species=$(grep -w "$core" "$data3")

    if [[ -n "$matching_species" ]]; then
      echo "$matching_species" >> "$output_file"
    fi
  done
 fi
done < "$data1"
