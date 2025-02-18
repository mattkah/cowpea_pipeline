#!/bin/zsh
 
# Data Inputs Vigna_unguiculata.ASM411807v1.60.chr.gff3 (filtered_assembly.gff3), gwas_regions.tsv (new.tsv), Compara_cowpea_homologs_ara.tsv (new2.tsv)
# Filtering Ensembl dataset for specific columns regarding chromosome number, start and end regions of the gene, and gene id

awk '$3 == "gene" {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $9}' Vigna_unguiculata.ASM411807v1.60.chr.gff3 | \
sed -e 's/ID=gene://g ; s/\;.*//g' > filtered_assembly.gff3

# Refining gwas_search_regions.tsv and Compara for future manipulations

sed -e 's/Vu0//g ; s/Vu//g' rda_data.tsv | tail -n +2 > new.tsv
       
cut -f1,3-8 Compara_cowpea_homologs_ara.tsv | tail -n +2 > new2.tsv
  
# Creating Variables

data1="filtered_assembly.gff3"
data2="new.tsv"
data3="new2.tsv"

# Joining filtered_assembly.gff3 and new2.tsv to return a combined file which includes gene location and gene descriptions/homolog info

join -1 5 -2 1 <(sort -k5 "$data1") <(sort -k1 "$data3" | sed 's/ /_/g') | sed 's/ /\t/g' > data4.tsv

# Creating a list of SNPs to match to corresponding regions 

cut -f22-24,26 "$data2" > snp_list.tsv

# Creating a While Loop to match corresponding SNPs to Gene IDs

output="combined_snp.tsv"

while IFS=$'\t' read -r gene_id chr gene start end homolog_id homolog_type gene_name species gene_desc protein_coding; do
  snp_value=""
   
    while IFS=$'\t' read -r chr_s start_s end_s snp_s; do
     if [[ "$chr" == "$chr_s" && "$start" -ge "$start_s" && "$end" -le "$end_s" ]]; then
      snp_value="$snp_s"
       break
     fi
    done < snp_list.tsv
  echo -e "$gene_id\t$chr\t$gene\t$start\t$end\t$homolog_id\t$homolog_type\t$gene_name\t$species\t$gene_desc\t$protein_coding\t$snp_value" >> "$output"
done < data4.tsv

# Removing records which do not contain a SNP match, since they are not of interest for GWAS analysis

awk -F '\t' '$12 != ""' combined_snp.tsv > tmp.tsv

# Joining corresponding SNP info via matching with SNP IDs

join -1 12 -2 26 <(sort -k12 tmp.tsv) <(sort -k26 "$data2") | sed 's/ /\t/g' > snp_new.tsv

# Trimming Dataset for only necessary fields & adding header

echo -e "snp_id\tgene_id\tchr\tgene_start\tgene_end\thomolog_id\thomolog_type\tortho_name\tspecies\tgene_desc\tprotein_coding\tsn\tsnp\tRDA1\tRDA2\tRDA3\tRDA4\tRDA5\taxis\
tloading\tvpdav01\ttmin_02\ttmax_08\tbio_3\tprec_05\taridity02\taridity08\tprec_12\tbio_19\trhumavsep5\tpredictor\tcorrelation\tsnp_start\tsnp_end\tpos\tallele0\tallele1" > rda_output.tsv
cut -f1-3,5-9,11-33,35-39 snp_new.tsv | sort -u >> rda_output.tsv

# Removing Unnecessary output files created by this script
rm -r snp_new.tsv tmp.tsv combined_snp.tsv snp_list.tsv final.tsv new*.tsv filtered_assembly.gff3 data4.tsv


