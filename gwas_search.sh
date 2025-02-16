#!/bin/zsh

# Data Inputs Vigna_unguiculata.ASM411807v1.60.chr.gff3 (filtered_assembly.gff3), gwas_regions.tsv (new.tsv), Compara_cowpea_homologs_ara.tsv (new2.tsv)
# Filtering Ensembl dataset for specific columns regarding chromosome number, start and end regions of the gene, and gene id

awk '$3 == "gene" {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $9}' Vigna_unguiculata.ASM411807v1.60.chr.gff3 | \
sed -e 's/ID=gene://g ; s/\;.*//g' > filtered_assembly.gff3

# Refining gwas_search_regions.tsv and Compara for future manipulations

sed -e 's/Vu0//g ; s/Vu//g' gwas_regions.tsv | tail -n +2 > new.tsv 

# Updating new.tsv to include only unique SNPs, storing output as variable, final 

final=$(
awk -F '\t' '{print $3}' new.tsv | sort -u | while read snp; do

max=1
low=""

low=$(awk -F '\t' -v snp="$snp" -v max="$max" '
  $3 == snp {
   if ($11 < max) {
    max = $11;
    low = $0
    }
  }
   END {print low}' new.tsv)
  echo "$low"
done
)

cut -f1,3-8 Compara_cowpea_homologs_ara.tsv | tail -n +2 > new2.tsv

# Creating Variables 

data1="filtered_assembly.gff3"
data3="new2.tsv"
echo "$final" | sort -k3 > final.tsv

# Joining filtered_assembly.gff3 and new2.tsv to return a combined file which includes gene location and gene descriptions/homolog info 

join -1 5 -2 1 <(sort -k5 "$data1") <(sort -k1 "$data3" | sed 's/ /_/g') | sed 's/ /\t/g' > data4.tsv

# Creating a list of SNPs to match to corresponding regions

cut -f2-5 final.tsv > snp_list.tsv

# Creating a While Loop to match corresponding SNPs to Gene IDs

output="combined_snp.tsv"

while IFS=$'\t' read -r gene_id chr gene start end homolog_id homolog_type gene_name species gene_desc protein_coding; do

  snp_value="" 

    while IFS=$'\t' read -r chr_s snp_s start_s end_s; do    
     if [[ "$chr" == "$chr_s" && "$start" -ge "$start_s" && "$end" -le "$end_s" ]]; then
      snp_value="$snp_s"
       break
     fi
    done < snp_list.tsv
  echo -e "$gene_id\t$chr\t$gene\t$start\t$end\t$homolog_id\t$homolog_type\t$gene_name\t$species\t$gene_desc\t$protein_coding\t$snp_value" >> "$output"
done < data4.tsv

# Removing records which do not contain a SNP match, since they are not of interest for GWAS analysis

awk -F '\t' '$12 != ""{print $0}' combined_snp.tsv > tmp.tsv

# Joining corresponding SNP info via matching with SNP IDs 

join -1 12 -2 3 <(sort -k12 tmp.tsv) <(sort -k3 final.tsv) | sed 's/ /\t/g' > snp_new.tsv

# Trimming Dataset for only necessary fields & adding header

echo -e "snp\tgene_id\tchr\tstart_gene\tend_gene\thomolog_id\thomolog_type\thomolog_name\tgene_desc\tprotein_coding\tid\tsnp_start\tsnp_end\tbp\tallele1\tallele0\taf\tbeta\tp_lrt" > snp_gwas.tsv
cut -f1-3,5-9,11-13,15-22 snp_new.tsv | sort -u >> snp_gwas.tsv

# Removing Unnecessary output files created by this script
rm -r snp_new.tsv tmp.tsv combined_snp.tsv snp_list.tsv final.tsv new*.tsv filtered_assembly.gff3 data4.tsv
