#!/bin/zsh

target="/Users/mattkahler/Desktop/cowpea_test/output_data/"

while true; do
  echo -e "Enter SNP: \c"
  read snp

cd /Users/mattkahler/Desktop/cowpea_test
if grep -wq "$snp" new.tsv; then
  cd output_data/ 
  file="$(find . -type f -iname "*.$snp.tsv")"
  open "$file"
else
  echo "Please try again."
fi
done


