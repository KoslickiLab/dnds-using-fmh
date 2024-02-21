#find does not include the entire pathway
sort gtdb_proteins_aa_rep_list.txt >> gtdb_proteins_aa_rep_list_sorted.txt
sort gtdb_proteins_nt_rep_list.txt >> gtdb_proteins_nt_rep_list_sorted.txt
cat gtdb_proteins_aa_rep_list_sorted.txt | cut -d'/' -f4 | sed "s/_protein.*//" >> gtdb_proteins_rep_list_sort_names.txt

awk 'NR==FNR{a[NR]=$0; next} {print a[FNR] ",", $0}' gtdb_proteins_nt_rep_list_sorted.txt gtdb_proteins_aa_rep_list_sorted.txt | sed "s/ //" >> gtdb_proteins_rep_list_concat.txt

awk 'NR==FNR{a[NR]=$0; next} {print a[FNR] ",", $0}' gtdb_proteins_rep_list_sort_names.txt gtdb_proteins_rep_list_concat.txt | sed "s/ //" >> dataset.csv

sed -i '1s/^/name,genome_filename,protein_filename\n/' dataset.csv
