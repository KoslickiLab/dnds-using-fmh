rm /home/grads/jzr5814/data/tests/prefetch.txt
rm /home/grads/jzr5814/data/tests/*.sig*
rm /home/grads/jzr5814/data/tests/results_kmers.csv
rm /home/grads/jzr5814/data/tests/SOURMASH-MANIFEST.csv
rm /home/grads/jzr5814/data/tests/prefetch_res_*.csv
rm test_sketch_ref.txt test_sketch_sample.txt test_prefetch*.txt

ref=../../data/tests/50_sequences.fna 
samples=../../data/tests/5_uniprot_seqs.fasta 
ref_output=../../data/tests/50_sequences.sig
samp_output=../../data/tests/5_uniprot_seqs.sig.zip
ref_scaled=1
samp_scaled=1
tbp=1
wd=../../data/tests/

nohup sourmash sketch protein -p k=7,k=14,k=21,k=28,scaled=$ref_scaled $ref -o $ref_output > test_sketch_ref.txt 2>&1 &

nohup sourmash sketch protein -p k=7,k=14,k=21,k=28,scaled=$samp_scaled $samples -o $samp_output > test_sketch_sample.txt 2>&1 &

for K in 7 14 21 28
do
	echo $K
	nohup sourmash prefetch $ref_output $samp_output --protein --o ${wd}prefetch_res_${K}.csv --threshold-bp $tbp --ksize $K > test_prefetch${K}.txt 2>&1 & #assign to diff log files

done

