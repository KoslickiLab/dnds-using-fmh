#### Gene Level Evaqluations

# dnds

describe

## to create ground truth file run from p rate mutations
python step0_produce_dNdS_ground_truth.py --reference_input /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.001_ksizes_5_30_1000_cdn/ground_truth_ref.fna --ground_truth_output dNdS_ground_truth.csv --mutation_rate_p 0.001 --ground_truth_ref_output ground_truth_ref_used_for_dNdS.fna --ground_truth_queries_output ground_truth_mutated_queries_used_for_dNdS.fna --iterations 100 --wd /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.001_ksizes_5_30_1000_cdn/

python step0_produce_dNdS_ground_truth_between_real_sequences.py --reference_input /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_K12567/ref_K12567.fna --query_input /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_K12567/K12567.fna --ground_truth_output dNdS_ground_truth.csv --wd /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_K12567/

## To run
bash step1_sourmash_compare_nt.sh
bash step1_sourmash_compare_prot.sh
bash step2_estimate_dnds.sh

# concatenate csv files keeping one header
awk 'NR == 1 || FNR > 1'

#### Genome Level Evaqluations

## Delete empty files in a directory

find . -type f -empty -print -delete
