# dnds
describe

# to create ground truth file run

python step0_produce_dNdS_ground_truth.py --reference_input /data/jzr5814/sourmash_dnds_estimation/tests/data/ground_truth_data/ground_truth_ref.fna --ground_truth_output /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/dNdS_ground_truth.csv --mutation_rate_p 0.1,0.15,0.2 --ground_truth_ref_output /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/ground_truth_ref_used_for_dNdS.fna --ground_truth_queries_output /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/ground_truth_mutated_queries_used_for_dNdS.fna --iterations 100

## To run
bash step1_sourmash_compare_nt.sh
bash step1_sourmash_compare_prot.sh
bash step2_estimate_dnds.sh

