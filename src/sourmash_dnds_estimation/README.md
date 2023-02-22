# dnds
describe

# to create ground truth file run
python step0_produce_dNdS_ground_truth.py --reference_input /data/jzr5814/sourmash_dnds_estimation/tests/data/ground_truth_data/ground_truth_ref.fna --ground_truth_output /data/jzr5814/sourmash_dnds_estimation/tests/data/ground_truth_data/dNdS_ground_truth.csv --mutation_rate_p 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1

## To run
bash step1_sourmash_compare_nt.sh
bash step1_sourmash_compare_prot.sh
bash step2_estimate_dnds.sh

