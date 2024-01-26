from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


p=0.01
dir='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01'
selection='positive'
file=f'{selection}_selection_queries_10002_0.01.fna'
#input_file = open(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/{select}_selection_queries_10002_{p}.fna")
input_file = open(f"{dir}/{file}")

my_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

ref_seq = my_dict['ref_gene'].seq
ref_translated = ref_seq.translate()

report_dict={}

for key in my_dict:
#    if my_dict[key].seq != ref_seq:
    query_seq = my_dict[key].seq
    query_translated = query_seq.translate()
    nt_diffs=0
    aa_diffs=0
    for nuc_ref_pos in range(len(ref_seq)):
        if ref_seq[nuc_ref_pos] != query_seq[nuc_ref_pos]:
            nt_diffs+=1
    for aa_ref_pos in range(len(ref_translated)):
        if ref_translated[aa_ref_pos] != query_translated[aa_ref_pos]:
            aa_diffs+=1
    ANI = (len(ref_seq)-nt_diffs)/len(ref_seq)
    AAI = (len(ref_translated)-aa_diffs)/len(ref_translated)
    if key not in report_dict:
        report_dict[key] = {}
        report_dict[key]['ANI'] = ANI
        report_dict[key]['AAI'] = AAI

        PdN = 1-AAI
        report_dict[key]['PdN'] = PdN

        PdS = AAI-(ANI)**3
        report_dict[key]['PdS'] = PdS

        try:
            dNdS = PdN/PdS*(0.77/2.23)
        except:
            dNdS = None

        report_dict[key]['dNdS'] = dNdS

df = pd.DataFrame(report_dict).T
#df.to_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/ground_truth/{select}_{p}.csv")
df.to_csv(f"{dir}/ground_truth/ground_truth_{p}_{selection}.csv")