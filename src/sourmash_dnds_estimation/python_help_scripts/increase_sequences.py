import os

#remove new lines in fasta file using bioawk unix command
input = "/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/test_METT_making_sequences_bigger/METT_queries.fna"
output = "out.fna"
cmd="bioawk -c fastx "+"'{ gsub(/\\n/,\"\",seq); print \">\"$name; print $seq }'"+f" {input} > {output}"
os.system(cmd)

#increase each sequence by a factor of 15
file=open(output, "r")
lines = file.readlines()
out=open("/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/test_METT_making_sequences_bigger/METT_queries_longer_seqs.fna",'w')
factor=15
for line in lines:
    new_line=line.rstrip()
    if new_line[0] == '>':
        out.write(new_line+"\n")
    elif new_line[0] != '>':
        for x in range(factor):
            out.write(new_line+"\n")

#remove intermediate files
os.system(f"rm {output}")
