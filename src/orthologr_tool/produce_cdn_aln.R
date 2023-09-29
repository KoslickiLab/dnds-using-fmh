library(orthologr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
    } 

## program...
protein_aln = args[1] # string of file name like '/data/jzr5814/orthologr_analysis/HIT000324409_pairwise/sequences_with_no_line_breaks.aln'
fna = args[2] # string of file name like '/data/jzr5814/orthologr_analysis/HIT000324409_pairwise/sequences_with_no_line_breaks_ref.fna'
output = args[3] #output dnds table

cdn_aln <- codon_aln(file_aln=protein_aln,
    file_nuc=fna,
    format="clustal",
    tool="pal2nal",
    get_aln=TRUE
)

cdn_aln

#dnds <- dNdS(query_file      = query,
#     subject_file    = ref,
#     delete_corrupt_cds = FALSE, # coding sequences that cannot be divided by 3 (triplets) will be removed
#     ortho_detection = "BRH", # perform BLAST best reciprocal hit orthology inference
#     aa_aln_type     = "pairwise", # perform pairwise global alignments of AA seqs 
#     aa_aln_tool     = "NW", # using Needleman-Wunsch
#     codon_aln_tool  = "pal2nal", # perform codon alignments using the tool Pal2Nal
#     dnds_est.method = "NG", # use NG method for dN/dS inference
#     comp_cores      = 1 )

#write.csv(data.frame    (dnds),output,row.names=FALSE,quote=FALSE)
