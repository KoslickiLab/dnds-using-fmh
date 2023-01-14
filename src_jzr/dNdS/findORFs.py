"""Module containing code for identifying open reading frames"""

from more_itertools import sliced

def reverse_complement(nt_seq):
    """Function reports the reverse complement of a DNA sequence"""
    complementaryDNA = {"A":"T", "C":"G", "G":"C", "T":"A"}
    reverse = nt_seq[::-1].upper()
    complement = ''
    for base in reverse:
        complement += complementaryDNA[base]
    return(complement)

def frame_cds(nt_seq,frame):
    """Function reports the coding sequence of a reading frame"""
    if frame < 4:
        coding_seq = list(sliced(nt_seq[frame-1:].upper(),3))    
    elif frame > 3:
        sequence = reverse_complement(nt_seq)
        coding_seq = list(sliced(sequence[frame-4:],3))
    return(coding_seq)

def listORFs(cds):
    """Function reports a list of ORFs found in a reading frame"""
    stop_codons = ['TAA', 'TAG', 'TGA']
    ORFs = []
    tempORFs = []
    for codon in cds:
        if codon == 'ATG':
            tempORFs.append(codon)
        elif codon in stop_codons:
            if tempORFs:
                ORFs.append(tempORFs)
            tempORFs = []
        elif codon != 'ATG' and codon not in stop_codons and len(codon) == 3:
            if len(codon) == 3:
                if tempORFs:
                    tempORFs.append(codon)
            else:
                tempORFs = []
    return(ORFs)

def translate_ORFs(cds_seq):
    """Function translates open reading frames using the coding sequence"""

    codontab = {
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',    # Serine
        'TTC': 'F', 'TTT': 'F',                                                    # Phenilalanine
        'TTA': 'L', 'TTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucine
        'TAC': 'Y', 'TAT': 'Y',                                                    # Tirosine
        'TAA': '*', 'TAG': '*', 'TGA': '*',                                        # Stop
        'TGC': 'C', 'TGT': 'C',                                                    # Cisteine
        'TGG': 'W',                                                                # Tryptophane
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',                            # Proline
        'CAC': 'H', 'CAT': 'H',                                                    # Histidine
        'CAA': 'Q', 'CAG': 'Q',                                                    # Glutamine
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',                            # Arginine
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',                                        # Isoleucine
        'ATG': 'M',                                                                # Methionine
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',                            # Threonine
        'AAC': 'N', 'AAT': 'N',                                                    # Asparagine
        'AAA': 'K', 'AAG': 'K',                                                    # Lysine
        'AGA': 'R', 'AGG': 'R',                                                    # Arginine
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',                            # Valine
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',                            # Alanine
        'GAC': 'D', 'GAT': 'D',                                                    # Aspartic Acid
        'GAA': 'E', 'GAG': 'E',                                                    # Glutamic Acid
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'                             # Glycine
    }

    translation=[]
    for codon in cds_seq:
        if codon in codontab:
            translation.append(codontab[codon])
    translation = ''.join(translation)
    return(translation)

def CDS_from_six_reading_frames_ORFs(nt_seq):
    """Function reports a dictionary for all ORFs found in the six reading frames of a DNA sequence"""
    six_reading_frames_ORFS = {1 : [], 2 : [], 3 : [], 4 : [], 5 : [], 6 : []}
    for frame in range(0,6):
        reading_frame_to_ORF_to_CDS = []
        cds_seq = frame_cds(nt_seq, frame+1)
        ORFs_list = listORFs(cds_seq)
        if ORFs_list:
            for ORF in ORFs_list:
                translation = translate_ORFs(ORF)
                reading_frame_to_ORF_to_CDS.append(translation)
            six_reading_frames_ORFS[frame+1] = reading_frame_to_ORF_to_CDS
    return(six_reading_frames_ORFS)

def ORFs_file(INFILE,OUTPUT_FILENAME='ORFs.faa'):
    """Function outputs a fasta file with ORFs found"""
    with open(INFILE, 'r', encoding="utf-8") as infile:
        lines = infile.readlines()
    with open(OUTPUT_FILENAME, 'w', encoding="utf-8") as output:
        for line in lines:
            if line[0] == '>':
                temp = line.strip('\n')
            else:
                sequence = line.strip('\n').upper()
                ORFs = CDS_from_six_reading_frames_ORFs(sequence)
                for frame in ORFs:
                    if ORFs[frame]:
                        for ORF in range(len(ORFs[frame])):
                            output.write(''.join(filter(str.isalnum, temp))+'_'+str(frame)+str(ORF)+'\n')
                            output.write(''.join(ORFs[frame][ORF])+'\n')
    output.close()
