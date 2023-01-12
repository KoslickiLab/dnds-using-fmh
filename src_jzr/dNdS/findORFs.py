"""Module containing code for identifying open reading frames"""

from more_itertools import sliced

complementaryDNA = {"A":"T", "C":"G", "G":"C", "T":"A"}

stop_codons = ['TAA', 'TAG', 'TGA']

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

def reverse_complement(nt_seq, **complementaryDNA):
    """Function reports the reverse complement of a DNA sequence"""
    reverse = nt_seq[::-1]
    complement = ''
    for base in reverse:
        complement += complementaryDNA[base]
    return(complement)

def frame_cds(nt_seq,frame):
    """Function reports the coding sequence of a reading frame"""
    if frame < 3:
        coding_seq = list(sliced(sequence[frame:],3))    
    elif frame > 2:
        sequence = reverse_complement(nt_seq)
        coding_seq = list(sliced(sequence[frame-3:],3))
    return(coding_seq)

def listORFs(cds, **stop_codons):
    """Function reports a list of ORFs found in a reading frame"""
    ORFs = []
    tempORFs = []
    for codon in cds:
        if codon == 'ATG':
            tempORFs.append(codon)
        elif codon in stop_codons:
            ORFs.append(tempORFs)
            tempORFs = []
        elif codon != 'ATG' and codon not in stop_codons:
            if tempORFs:
                tempORFs.append(codon)
    return(ORFs)

def translate_ORFs(cds_seq, **codontab):
    """Function translates open reading frames using the coding sequence"""
    translation=[]
    for codon in cds_seq:
        if codon in codontab:
            translation.append(codontab[codon])
    translation = ''.join(translation)
    return(translation)

def CDS_from_six_reading_frames_ORFs(nt_seq):
    """Function reports a dictionary for all ORFs found in the six reading frames of a DNA sequence"""
    six_reading_frames_ORFS = {1 : None, 2 : None, 3 : None, 4 : None, 5 : None, 6 : None}
    reading_frame_to_ORF_to_CDS = []
    for frame in range(0,6):
        cds_seq = frame_cds(nt_seq, frame)
        ORFs_list = istORFs(cds_seq)
        for ORF in ORFs_list:
            translation = translate_cds(ORF)
            reading_frame_to_ORF_to_CDS.append(translation)
        six_reading_frames_ORFS[frame+1] = reading_frame_to_ORF_to_CDS
    return(six_reading_frames_ORFS)

def translation_file(INFILE,OUTPUT_FILENAME='ORFs.faa'):
    """Function outputs a fasta file with ORFs found"""
    with open(OUTPUT_FILENAME, 'w', encoding="utf-8") as output:
        for line in lines:
            if line[0] == '>':
                temp = line.strip('\n')
            else:
                sequence = line.strip('\n').upper()
                ORFs = CDS_from_six_reading_frames_ORFs(sequence)
                for frame in ORFs:
                    for ORF in range(len(ORFs[frame])):
                            output.write(''.join(filter(str.isalnum, temp))+'_'+frame+ORF+'\n')
                            OUTPUT_FILENAME.write(''.join(ORFs[ORF])+'\n')
    OUTPUT_FILENAME.close()
