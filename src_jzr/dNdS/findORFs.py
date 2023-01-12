"""Module containing code for identifying open reading frames"""

from more_itertools import sliced

# class for translating 

complementaryDNA = {"A":"T", "C":"G", "G":"C", "T":"A"}

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
    reverse = nt_seq[::-1]
    complement = ''
    for base in reverse:
        complement += complementaryDNA[base]
    return(complement)

def obtain_frame(nt_seq,frame):
    """Function to obtain a specific reading frame"""
    if frame > 2:
        sequence = reverse_complement(nt_seq)
    nt_frame = list(sliced(sequence[frame:],3))
    return(nt_frame)

def translate_sequence(nt_seq, **codontab):
    translation=[]
    for codon in nt_seq:
        if codon in codontab:
            translation.append(codontab[codon])
    return(translation)

def translate_file(INFILE,OUTPUT_FILENAME='ORFs.faa'):
    with open(OUTPUT_FILENAME, 'w', encoding="utf-8") as output:``
        for line in lines:
            if line[0] == '>':
                temp = line.strip('\n')
            else:
                sequence = line.strip('\n').upper()
                for frame in range(0,6):
                    output.write(str(temp)+'frame_'+str(frame+1)+'\n')
                    nt_frame = obtain_frame(sequence, frame)
                    translation = translate_sequence(nt_frame)
                    OUTPUT_FILENAME.write(''.join(translation)+'\n')
    OUTPUT_FILENAME.close()
