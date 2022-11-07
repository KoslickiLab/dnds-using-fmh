"""Module providing sliced Function to obtain codon sequence list."""
from more_itertools import sliced

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

INFILE='../data/500_sequences.fna'
with open(INFILE, 'r', encoding="utf-8") as infile:
    lines = infile.readlines()

OUTPUT='../data/frames_500_sequences.faa'
#output = open(OUTPUT, 'w', encoding="utf-8")
with open(OUTPUT, 'w', encoding="utf-8") as output:
    for line in lines:
        if line[0] == '>':
            temp = line.strip('\n')[]
        else:
            sequence = line.strip('\n').upper()
            for frame in range(0,6):
                output.write(str(temp)+'frame_'+str(frame+1)+'\n')
                translation=[]
                if frame < 3:
                    nt = list(sliced(sequence[frame:],3))
                elif frame > 2:
                    nt = list(sliced(sequence[::-1][frame-3:],3))
                for codon in nt:
                    if codon in codontab:
                        translation.append(codontab[codon])
                output.write(''.join(translation)+'\n')

output.close()
