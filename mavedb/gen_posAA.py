""" generates the possible amino acids given the sequence files from mavedb """

import os

# ------------- #
# required data #
# ------------- #

dna_to_aa = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
            'TAT':'Y','TAC':'Y','TAA':'X','TAG':'X',
            'TGT':'C','TGC':'C','TGA':'X','TGG':'W',
            'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
            'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
            'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

aa_abrev = {'Ala': 'A', 'Cys': 'C', 'Asp':'D', 'Glu':'E','Phe':'F','Gly':'G',
'His':'H','Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N','Pro':'P',
'Gln':'Q','Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W','Tyr':'Y',
'Ter':'X','del':'del','ins':'ins','dup':'dup','Del':'del','Dup':'dup','Ter':'X'}
aa_unabrev = {v: k for k, v in aa_abrev.items()}

basepairs = ['T','C','A','G']

# ----------- #
#  functions  #
# ----------- #

def codonconv(codon):
    # find the amino acid for a given codon using the dictionary dna_to_aa
    try:
        aa = dna_to_aa[codon]
        return(aa)
    except:
        return(None)

def allposs(codon, orig):
    # find all posible protein consequences given a single codon
    aalist = []
    for pos, letter in enumerate(codon):
        for mut in basepairs:
            if pos == 0:
                mutcodon = mut + codon[1:] 
            if pos == 1:
                mutcodon = codon[0] + mut + codon[2]
            if pos == 2:
                mutcodon = codon[:2] + mut
            aa = aa_unabrev[codonconv(mutcodon)]
            if aa == None:
                continue
            aalist.append(aa)
    aaset = set(aalist)
    if orig in aaset:
        aaset.remove(orig)
    codonlist = []
    for pos, letter in enumerate(codon):
        for mut in basepairs:
            if pos == 0:
                mutcodon = mut + codon[1:] 
            if pos == 1:
                mutcodon = codon[0] + mut + codon[2]
            if pos == 2:
                mutcodon = codon[:2] + mut
            aa = aa_unabrev[codonconv(mutcodon)]
            if aa in aaset:
                codonlist.append(mutcodon)
                aaset.remove(aa)
    codonset = set(codonlist)
    return(codonset)

def get_posAA(orf):
    """
    figures out all possible snp mutation of the given transcript
    """
    codonlist = []
    codonorf = [orf[i:i+3] for i in range(0, len(orf), 3)]
    for i, codon in enumerate(codonorf):
        orig_aa = aa_unabrev[dna_to_aa[codon]]
        allposs_codon = allposs(codon, orig_aa)
        for mutcodon in allposs_codon:
            mutaa = aa_unabrev[codonconv(mutcodon)]
            mut_abrev = aa_abrev[mutaa]
            orig_abrev = aa_abrev[orig_aa]
            name = orig_abrev + str(i+1) + mut_abrev
            hgvs_pro = 'p.' + orig_aa + str(i+1) + mutaa
            #primer = primermaker("","",codonorf, i, 6, mutcodon)
            #codonrow = [orig_aa, str(i+1), mutaa, codon, mutcodon, primer,
            #        reversecomp(primer)]
            codonrow = [hgvs_pro, orig_aa, str(i+1), mutaa, codon, mutcodon]
            codonlist.append(codonrow)
    return codonlist

# ------------ #
# Start script #
# ------------ #

# iterate through targets and files
for target in os.listdir('scoresets'):
    for filename in os.listdir('scoresets/' + target):
        #print(filename)
        if filename.endswith('_seq.fasta'):
            #print(filename)
            # this is the the sequence we want
            with open('scoresets/' + target + '/' + filename, 'r') as f:
                seq = f.readlines()[0]
            posAA_filename = filename[:-10] + '_posAA.csv'
            posAA = get_posAA(seq)
            with open('scoresets/' + target + '/' + posAA_filename, 'w') as f:
                for row in posAA:
                    f.write(','.join(row)+'\n')
    print(target)