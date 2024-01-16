#### MAIN FILE ####
# Add a new file with 2-3 of the functions created in this course
# Saved them as Fuctions.py.. the module name is Functions
 
# METHOD 1: IMPORTING MODULES 
# a) Importing whole module 
from Sequenceclass import DNASequence

dna_sequence = DNASequence('seq1','ATGGTCCGATTGCAATCGATCGATTAGC')
print(f'The GC content of the sequence is {dna_sequence.calc_gc()}%')


from Sequenceclass import Proteins
pro_sequence1 = Proteins('aa1', 'MDKPLFATH*', 'sample_protein_1')
pro_sequence2 = Proteins('aa2', 'MKSDEHAATPLZR*', 'sample_protein_2')

print(f'The percentage of hydrophobic residues in the first sequence is {pro_sequence1.hydrophobic_percent()}\n' \
      f'while that of the second protein sequence is {pro_sequence2.hydrophobic_percent()}')

# METHOD 2: IMPORTING THE ENTIRE PACKAGE 
# Problem: Using a genome fasta file, write out a protein translation in fasta format
from Sequenceclass import * 

fasta_file = 'sample_seq.txt'

with open('test.out', 'w') as output:
    # make a FastaFile object
    fasta = Fastafile(fasta_file)
    # the get_seq_record method yields, so I can use it in a for loop
    for record in fasta.get_seq_record(DNASequence):
        if len(record)%3 == 0:
            # make a Protein Sequence record with the id and the translated sequence
            protein = Proteins(record.id, record.translate())
            # write this to a fasta file using our create_fasta function
            output.write(protein.get_fasta())

