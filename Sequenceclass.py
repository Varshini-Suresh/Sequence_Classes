### SEQUENCES SUPERCLASS AND SUBCLASSES ###

## CREATING THE CLASSES ##  
## Creating a Sequence superclass ##
class Sequence: 

    # Adding constrictor to set shared attributes 
    def __init__(self, id, seq):
        self._id = id
        self._seq = seq.upper()

    # Specifying the getters 
    @property
    def id(self):
        return self._id
    
    @property
    def seq(self):
        return self._seq
    
    # Adding similar methods (shared content)
    # Print the class attributes 
    def __str__(self):
        return f'Sequence Object:\nID: {self.id}; Seq: {self.seq};'
    
    # Returns length when len is called on the object 
    def __len__(self):
        return len(self.seq)
    
    # Iterates over sequence when a for loop is called 
    def __iter__(self):
        return iter(self.seq)

    # Returns string in a fasta format  
    def get_fasta(self):
        return f'>{self.id}\n{self.seq}\n'

# Next the subclasses for DNA and Protein sequences (with methods relevant only to those)
# Creating the DNA subclass 
class DNASequence(Sequence): # inherits from superclass 

    # Returns the GC content of the DNA sequence
    def calc_gc(self, dp=2):
        C_count = self.seq.count('C')
        G_count = self.seq.count('G')
        GC_count = round((C_count + G_count)/len(self.seq) *100, dp)
        return f'{GC_count}%'


    # Translates the sequences 
    def translate(self):
        # Reference codon table 
        bases = "tcag".upper()
        codons = [a + b + c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))

        # Translation
        translation = ''
        for i in range(0, len(self.seq)-2, 3): # 
            codon = self.seq[i:i+3]
            aa_seq = codon_table[codon]
            translation += aa_seq
        return translation
    
    # Gets the length of its protein sequence 
    def get_protein_length(self):
        return len(self.seq)//3

# Creating the protein subclass 
class Proteins(Sequence):

    # Overriding the constructor, to have additional description attribute
    def __init__(self, id, seq, descr=''):
        super().__init__(id, seq)
        self._descr = descr

    # Adding a getter for description
    @property
    def descr(self):
        return self._descr
    
    # Methods unique to protein
    # Returns percentage of hydrophobic residues in the sequence
    def hydrophobic_percent(self, dp=2):
        hydro_residues = ['A','I','L','M','F','W','Y','V']
        count = 0
        for i in self.seq:
            if i in hydro_residues:
                count += 1
        percentage = round((count/len(self.seq))*100, dp)
        return (f'{percentage}%')
    
    # Returns length of protein sequence 
    def get_protein_len(self):
        return self._file
    
## Creating the fasta class ## 
class Fastafile:

    # Constrictor - only attribute is string containing a file name
    def __init__(self, file):
        self._file = file
    
    # adding getter for filename
    @property
    def file(self):
        return self._file
    
    # Method to retrieve DNASequence records out of the fasta file
    def get_seq_record(self, sequence_class):
        """
        yields so we can iterate all the way over the file
        and return a DNASequence for every pair of lines
        this function now takes the sequence class as a parameter
        """
        with open(self.file) as fasta:
            for line in fasta:
                if line.startswith('>'):
                    id = line.rstrip().lstrip('>')
                    seq = next(fasta).rstrip()
                yield sequence_class(id, seq)


## SOME TOP LEVEL CODE ## 
           
if __name__ == '__main__': # so it is not run when imported 

    # Creating instances
    dna = DNASequence('geneA', 'ATGCCAGATTAG')   
    protein = Proteins('proteinA', 'MIDAKSEKH*', "sample_protein")

    # testing that my __str__ and __len__methods work
    print(f'{dna}')
    assert len(dna) == 12
    # testing that my __iter__ works 
    for base in dna:
        print(base)
    # testing the translation methods work
    assert dna.translate() == 'MPD*'
    assert dna.get_protein_length() == 4

    # check that my override works and I can add a description to a Protein Sequence
    # if I wanted to print my description here, I would also need to override the __str__ method
    print(f'{protein}')
    # but I can see that the description is there and the getter works
    print(f'The descriptions is: {protein.descr}')
    # testing Proteins class methods 
    print(f'The fasta format of the protein is\n{protein.get_fasta()}')
