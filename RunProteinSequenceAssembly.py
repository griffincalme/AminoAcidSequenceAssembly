#This is the user-interface for the protein sequence assembly library

from ProteinSequenceAssemblyLib import open_fragment_library, sequence_assembler

filepath = 'YOUR_LINE_SEPARATED_PEPTIDES.txt'   # Put your filepath here (line-separated sequences)
my_library = open_fragment_library(filepath)

# sequence_assembler returns two lists, one with the assembled peptides and one with the leftover non-matching peptides
my_assembled_peptide_library, my_unmatched_peptides = sequence_assembler(my_library, min_overlap=5, min_match_score=40)


print('Assembled peptides: ')
print(my_assembled_peptide_library)

print('\nLeftover unmatched peptides: ')
print(my_unmatched_peptides)

with open(filepath[:-4] + '_assembled_peptides.txt', 'a') as f1:
    assembled_peptides_string = '\n'.join(my_assembled_peptide_library)
    f1.write(assembled_peptides_string)

with open(filepath[:-4] + '_unmatched_peptides.txt', 'a') as f2:
    unmatched_peptides_string = '\n'.join(my_unmatched_peptides)     # Join subsequences by a newline, then save as txt
    f2.write(unmatched_peptides_string)
