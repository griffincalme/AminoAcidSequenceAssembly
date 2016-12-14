#Must install more_itertools and scikit-bio(need to run "python setup.py install" from source for skbio)

from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from more_itertools import unique_everseen
import random


def contig_merger(first_seq, second_seq):
    merged_contig = []  # List for containing each AA of merged sequence
    counter = 0

    for AA in first_seq:
        if AA == '-':             # if the item is hyphen, take letter from other sequence
            merged_contig.append(second_seq[counter])
        elif AA.isalpha():        # else if the item is a letter, use that letter
            merged_contig.append(AA)
        else:                     # otherwise there must be an error, delete both sequences
            merged_contig = []
            break
        counter += 1
    merged_contig = ''.join(merged_contig)

    return merged_contig



min_overlap = 6

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
sequence_length = 200
max_subsequence_length = 8
min_subsequence_length = 7
Validation_Runs = 5

correct_list = []

for i in range(0,Validation_Runs):
    sequence = ''.join(random.choice(amino_acids) for i in range(sequence_length))


    subsequence_number = 1000
    subsequence_library = [sequence]


    for i in range(0, subsequence_number):
        subsequence_length = random.randrange(min_subsequence_length, max_subsequence_length + 1)

        start = random.randint(0, sequence_length - subsequence_length)  # Random starting point for slice
        end = start + subsequence_length  # Ending point for slice

        subsequence = sequence[start:end]
        subsequence_library.append(subsequence)



    original_supersequence = subsequence_library[0] #DO NOT USE IN PRODUCTION, for validating original sequence
    subsequence_library = subsequence_library[1:]   #DO NOT USE IN PRODUCTION, for validating original sequence

    no_duplicates_library = list(unique_everseen(subsequence_library))  # Remove duplicate subsequences
    no_duplicates_library = [x.upper() for x in no_duplicates_library]  # Force uppercase amino acids

    growing_sequence = no_duplicates_library[0]  # Assign the first subsequence to be the seed
    no_duplicates_library.remove(growing_sequence)  # Remove the seed from the subsequence library

    print('\nTotal number of amino acid subsequences: ' + str(len(subsequence_library)))
    print('Unique amino acid subsequences: ' + str(len(no_duplicates_library)))
    print('\nSeed sequence is ' + growing_sequence)
    print('\nSubsequence library with seed and any duplicates removed:')
    print(no_duplicates_library)
    print('\n')

    working_library = no_duplicates_library

    #this part runs the assembly
    for i in range(0,100):
        for j in working_library:
            if j in growing_sequence:
                working_library.remove(j)

            else:
                aln, _, _ = global_pairwise_align_protein(Protein(growing_sequence), Protein(j), penalize_terminal_gaps=False)
                seq1, seq2 = aln

                match = seq1.match_frequency(seq2, relative=False)
                if match >= min_overlap:
                    merged_contig = contig_merger(str(seq1), str(seq2))
                    growing_sequence = merged_contig


                    print(growing_sequence)



    print('\nLeftover unmatched subsequences: ')
    if len(working_library) == 0:
        print('None!')
    else:
        print(working_library)

    print('\nYour original supersequence is: ')
    print(original_supersequence)

    print('\nYour guessed supersequence is: ')
    print(growing_sequence)

    if original_supersequence == growing_sequence:
        print('\nCorrect!')
        correct_list.append(True)
    else:
        correct_list.append(False)


print(correct_list)

print(str(sum(correct_list)/len(correct_list) * 100) + "%")
