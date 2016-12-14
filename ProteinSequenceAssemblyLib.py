# Must install more_itertools and scikit-bio(need to run "python setup.py install" from source for skbio)
# Written by Griffin Calme (2016)
# Runs a Needleman-Wunsch global sequence alignment and iteratively merges the overlapping sequences
# Takes a line-separated text file of overlapping peptide subsequences and outputs the original supersequence
# This is similar to shotgun sequencing of DNA, but for peptides


from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from more_itertools import unique_everseen
import sys


def open_fragment_library(filename):  # Not a dependency
    with open(filename) as f:
        subsequence_library = f.read().splitlines()  # Open file and save to list, split items by line

    no_duplicates_library = list(unique_everseen(subsequence_library))  # Remove duplicate subsequences
    no_duplicates_library = [x.upper() for x in no_duplicates_library]  # Force uppercase amino acids

    print('\nFilename: ' + filename)
    print('\nTotal number of amino acid subsequences: ' + str(len(subsequence_library)))
    print('Unique amino acid subsequences: ' + str(len(no_duplicates_library)))

    print('\nSubsequence library with duplicates removed:')
    print(no_duplicates_library)
    print('\n')

    return no_duplicates_library


def contig_merger(growing_seq, compared_seq, original_growing_sequence):
    merged_contig = []  # List for containing each AA of merged sequence

    for index, letter1 in enumerate(growing_seq):
        letter2 = compared_seq[index]

        if letter1 == '-':             # If the letter in seq1 is hyphen, take letter from seq2
            merged_contig.append(letter2)

        elif letter2 == '-':
            merged_contig.append(letter1)  # If the letter in seq2 is a hyphen, take letter from seq1

        elif letter1 != letter2:           # If the letters do not match anywhere in the sequence, abort merging
            merged_contig = original_growing_sequence    # But return the growing sequence
            break

        elif letter1 == letter2:        # If the letters match, take the letter from the growing seq
            merged_contig.append(letter1)

    # In the case that merging was aborted, this will do nothing to a string type
    merged_contig = ''.join(merged_contig)  # Takes the list of letters and merges back to a sequence string

    return merged_contig


def sequence_assembler(fragment_library, min_overlap, min_match_score=40):
    working_library = fragment_library
    assembled_peptide_library = []

    while len(working_library) != 0:
        # try first:
        growing_sequence = working_library[0]
        working_library.remove(growing_sequence)

        # This part runs the assembly
        for sequence in working_library:

            aln, score, position_list = global_pairwise_align_protein(Protein(growing_sequence), Protein(sequence),
                                                      gap_open_penalty=10000, penalize_terminal_gaps=False)
            grow_seq, compare_seq = aln
            match = grow_seq.match_frequency(compare_seq, relative=False)

            if match >= min_overlap and score > min_match_score:
                my_merged_contig = contig_merger(str(grow_seq), str(compare_seq), growing_sequence)
                growing_sequence = my_merged_contig

                sys.stdout.write('\r')
                sys.stdout.write(growing_sequence)
                sys.stdout.flush()

                # removes used fragments from working library
                working_library = [fragment for fragment in working_library if fragment not in growing_sequence]

        print('')
        assembled_peptide_library.append(growing_sequence)

    # Remove peptides that have not found any matches
    unmatched_peptides = [peptide for peptide in assembled_peptide_library if peptide in fragment_library]
    assembled_peptide_library = [peptide for peptide in assembled_peptide_library if peptide not in fragment_library]

    return assembled_peptide_library, unmatched_peptides
