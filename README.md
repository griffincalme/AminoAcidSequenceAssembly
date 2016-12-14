# AminoAcidSequenceAssembly
Sequence assembler for amino acid subsequences. Takes a line-separated text file and returns the assembled/merged peptides (using Needleman-Wunsch global alignment) as a text file.


ProteinSequenceAssemblyLib.py contains functions for the "RunProteinSequenceAssembly.py" script.

Edit the RunProteinSequenceAssembly.py with your line-separated text file of peptide subsequences and run it.


In the "Validation" folder, "ProteinSequenceAssembly_NW_Validation.py" validates the algorithm by generating random polypeptides, randomly fragmenting them, and running the assembly on the resulting subsequences. Probably needs more improvement.
