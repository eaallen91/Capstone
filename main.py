'''
Eric Allen

main program to control building and testing of 3 seq alignments
'''

from Memo import align3D
from Semi_Memo import alignSemi3D
from LCS import LCS3D
from Local import local3D

### Executes the alignmets allowing for testing all of the 
### different tests
def executeAlignment(seq1, seq2, seq3):
    
    #seq1 = input("Please input the first string: ")
    #seq2 = input("Please input the second string: ")
    #seq3 = input("Please input the third string: ")
    alignedSeqs = []
    alignment = input("Would you like \"Global\", \"Semi\", \"LCS\" \"Local\": ")
    if alignment == "LCS":
        alignedSeqs = LCS3D(seq1, seq2, seq3)
    else:
        gap = eval(input("Choose a gap penalty: "))
        mismatch = eval(input("Choose a mismatch penalty: "))
        match = eval(input("Choose a match penalty: "))
        minMax = input("Would you like \"Min\" or \"Max\": ")
        if alignment == "Global":
            alignedSeqs = align3D(seq1, seq2, seq3, gap, mismatch, match, minMax)
        elif alignment == "Semi":
            alignedSeqs = alignSemi3D(seq1, seq2, seq3, gap, mismatch, match, minMax)
        elif alignment == "Local":
            alignedSeqs = local3D(seq1, seq2, seq3, gap, mismatch, match)
        
    print(seq1)
    print(seq2)
    print(seq3)
    print("")
    for k in range(len(alignedSeqs)):
        print(alignedSeqs[k])
    
def main():

    #### Test 1: From class slides
    seq1 = "AGT"
    seq2 = "CGT"
    seq3 = "ACT"
    print("TEST 1")
    #executeAlignment(seq1, seq2, seq3)

    #### Test 2: Made up sequences
    seq1 = "GAAT"
    seq2 = "GAT"
    seq3 = "AATG"
    print("TEST 2")
    executeAlignment(seq1, seq2, seq3)

    #### Test 3: First 140 characters from 3 seqs in mtDNA.fasta
    seq1 = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTC"
    seq2 = "GTTTATGTAGCTTACCCCCTCAAAGCAATACACTGAAAATGTTTCGACGGGTTTACATCACCCCATAAACAAACAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGCCCCGTG"
    seq3 = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTC"
    print("TEST 3")
    #executeAlignment(seq1, seq2, seq3)

    #### Test 4: First 20 characters from test 3
    seq1 = "GATCACAGGTCTATCACCCT"
    seq2 = "GTTTATGTAGCTTACCCCCT"
    seq3 = "GATCACAGGTCTATCACCCT"
    print("TEST 4")
    #executeAlignment(seq1, seq2, seq3)

    #### Test 5: First 320 character from 3 seqs in litmus.fasta
    seq1 = "GTTAACTACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCC"
    seq2 = "GTTAACTACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCC"
    seq3 = "GTTAACTACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCC"
    print("TEST 5")
    #executeAlignment(seq1, seq2, seq3)

    #### Test 6: For semi global cases
    seq1 = "GAAATAATA"
    seq2 = "GAATTAA"
    seq3 = "ATTA"
    print("TEST 6")
    executeAlignment(seq1, seq2, seq3)

    #### Test 7: Actual words
    seq1 = "Hello"
    seq2 = "ello"
    seq3 = "Beheld"
    print("TEST 7")
    #executeAlignment(seq1, seq2, seq3)

    ### Test 8
    seq1 = "SHAKE"
    seq2 = "SPEARE"
    seq3 = "THEATRE"
    print("TEST 8")
    executeAlignment(seq1, seq2, seq3)

    ### Test 9
    seq1 = "POTATOES"
    seq2 = "TOMATOES"
    seq3 = "PATS"
    print("TEST 9")
    executeAlignment(seq1, seq2, seq3)

main()
