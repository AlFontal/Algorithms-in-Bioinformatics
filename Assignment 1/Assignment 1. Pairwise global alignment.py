
#!/usr/bin/env python
"""
Author: Alejandro Fontal
Pairwise global alignment with Dynamic programming
"""
from __future__ import division
import time
import operator

blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 23:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            blosum_matrix.append(map(int,parts[1:]))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

# write your own functions below here

def arrayAsMatrix(array):
    """
    This function takes an array and prints it in a way that is visualized as
    a matrix. It doesn't return any value, it's just for visualization purposes.
    """
    print('\n'.join([''.join(['{:6}'.format(item) for item in row])for row in array]))

def globalAlignment(seq1,seq2, gap_open_pen = -4, gap_end_pen = 0,
                    dictionary = BLOSUM62_ORDER,
                    scoring_matrix = BLOSUM62_MATRIX,
                    printInitMatrix = False,
                    printFilledMatrix = True):
    """
    :param seq1: string object containing an aminoacidic sequence to be aligned
    :param seq2: string object containing the other aminoacidic sequence to be aligned
    :param gap_open_pen: Internal gaps penalty (default = -4)
    :param gap_end_pen: End gaps penalty (default = 0)
    :param dictionary: Dictionary object containing the aminoacidic residues and their index in the scoring matrix
    :param scoring_matrix: Scoring matrix used to obtain the scores for the matches (Default = BLOSUM62)
    :param printInitMatrix: To decide whether to print or not the filled matrix (Default = True)
    :param printFilledMatrix: To decide whether to print or not the initial matrix (Default = True)
    :return:Score for the optimal alignment, filled matrix with all the scores per cell.
    """

    sorted_key = sorted(dictionary.items(), key=operator.itemgetter(1))
    alphabet = []
    for i, _ in sorted_key:
        alphabet.append(i)

    #Create matrix object to be filled
    A = []
    init_matrix = []
    endrow = len(seq1) 
    endcol = len(seq2) 
    
    #Define the matrix dimensions
    for row in range(endrow+1):
        A.append([0] * (endcol+1))

    #Fill first column
    for row in range(1, endrow+1):
        A[row][0] = A[row-1][0] + gap_end_pen
    
    #Fill first row
    for col in range(1, endcol+1):
        A[0][col] = A[0][col-1] + gap_end_pen

    if printInitMatrix == True:

        init_matrix = A
        print("The initial matrix, once first row and column have been filled is:")
        arrayAsMatrix(init_matrix)
    """
    Calculate score for every possible movement along the matrix, add gap penalty for Hor/Ver movements and check
    the scoring matrix for diagonal movements. Calculate the maximum of the three values and store it in the cell.
    """

    #Filling the interior values of the matrix (So gap penalty = gap_open_pen)
    for row in range(1, endrow):
        for col in range(1, endcol):

            moveVer = A[row-1][col] + gap_open_pen
            moveHor = A[row][col-1] + gap_open_pen
            moveDiag = A[row-1][col-1] + scoring_matrix[alphabet.index(seq1[row-1])][alphabet.index(seq2[col-1])]

            #Check which of the 3 possibilities would yield the max score and store it as that.
            A[row][col] = max(moveHor, moveVer, moveDiag)

    """
    Now we have to define the gap penalty on the end row and column so that gap penalty = gap_end_pen for the
    vertical movements in the last column and the horizontal movements in the last row. The cell in the last
    row and the last column will have for both horizontal and vertical movements gap penalty = gap_end_pen
    """

    #Filling the last row
    for col in range(1, endcol):
            row = endrow
            moveVer = A[row-1][col] + gap_open_pen
            moveHor = A[row][col-1] + gap_end_pen
            moveDiag = A[row-1][col-1] + scoring_matrix[alphabet.index(seq1[row-1])][alphabet.index(seq2[col-1])]

            #Check which of the 3 possibilities would yield the max score and store it as that.
            A[row][col] = max(moveHor, moveVer, moveDiag)
            
    #Filling the last column
    for row in range(1, endrow):
            col = endcol
            moveVer = A[row-1][col] + gap_end_pen
            moveHor = A[row][col-1] + gap_open_pen
            moveDiag = A[row-1][col-1] + scoring_matrix[alphabet.index(seq1[row-1])][alphabet.index(seq2[col-1])]

            #Check which of the 3 possibilities would yield the max score and store it as that.
            A[row][col] = max(moveHor, moveVer, moveDiag)

    #Filling the last value
    row = endrow
    col = endcol
    moveVer = A[row-1][col] + gap_end_pen
    moveHor = A[row][col-1] + gap_end_pen
    moveDiag = A[row-1][col-1] + scoring_matrix[alphabet.index(seq1[row-1])][alphabet.index(seq2[col-1])]

    #Check which of the 3 possibilities would yield the max score and store it as that.
    A[row][col] = max(moveHor, moveVer, moveDiag)    

    if printFilledMatrix == True:
        print("The filled matrix, once the scores for all the cells have been calculated is:")
        arrayAsMatrix(A)

    print("The score of the best alignment is " + str(A[-1][-1]))
    return A[-1][-1], A


"""
For the traceback strategy, the algorithm will start in the last position of the matrix (the actual score) and then,
depending on if the difference between the cell and the cell besides or above is equal to the gap penalty the first
movement will be stored as an UP("U") or LEFT("L") movement and the next loop will evaluate the cell above or left,
respectively. In case not the cell above nor the one besides are possible steps, then it will be assumed that the
movement was a diagonal movement and the value "D" will be stored. The fact that the algorithm has been constructed
like this means that the traceback will yield just one of the possible optimal solutions to the alignment, since it
first evaluates if there is a possible path with a gap and if there is not it assumes it had to be an alignment.
"""

def getTraceback(seq1, seq2, filled_matrix, gap_open_pen = -4, gap_end_pen = 0,
                 printAlignment = True):
    """
    :param seq1: string object containing an aminoacidic sequence to be aligned
    :param seq2: string object containing the other aminoacidic sequence to be aligned
    :param filled_matrix: array object containing the filled matrix of a global alignment
    :param gap_open_pen: Internal gaps penalty (default = -4)
    :param gap_end_pen: End gaps penalty (default = 0)
    :param printAlignment: decides whether the alignment is printed or not
    :return:
    """

    A = filled_matrix
    traceback = []
    endrow = len(seq1)
    endcol = len(seq2)
    row = endrow
    col = endcol
    rows = []
    cols = []

    while col > 0 or row > 0:
        if col == 0:
            traceback.append("U")
            cols.append(col)
            rows.append(row)
            row = row  -1

        elif row == 0:
            traceback.append("L")
            cols.append(col)
            rows.append(row)
            col = col -1


        elif (col != endcol and A[row][col] - A[row-1][col] == gap_open_pen)\
            or (col == endcol and (A[row][col] - A[row-1][col]) == gap_end_pen):
            traceback.append("U")
            cols.append(col)
            rows.append(row)
            row = row -1

        elif (row != endrow and A[row][col] - A[row][col-1] == gap_open_pen)\
                or row == endrow and A[row][col] - A[row][col-1] == gap_end_pen:
            traceback.append("L")
            cols.append(col)
            rows.append(row)
            col = col -1
        else:
            traceback.append("D")
            cols.append(col)
            rows.append(row)
            row = row -1
            col = col -1


    """
    Until now, we get a traceback list that shows what path was followed in the matrix to produce the best scoring
    alignment. This list of "L", "U" and "D" movements, however, isn't as readable as the typical output
    of an alignment that shows one sequence in front of the other with "-" representing the gaps. In order to
    produce a more readable output, and based on the movements of the traceback, I generate two lists that show
    both sequences in a way that allows to compare them visually.
    I did this, first reversing the traceback and then, step by step, adding gaps to the 2nd sequence when the
    movement is vertical ("U") and adding gaps to the 1st sequences when the movement is horizontal ("L").
    The length of the alignment is longer than either of the sequences, so I solved this by generating two
    variables "l" and "m" that increase whenever a gap is found in the 1st or 2nd sequence respectively. The
    index where to find the character of the sequence is therefore n-l-1 or n-m-1. (n is the index of the alignment).
    Finally, I also create an "alignment joint" where, in order to visualize the alignment, identical residues are
    united by a symbol "|" and alignments of non identical residues are united by ":".
    """

    trace = []
    for i in reversed(traceback):
        trace.append(i)
    alignmentseq1 = []
    alignmentseq2 = []
    alignment = []
    n = 0
    m = 0
    l = 0
    for move in trace:
        n += 1
        if move == "U":
            m += 1

            alignmentseq1.append(seq1[n-l-1])
            alignmentseq2.append("-")
        elif move == "L":
            l += 1

            alignmentseq1.append("-")
            alignmentseq2.append(seq2[n-m-1])

        else:
            alignmentseq2.append(seq2[n-m-1])
            alignmentseq1.append(seq1[n-l-1])
    for idx in range(len(traceback)):

        if alignmentseq1[idx] == alignmentseq2[idx]:
            alignment.append("|")

        elif alignmentseq1[idx] == "-" or alignmentseq2[idx] == "-":
            alignment.append(" ")

        else:
            alignment.append(":")

    if printAlignment == True:
        print("Seq 1:  "+''.join(alignmentseq1)+ "\n        " + ''.join(alignment) + "\nSeq 2:  "+\
              ''.join(alignmentseq2))


    return alignment

"""
Having the traceback done, getting the % of identical residues in the alignment is rather easy.
The amount of identical residues can be obtained by counting the number of "|" symbols in the
alignment list generated by the getTraceback() function. The proportion of identity is therefore
just this number of identical residues divided by the length of the alignment.
Multiplying this by 100 we get the percentage.
"""
def getIdentity(alignment):
    """
    :param alignment: alignment list generated by traceback function
    :return: percentage of identity between the residues of the two sequences
    """

    num_identical_residues = alignment.count("|")
    aln_len = len(alignment)
    perc_identity = (num_identical_residues/aln_len)*100
    number_gaps = alignment.count(" ")
    perc_gaps = (number_gaps/aln_len)*100
    print("There is a " + str(round(float(perc_identity),2))) +\
    "% of identical residues (" + str(num_identical_residues) + ")"
    print("Number of gaps: " + str(number_gaps) + " (" + str(round(float(perc_gaps),2)) + "%)")
    return perc_identity

    

    
if __name__=="__main__":


    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    # seq3: GPA1_ARATH
    seq3 = "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYM" \
           "LSSESIAIGEKLSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVY" \
           "RLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIE" \
           "HAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
    # seq4: GPA1_ORYSI
    seq4 = "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTGFDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKY" \
           "VISPDNQEIGEKLSDIDGRLDYPLLNKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVVQIQFSPVGENKRGGEV" \
           "YRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMMETKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEV" \
           "EHAYEFVKKKFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"


    # input for the first sequence
    seq1_input = raw_input("Put here the first sequence to be aligned:\n").upper()
    # input for the second sequence
    seq2_input = raw_input("Put here the second sequence to be aligned:\n").upper()
    # input for the open gap penalty
    open_penalty = int(raw_input("What should be the open gap penalty?\n"))
    if raw_input("Should there be a separate end gap penalty? (Y/N)\n").upper() == "N":
        end_penalty = open_penalty
    else:
        end_penalty = int(raw_input("What should be the end gap penalty?\n"))

    printMatrix = (raw_input("Should the filled matrix be printed? (Y/N)\n").upper() == "Y")
    printAln = (raw_input ("Should the alignment of the sequences be printed? (Y/N)\n").upper() == "Y")


    #Calling the functions:

    score, filled_matrix = globalAlignment(seq1_input,seq2_input, printFilledMatrix= printMatrix,
                                           gap_open_pen = -open_penalty,
                                           gap_end_pen = -end_penalty)
    alignment_list = getTraceback(seq1_input, seq2_input, filled_matrix, gap_open_pen = -open_penalty,
                                  gap_end_pen = -end_penalty, printAlignment = printAln)
    identity = getIdentity(alignment_list)

"""
    #Question 3:
    for i in range (1,21):
        print(i)
        score, filled_matrix = globalAlignment(seq1, seq2, gap_open_pen = -i,
                                               gap_end_pen = -i,
                                               printFilledMatrix = False)
        alignment_list = getTraceback(seq1, seq2, filled_matrix,
                                      gap_open_pen = -i, gap_end_pen = -i)
"""




