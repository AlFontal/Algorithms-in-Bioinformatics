
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


def globalAlignment(x,y, gap_pen = -4, gap_pen_end = 0, dictionary = BLOSUM62_ORDER, scoring_matrix = BLOSUM62_MATRIX ):
    """
    :param x: Peptide sequence 1 to be aligned
    :param y: Peptide sequence 2 to be aligned
    :param gap_pen: Internal gaps penalty (default = -4)
    :param gap_pen_end: End gaps penalty (default = 0)
    :param dictionary: Dictionary object containing the aminoacidic residues and their index in the scoring matrix
    :param scoring_matrix: Scoring matrix used to obtain the scores for the matches (Default = BLOSUM62)
    :return: Score for the optimal alignment
    """

    sorted_key = sorted(dictionary.items(), key=operator.itemgetter(1))
    alphabet = []
    for i, _ in sorted_key:
        alphabet.append(i)

    #Create matrix object to be filled
    D = []
    endrow = len(x)
    endcol = len(y)

    #Define the matrix dimensions
    for i in range(len(x)+1):
        D.append([0] * (len(y)+1))

    #Fill first column
    for i in range(1, len(x)+1):
        if i == 1:
            D[i][0] = gap_pen_end
        else:
            D[i][0] = D[i-1][0] + gap_pen_end

    #Fill first row
    for j in range(1, len(y)+1):
        if j == 1:
            D[0][j] = gap_pen_end
        else:
            D[0][j] = D[0][j-1] + gap_pen_end



    #Calculate score for every possible movement along the matrix, add gap penalty for Hor/Ver movements and
    #Check scoring matrix for diagonal movements

    #This refers only to the interior values of the matrix (So gap penalty = gap_pen)
    for i in range(1, len(x)):
        for j in range(1, len(y)):

            distHor = D[i-1][j] + gap_pen
            distVer = D[i][j-1] + gap_pen
            distDiag = D[i-1][j-1] + scoring_matrix[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            #Check which of the 3 possibilities would yield the max score and store it as that.
            D[i][j] = max(distHor, distVer, distDiag)

    #Now we have to define the gap penalty on the end row and column so tat gap penalty = gap_pen:end
    for i in range(len(x), len(x)+1):
        for j in range(1, len(y)+1):

            distHor = D[i-1][j] + gap_pen_end
            distVer = D[i][j-1] + gap_pen_end
            distDiag = D[i-1][j-1] + scoring_matrix[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            #Check which of the 3 possibilities would yield the max score and store it as that.
            D[i][j] = max(distHor, distVer, distDiag)

    for i in range(1, len(x)+1):
        for j in range(len(y), len(y)+1):

            distHor = D[i-1][j] + gap_pen_end
            distVer = D[i][j-1] + gap_pen_end
            distDiag = D[i-1][j-1] + scoring_matrix[alphabet.index(x[i-1])][alphabet.index(y[j-1])]

            #Check which of the 3 possibilities would yield the max score and store it as that.
            D[i][j] = max(distHor, distVer, distDiag)

    #Print matrix once filled (helps visibility)
    print('\n'.join([''.join(['{:6}'.format(item) for item in row])for row in D]))



    """
    For the traceback strategy, we'll start in the last position of the matrix (the actual score) and then,
    depending on if the difference between the cell and the cell besides or above is equal to the gap penalty
    the first movement will be stored as an "UP" or "LEFT" movement and the next loop will evaluate the cell
    above or left, respectively. In case not the cell above nor the one besides are possible steps, then it
    will be assumed that the movement was a diagonal movement (a match). The fact that the algorithm has been
    constructed like this means that the traceback will yield just one of the possible optimal solutions to the
    alignment, since it first evaluates if there is a possible path with a gap and if there is not it assumes it
    had to be a match.
    """

    traceback = []
    row = endrow
    col = endcol
    cols = []
    rows = []
    while col > 0 or row > 0:
        if col == 0:
            traceback.append("UP")
            cols.append(col)
            rows.append(row)
            row = row  -1

        elif row == 0:
            traceback.append("LEFT")
            cols.append(col)
            rows.append(row)
            col = col -1


        elif (col != endcol and D[row][col] - D[row-1][col] == gap_pen)\
                or col == endcol and D[row][col] - D[row-1][col] == gap_pen_end:
            traceback.append("UP")
            cols.append(col)
            rows.append(row)
            row = row -1
        elif (row != endrow and D[row][col] - D[row][col-1] == gap_pen)\
                or row == endrow and D[row][col] - D[row][col-1] == gap_pen_end:
            traceback.append("LEFT")
            cols.append(col)
            rows.append(row)
            col = col -1
        else:
            traceback.append("DIAG")
            cols.append(col)
            rows.append(row)
            row = row -1
            col = col -1
    print(traceback)


    """
    Until now, we get a traceback list that shows what path was followed in the matrix to produce the best scoring
    alignment. This list of "LEFT", "UP" and "DIAG" movements, however, isn't as readable as the typical output
    of an alignment that shows one sequence in front of the other with "-" representing the gaps. In order to
    produce a more readable output, and based on the movements of the traceback, I generate two lists that show
    both sequences in a way that allows to compare them visually.
    I did this, first reversing the traceback and then, step by step, adding gaps to the 2nd sequence when the
    movement is vertical ("UP") and adding gaps to the 1st sequences when the movement is horizontal ("LEFT").
    The length of the alignment is longer than either of the sequences, so I solved this by generating two
    variables "l" and "m" that increase whenever a gap is found in the 1st or 2nd sequence respectively. The
    index where to find the character of the sequence is therefore n-l-1 or n-m-1. (n is the index of the alignment).
    """

    trace = []
    for i in reversed(traceback):
        trace.append(i)

    alignment1 = []
    alignment2 = []
    n = 0
    m = 0
    l = 0
    for move in trace:
        n += 1
        if move == "UP":
            m += 1

            alignment1.append(x[n-l-1])
            alignment2.append("-")
        elif move == "LEFT":
            l += 1

            alignment1.append("-")
            alignment2.append(y[n-m-1])

        else:
            alignment1.append(y[n-m-1])
            alignment2.append(x[n-l-1])

    """
    Having the traceback done, getting the % of identical residues in the alignment is rather easy.
    What I did is, using a for loop to check every element of the alignment and decide whether it's
    the same or not. If it's the same, a variable (k) gets increased by one. At the end of the loop
    k will be equal to the number of times two elements with the same index in the alignment were found
    to be the same. The proportion of identity is therefore just k divided by the length of the alignment
    and multiplying this by 100 we get the percentage.
    """

    k = 0
    length = len(traceback)
    for idx in range(length):
        if alignment1[idx] == alignment2[idx]:
            k += 1

    perc_identity = (k/len(alignment1))*100
    print("Identical residues percentage is " + str(perc_identity))

    print(alignment1)
    print(alignment2)


    return D[-1][-1], traceback, alignment1, alignment2, perc_identity


if __name__=="__main__":


    start_time = time.time()


    x = 'ALEGANDR'
    y = 'ILDEFNDR'
    print("--- %s seconds ---" % (time.time() - start_time))



    #x = "THISLINE"
    #y = "ISALIGNED"
    print ('X: '+x)
    print ('Y: '+y)
    a = (globalAlignment(x,y))
    print(a[1])



