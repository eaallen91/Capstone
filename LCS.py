'''
Eric Allen

Optimal alignment of three sequences

The most current testing version. Each position contains an array of [val, case].
Once the best val has been determined then the case that generated the val is stored.
The path back through the matrix is then known.
'''

'''
These are the cases that I refer to as case 1, case 2, . . . case 7

Case 1   Case 2   Case 3   Case 4   Case 5   Case 6   Case 7
  X        -       -         X        X        -        X   
  -        X       -         X        -        X        X
  -        -       X         -        X        X        X
'''

from CommonFunctions import findMax

'''
3D Sequence alignment that will fill a 3D
matrix and with the help of a helper variable
backtrace through the matrix and output the
optimal solution
'''
def LCS3D(dnaSeq1, dnaSeq2, dnaSeq3):

    dnaMatrix1 = []
    dnaMatrix2 = []
    dnaMatrix3 = []

    for i in dnaSeq1:
        dnaMatrix1.append(i)

    for j in dnaSeq2:
        dnaMatrix2.append(j)

    for k in dnaSeq3:
        dnaMatrix3.append(k)

    ## Add the empty spot onto the beginning of the strings. So that they may
    ## match up with the matrix
    dnaMatrix1.insert(0, "_")
    dnaMatrix2.insert(0, "_")
    dnaMatrix3.insert(0, "_")

    ### Create a 3D Scoring matrix and a 3D memo matrix
    scoringMatrix = [ [ [ 0 for i in range(len(dnaMatrix1)) ] for j in range(len(dnaMatrix2)) ] for k in range(len(dnaMatrix3)) ]

    ### Fills the entire matrix and handles all special cases
    for level in range(0, len(dnaMatrix3)):
        for col in range(0, len(dnaMatrix1)):
            for row in range(0, len(dnaMatrix2)):

                oneToTwo = bool(dnaMatrix1[col] == dnaMatrix2[row])
                oneToThree = bool(dnaMatrix1[col] == dnaMatrix3[level])
                twoToThree = bool(dnaMatrix2[row] == dnaMatrix3[level])

                ### Special cases exist for the first level of the 3D matrix
                ### it is handled as a 2D matrix
                if level == 0:
                    if (row == 0) and (col == 0):
                        scoringMatrix[level][row][col] = 0
                    elif row == 0:
                        scoringMatrix[level][row][col] = 0
                    elif col == 0:
                        scoringMatrix[level][row][col] = 0
                    else:
                        scoringMatrix[level][row][col] = 0

                ### There are special cases for edges throughout the rest of the matrix
                ### but other than that there are no special cases. The majority of the 
                ### the cells in the matricies will be the final else.
                else:
                    if (row == 0) and (col == 0):
                        scoringMatrix[level][row][col] = 0
                    elif row == 0:
                        scoringMatrix[level][row][col] = 0
                    elif col == 0:
                        scoringMatrix[level][row][col] = 0
                    else:
                        if oneToTwo and twoToThree:
                            scoringMatrix[level][row][col] = scoringMatrix[level-1][row-1][col-1] + 1     ### Case 7
                        else:
                            scoringMatrix[level][row][col] = max([scoringMatrix[level][row][col-1],       ### Case 1
                                                                      scoringMatrix[level][row-1][col],   ### Case 2
                                                                      scoringMatrix[level-1][row][col],   ### Case 3
                                                                      scoringMatrix[level][row-1][col-1], ### Case 4
                                                                      scoringMatrix[level-1][row][col-1], ### Case 5
                                                                      scoringMatrix[level-1][row-1][col]])### Case 6
                            
    #################
    #for i in range(len(scoringMatrix)):
    #    for j in range(len(scoringMatrix[i])):
    #        print(scoringMatrix[i][j],"          ", memoMatrix[i][j])
    #    print("")
    #################
    
    ## Do backtracing starting at the bottom right corner.
    resultSeq1 = []
    resultSeq2 = []
    resultSeq3 = []

    col = len(dnaMatrix1)-1
    row = len(dnaMatrix2)-1
    level = len(dnaMatrix3)-1
    done = False

    currentScore = scoringMatrix[level][row][col]
    finalScore = scoringMatrix[level][row][col]
    if currentScore == 0:
        done = True

    while (done == False):
        
        caseOne = scoringMatrix[level][row][col-1]
        caseTwo = scoringMatrix[level][row-1][col]
        caseThree = scoringMatrix[level-1][row][col]
        caseFour = scoringMatrix[level][row-1][col-1]
        caseFive = scoringMatrix[level-1][row][col-1]
        caseSix = scoringMatrix[level-1][row-1][col]

        if currentScore == caseOne:
            col -= 1
        elif currentScore == caseTwo:
            row -= 1
        elif currentScore == caseThree:
            level -= 1
        elif currentScore == caseFour:
            col -= 1
            row -= 1
        elif currentScore == caseFive:
            col -= 1
            level -= 1
        elif currentScore == caseSix:
            row -= 1
            level -= 1
        else:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, dnaMatrix2[row])
            resultSeq3.insert(0, dnaMatrix3[level])
            col -= 1
            row -= 1
            level -= 1

        currentScore = scoringMatrix[level][row][col]
        if (currentScore == 0):
            done = True

    return [resultSeq1, resultSeq2, resultSeq3], finalScore
