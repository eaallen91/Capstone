'''
Eric Allen

Optimal alignment of three sequences

The most current testing version. Each position contains an array of [val, case].
Once the best val has been determined then the case that generated the val is stored.
The path back through the matrix is then known.
'''

from CommonFunctions import chooseMinMax

'''
These are the cases that I refer to as case 1, case 2, . . . case 7

Case 1   Case 2   Case 3   Case 4   Case 5   Case 6   Case 7
  X        -       -         X        X        -        X   
  -        X       -         X        -        X        X
  -        -       X         -        X        X        X
'''

'''
3D Sequence alignment that will fill a 3D
matrix and with the help of a helper variable
backtrace through the matrix and output the
optimal solution
'''
def newSemi3D(dnaSeq1, dnaSeq2, dnaSeq3, gap, mismatch, match, minMax):
    
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
    dnaMatrix1.insert(0, "-")
    dnaMatrix2.insert(0, "-")
    dnaMatrix3.insert(0, "-")

    ### Create a 3D Scoring matrix and a 3D memo matrix
    scoringMatrix = [ [ [ 0 for i in range(len(dnaMatrix1)) ] for j in range(len(dnaMatrix2)) ] for k in range(len(dnaMatrix3)) ]
    memoMatrix =    [ [ [ 0 for i in range(len(dnaMatrix1)) ] for j in range(len(dnaMatrix2)) ] for k in range(len(dnaMatrix3)) ]

    ### Fills the entire matrix and handles all special cases
    for level in range(0, len(dnaMatrix3)):
        for col in range(0, len(dnaMatrix1)):
            for row in range(0, len(dnaMatrix2)):

                oneToTwo = bool(dnaMatrix1[col] == dnaMatrix2[row])
                oneToThree = bool(dnaMatrix1[col] == dnaMatrix3[level])
                twoToThree = bool(dnaMatrix2[row] == dnaMatrix3[level])
                
                ##### Case 4
                caseFour = 0
                if oneToTwo:
                    caseFour =  match
                else:
                    caseFour =  mismatch

                #### Case 5
                caseFive = 0
                if oneToThree:
                    caseFive =  match
                else:
                    caseFive =  mismatch

                #### Case 6
                caseSix = 0
                if twoToThree:
                    caseSix =  match
                else:
                    caseSix =  mismatch

                #### Case 7
                caseSeven = 0
                if oneToTwo:
                    if oneToThree:
                        caseSeven = match * 2
                    else:
                        caseSeven = match + mismatch
                else:
                    if oneToThree:
                        caseSeven = match + mismatch
                    else:
                        caseSeven = mismatch * 2

                ### Special cases exist for the first level of the 3D matrix
                ### it is handled as a 2D matrix
                if level == 0:
                    if (row == 0) and (col == 0):
                        memoMatrix[level][row][col] = 0
                    elif row == 0:
                        memoMatrix[level][row][col] = 1
                    elif col == 0:
                        memoMatrix[level][row][col] = 2
                    else:
                        scoringMatrix[level][row][col], memoMatrix[level][row][col] = chooseMinMax([scoringMatrix[level][row][col-1] + (gap*2), 1,
                                                                                                    scoringMatrix[level][row-1][col] + (gap*2), 2,
                                                                                                    scoringMatrix[level][row-1][col-1] + gap + caseFour, 4],
                                                                                                   minMax)
                ### There are special cases for edges throughout the rest of the matrix
                ### but other than that there are no special cases. The majority of the 
                ### the cells in the matricies will be the final else.
                else:
                    if (row == 0) and (col == 0):
                        memoMatrix[level][row][col] = 3
                    elif row == 0:
                        scoringMatrix[level][row][col], memoMatrix[level][row][col] = chooseMinMax([scoringMatrix[level][row][col-1] + (gap*2), 1,
                                                                                               scoringMatrix[level-1][row][col] + (gap*2), 3,
                                                                                               scoringMatrix[level-1][row][col-1] + gap + caseFive, 5],
                                                                                              minMax)
                    elif col == 0:
                        scoringMatrix[level][row][col], memoMatrix[level][row][col] = chooseMinMax([scoringMatrix[level][row-1][col] + (gap*2), 2,
                                                                                               scoringMatrix[level-1][row][col] + (gap*2), 3,
                                                                                               scoringMatrix[level-1][row-1][col] + gap + caseSix, 6],
                                                                                              minMax)
                    else:
                        scoringMatrix[level][row][col], memoMatrix[level][row][col] = chooseMinMax([scoringMatrix[level][row][col-1] + (gap*2), 1,     ### Case 1
                                                                                               scoringMatrix[level][row-1][col] + (gap*2), 2,          ### Case 2
                                                                                               scoringMatrix[level-1][row][col] + (gap*2), 3,          ### Case 3
                                                                                               scoringMatrix[level][row-1][col-1] + gap + caseFour, 4, ### Case 4
                                                                                               scoringMatrix[level-1][row][col-1] + gap + caseFive, 5, ### Case 5
                                                                                               scoringMatrix[level-1][row-1][col] + gap + caseSix, 6,  ### Case 6
                                                                                               scoringMatrix[level-1][row-1][col-1] + caseSeven, 7],   ### Case 7
                                                                                              minMax)
    #################
    #for i in range(len(scoringMatrix)):
    #    for j in range(len(scoringMatrix[i])):
    #        print(scoringMatrix[i][j], "            ", memoMatrix[i][j])
    #    print("")
    #################
    
    ## Do backtracing starting at the bottom right corner.
    resultSeq1 = []
    resultSeq2 = []
    resultSeq3 = []

    col = len(dnaMatrix1)-1
    row = len(dnaMatrix2)-1
    level = len(dnaMatrix3)-1

    finalScore = scoringMatrix[level][row][col]

    done = False
    while(done == False):
        
        currentCase = memoMatrix[level][row][col]
        
        ### Case 1
        if currentCase == 1:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, "-")
            resultSeq3.insert(0, "-")
            col = col - 1
        ### Case 2
        elif currentCase == 2:
            resultSeq1.insert(0, "-")
            resultSeq2.insert(0, dnaMatrix2[row])
            resultSeq3.insert(0, "-")
            row = row - 1
        ### Case 3
        elif currentCase == 3:
            resultSeq1.insert(0, "-")
            resultSeq2.insert(0, "-")
            resultSeq3.insert(0, dnaMatrix3[level])
            level = level - 1
        ### Case 4
        elif currentCase == 4:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, dnaMatrix2[row])
            resultSeq3.insert(0, "-")
            col = col - 1
            row = row - 1
        ### Case 5
        elif currentCase == 5:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, "-")
            resultSeq3.insert(0, dnaMatrix3[level])
            col = col - 1
            level = level - 1
        ### Case 6
        elif currentCase == 6:
            resultSeq1.insert(0, "-")
            resultSeq2.insert(0, dnaMatrix2[row])
            resultSeq3.insert(0, dnaMatrix3[level])
            row = row - 1
            level = level - 1
        ### Case 7
        else:
            resultSeq1.insert(0, dnaMatrix1[col])
            resultSeq2.insert(0, dnaMatrix2[row])
            resultSeq3.insert(0, dnaMatrix3[level])
            col = col - 1
            row = row - 1
            level = level - 1

        if (row == 0) and (col == 0) and (level == 0):
            done = True

    return [resultSeq1, resultSeq2, resultSeq3]
