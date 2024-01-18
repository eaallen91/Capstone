'''
Eric Allen

Optimal alignment of three sequences

This is the old version that calculates all numbers and then does the back tracing by looking at 7 different cases.
Returns mostly perfect scores but not always the optimal. It also has half the run time of the newer current version
'''
from time import time

'''
Manually finds the smallest number between
all possible cases. It will return the last
smallest number if there are multiple cases
with the same value
'''
def findMin(vals):
    smallestNum = vals[0]

    for i in range(len(vals)):
        if vals[i] <= smallestNum:
            smallestNum = vals[i]

    return smallestNum

def align3D(dnaSeq1, dnaSeq2, dnaSeq3, gap, mismatch, match):

    compStartTime = time()
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

    scoringMatrix = [ [ [0 for i in range(len(dnaMatrix1)) ] for j in range(len(dnaMatrix2)) ] for k in range(len(dnaMatrix3)) ]

    ############################### INITIALIZE THE EDGES FOR THE TOP LEVEL ######################################
    ## Fill the empty string col in the first level
    for col in range(1, len(dnaMatrix1)):
        scoringMatrix[0][0][col] = scoringMatrix[0][0][col-1] + (gap*2)
        
    ## Fill the empty string row in the first level
    for row in range(1, len(dnaMatrix2)):
        scoringMatrix[0][row][0] = scoringMatrix[0][row-1][0] + (gap*2)
    #############################################################################################################


    ############################## INITIALIZE THE EDGES FOR THE MIDDLE LEVELS ###################################
    ## Computes the edges for the interior levels
    for level in range(1, len(dnaMatrix3)):
        scoringMatrix[level][0][0] = scoringMatrix[level-1][0][0] + (gap * 2)
        
        ## Fill the empty string col in the middle levels
        for col in range(1, len(dnaMatrix1)):
            singleGap = scoringMatrix[level-1][0][col-1] + gap
            if dnaMatrix1[col] == dnaMatrix3[level]:
                singleGap = singleGap + match
            else:
                singleGap = singleGap + mismatch

            scoringMatrix[level][0][col] = findMin([scoringMatrix[level][0][col-1] + (gap*2),
                                                    scoringMatrix[level-1][0][col] + (gap*2),
                                                    singleGap])

        ## Fill the empty string row in the middle levels
        for row in range(1, len(dnaMatrix2)):
            singleGap = scoringMatrix[level-1][row-1][0] + gap
            if dnaMatrix2[row] == dnaMatrix3[level]:
                singleGap = singleGap + match
            else:
                singleGap = singleGap + mismatch
            scoringMatrix[level][row][0] = findMin([scoringMatrix[level][row-1][0] + (gap*2),
                                                   scoringMatrix[level-1][row][0] + (gap*2),
                                                   singleGap])
    #############################################################################################################



    ########################## Fill the middle of the matrix for the first level ################################
    for col in range(1, len(dnaMatrix1)):
        for row in range(1, len(dnaMatrix2)):
            singleGap = scoringMatrix[0][row-1][col-1] + gap
            if dnaMatrix1[col] == dnaMatrix2[row]:
                singleGap = singleGap + match
            else:
                singleGap = singleGap + mismatch
            
            scoringMatrix[0][row][col] = findMin([scoringMatrix[0][row][col-1] + (gap*2),
                                                  scoringMatrix[0][row-1][col] + (gap*2),
                                                  singleGap])
    #############################################################################################################

    middleStartTime = time()
    for level in range(1, len(dnaMatrix3)):
        for col in range(1, len(dnaMatrix1)):
            for row in range(1, len(dnaMatrix2)):

                ##### Case 4
                singleGap1 = scoringMatrix[level][row-1][col-1] + gap
                if dnaMatrix1[col] == dnaMatrix2[row]:
                    singleGap1 = singleGap1 + match
                else:
                    singleGap1 = singleGap1 + mismatch

                #### Case 5
                singleGap2 = scoringMatrix[level-1][row][col-1] + gap
                if dnaMatrix1[col] == dnaMatrix3[level]:
                    singleGap2 = singleGap2 + match
                else:
                    singleGap2 = singleGap2 + mismatch

                #### Case 6
                singleGap3 = scoringMatrix[level-1][row-1][col] + gap
                if dnaMatrix2[row] == dnaMatrix3[level]:
                    singleGap3 = singleGap3 + match
                else:
                    singleGap3 = singleGap3 + mismatch

                #### Case 7
                singleGap4 = scoringMatrix[level-1][row-1][col-1]
                if dnaMatrix1[col] == dnaMatrix2[row]:
                    singleGap4 =  singleGap4 + match
                    if dnaMatrix1[col] == dnaMatrix3[level]:
                        singleGap4 = singleGap4 + match
                    else:
                        singleGap4 = singleGap4 + mismatch
                else:
                    singleGap4 = singleGap4 + mismatch
                    if (dnaMatrix3[level] == dnaMatrix1[col]) or (dnaMatrix3 == dnaMatrix2[row]):
                        singleGap4 = singleGap4 + match
                    else:
                        singleGap4 = singleGap4 + mismatch
                        
                scoringMatrix[level][row][col] = findMin([scoringMatrix[level][row][col-1] + (gap*2),    ## Case 1
                                                          scoringMatrix[level][row-1][col] + (gap*2),    ## Case 2
                                                          scoringMatrix[level-1][row][col] + (gap*2),    ## Case 3
                                                          singleGap1,                                    ## Case 4
                                                          singleGap2,                                    ## Case 5
                                                          singleGap3,                                    ## Case 6
                                                          singleGap4])                                    ## Case 7
    
    endMiddleTime = time()
    print("Filling middle of array time = ", endMiddleTime - middleStartTime)
    #################
    #for i in range(len(scoringMatrix)):
    #    for j in range(len(scoringMatrix[i])):
    #        print(scoringMatrix[i][j])
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

    startTime = time()
    done = False
    while(done == False):
        ##### The current score of the scoring matrix and the
        ##### seven possible cases that resulted in the current
        ##### score
        

        currentScore = scoringMatrix[level][row][col]
        caseOne = scoringMatrix[level][row][col-1]
        caseTwo = scoringMatrix[level][row-1][col]
        caseThree = scoringMatrix[level-1][row][col]
        caseFour = scoringMatrix[level][row-1][col-1]
        caseFive = scoringMatrix[level-1][row][col-1]
        caseSix = scoringMatrix[level-1][row-1][col]
        caseSeven = scoringMatrix[level-1][row-1][col-1]
        
        #startTwo = time()
        #### If level = 0 then it is just a 2D matrix on the top of the cube
        if level == 0:
            ######## Edge Cases #######
            ### Case 1
            if row == 0:
                resultSeq1.insert(0, dnaMatrix1[col])
                resultSeq2.insert(0, "-")
                resultSeq3.insert(0, "-")
                col = col - 1
            ### Case 2
            elif col == 0:
                resultSeq1.insert(0, "-")
                resultSeq2.insert(0, dnaMatrix2[row])
                resultSeq3.insert(0, "-")
                row = row - 1
            ###########################
            ### Case 1
            elif currentScore == (caseOne + (gap*2)):
                resultSeq1.insert(0, dnaMatrix1[col])
                resultSeq2.insert(0, "-")
                resultSeq3.insert(0, "-")
                col = col -1
            ### Case 2
            elif currentScore == (caseTwo + (gap*2)):
                resultSeq1.insert(0, "-")
                resultSeq2.insert(0, dnaMatrix2[row])
                resultSeq3.insert(0, "-")
                row = row - 1
            ### Case 4
            else:
                resultSeq1.insert(0, dnaMatrix1[col])
                resultSeq2.insert(0, dnaMatrix2[row])
                resultSeq3.insert(0, "-")
                col = col - 1
                row = row - 1
        else:
            ### Case 3
            if (col == 0) and (row == 0):
                resultSeq1.insert(0, "-")
                resultSeq2.insert(0, "-")
                resultSeq3.insert(0, dnaMatrix3[level])
                level = level - 1
            elif col == 0:
                ### Case 2
                if currentScore == (caseTwo + (gap*2)):
                    resultSeq1.insert(0, "-")
                    resultSeq2.insert(0, dnaMatrix2[row])
                    resultSeq3.insert(0, "-")
                    row = row - 1
                ### Case 3
                elif currentScore == (caseThree + (gap*2)):
                    resultSeq1.insert(0, "-")
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, dnaMatrix3[level])
                    level = level - 1
                ### Case 5
                else:
                    resultSeq1.insert(0, "-")
                    resultSeq2.insert(0, dnaMatrix2[row])
                    resultSeq3.insert(0, dnaMatrix3[level])
                    row = row - 1
                    level = level -1
            elif row == 0:
                ### Case 1
                if currentScore == (caseOne + (gap*2)):
                    resultSeq1.insert(0, dnaMatrix1[col])
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, "-")
                    col = col - 1
                ### Case 3
                elif currentScore == (caseThree + (gap*2)):
                    resultSeq1.insert(0, "-")
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, dnaMatrix3[level])
                    level = level - 1
                ### Case 4
                else:
                    resultSeq1.insert(0, dnaMatrix1[col])
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, dnaMatrix3[level])
                    col = col - 1
                    level = level - 1
            else:
                ### Case 7
                #if (currentScore == (caseSeven + (match*2))) or (currentScore == (caseSeven +(mismatch*2))) or (currentScore == (caseSeven + match + mismatch)):
                #    resultSeq1.insert(0, dnaMatrix1[col])
                #    resultSeq2.insert(0, dnaMatrix2[row])
                #    resultSeq3.insert(0, dnaMatrix3[level])
                #    col = col - 1
                #    row = row - 1
                #    level = level - 1
                ### Case 1
                if currentScore == (caseOne + (gap*2)):
                    resultSeq1.insert(0, dnaMatrix1[col])
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, "-")
                    col = col - 1
                ### Case 2
                elif currentScore == (caseTwo + (gap*2)):
                    resultSeq1.insert(0, "-")
                    resultSeq2.insert(0, dnaMatrix2[row])
                    resultSeq3.insert(0, "-")
                    row = row - 1
                ### Case 3
                elif currentScore == (caseThree + (gap*2)):
                    resultSeq1.insert(0, "-")
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, dnaMatrix3[level])
                    level = level - 1
                ### Case 4
                elif (currentScore == (caseFour + gap + match)) or (currentScore == (caseFour + gap + mismatch)):
                    resultSeq1.insert(0, dnaMatrix1[col])
                    resultSeq2.insert(0, dnaMatrix2[row])
                    resultSeq3.insert(0, "-")
                    col = col - 1
                    row = row - 1
                ### Case 5
                elif (currentScore == (caseFive + gap + match)) or (currentScore == (caseFive + gap + mismatch)):
                    resultSeq1.insert(0, dnaMatrix1[col])
                    resultSeq2.insert(0, "-")
                    resultSeq3.insert(0, dnaMatrix3[level])
                    col = col - 1
                    level = level - 1
                ### Case 6
                elif (currentScore == (caseSix + gap + match)) or (currentScore == (caseSix + gap + mismatch)):
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
        #endTwo = time()
        #print("Previous case time = ", endTwo - startTwo)
        if (row == 0) and (col == 0) and (level == 0):
            done = True
    
    endTime = time()
    print("Backtrace time = ", endTime - startTime)                                                    

    finalSeqs = [resultSeq1, resultSeq2, resultSeq3]
    compEndTime = time()
    print("Comparison Run Time = ", compEndTime - compStartTime)
    return finalSeqs
