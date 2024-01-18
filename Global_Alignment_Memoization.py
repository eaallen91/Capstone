'''
Eric Allen

Optimal alignment of three sequences

The most current testing version. Each position contains an array of [val, case].
Once the best val has been determined then the case that generated the val is stored.
The path back through the matrix is then known.
'''

'''
Finds the smallest val between the seven cases
It returns the FIRST smallestNum and its case if
ther are multiple with the same val
'''
def findMin(vals):

    nums = []
    cases = []
    ### Seperates the vals and the cases
    for i in range(len(vals)):
        nums.append(vals[i][0])
        cases.append(vals[i][1])

    ### Looks for the smallest val and its
    ### corresponding case
    smallestNum = nums[0]
    case = cases[0]
    for i in (range(len(nums))):
        if nums[i] <= smallestNum:
            smallestNum = nums[i]
            case = cases[i]

    return(smallestNum, case)

'''
3D Sequence alignment that will fill a 3D
matrix and with the help of a helper variable
backtrace through the matrix and output the
optimal solution
'''
def align3D(dnaSeq1, dnaSeq2, dnaSeq3, gap, mismatch, match):

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

    ### scoringMatix[level][row][col][data]
    ### [data] = [scoreVal, case]
    scoringMatrix = [ [ [ [0,0] for i in range(len(dnaMatrix1)) ] for j in range(len(dnaMatrix2)) ] for k in range(len(dnaMatrix3)) ]

    ############################### INITIALIZE THE EDGES FOR THE TOP LEVEL ######################################
    ## Fill the empty string col in the first level
    for col in range(1, len(dnaMatrix1)):
        scoringMatrix[0][0][col][0] = scoringMatrix[0][0][col-1][0] + (gap*2)
        scoringMatrix[0][0][col][1] = 1
        
    ## Fill the empty string row in the first level
    for row in range(1, len(dnaMatrix2)):
        scoringMatrix[0][row][0][0] = scoringMatrix[0][row-1][0][0] + (gap*2)
        scoringMatrix[0][row][0][1] = 2
    #############################################################################################################


    ############################## INITIALIZE THE EDGES FOR THE MIDDLE LEVELS ###################################
    ## Computes the edges for the interior levels
    for level in range(1, len(dnaMatrix3)):
        scoringMatrix[level][0][0][0] = scoringMatrix[level-1][0][0][0] + (gap * 2)
        scoringMatrix[level][0][0][1] = 3
        
        ## Fill the empty string col in the middle levels
        for col in range(1, len(dnaMatrix1)):
            caseOne = scoringMatrix[level][0][col-1][0] + (gap*2)
            caseThree = scoringMatrix[level-1][0][col][0] + (gap*2)
            singleGap = scoringMatrix[level-1][0][col-1][0] + gap
            if dnaMatrix1[col] == dnaMatrix3[level]:
                singleGap = singleGap + match
            else:
                singleGap = singleGap + mismatch

            cases = [[caseOne, 1], [caseThree, 3], [singleGap, 5]]
            scoringData = findMin(cases)
            scoringMatrix[level][0][col][0] = scoringData[0]
            scoringMatrix[level][0][col][1] = scoringData[1]


        ## Fill the empty string row in the middle levels
        for row in range(1, len(dnaMatrix2)):
            caseTwo = scoringMatrix[level][row-1][0][0] + (gap*2)
            caseThree = scoringMatrix[level-1][row][0][0] + (gap*2)
            singleGap = scoringMatrix[level-1][row-1][0][0] + gap
            if dnaMatrix2[row] == dnaMatrix3[level]:
                singleGap = singleGap + match
            else:
                singleGap = singleGap + mismatch
            cases = [[caseTwo, 2], [caseThree, 3], [singleGap, 6]]
            scoringData = findMin(cases)
            scoringMatrix[level][row][0][0] = scoringData[0]
            scoringMatrix[level][row][0][1] = scoringData[1]
    #############################################################################################################



    ########################## Fill the middle of the matrix for the first level ################################
    for col in range(1, len(dnaMatrix1)):
        for row in range(1, len(dnaMatrix2)):
            caseOne = scoringMatrix[0][row][col-1][0] + (gap*2)
            caseTwo = scoringMatrix[0][row-1][col][0] + (gap*2)
            singleGap = scoringMatrix[0][row-1][col-1][0] + gap
            if dnaMatrix1[col] == dnaMatrix2[row]:
                singleGap = singleGap + match
            else:
                singleGap = singleGap + mismatch

            cases = [[caseOne, 1], [caseTwo, 2], [singleGap, 4]]
            scoringData = findMin(cases)
            
            scoringMatrix[0][row][col][0] = scoringData[0]
            scoringMatrix[0][row][col][1] = scoringData[1]
            
            
    #############################################################################################################

    for level in range(1, len(dnaMatrix3)):
        for col in range(1, len(dnaMatrix1)):
            for row in range(1, len(dnaMatrix2)):
                caseOne = scoringMatrix[level][row][col-1][0] + (gap*2)   ## Case 1
                caseTwo = scoringMatrix[level][row-1][col][0] + (gap*2)   ## Case 2
                caseThree = scoringMatrix[level-1][row][col][0] + (gap*2) ## Case 3

                ##### Case 4
                singleGap1 = scoringMatrix[level][row-1][col-1][0] + gap
                if dnaMatrix1[col] == dnaMatrix2[row]:
                    singleGap1 = singleGap1 + match
                else:
                    singleGap1 = singleGap1 + mismatch

                #### Case 5
                singleGap2 = scoringMatrix[level-1][row][col-1][0] + gap
                if dnaMatrix1[col] == dnaMatrix3[level]:
                    singleGap2 = singleGap2 + match
                else:
                    singleGap2 = singleGap2 + mismatch

                #### Case 6
                singleGap3 = scoringMatrix[level-1][row-1][col][0] + gap
                if dnaMatrix2[row] == dnaMatrix3[level]:
                    singleGap3 = singleGap3 + match
                else:
                    singleGap3 = singleGap3 + mismatch

                #### Case 7
                singleGap4 = scoringMatrix[level-1][row-1][col-1][0]
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

                cases = [[caseOne, 1], [caseTwo, 2], [caseThree, 3], [singleGap1, 4], [singleGap2, 5], [singleGap3, 6], [singleGap4, 7]]
                scoringData = findMin(cases)
                scoringMatrix[level][row][col][0] = scoringData[0]
                scoringMatrix[level][row][col][1] = scoringData[1]
                
                
    #################
    for i in range(len(scoringMatrix)):
        for j in range(len(scoringMatrix[i])):
            print(scoringMatrix[i][j])
        print("")
    #################
    
    ## Do backtracing starting at the bottom right corner.
    resultSeq1 = []
    resultSeq2 = []
    resultSeq3 = []

    col = len(dnaMatrix1)-1
    row = len(dnaMatrix2)-1
    level = len(dnaMatrix3)-1

    finalScore = scoringMatrix[level][row][col][0]

    done = False
    while(done == False):
        print(level, row, col, scoringMatrix[level][row][col][0], scoringMatrix[level][row][col][1])
        print("")
        currentCase = scoringMatrix[level][row][col][1]

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
        
                                                        

    finalSeqs = [resultSeq1, resultSeq2, resultSeq3]
    return finalSeqs
