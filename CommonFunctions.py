'''
Finds the smallest val between the cases. It returns
the LAST smallestNum and its case if there are
multiple with the same val
'''
def findMin(vals):

    smallestNum = vals[0]
    case = vals[1]

    ### Looks for the smallest val and its 
    ### corresponding case
    for i in range(2, len(vals), 2):
        if vals[i] <= smallestNum:
            smallestNum = vals[i]
            case = vals[i+1]

    return(smallestNum, case)

'''
Finds the largest val between the cases. It returns
the LAST largestNum and its case if there are
multiple with the same val
'''
def findMax(vals):
    
    largestNum = vals[0]
    case = vals[1]
    
    ### Looks for the largest val and its
    ### corresponding case
    for i in range(2, len(vals), 2):
        if vals[i] >= largestNum:
            largestNum = vals[i]
            case = vals[i + 1]

    return(largestNum, case)

'''
Used to switch between min and max depending on the
user's preference so that there is not execessive
amounts of duplicated code
'''
def chooseMinMax(vals, minMax):
    
    resultNum = 0
    resultCase = 0
    if (minMax == "Min"):
        resultNum, resultCase = findMin(vals)
    else:
        resultNum, resultCase = findMax(vals)

    return(resultNum, resultCase)
