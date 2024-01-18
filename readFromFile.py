'''
Reads a file that is a simplified fasta format.
The labels for each seq is on its own line that
is indicated by a ">". There should only be 3 seqs
per file. Any more and they will not get processed,
any less and the program will return an error message
to that there is not enough seqs
'''

'''
Opens a file and reads line by line,
storing the contents into the proper
arrays for later use
'''
def readFromFile(filename):
    resultList = []

    ### Test to make sure the file exists before opening it
    try:
        with open(filename): pass
    except IOError:
        return False

    ### Open the file for reading format
    infile = open(filename, 'r')

    ### Read in the first line. It should be a header line
    line = infile.readline()
    
    ### If the first character of the file is not a ">"
    ### then return false
    if line[0] != ">":
        return False

    ### Retrieve the header for the first seq
    label = line[1:]

    seq = ''
    
    ### Loop through the file and obtain the seqs that are
    ### contained within them
    for line in infile:
        line = line.rstrip()

        ### Skip empty lines
        if line == '':
            continue
        ### Once a new seq label has been encounted, store the
        ### old seq and label, and get the new label
        elif line[0] == '>':
            resultList.append([label, seq])
            label = line[1:]

            seq = ''
        ### Otherwise add the next line on to seq
        else:
            seq += line

    ### Close the file, append the last label and seq, and return
    infile.close()
    resultList.append([label, seq])
    return resultList
