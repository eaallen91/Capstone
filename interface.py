'''
Eric Allen

The GUI interface for the program that allows the user 
to specify the type of alignment, reading in from a file
or through manual input, writing to a file, printing 
on screen, specifying the gap, mismatch, and match penalties
and the ablitity to choose between max or min alignments
'''
'''
Import the tkinter packages needed to create the 
GUI
'''
from tkinter import *
from tkinter import ttk

'''
Import the different alignment files, so that the user
may align global, semi-global, Longest common subsequence,
or local
'''
from Memo import align3D
from Semi_Memo import alignSemi3D
from LCS import LCS3D
from Local import local3D

''' 
Adds the ability to read the strings in from a file
'''
from readFromFile import readFromFile

'''
If the read from file checkbox is selected/unselected
the states of the the entry fields need to change
'''
def flipSeqStates():
    # Flip the states of the seq entries
    if inputSeq1Entry.instate(['!disabled']):
        inputSeq1Entry.state(['disabled'])
        inputSeq2Entry.state(['disabled'])
        inputSeq3Entry.state(['disabled'])
        fileEntry.state(['!disabled'])
    else:
        inputSeq1Entry.state(['!disabled'])
        inputSeq2Entry.state(['!disabled'])
        inputSeq3Entry.state(['!disabled'])
        fileEntry.state(['disabled'])

'''
If the checkbox for writing to a file is clicked then the
state of the checkbox needs to switch accordingly
'''
def flipWrite():
    if outputEntry.instate(['!disabled']):
        outputEntry.state(['disabled'])
    else:
        outputEntry.state(['!disabled'])


'''
Depending on which form of alignment is chosen, flip on or off
the proper entry boxes.

Longest Common Subsequence does not require gap, mismatch, or
match penalties, it also does not require min or max distinction

Local alignment does not require min max distinction

Semi and Global will make sure that all fields are turned back on
'''
def flipPenalties():
    entries = [gapEntry, mismatchEntry, matchEntry, maxBox]
    
    ### If Longest Common Subsequence was selected, disable
    ### all penalty boxes, and the maxBox option
    if alignment.get() == 'LCS':
        for item in entries:
            if item.instate(['!disabled']):
                item.state(['disabled'])
    ### If Local was selected, turn off maxBox, and turn on
    ### anything else that is currently disabled
    elif alignment.get() == 'Local':
        for item in entries:
            if item.instate(['disabled']):
                if item != maxBox:
                    item.state(['!disabled'])
            else:
                if item == maxBox:
                    item.state(['disabled'])
    ### For Semi and Global turn all options back on
    else:
        for item in entries:
            if item.instate(['disabled']):
                item.state(['!disabled'])

'''
Blank the subframe so that it may be used again for
the next round of displays
'''     
def blankSubframe():
    del Seq1Str[:]
    del Seq2Str[:]
    del Seq3Str[:]

    Seq1Str.append(StringVar())
    Seq2Str.append(StringVar())
    Seq3Str.append(StringVar())

'''
Displays any messages that need to inform the user
'''
def displayMessage(messages):
    for i in range(len(messages)):
        seqStrs[i][0].set(messages[i])
        printLabel = ttk.Label(subframe, textvariable=seqStrs[i][0])
        printLabel.grid(column=0, row=i)

'''
Retrieves all of the information that the user provided
and passes it into the proper sequence alignment. It
then either displays the aligned seqs or it writes it
to a user specified file
'''
def align(*args):
    try:
        ### If the file button is not selected then get the 
        ### seq entries from the text fields, otherwise read
        ### in the file specified by the user
        if fileButton.get() == 0:
            alignedSeqs = []
            seq1 = inputSeq1.get()
            seq2 = inputSeq2.get()
            seq3 = inputSeq3.get()
        else:
            seqs = readFromFile(fileName.get())
            if (not(seqs)):
                blankSubframe()
                messages = []
                message = "There was an error with you input file"
                messages.append(message)
                messages.append("Please review the correct format and retry")
                displayMessage(messages)
                return
            elif(len(seqs) < 3):
                blankSubframe()
                messages = []
                messages.append("There are too few seqs in your file")
                messages.append("Please review, and correct mistakes")
                displayMessage(messages)
                return
            seq1 = seqs[0][1]
            seq2 = seqs[1][1]
            seq3 = seqs[2][1]
            
        ### Retrive the gap, mismatch, and match values
        gapVal = gap.get()
        mismatchVal = mismatch.get()
        matchVal = match.get()
        
        ### Get the minMax value, if it is blank, as StringVal defaults to,
        ### set the value to "Min"
        if (minMax.get() == ''):
            minMax.set("Min")
        
        # Gets the type of alignment and then properly chooses the 
        # the corresponding alignment method
        alignType = alignment.get()
        
        if alignType == "LCS":
            alignedSeqs,finalScore = LCS3D(seq1, seq2, seq3)
        elif alignType == "Semi":
            alignedSeqs,finalScore = alignSemi3D(seq1, seq2, seq3, gapVal, mismatchVal, matchVal, minMax.get())
        elif alignType == "Local":
            alignedSeqs,finalScore = local3D(seq1, seq2, seq3, gapVal, mismatchVal, matchVal)
        else:
            alignedSeqs,finalScore = align3D(seq1, seq2, seq3, gapVal, mismatchVal, matchVal, minMax.get())

        # If the output file checkbox was marked then we open the file
        # the user specified, if it does not exist it is created,
        # and write the aligned seqs to the file
        if outputButton.get() == 1:
            try:
                with open(outputFile.get()): pass
            except IOError:
                blankSubframe()
                messages = []
                messages.append("The file you are trying to write to is invalid")
                displayMessage(messages)
                return False
            outFile = open(outputFile.get(), "a")
            # Format and provide information regarding the test parameters
            outFile.write(">" + alignType +" Final Score: "+ str(finalScore) + ", Gap: "+ str(gap.get()) + ", Mismatch: " + str(mismatch.get())
                          + ", Match: " + str(match.get()) + ", " + minMax.get())
            if fileButton.get() == 1:
                outFile.write(", from file " + fileName.get() + "\n")
            else:
                outFile.write(", from user input\n")
            for k in range(len(alignedSeqs)):
                letterCounter = 0
                ## Write to the file, a letter at a time, and insert a newline
                ## every 80 characters
                for letter in alignedSeqs[k]:
                    outFile.write(letter)
                    letterCounter += 1
                    if letterCounter == 80:
                        outFile.write("\n");
                        letterCounter = 0
                outFile.write("\n\n");
                
            ### Once the file is written we blank the subframe and put the result label
            ### on it
            blankSubframe()
            messages = []
            messages.append("Seqs aligned and printed to file " + outputFile.get())
            displayMessage(messages)
            outFile.close()
        ### If the user did not specify an output file, then print on the GUI
        else:
            finalScoreStr.set(finalScore)
            blankSubframe()
            for i in range(len(alignedSeqs)):
                for k in range(len(alignedSeqs[i])):
                    seqStrs[i].append(StringVar())
                    seqStrs[i][k].set(alignedSeqs[i][k])
            for i in range(len(Seq1Str)):
                printLabel = ttk.Label(subframe, textvariable=Seq1Str[i])
                printLabel.grid(column=i, row=0)

            for i in range(len(Seq2Str)):
                printLabel = ttk.Label(subframe, textvariable=Seq2Str[i])
                printLabel.grid(column=i, row=1)

            for i in range(len(Seq3Str)):
                printLabel = ttk.Label(subframe, textvariable=Seq3Str[i])
                printLabel.grid(column=i, row=2)
    except ValueError:
        blankSubframe()
        messages = []
        messages.append("There was an error with your input")
        messages.append("Please review your values")
        displayMessage(messages)
        pass


### The master frame that all widgets reside within
root = Tk()
root.title("3D Sequence Alignment Tool")


Seq1Str = [StringVar()]
Seq2Str = [StringVar()]
Seq3Str = [StringVar()]
seqStrs = [Seq1Str, Seq2Str, Seq3Str]
finalScoreStr = StringVar()
### All of the variables used by the various parts of
### the interface.

### The seqs input by the user
inputSeq1 = StringVar()
inputSeq2 = StringVar()
inputSeq3 = StringVar()

### Gap, mismatch, match values, specified by the
### user
gap = IntVar()
mismatch = IntVar()
match = IntVar()

### Checkboxes and radiobutton variables,
minMax = StringVar()
alignment = StringVar()
fileButton = IntVar()
outputButton = IntVar()

### The file names for input and output
fileName = StringVar()
outputFile = StringVar()

# Create the root frame that everything will go on
mainframe = ttk.Frame(root, padding=(10, 10, 10, 10))
mainframe.pack(fill="x")

# Create the radio buttons that allow switching between alignment types
globalAlign = ttk.Radiobutton(mainframe, text="Global (default choice)", variable=alignment, command=flipPenalties, value="Global")
semiAlign = ttk.Radiobutton(mainframe, text="Semi-global", variable=alignment, command=flipPenalties, value="Semi")
lcsAlign = ttk.Radiobutton(mainframe, text="Longest Common Subsequence", variable=alignment, command=flipPenalties, value="LCS")
localAlign = ttk.Radiobutton(mainframe, text="Local", variable=alignment, command=flipPenalties, value="Local")

# Put the radio buttons onto the frame
globalAlign.grid(column=1, row=1, columnspan=2, sticky=W)
semiAlign.grid(column=1, row=2, sticky=W)
lcsAlign.grid(column=1, row=3, columnspan=2, sticky=W)
localAlign.grid(column=1, row=4, sticky=W)

# Create the readfile checkbox that allows the program to read from a file
readFile = ttk.Checkbutton(mainframe, text='Read from a file', command=flipSeqStates, 
                           variable=fileButton, onvalue='1', offvalue='0')
readFile.grid(column=1, row=5, sticky=W)

# Create the file name entry field
ttk.Label(mainframe, text="Input file").grid(column=2, row=5, sticky=(W, E))
fileEntry = ttk.Entry(mainframe, textvariable=fileName)
fileEntry.state(['disabled'])
fileEntry.grid(column=3, row=5, columnspan=3)

# Create the writefile checkbox that allows the program to write the results
# to a file
writeFile = ttk.Checkbutton(mainframe, text='Write to a file', command=flipWrite,
                            variable=outputButton, onvalue='1', offvalue='0')
writeFile.grid(column=1, row=6, sticky=W)

# Create the output file entry field
ttk.Label(mainframe, text="Output file").grid(column=2, row=6, sticky=(W, E)) 
outputEntry = ttk.Entry(mainframe, textvariable=outputFile)
outputEntry.state(['disabled'])
outputEntry.grid(column=3, row=6, columnspan=3)

# Create the entry fields for the 3 seqs, if not reading from a file
inputSeq1Entry = ttk.Entry(mainframe, textvariable=inputSeq1)
inputSeq2Entry = ttk.Entry(mainframe, textvariable=inputSeq2)
inputSeq3Entry = ttk.Entry(mainframe, textvariable=inputSeq3)

# Put the string 1 entry onto the frame
inputSeq1Entry.grid(column=2, row=7, columnspan=10, sticky=(W, E))
ttk.Label(mainframe, text="Sequence 1").grid(column=1, row=7, sticky=W)

# Put the string 2 entry onto the frame
inputSeq2Entry.grid(column=2, row=8, columnspan=10, sticky=(W, E))
ttk.Label(mainframe, text="Sequence 2").grid(column=1, row=8, sticky=W)

# Put the string3 entry onto the frame
inputSeq3Entry.grid(column=2, row=9, columnspan=10, sticky=(W, E))
ttk.Label(mainframe, text="Sequence 3").grid(column=1, row=9, sticky=W)

# Create the entry fields for gap, mismatch, and match penalties
gapEntry = ttk.Entry(mainframe, width=7, textvariable=gap)
mismatchEntry = ttk.Entry(mainframe, width=7, textvariable=mismatch)
matchEntry = ttk.Entry(mainframe, width=7, textvariable=match)

# Place the gap entry field onto the frame
ttk.Label(mainframe, text="Gap").grid(column=1, row=10, sticky=W)
gapEntry.grid(column=1, row=10)

# Place the mismatch entry field onto the frame
ttk.Label(mainframe, text="Mismatch").grid(column=2, row=10, sticky=W)
mismatchEntry.grid(column=3, row=10)

# Place the match entry field onto the frame
ttk.Label(mainframe, text="Match").grid(column=5, row=10, sticky=W)
matchEntry.grid(column=6, row=10)

# Place the min max checkbox on the frame
maxBox = ttk.Checkbutton(mainframe, text="Max(default is min)", variable=minMax, onvalue="Max", offvalue="Min")
maxBox.grid(column=1, row=13, columnspan=15)

### The command button that causes everything to process and happen
ttk.Button(mainframe, text="Align", command=align).grid(column=1, row=14, columnspan=15)

finalScoreDisplay = ttk.Label(mainframe, textvariable=finalScoreStr).grid(column=1, row=15, columnspan=15)
ttk.Label(mainframe, text="Final Score").grid(column=0, row=15, columnspan=3, sticky=E)
### Subframe is located below the align button and it holds the results of aligned seqs
### or of any notification messages
subframe = ttk.Frame(root, borderwidth=5, relief="sunken")
subframe.pack(expand="yes", fill="both")

### Loop through all the elements in the window and add padding to them
for child in mainframe.winfo_children(): 
    child.grid_configure(padx=5, pady=5)

### Used to put blank labels onto the bottom frame.
### By starting it off this way I am able to update the textvariables
### and dynamically change what is displayed
for i in range(len(Seq1Str)):
    printLabel = ttk.Label(subframe, textvariable=Seq1Str[i])
    printLabel.grid(column=i, row=0)

for i in range(len(Seq2Str)):
    printLabel = ttk.Label(subframe, textvariable=Seq2Str[i])
    printLabel.grid(column=i, row=1)

for i in range(len(Seq3Str)):
    printLabel = ttk.Label(subframe, textvariable=Seq3Str[i])
    printLabel.grid(column=i, row=2)

root.mainloop()
