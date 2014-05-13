# Import Tkinter to be able to graphically retrieve the directory
from Tkinter import Tk
from tkFileDialog import askdirectory
from os import listdir, mkdir
from os.path import exists
from sets import Set
import re
import collections
import csv
import shutil

# Get the directory
directory = askdirectory()
files = listdir(directory)
txtFiles = []

# Get the txt files only
for f in files:
    if ".txt" in f:
        txtFiles.append(f)

# Get the plate numbers
numbers = []
for f in txtFiles:
    numbers.append(f.split("_")[0])

# Function to get the drug name for file naming later on


def nameFile(f):
    number = f.split("_")[0]
    drug = f.split("_")[1].split(".")[0]
    return directory + "/" + number + "_" + drug + "_stitched.txt"


def nameFile2(f):
    number = f.split("_")[0]
    drug = f.split("_")[1].split(".")[0]
    return directory + "/" + number + "_" + drug + "_complete.txt"


# Get the duplicate plate numbers
duplicates = [[n, times] for n, times in collections.Counter(numbers).items() if times > 1]

if len(duplicates) != 0:
    folderName = directory + "/UnstitchedData"
    if not exists(folderName):
        mkdir(folderName)

incFolderName = directory + "/IncompleteData"
if not exists(incFolderName):
    mkdir(incFolderName)

allRows = ["A", "B", "C", "D", "E", "F", "G", "H"]
allCols = range(1, 13)
allWells = Set()
for n in range(0, len(allRows)):
    for p in range(0, len(allCols)):
        allWells.add((allRows[n], str(allCols[p])))

# Handle the duplicate files
for i in range(0, len(duplicates)):
    skip = False
    drugNames = Set()
    answer = "y"
    for m in range(0, len(txtFiles)):
        if duplicates[i][0] in txtFiles[m] and "stitched" in txtFiles[m]:
            skip = True
        if duplicates[i][0] in txtFiles[m]:
            drugNames.add(txtFiles[m].split("_")[1].split(".")[0])
    if len(drugNames) != 1:
        string = ("All the drug names for " + str(duplicates[i][0]) +
                  " do not match, they include " + str(drugNames) +
                  ". Do you still want to stitch these files? (y/n): ")
        answer = raw_input(string)
    if answer == "y" or answer == "yes":
        if not skip:
            fileCount = 0
            allData = []
            wells = Set()
            for j in range(0, len(txtFiles)):
                if duplicates[i][0] in txtFiles[j]:
                    fileCount += 1
                    fileName = directory + "/" + txtFiles[j]
                    with open(fileName, "r") as f:
                        reader = csv.reader(f, delimiter="\t")
                        if fileCount == 1:
                            name = nameFile(txtFiles[j])
                            header = reader.next()
                            allData.append(header)
                        for row in reader:
                            if len(row) == len(header):
                                allData.append(row)
                    shutil.move(fileName, folderName)

            for k in allData:
                if k == header:
                    allData.remove(k)
                else:
                    well = k[2], k[3]
                    wells.add(well)

            missingWells = allWells.difference(wells)

            for q in missingWells:
                newRow = ["-1", "-1", str(q[0]), str(q[1]), "N",
                          "-1", "0", "-1", "-1", "-1", "-1", "-1",
                          "-1", "-1", "-1", "-1", "-1", "-1", "-1",
                          "-1", "-1", "-1", "-1", "-1", "-1", "-1"]
                allData.append(newRow)

            newFile = open(name, "w+")
            newFile.write("\t".join(header) + "\n")

            for l in allData:
                newFile.write("\t".join(l) + "\n")

            newFile.close()

            print ("Stitched " + str(duplicates[i][1]) + " " +
                   str(duplicates[i][1]) + " entries together and added dummy data for " +
                   str(len(missingWells)) + " missing wells. Original files moved to " +
                   folderName + ".\n")

# Check for the presence of every well in all of the files
for i in range(0, len(txtFiles)):
    fileName = directory + "/" + txtFiles[i]
    allData = []
    wells = Set()
    with open(fileName, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = reader.next()
        allData.append(header)
        for row in reader:
            if len(row) == len(header):
                allData.append(row)
        for k in allData:
            if k != header:
                well = k[2], k[3]
                wells.add(well)

        missingWells = allWells.difference(wells)

        if len(missingWells) == 0:
            shutil.move(fileName, incFolderName)

            for q in missingWells:
                newRow = ["-1", "-1", str(q[0]), str(q[1]), "N",
                          "-1", "0", "-1", "-1", "-1", "-1", "-1",
                          "-1", "-1", "-1", "-1", "-1", "-1", "-1",
                          "-1", "-1", "-1", "-1", "-1", "-1", "-1"]
                allData.append(newRow)

            name = nameFile2(txtFiles[i])

            newFile = open(name, "w+")
            newFile.write("\t".join(header) + "\n")

            for l in allData:
                newFile.write("\t".join(l) + "\n")

            newFile.close()

        print ("Added " + str(len(missingWells)) + " lines to " + txtFiles[i] +
               " to create " + name + ". Original file hase been moved to " +
               incFolderName + ".\n")
