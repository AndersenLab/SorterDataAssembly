#Import Tkinter to be able to graphically retrieve the directory
from Tkinter import Tk
from tkFileDialog import askdirectory
from os import listdir, mkdir
from os.path import exists
import re
import collections
import csv
import shutil

#Get the directory
directory = askdirectory()
files = listdir(directory)
txtFiles = []

#Get the txt files only
for f in files:
	if ".txt" in f:
		txtFiles.append(f)

#Get the plate numbers
numbers = []
for f in txtFiles:
	numbers.append(f[0:3])

#Function to get the drug name for file naming later on
def nameFile(f):
	number = f.split("_")[0]
	drug = f.split("_")[1].split(".")[0]
	return directory + "/" + number + "_" + drug + "_stitched.txt" 

#Get the duplicate plate numbers
duplicates = [[n, times] for n, times in collections.Counter(numbers).items() if times > 1]

if len(duplicates) != 0:
	folderName = directory + "/UnstitchedData"
	if not exists(folderName):
		mkdir(folderName)

for i in range(0,len(duplicates)):
	fileCount = 0
	allData = []
	tempData = []
	for j in range(0,len(txtFiles)):
		if duplicates[i][0] in txtFiles[j]:
			fileCount += 1
			fileName = directory + "/" + txtFiles[j]
			with open(fileName, "r") as f:
				reader = csv.reader(f, delimiter = "\t")
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
	
	newFile = open(name, "w+")
	newFile.write("\t".join(header) + "\n")

	for l in allData:
		newFile.write("\t".join(l) + "\n")

	newFile.close()








