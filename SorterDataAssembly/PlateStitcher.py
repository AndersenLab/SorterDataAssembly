#Import Tkinter to be able to graphically retrieve the directory
from Tkinter import Tk
from tkFileDialog import askdirectory
from os import listdir
import re
import collections

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

#Get the duplicate plate numbers
duplicates = [n for n, times in collections.Counter(numbers).items() if times > 1]

