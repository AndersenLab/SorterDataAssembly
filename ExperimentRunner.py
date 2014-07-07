import rpy2
from PlateStitcher import *
from os.path import relpath, expanduser
import subprocess

dirList = sys.argv[1:]

stitchAll(dirList)

print (colors.FINISH + "\nPlateSticher finished running, now starting " +
       "data assembly.\n" + colors.DEFAULT)

string = (colors.PROMPT + "\nDo you want to generate all of the " +
          "setup and score reports? (y/n, answering no will generate only " +
          "the final data frame .csv file): " + colors.DEFAULT)
answer = raw_input(string)

if answer == "y" or answer == "yes":
    generateReports = "TRUE"
else:
    generateReports = "FALSE"

# Uncomment if running old version of assembly code
# for i in range(0, len(dirList)):
#     dirList[i] = relpath(dirList[i], expanduser("~"))

directories = " ".join(dirList)

command = " ".join(["Rscript SimpleDataProcess.R", generateReports,
                    directories])

subprocess.call(command, shell=True)

print (colors.FINISH + "\nScript Complete\n" + colors.DEFAULT)
