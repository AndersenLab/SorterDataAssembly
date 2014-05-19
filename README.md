Sorter Data Assembly and Report Generation
==========================================

PlateStitcher.py
----------------
This program is run on a directory containing either the raw setup or raw score files for the high throughput assays. It is designed such that it can, and should, be run from a directory external to that in which the data files are stored.

####Naming Conventions
Naming conventions are critical to the correct execution of this code. The convention is outlined below:

+All plates should be named in the following manner: **plateNumber_drug_notes.txt**
	+Plate number should always be the letter "p" followed immediately by the two digit representation of the number (i.e. "08" or "15")
	+Be careful of spelling errors in the drug name, though this is not critical as the program will prompt you before stitching together plates with matching numbers but different drug names
	+Anything in the notes section of the name will be ignored by the program
	+The file ***must*** be saved as .txt

####Usage
1. cd into directory containing the PlateStitcher.py script
2. Execute $python PlateStitcher.py
3. Select the desired directory from the graphical interface
	+Note: The graphical interface does not automatically quit until the end of the script. You need to reselect your terminal window in order to answer any prompts that arise during the script's execution.
4. Program will run, printing any file modifications, movements, etc. User may be prompted for input if necessary (i.e. two or more plate numbers match, but drug names do not).