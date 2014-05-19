Sorter Data Assembly and Report Generation
==========================================

PlateStitcher.py
----------------
This program is run on a directory containing either the raw setup or raw score files for the high throughput assays. It is designed such that it can, and should, be run from a directory external to that in which the data files are stored. *This program should be run on both the setup and score files for each assay before the data are run through the report generator R script.*

####Naming Conventions
Naming conventions are critical to the correct execution of this code. The convention is outlined below:

+ All plates should be named in the following manner: **plateNumber_drug_notes.txt**
	+ Plate number should always be the letter "p" followed immediately by the two digit representation of the number (i.e. "08" or "15")
	+ Be careful of spelling errors in the drug name, though this is not critical as the program will prompt you before stitching together plates with matching numbers but different drug names
	+ Anything in the notes section of the name will be ignored by the program
	+ The file ***must*** be saved as .
	+ Plate number and drug name ***must*** be separated by an underscore

####Usage
1. cd into directory containing the PlateStitcher.py script
2. Execute $python PlateStitcher.py
3. Select the desired directory from the graphical interface
	+ Note: The graphical interface does not automatically quit until the end of the script. You need to reselect your terminal window in order to answer any prompts that arise during the script's execution.
4. Program will run, printing any file modifications, movements, etc. User may be prompted for input if necessary (i.e. two or more plate numbers match, but drug names do not).


SetupAndScoreReportGenerator.R
------------------------------
This program takes, as input in the first code block, a list of directories that collectively make up one experiment (i.e. GWAS1a and GWAS1b). This code outputs setup and score reports for each individual plate as well as a master dataframe for the complete experiment, with all variables.

####Naming Conventions
Data on the date, round, description, assay, plate number, and drug of each experiment are grabbed from the directory and file names. The file names for the setup and score data files should follow the same conventions as those described above for the PlateStitcher.py program. The directory name conventions are described below:

+ The directory name should follow the format: **YearMonthDay_DescriptionRoundAssay**
	+ Year should be the four digit representation of the year
	+ Month should be the two digit representation of the month
	+ Day should be the two digit representation of the day
	+ An underscore should separate the date from the rest of the info, not a dash as was used previously
	+Round must be a number, though the number digits associated with that number do not matter ("01" is the same as "1")
	+ Assay should be a letter for the assay round that is being completed ("a" for 1st, "b" for secound, etc.)
	+ Example: **20140317_GWAS1a**