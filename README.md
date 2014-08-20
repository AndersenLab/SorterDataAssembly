Sorter Data Assembly and Report Generation
==========================================

ExperimentRunner.py
-------------------
Use this program to run both PlateStitcher.py and SimpleDataProcess.R on experiment directories.

### Usage

1. Run from the command line `$ python ExperimentRunner.py [directory list]`
	+ Example `

PlateStitcher.py
----------------
This program is run on a directory containing either the raw setup or raw score files for the high throughput assays. It is designed such that it can, and should, be run from a directory external to that in which the data files are stored.

### Naming Conventions

Naming conventions are critical to the correct execution of this code. The convention is outlined below:

+ All plates should be named in the following manner: **plateNumber_condition-notes.txt**
	+ Plate number should always be the letter "p" followed immediately by the two digit representation of the number (i.e. "08" or "15")
	+ Be careful of spelling errors in the condition name, though this is not critical as the program will prompt you before stitching together plates with matching numbers and different condition names. Unfixed Spelling errors will cause issues later on in the processing and mapping code.
	+ Anything in the notes section of the name will be ignored by the program
	+ The file ***must*** be saved as .txt
	+ Plate number and condition name ***must*** be separated by an underscore
		+ Underscore are reserved characters, **DO NOT** use them to separate portions of the condition name, use dashes
		+ **Example:** 20140505_abamectin_25C.txt should be 20140505_abamectin-25C.txt if you would like to differentiate it from abamectin at 20 degrees C, etc.
		+ The notes section of the name is reserved for name appendices from the PlateStitcher.py script
		+ **File names should only contain one underscore (to separate the date and condition name) when originally saved**

### Usage

1. **Only run PlateStitcher through ExperimentRunner.py, it cannot be used on its own**


SimpleDataProcess.R
-------------------

#### Naming Conventions

Data on the date, round, description, assay, plate number, and condition of each experiment are grabbed from the directory and file names. The file names for the setup and score data files should follow the same conventions as those described above for the PlateStitcher.py program. The directory name conventions are described below:

+ The directory name should follow the format: **YearMonthDay_DescriptionRoundAssay**
	+ Year should be the four digit representation of the year
	+ Month should be the two digit representation of the month
	+ Day should be the two digit representation of the day
	+ An underscore should separate the date from the rest of the info, not a dash as was used previously
	+ Round must be a number, though the number digits associated with that number do not matter ("01" is the same as "1")
	+ Assay should be a letter for the assay round that is being completed ("a" for first, "b" for second, etc.)
	+ Example: **20140317_GWAS1a**

### Dependencies:

All of the following files should be in the same directory as the SetupAndScoreReportGenerator.R file:

+ MasterSetupReport2.Rmd - Markdown template for the setup reports
+ MasterScoreReport2.Rmd - Markdown template for the score reports
+ PresentationStyle.RData - Presentation style specs for plots in the reports
+ SimpleDataProcessFxns.R - File containing most of the functions to process the raw sorter data

All of the following files and directories should be in each experiment directory:

+ contamination.R - R file listing each plate and the associated contaminated wells
+ controls.R - R file listing each control plate and the associated condition plates it serves as the control for
+ strains.R - R file with a vector containing the names of the strains on each plate in order by row
+ setup - Directory containing all of the raw setup files
+ score - Directory containing all of the raw score files

### Human Input Files:
The following describes the setup for the three human input files. All three file names are case sensitive, including the file extension (must be ".R" not ".r").

#### contamination.R

This file contains a line for each plate following the above plate numbering convention, followed by the assignment of a vector with the well numbers that had contamination. An example of the head of one of these files is below:

```r
p01 <- c("D9", "D11")
p02 <- c("A9", "B5", "B11")
p03 <- c("")
p04 <- c("E7", "F7", "F9", "H9")
p05 <- c("G7") #same strains paralyzed
p06 <- c("F9", "H3", "H5", "H7")
p07 <- c("C7", "C9", "D7", "D9", "D11", "E3", "E7")
p08 <- c("C3", "E11", "G5")
p09 <- c("A9", "G1", "G3", "H1")
p10 <- c("")
```

#### controls.R

This file contains two lists, controlPlates and testPlates. In controlPlates, the control plates are listed in groups based on which plates they control for. In testPlates, the test plates are listed in the same order as the corresponding control plates.

For example, if plate 1 is the control for plate 2, the file would look like this:

```r
controlPlates <- list(1)
testPlates <- list(2)
```

If plate 1 is a control plate for a sequential list of plates, the colon character can be used along with the two endpoints to indicate this in the test plates list:

```r
controlPlates <- list(1)
testPlates <- list(2:10)
```

If multiple controls are present and they control for different plates, the indices in each list must line up. In this example, plate 1 controls for plates 2 through 10 and plate 11 controls for plates 12 through 20:

```r
controlPlates <- list(1, 11)
testPlates <- list(2:10, 12:20)
```

If the test plates or control plates are not sequential, they must be inserted into the list by vectors. In this example, plates 1 and 3 control for plates 4, 6, and 8 while plates 2 and 10 control for plates 5, 7, and 9:

```r
controlPlates <- list(c(1, 3), c(2, 10))
testPlates <- list(c(4, 6, 8), c(5, 7, 9))
```

#### strains.R

This file consists of a row-wise vector of all of the strains on the plate. The list is created right to left then down the plate, as if reading a book from A1 to H12. An example is given below:

```r
strains <- c("CX11285", NA, "ED3048", NA, "JU1200", NA, "CX11315", NA, "JU1440", NA, "CX11314", NA,
           "JU1242", NA, "KR314", NA, "LSJ1", NA, "CB4852", NA, "DL200", NA, "CB4856", NA,
           "ED3011", NA, "JU258", NA, "CX11307", NA, "JU1400", NA, "PB306", NA, "ED3052", NA,
           "WN2002", NA, "PS2025", NA, "CB4857", NA, "CB4858", NA, "AB4", NA, "ED3017", NA,
           "MY23", NA, "JU1246", NA, "JU1395", NA, "CB4854", NA, "JU1213", NA, "EG4724", NA,
           "EG4725", NA, "JU311", NA, "JU397", NA, "JU1581", NA, "JU1568", NA, "MY1", NA,
           "JU323", NA, "JU792", NA, "JU775", NA, "JU310", NA, "JU393", NA, "JU847", NA,
           "CX11271", NA, "CB4853", NA, "JU1586", NA, "ED3012", NA, "JU440", NA, "CX11276", NA)
```

### Usage

1. **Only run PlateStitcher through ExperimentRunner.py, it cannot be used on its own**











