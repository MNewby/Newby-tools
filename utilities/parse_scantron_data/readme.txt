Usage:

Simple way:

  Just run this file in the same directory as the .dta files of interest.  It will search the directory for each .dta file and output a spreadsheet-readable (i.e., Excel) .csv file for each.  The .csv file will have each question separated (by a comma).

More complicated:

  Command line arguments with a "." in them will be interpreted as filenames.  ".dta" suffixes will be interpreted as inputs to the script, all other ".xxx" (where xxx is any three-character pattern, but not "dta") inputs will be interpreted as output filenames. The first *.dta will be paired with the first *.xxx, the second *.dta will be paired with the second *.xxx, etc.  If there is an uneven number of .dta and .xxx inputs, generic filenames will be used for the excess .dta inputs, and extra .xxx inputs will be discarded.

Example:

python parse_scantron_dta.py exam1.dta exam2.dta exam1.csv exam1.csv

Will turn "exam1.dta" into a "exam1.csv" file, and "exam2.dta" into an "exam2.csv" file.


VERY USEFUL:

  Adding "letters", "alph", or "-l" to the command line will convert numerical answers into alphabetical answers.  Oddly, the pattern appear to be "01234" = "EDCBA".

Example:

python parse_scantron_dta.py exam1.dta exam1.csv -l

This is probably the most useful way to use this script.

Other arguments:

"delimiter", "delim", or "-d", follow by a character:  use the character as the delimiter in the generated file.  Example:  "python parse_scantron_dta.py exam1.dta exam1.csv -d \t" will separate each item by tab (\t) spaces.
