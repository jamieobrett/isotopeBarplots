# isotopeBarplots
8-30-24

Short script to take a natural abundance-corrected peak areas file (direct output from IsoCorrectoR), a compound list file, and a metadata file and plot stacked barplots (absolute and percentage) with errorbars.

Metadata file format: tab-delimited .txt file
- first column = Name (exact match to the column name of the abundances file)
- second column = Timepoint #any format, using a consistent unit (don't mix h and d)
- third column = CellLine
- fourth column = Treatment #could be BSA, Stearate, etc.
- fifth column = RealSample #TRUE or FALSE -- exclude blanks, QCpools, etc.
The first row contains these headers names exactly as written

 Compound file format: .txt file with one compound name per line (names must match the IsoCorrectoR file), no header

 Corrected csv file format: exact ouptut from IsoCorrectoR

 Sample names syntax (you will need to carefully rename your IsoCorrectoR output columns and your metadata Name entries if needed):
 - No spaces or symbols other than "_" or "."
 - Start with a letter, not a number or symbol
