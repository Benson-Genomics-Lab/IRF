
# IRF
**Version 3.08 Dec 2022**

## Table of Contents ##
- [Purpose](#purpose)   
- [Reference](#reference)
- [Authors](#authors)
- [License](#license)
- [Pre-compiled Versions](#pre-compiled-versions)
- [Instructions for Compiling](#instructions-for-compiling) 
- [Quick Start](#quick-start)
- [Using Command Line Version of Inverted Repeats Finder](#using-command-line-version-of-inverted-repeats-finder)  
- [IRF Online](#irf-online)  
- [What's New](#whats-new)


## Purpose ## 
An inverted repeat (IR) consists of two arms (or stems) of similar DNA or RNA—with one inverted and complemented relative to the other—adjacent to a central, usually nonhomologous, spacer sequence.  Smaller IRs stems are typically found in mRNA sequences where they create stem-loop or hairpin structures as the complementary stems form base pairs.  Larger IRs are found in genomic sequence as in the human Y-chromosome, as detailed in the reference below.  IRs are also found at the boundaries of transposable elements.  

This program detects approximate inverted repeats.  It can be tuned, using the parameters, to target different stem lengths, different spacer lengths, and different degrees of approximation.  It allows GT pairing, which occurs in RNA stem-loops.  It also allows detection of mirror repeats (where the stems are inverted, but not complemented).  

In order to use the program, the user submits a sequence in FASTA format. The output consists of two files: a repeat table file and an alignment file. The repeat table, viewable in a web browser, contains information about each repeat, including its location, size, and nucleotide content. Clicking on the location indices in one of the table entries opens a second browser page that shows an alignment of the stems to each other. 

This material is based upon work supported by the National Science Foundation under Grant No. DBI-0413462.

## Reference ##
P. E. Warburton, J. Giordano, F. Cheung, Y. Gelfand and G. Benson.
**Inverted Repeat Structure of the Human Genome:
The X-Chromosome Contains a Preponderance of Large,
Highly Homologous Inverted Repeats That Contain Testes Genes**,
Genome Research, 14:1861-1869, 2004. [10.1101/gr.2542904](https://dx.doi.org/10.1101%2Fgr.2542904)

## Authors ##
Gary Benson, Yevgeniy Gelfand

## License ##
Inverted Repeats Finder Copyright (C) 1999-2022 Gary Benson

IRF is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

IRF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with IRF. If not, see https://www.gnu.org/licenses/.

## Pre-compiled Versions ##
To obtain current and/or earlier, pre-compiled versions of IRF:

- locate and click on the releases link on the right of the main github page for this repository
- check in each of the listed releases for executable files

## Instructions for Compiling ##
Change to the `src` directory and type make: 
```
cd src
make
```  
Ignore the warnings.

Change back to the outer directory and test the compiled version with the following command to process the yeast1.fa file included with this distribution: 
```
cd ..
bin/irf.exe yeast1.fa 2 3 5 80 10 20 100000 1000000 -d -ngs
```
This should produce 4 files as well as full output on the command line with 198 repeats:
```
yeast1.fa.2.3.5.80.10.20.100000.1000000.1.txt.html
yeast1.fa.2.3.5.80.10.20.100000.1000000.2.txt.html
yeast1.fa.2.3.5.80.10.20.100000.1000000.1.html
yeast1.fa.2.3.5.80.10.20.100000.1000000.2.html
```

Open the `yeast1.fa.2.3.5.80.10.20.100000.1000000.1.html` file in your browser and you should see a table with the repeats and links to the other files.

## Quick Start ##

The following is a recommended command line to run IRF. Parameters are explained further below. This assumes the executable has been renamed `irf.exe`.

```
irf.exe yourfile.fa 2 3 5 80 10 20 100000 1000000
```

## Using Command Line Version of Inverted Repeats Finder ##
To see the full instructions, run the executable file without any parameters:

```
bin/irf.exe 

Inverted Repeats Finder, Version 3.08
Copyright (C) Dr. Gary Benson 2002-2022. All rights reserved.


Please use: bin/irf.exe File Match Mismatch Delta PM PI Minscore Maxlength MaxLoop [options]

Where: (all weights, penalties, and scores are positive)
  File = sequences input file
  Match  = matching weight
  Mismatch  = mismatching penalty
  Delta = indel penalty
  PM = match probability (whole number)
  PI = indel probability (whole number)
  Minscore = minimum alignment score to report
  MaxLength = maximum stem length to report (10,000 minimum and no upper limit, but system will run out memory if this is too large)
  MaxLoop = filters results to have loop less than this value (will not give you more results unless you increase -t4,-t4,-t7 as well)
  [options] = one or more of the following :
               -m    masked sequence file
               -f    flanking sequence
               -d    data file
               -h    suppress HTML output

               -l    lowercase letters do not participate in a k-tuple match, but can be part of an alignment
               -gt   allow the GT match (gt matching weight must follow immediately after the switch)
               -mr   target is mirror repeats
               -r    set the identity value of the redundancy algorithm (value 60 to 100 must follow immediately after the switch)

               -la   lookahead test enabled. Results are slightly different as a repeat might be found at a different interval. Faster.
               -a3   perform a third alignment going inward. Produces longer or better alignments. Slower.
               -a4   same as a3 but alignment is of maximum narrowband width. Slightly better results than a3. Much slower.
               -i1   Do not stop once a repeat is found at a certain interval and try larger intervals at nearby centers. Better(?) results. Slower.
               -i2   Do not stop once a repeat is found at a certain interval and try all intervals at same and nearby centers. Better(?) results. Much slower.
               -r0   do not eliminate redundancy from the output
               -r2   modified redundancy algorithm, does not remove stuff which is redundant to redundant. Slower and not good for TA repeat regions, would not leave the largest, but a whole bunch.

               -t4   set the maximum loop separation for tuple of length4 (default 154, separation <=1,000 must follow)
               -t5   set the maximum loop separation for tuple of length5 (default 813, separation <=10,000 must follow)
               -t7   set the maximum loop separation for tuple of length7 (default 14800, limited by your system's memory, make sure you increase maxloop to the same value)
               -ngs  more compact .dat output on multisequence files, returns 0 on success. 

Note the sequence file should be in FASTA format:

>Name of sequence
   aggaaacctg ccatggcctc ctggtgagct gtcctcatcc actgctcgct gcctctccag
   atactctgac ccatggatcc cctgggtgca gccaagccac aatggccatg gcgccgctgt
   actcccaccc gccccaccct cctgatcctg ctatggacat ggcctttcca catccctgtg
```
## IRF Online ##

IRF is available for online submission of sequences at:
https:\\tandem.bu.edu\irf

More information on IRF is available at that website.

## What's New ##

### V3.08 ###
- fixed bug in irf.3.c in memory allocation of RCcodesSimilar


### V3.07 ###
- added -NGS flag (unix only)


### V3.06 ###
- compiles for linux now, uncomment UNIXCONSOLE 
- not sure if a big deal, but criticalquit_error function wasn’t actually exiting,
changed it to exit with exit code 1

### V3.05 ###
- Found a bug when calling search_for_centerseen_exact, wrong indices being sent in. Change does not effect results, but makes the probram much faster because the list no longer has duplicates

### V3.04 ###
- Now search distances are unlimited. You are only limited by the system memory.

### V3.03 ###
- Fixed a bug in the loading routine (on old bug that migrated from trf). When reading a FASTA file on windows which was transferred to windows from unix as binary, we must read in binary mode or we lose 1 or more characters at the end of the sequence.
- To do: Lookahead can be much better for GT if we read the next character and  checke if it is either a G or T combination too. Right now I don’t do it because it is a lot of work. An easy way is too use the match() macro, but that macro does not take into account case of letters (in fact it ignores it). Remember: tuples only referr to uppercase letters in masked files.
- Commented out all memory checks because scalloc already does that. No need to do the same thing twice.
- Added this change to alignment to make sure alignment jumps correctly when a run of GT matches is found: if(Gtmatchyesno&&(*pcurr==(*pdiag+GTAlpha))&&(match_yes_no==GTAlpha)) matchatmax_col=c;
- Made parameter string to be a string, so all parameters would are displayed. IRDB CHANGE!!!
- Added a third alignment going inward to get better inner endpoints. Program is slower but produces more accurate results. Must think more about the bandwidth for this alignment. Added a switch for it.
- Found a bug where RCCodes was sometimes overwriting a value just inserted in the TUPLEHASH (when code and rccodes were the same and distance exceeded on the first try). After fixing, I noticed the program runs a little slower but that is to be expected: before it was throwing out a lot of matches. It finds a few more repeats too.
- Added Gtpairs to output. Added center, avg center to output. IRDB CHANGE!!! Print a line of field descriptors in the .dat file?
- Added memory check library from EasyLife. Remove memory errors, as EasyLife already checks, even in release mode 
- Tried: “center waiting to align” modification to lookahead to make sure it will be tried, can’t do it because the lookeahead test is before the testcriteria function.
- Add a flag to cancel HTML output
- Not active: Added interval check for centerseen list, to do either larger intervals, or all other intervals. Both seem to be advantageous, but few repeats were observed to be worse (redundancy?) Added as a switch.
- changed the –t3 to –t4 IRDB CHANGE!!!
- Added another redundancy function Added flags to use it on no redundancy whatsoever, r2 and r0, r1 is the default one.
- Added tuple detection for GT. Some repeats also disappear. I am pretty sure this is because due to other matches, other centers appear to have more matches and testcriteria function does not pass. It is also possible, that because these other centers are due to GT matches and GT align scores are lower then regular matches. 
- See above problem: Changed gtmatch #matches to a double to reflect GT match Now partial matches are counted for GT matches based on score. I use GTDetectMatch = ((double)(GTAlpha - Beta))  /  ((double)(Alpha - Beta));
- Found 2 very serious bugs while doing above with my GT scheme. 
a.	My gtSimilarRrcodes were not correctly computed. 
b.	for (nsc=-1; simCode!=-1; simCode=RCcodesSimilar[g][Tuplerccode[g]][nsc]) was written as (nsc=-1; simCode!=-1; simCode=RCcodesSimilar[g][simCode][nsc]) Means instead of reading the same array, it was jumping around

### V3.02 ###
- Added a fix to NOT output the repeat to disk if it is in the center_seen_list. That keeps the center_seen_list shorter and the output file smaller.
- Not active: Added code to extend the alignment through the original point where it was found. This has a weird effect that two repeats would be joined with a lower score than an individual (outer) repeat. This is because the score is higher than the inner repeat’s score and outer repeat never gets a chance to get tested. I think the test is counterproductive. ???
- Not active: Added code to skip alignment if determined that the next character matches as well (to remove redundant alignments) It significantly speeds up the algorithm, but mysteriously some repeats disappear. (see “Look Ahead test differences, possible explanation” document)
- Added option for GT match with gtweight
- Added option for Mirror Repeat Target
- Added options for tuple max loop user settings

### v3.01 ###
- A new compile time flag is added: IRF_MEGA. When defined, maxwraplength is increased to 500,000 (from 100,000) and MINBANDRADIUS_OUTER is increased to 100.
- Fixed a bug that caused program to crash when running on a multiple FASTA file. The problem was that the new files have loopsize postfix and there was a place which did not check if the file was opened successfully.
- Added a fix that increases user set loopsize by 2*MAXINTERVALSIZE. This is done because apparently the program does not find repeats with separation close but still less than the loopsize specified. Because of this, some repeats are found with separation greater than loopsize. These repeats are removed at the end by the RemoveBySizeAndLoop function.

### V3.00 ###
- This version is significantly different from the previous as the match arrays are now dynamic (they were arrays before). A few bugs were found mostly concerning matches on the sides of the intervals.
- Added a variable which lets the user set the loopsize.

### V2.0N ###
- A change to ignore lowercase letters made.
- Some statistical work done on the RNKP values.
- Versioning not yet kept so a lot of other stuff was changed undocumented.


