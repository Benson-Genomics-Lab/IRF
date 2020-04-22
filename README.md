# IRF
**Version 3.07 Feb 2015**

## Table of Contents ##
- [Purpose](#purpose)   
- [Reference](#reference)
- [Authors](#authors)
- [License](#license)
- [Pre-compiled Versions](#pre-compiled-versions)
- [Instructions for Compiling](#instructions-for-compiling) 
- [Quick Start](#quick-start)
- [Using Command Line Version of Inverted Repeats Finder](#using-command-line-version-of-inverted-repeats-finder)  
- [IRF Definitions](#trf-definitions)  
- [FASTA Format](#fasta-format)
- [Table Explanation](#table-explanation)
- [Alignment Explanation](#alignment-explanation)
- [How does Inverted Repeats Finder work?](#how-does-inverted-repeats-finder-work)  
- [What's New](#whats-new)

## What's New ##

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
