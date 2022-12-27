#makefile for program IRF

#run by typing 'make'
#cleanup by typing 'make clean'

#ideas from
#http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

#CC is the compiler
CC=gcc

#CFLAGS are options sent to the compiler
# -I means include files directory
#-O2 is an optimizer code
CFLAGS=-O2

#LFLAGS are options sent to the linker
# -lm is the math library to include
LFLAGS=-lm

#Makefile is stored in the src directory
#output directory is bin, a sibling directory of src
ODIR=../bin/

#IRF code is stored in the src directory
IRFDIR=src/

#dependent files
IRFDEPS=irf.h easylife.h memorycheck.h test.interval.sums.h centerlist.h krunsums2Anewt.h intervals.h centerlist.11.03.h centertag.h new.centerlists.h irfrun.h print.interval.sums.h profiler.h align.h

EASYLIFEDEPS=irf.h

IRFOBJS=irf.o easylife.o

all: $(ODIR)irf.exe

#target: dependencies
#$(x) is replaced by macro x
#-c means make an object file (.o file)
#-o is the name of the output file
#$@ is the name on the left side of the colon (:)
#$< is the file to compile, it stands for the first file in the dependencies list on the right side of the coloon
irf.o: irf.3.c $(IRFDEPS)
	$(CC) -c -o $@  $< $(CFLAGS)

easylife.o: easylife.c $(EASYLIFEDEPS)
	$(CC) -c -o $@  $< $(CFLAGS)

#final target irf.exe
#depends on the object files
#note there is no -c here because this is making the executable file
#-o is the name of the output file
#$@ is the name on the left side of the colon (:)
#$(ODIR)$@ decodes to ../bin/aln
$(ODIR)irf.exe: $(IRFOBJS)
	$(CC) -o $@ $(IRFOBJS) $(LFLAGS)
	

#clean up by removing all the object files (*.o = anything.o)
#clean will run only when typing make clean after producing the executable by typing make
#don't use until the executable is correct because it removes all the intermediate .o files
#which are required when debugging so they don't have to be recompiled if not changed
clean:
	rm -f *.o 
