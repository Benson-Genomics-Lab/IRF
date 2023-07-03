/**************************************************************
*   trfrun.h :  This file contains the code that calls the IRF
*               algorithm.
*               It declares all the variables that control how
*               the program runs and where the output is saved.
*               If the input file contains more than one
*               sequence then the input in broken into single
*               sequences by the control routine and fed to the
*               algorithm one sequence at a time.
*               The output is assembled by the control routine
*               as described in the file readme.txt
*
*												  July 16, 2002
*
***************************************************************/


#ifndef IRFRUN_H
#define IRFRUN_H

#include <stdio.h>
#include <stdlib.h>


#ifndef _MAX_PATH
#define _MAX_PATH 260
#endif

#ifndef MAX_TEXT_FILE_LINE
#define MAX_TEXT_FILE_LINE 500
#endif

#ifndef MAX_RECORD_LABEL
#define MAX_RECORD_LABEL 60
#endif

/* Global strings to store non-tabulated information in html file */
char hsequence[_MAX_PATH+1];
char hparameters[_MAX_PATH+1];
char hlength[_MAX_PATH+1];

/* max # of items in tables for extended output format */
#define EO_MAX_TBL 120

struct index_list
{
    char    ref[MAX_RECORD_LABEL+1];    /* records label for linking */
    int     count;      /* indicates order in original file */
    int     leftfirst;      /* first index */
    int     leftlast;       /* last index */
    int     rightfirst;      /* first index */
    int     rightlast;       /* last index */
	int		leftlength;
	int		rightlength;
	int		separation;
    int     score;
    int     center, avecenter;
    float   matches;
    float   indels;
    float   atcount;
    float   gccount;
    float   atpairs;
    float   gcpairs;
    float   gtpairs;
    /* int     redundant; */
    /* G. Benson 4/10/23 replaced above to fix compiler warning at line 844, int to pointer cast */
    long long int     redundant;
    int     markedfordeletion;
    struct index_list* next;
};
typedef struct index_list IL;


int		LoadSequenceFromFile(FASTASEQUENCE *pseq,FILE* fp);
void    OutputHTML(IL * headptr, char * tablefile, char *alignmentfile);

IL*     GetList(char * datafile);
void    FreeList(IL * headptr);
void    MakeFileName(char* newname, char* oldname, int tag);
void    OutputHeading(FILE* fp, char* tablefile, char* alignmentfile);
IL*     RemoveBySizeAndLoop(IL * headptr, int maxsize, int maxloop);
IL*     SortByIndex(IL * headptr);

IL*     SortByCount(IL * headptr);
void    CleanAlignments(IL * headptr, char * alignmentfile);
void    BreakAlignments(IL * headptr, char * alignmentfile);
void    MakeDataFile(IL * headptr,char * datafile,int data);
void    MakeMaskedFile(IL* headptr,int masked,char*  Sequence,char* maskfile);

IL*     RemoveRedundancy(IL * headptr);
IL*     RemoveRedundancy2(IL * headptr);
int     IntervalOverlap(IL* iptr, IL* jptr);
int     IsRedundant(IL* iptr, IL* jptr);
/***********************************
*   Chief Procedure Definition
***********************************/

void IRFClean( char * datafile, char * tablefile, char * alignmentfile, int data, int maxsize, int maxloop, int masked, char* maskfile )
{
	int i;
    IL *headptr=NULL,*currptr;
 
    headptr = GetList(datafile);

    headptr = RemoveBySizeAndLoop(headptr, maxsize, maxloop);

    headptr = SortByIndex(headptr);

    if (paramset.redalg==1)
        headptr = RemoveRedundancy(headptr);
    else if (paramset.redalg==2)
        headptr = RemoveRedundancy2(headptr);

    headptr = SortByCount(headptr);

if (!paramset.HTMLoff) {
/*	fprintf(stderr,"press any key to clean alignments....\n"); fflush(stderr);
	getch();
	getch(); */
    CleanAlignments(headptr, alignmentfile);
	//fprintf(stderr,"ok..\n"); fflush(stder);

/*	fprintf(stderr,"press any key to break alignments...\n"); fflush(stderr);
	getch();
	getch(); */
    BreakAlignments(headptr, alignmentfile);
	//fprintf(stderr,"ok..\n"); fflush(stderr);

/*	fprintf(stderr,"press any key to output html and finish other stuff....\n"); fflush(stderr);
	getch();
	getch(); */
    OutputHTML(headptr, tablefile, alignmentfile);
	//fprintf(stderr,"ok..\n"); fflush(stderr);
}    

    MakeDataFile(headptr,datafile,data);

    MakeMaskedFile(headptr, masked, Sequence, maskfile);

    /* update the global result */
    for(i=0,currptr=headptr;currptr!=NULL;i++,currptr=currptr->next);
    paramset.outputcount= i;
	

    FreeList(headptr);

}

/*******************************************************
*   loads a sequence from an open file. If the routine
*   comes to another sequence header while reading a
*   sequence the sequence returns 1 to indicate that
*   more sequences are present. Otherwise the routine
*   returns 0 to indicate EOF or -1 to indicate error.
*   The member sequence must be NULL before the routine
*   is called. The calling function must free the allo-
*   cated memory after use.
********************************************************/

int LoadSequenceFromFile(FASTASEQUENCE *pseq,FILE* fp)
{
    int letter,i,pos1,length,next;
    char *ptext,*ptextcomp,ptextcopy;

    /* read the FASTA '>' symbol */
    letter = getc(fp);
    if(letter!='>') return -1; /* invalid format */

    /* read name and description text */
    for(i=0,ptext=pseq->name;i<(MAXSEQNAMELEN-1);i++,ptext++)
    {
        *ptext = getc(fp);
        if(*ptext==10||*ptext==13)
        {
            break;
        }
    }
    *ptext='\0';
    /* if line was not read completely flush the rest */
    if(i==(MAXSEQNAMELEN-1))
    {
        letter = 0;
        while(letter!=13&&letter!=10&&letter!=EOF) letter = getc(fp);
    }

    /* get the length of the sequence */
    pos1 = ftell(fp);
    length=0;
    for(;;)
    {
        letter = getc(fp);
        if(letter==EOF||letter=='>') break;
        length++;
    }

    /* allocate memory including white space for both sequence and sequencecomplement */
    memory_stats_print("\n LoadSequenceFromFile(sequence): requesting", (length+1)*sizeof(char) );
    pseq->sequence = (char*) scalloc(length+1,sizeof(char));
    if(pseq->sequence==NULL) return -2;
	pseq->sequencecomplement=(char *) scalloc(length+1,sizeof(char));
    memory_stats_print("\n LoadSequenceFromFile(sequencecomplement): requesting", (length+1)*sizeof(char) );
	if(pseq->sequencecomplement==NULL) return -2;



    /* read sequence into buffer */
    fseek(fp,pos1,SEEK_SET);
    pseq->length=0;
    ptext = pseq->sequence;
    ptextcomp = pseq->sequencecomplement;
    for(i=0;i<26;i++) pseq->composition[i]=0;
    for(;;)
    {
        /* get a character from file */
        *ptext = getc(fp);

        /* break if end of file */
        if(*ptext==EOF) { next = 0; break; }

        /* break if another sequence found */
        if(*ptext=='>') { next = 1; ungetc('>',fp); break; }

        /* if character is in range of alpha characters */
        if(*ptext>='A'&&*ptext<='z') /* in alpha range */
        {

			//make a copy of *ptext;
			ptextcopy=*ptext;

            if(ptextcopy<='Z') /* in upper case range */
            {
				*ptextcomp=Complementascii[*ptext];


				/* The following line does the following:
				   Because we don't want weird characters to match themselves, 
				   we change them to lower case in the reversed string.
				   Once we are printing the results, we will have to change them back. */
				if (*ptextcomp!='A' && *ptextcomp!='C' && *ptextcomp!='G' && *ptextcomp!='T') (*ptextcomp)-=('A'-'a'); 

                pseq->composition[toupper(*ptext)-'A']++;
                ptext++;
				ptextcomp++;
                pseq->length++;
            }
            else if(ptextcopy>='a') /* in lower case range */
            {

				/* to accomodate the new scheme of k-tuple exclusion of lowercases and alignment matching of them */
				if (paramset.lowercase==1) {
					/* make upper case only if it is not a good chrachter which was specifically labeled lowercase */
					if (*ptext!='a' && *ptext!='c' && *ptext!='g' && *ptext!='t') (*ptext)+=('A'-'a');
				} else {
					/* make upper case always */
					(*ptext)+=('A'-'a');
				}
                pseq->composition[toupper(*ptext)-'A']++; 

                /* we just make a regular copy */
				*ptextcomp=Complementascii[ptextcopy];

				/* The following line does the following:
				   Because we don't want weird characters to match themselves, 
				   we change them to lower case in the reversed string.
				   Once we are printing the results, we will have to change them back. 
				
				   NOTE: we still want masked out lowercase characters to be lowercase
				
				*/

				if (*ptextcomp=='a' || *ptextcomp=='c' || *ptextcomp=='g' || *ptextcomp=='t') (*ptextcomp)+=('A'-'a'); 


                ptext++;
				ptextcomp++;
                pseq->length++;
            }

        }
    }

    /* terminate sequence text as a string */
    *ptext='\0';

    /* compute notacgt member */
    pseq->nucleotides =
		pseq->composition['A'-'A'] + pseq->composition['C'-'A'] +
		pseq->composition['G'-'A'] + pseq->composition['T'-'A'];
	
    return next;
}

/***********************  OutputHTML *********************************************************/
void    OutputHTML(IL * headptr, char * tablefile, char *alignmentfile)
{
    FILE* fp;
    IL* currptr;
    int i,j;
    int alignments;
    int nfiles;
    char outfile[_MAX_PATH+1];
    char linkfile[_MAX_PATH+1];
    char namebuffer[_MAX_PATH+1];


    /* find out how many elements there are and how
       many files will be needed */
    for(alignments=0,currptr=headptr;
        currptr!=NULL;
        currptr=currptr->next,alignments++);
    nfiles = alignments/EO_MAX_TBL;
    if((alignments%EO_MAX_TBL)>0)nfiles++;
    if(nfiles==0) nfiles=1; /* make sure at least one file is generated */

    /* loop creating files */
    currptr=headptr;
    for(i=1;i<=nfiles;i++)
    {
        /* create name of ith file */
        MakeFileName( outfile, tablefile, i);

        /* create name of file where the alignments are found */
        MakeFileName( linkfile, alignmentfile, i);
			
        /* open the file for writing */
        fp = fopen(outfile,"w");

        /* output heading */
        OutputHeading(fp, outfile,alignmentfile);

        /* print links to other tables */
        fprintf(fp,"\n<P><PRE>Tables:   ");
        for(j=1;j<=nfiles;j++)
        {
            /* output a link for all but the current page */ 
            if(j!=i)
            {
                MakeFileName(namebuffer,tablefile,j);
                fprintf(fp,"<A HREF=\"%s\" target=\"_self\">%d</A>   ", namebuffer, j);
            }
            else
            {
                fprintf(fp,"%d   ", j);
            }
            if(j%16==0&&j<nfiles) fprintf(fp,"\n          ");
        }

        /* printf "Table n of N" line */
        fprintf(fp,"\n\nThis is table  %d  of  %d  ( %d repeats found )\n</PRE>",i,nfiles,alignments);

        /*print help lines */
        fprintf(fp,"<PRE>\nClick on indices to view alignment\n");

        #if defined(WINDOWSGUI)
            fprintf(fp,"</PRE><P>See <FONT COLOR=\"#0000FF\">Table Explanation</FONT> in Inverted Repeats Finder Help</P><BR>\n");
        #elif defined(WINDOWSCONSOLE)
            fprintf(fp,"</PRE><A HREF=\"http://tandem.bu.edu/irf/irf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");
        #elif defined(UNIXGUI)
            #error: Unix GUI code not implemented
        #elif defined(UNIXCONSOLE)
            fprintf(fp,"</PRE><A HREF=\"http://tandem.bu.edu/irf/irf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n");
        #endif



        /* print beginning of table */
        fprintf(fp,"<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");

        /* output rows of data */
        for(j=0;j<EO_MAX_TBL && currptr!=NULL;currptr=currptr->next,j++)
        {
            
            if(j%22==0) fprintf(fp,"<TR><TD WIDTH=140><CENTER>Left Indices</CENTER></TD><TD WIDTH=80><CENTER>Length</CENTER></TD><TD WIDTH=70><CENTER>Right Indices</CENTER></TD><TD WIDTH=70><CENTER>Length</CENTER></TD><TD WIDTH=70><CENTER>Loop</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=70><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>A/T</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>G/C</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>A-T pairs</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>G-C pairs</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>G-T pairs</CENTER></TD><TD WIDTH=40><CENTER>Center times 2</CENTER></TD><TD WIDTH=40><CENTER>Average Center times 2</CENTER></TD></TR>\n");
            fprintf(fp,"<TR><TD><CENTER><A HREF=\"%s#%s\">%d--%d</A></CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d--%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%.2f</CENTER></TD><TD><CENTER>%d</CENTER></TD><TD><CENTER>%d</CENTER></TD></TR>\n",
                    linkfile,currptr->ref,currptr->leftfirst,currptr->leftlast,currptr->leftlength,currptr->rightfirst,currptr->rightlast,currptr->rightlength,currptr->separation,currptr->matches,currptr->indels,currptr->score,currptr->atcount,currptr->gccount,currptr->atpairs,currptr->gcpairs,currptr->gtpairs, currptr->center, currptr->avecenter);
			//printf("%s\n",currptr->ref);
        }

        /* if empty list print at least one heading */
        if(headptr==NULL) fprintf(fp,"<TR><TD WIDTH=140><CENTER>Left Indices</CENTER></TD><TD WIDTH=80><CENTER>Length</CENTER></TD><TD WIDTH=70><CENTER>Right Indices</CENTER></TD><TD WIDTH=70><CENTER>Length</CENTER></TD><TD WIDTH=70><CENTER>Loop</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=70><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>A/T</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>G/C</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>A-T pairs</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>G-C pairs</CENTER></TD><TD WIDTH=40><CENTER>Percent<BR>G-T pairs</CENTER></TD><TD WIDTH=40><CENTER>Center times 2</CENTER></TD><TD WIDTH=40><CENTER>Average Center times 2</CENTER></TD></TR>\n");

        /* close table */
        fprintf(fp,"\n</TABLE>\n");

        /* if no repeats print message */
        if(headptr==NULL) fprintf(fp, "\nNo Repeats Found!<BR>");/*for empty list*/

        /* print links to other tables (again)*/
        fprintf(fp,"\n<P><PRE>Tables:   ");
        for(j=1;j<=nfiles;j++)
        {
            /* output a link for all but the current page */ 
            if(j!=i)
            {
                MakeFileName(namebuffer,tablefile,j);
                fprintf(fp,"<A HREF=\"%s\" target=\"_self\">%d</A>   ", namebuffer, j);
            }
            else
            {
                fprintf(fp,"%d   ", j);
            }
            if(j%16==0&&j<nfiles) fprintf(fp,"\n          ");
        }
        fprintf(fp,"\n</PRE>");

        if(i==nfiles) fprintf(fp,"<P>The End!\n");

        fprintf(fp,"\n</BODY></HTML>\n");
        fclose(fp);
    }

    return;
}
/***********************  MakeFileName *********************************************************/
void    MakeFileName(char* newname, char* oldname, int tag)
{
    char newext[20];
    char oldext[10];
    int numindex;
    int count;

    /*copy oldname to newname buffer*/
    strcpy(newname,oldname);

    /*find position of last number in name*/
    for(count=0, numindex=0; newname[count]!='\0';count++)
    {
        if(isdigit((int) newname[count])) numindex = count;
    }


    /*get old extension */
    strcpy(oldext,&newname[numindex+1]);

    /* create new extension based on tag*/
    sprintf(newext,".%d%s",tag,oldext);

    /* copy newext in place over old newname */
    strcpy(&newname[numindex+1], newext);

    return;
}

/***********************  OutputHeading *********************************************************/
void    OutputHeading(FILE* fp, char * tablefile, char *alignmentfile)
{
    /* output fixed (old) heading */
    fprintf(fp,"<HTML><HEAD><TITLE>%s</TITLE><BASE TARGET=\"%s\"></HEAD><BODY bgcolor=\"#FBF8BC\"><BR><PRE>Inverted Repeats Finder Program written by:</PRE><PRE><CENTER>Gary Benson<BR>Bioinformatics Program<BR>Boston University<BR>Version %s<BR></CENTER>",tablefile, alignmentfile, versionstring);
    fprintf(fp,"\nPlease cite:\nP. E. Warburton, J. Giordano, F. Cheung, Y. Gelfand and G. Benson. \n\"Inverted Repeat Structure of the Human Genome: The X-Chromosome \nContains a Preponderance of Large, Highly Homologous Inverted \nRepeats That Contain Testes Genes\", \nGenome Research, 14:1861-1869, 2004. 10.1101/gr.2542904.\n");
    
    fprintf(fp,"\n%s",hsequence);
    fprintf(fp,"%s",hparameters);
    fprintf(fp,"%s</PRE>\n",hlength);

    return;
}


/***********************  GetList *********************************************************/
IL*     GetList(char * datafile)
{
    FILE    *fp;
    IL      *headptr, *newptr, *lastptr;
    int     counter, i;

    headptr=newptr=lastptr=NULL;

    /* open file */
    fp = fopen(datafile, "r");
    if ( NULL==fp ) return NULL;

    /* get hsequence line on ninth line of data file*/
    for (counter =0; counter<9; counter++) fgets(hsequence, _MAX_PATH, fp);

    /* get hparameters line on twelvth line of data file*/
    for (counter =0; counter<3; counter++) fgets(hparameters, _MAX_PATH, fp);

    /* get hlength from another global variable (bad practice)*/
    sprintf(hlength, "Length:  %d", Length);


    /* loop to fill out list from buffer */
    counter = 1; /* keeps track of order they are found */
    while(1)
    {
        /* create new index list element */
        newptr = (IL*)scalloc( 1, sizeof(IL));
        memory_stats_print("\n GetList(newptr): requesting", sizeof(IL) );
        /*if(newptr==NULL)
        {
            {debugerror("\nUnable to allocate memory to newptr. Aborting!"); exit(0); }
            FreeList(headptr);
            return NULL;
        }*/
        newptr->count = counter++;
		newptr->redundant=0;
        newptr->markedfordeletion=0;

        /* get data from file */
        i = fscanf(fp, "%s %d %d %d %d %d %d %d %f %f %d %f %f %f %f %f %d %d",
                   newptr->ref, 
				   &newptr->leftfirst, &newptr->leftlast, &newptr->leftlength,
				   &newptr->rightfirst, &newptr->rightlast, &newptr->rightlength, &newptr->separation,
				   &newptr->matches,&newptr->indels,&newptr->score,
                   &newptr->atcount, &newptr->gccount, &newptr->atpairs, &newptr->gcpairs, &newptr->gtpairs, &newptr->center, &newptr->avecenter);

        if(i==EOF)
        {
            sfree(newptr);
            break;
        }


        if (headptr == NULL) /* first element */
        {
            headptr = lastptr = newptr;
            lastptr->next = NULL;
        }
        else /* add new element to end of list */
        {
        lastptr->next = newptr;
        lastptr = newptr;
        lastptr->next = NULL;
        }
    }

    fclose(fp);
    return  headptr;
}

/***********************  FreeList *********************************************************/
void   FreeList(IL * headptr)
{
    IL* nextptr;
    IL* holdptr;

    nextptr = headptr;
    while(nextptr!=NULL)
    {
        holdptr = nextptr;
        nextptr = nextptr->next;
        sfree(holdptr);
    }
    
    return;
}

/************ RemoveBySizeAndLoop() ********************************************************/

IL*     RemoveBySizeAndLoop(IL * headptr, int maxsize, int maxloop)
{
    IL* currptr;
    IL* prevptr;
	int averagelength;

    /* loop thru list removing all elements with period > maxsize */
    for(currptr=headptr; currptr!=NULL;)
    {
		averagelength=(int)ceil((currptr->leftlength+currptr->rightlength)/2.0);
        if(averagelength>maxsize || currptr->separation>maxloop) /* remove */
        {
            if(currptr==headptr)
            {
                headptr=headptr->next;
                sfree(currptr);
                currptr=headptr; 
            }
            else
            {
                prevptr->next = currptr->next;
                sfree(currptr);
                currptr = prevptr->next;
            }

        }
        else
        {
            prevptr= currptr;
            currptr= currptr->next;
        }
    }

    return headptr;

}


IL*     SortByIndex(IL * headptr)
{
    IL* currptr;
    IL* holdptr;
    IL* prevptr;
    int dif; /* flags when changes occur in one pass */

    if(headptr==NULL) return headptr; /* return if no elements */
    if(headptr->next==NULL) return headptr; /* return if one element only */

    dif=1;
    currptr=headptr;

    while(dif)
    {
        dif = 0;
        /* repeat inner loop until end is reached */
        while(currptr->next!=NULL)
        {
            if(currptr->leftfirst > currptr->next->leftfirst) /* swap */
            {
                if(currptr==headptr)
                {
                    holdptr=currptr->next->next;
                    headptr=currptr->next;
                    headptr->next=currptr;
                    currptr->next=holdptr;
                    prevptr=headptr;

                }
                else
                {
                    prevptr->next=currptr->next;
                    holdptr=currptr->next->next;
                    prevptr->next->next=currptr;
                    currptr->next=holdptr;
                    prevptr=prevptr->next;
                }

                dif=1; /* mark as changed */
            }
            else
            {
                prevptr=currptr;
                currptr=currptr->next;
            }

        }
        currptr=headptr; /* restart from begining */
    }
    
    return headptr;
}


IL*     SortByCount(IL * headptr)
{
    IL* currptr;
    IL* holdptr;
    IL* prevptr;
    int dif; /* flags when changes occur in one pass */

    if(headptr==NULL) return headptr; /* return if no elements */
    if(headptr->next==NULL) return headptr; /* return if one element only */

    dif=1;
    currptr=headptr;

    while(dif)
    {
        dif = 0;
        /* repeat inner loop until end is reached */
        while(currptr->next!=NULL)
        {
            if(currptr->count > currptr->next->count) /* swap */
            {
                if(currptr==headptr)
                {
                    holdptr=currptr->next->next;
                    headptr=currptr->next;
                    headptr->next=currptr;
                    currptr->next=holdptr;
                    prevptr=headptr;

                }
                else
                {
                    prevptr->next=currptr->next;
                    holdptr=currptr->next->next;
                    prevptr->next->next=currptr;
                    currptr->next=holdptr;
                    prevptr=prevptr->next;
                }

                dif=1; /* mark as changed */
            }
            else
            {
                prevptr=currptr;
                currptr=currptr->next;
            }

        }
        currptr=headptr; /* restart from begining */
    }
    
    return headptr;
}

IL*     RemoveRedundancy(IL * headptr)
{
    IL* iptr;       /* first pointer of a pair being examined */
    IL* jptr;       /* second pointer of a pair being examined */
    IL* previptr;   /* points to element before iptr */
    IL* prevjptr;   /* points to element before jptr */
    int overlap;    /* overlap of two intervals */
    int iinterval;  /* first pointer index interval */
    int jinterval;  /* second pointer index interval */

    if(headptr==NULL) return headptr; /* return if no elements */
    if(headptr->next==NULL) return headptr; /* return if one element only */

    iptr = headptr;       /* initialize to start at head of list */

    while(iptr!=NULL) /* loop from start to end of list*/
    {
		iinterval = iptr->leftlast - iptr->leftfirst + 1 + iptr->rightlast - iptr->rightfirst + 1;

        jptr=iptr->next;
        prevjptr=iptr;
        while(jptr!=NULL) /* loop until not enough overlap */
        {
            jinterval = jptr->leftlast - jptr->leftfirst + 1 + jptr->rightlast - jptr->rightfirst + 1;
            overlap = IntervalOverlap(iptr, jptr);

            /* if one is  outside of bounds, so will be all the next ones, so don't look at them */
            if(jptr->leftfirst>iptr->leftlast) { 
				break;
			} 

            /* overlap in relation to iinterval */
            if( !(overlap/(double)iinterval<(paramset.redident)) && IsRedundant(iptr,jptr))
            {
				iptr->redundant=1;
            } else

            /* overlap in relation to jinterval */
            if( !(overlap/(double)jinterval<(paramset.redident)) && IsRedundant(jptr,iptr))
            {
				jptr->redundant=1;

            }

            /* update and continue to next iteration of inner loop */
            prevjptr = jptr;
            jptr = jptr->next;

        }
        previptr = iptr;
        iptr = iptr->next;
    }

	/* now go though the list and remove all the redundant ones */
    for(iptr = headptr;iptr!=NULL; ) {

		if (iptr->redundant==1) { 
                if(iptr==headptr)
                {
                    headptr=headptr->next;
                    sfree(iptr);
                    iptr=headptr;
                    continue;
                }
                else
                {
                    previptr->next = iptr->next;
                    sfree(iptr);
					iptr = previptr->next;
                    continue;
                }
		}
	
		previptr = iptr;
		iptr=iptr->next;

	}
 


    return headptr;
}

/* Redundant to redundant, should not be considered redundant. This is what this function adresses. 
   Unfortunately, it is slow. (RemoveRedundancy2=RemoveRedundancy*LengthOfTheLargestPyramid) Gelfand, July 14, 2004 */
IL*     RemoveRedundancy2(IL * headptr)
{
    IL* iptr;       /* first pointer of a pair being examined */
    IL* jptr;       /* second pointer of a pair being examined */
    IL* previptr;   /* points to element before iptr */
    IL* prevjptr;   /* points to element before jptr */
    int overlap;    /* overlap of two intervals */
    int iinterval;  /* first pointer index interval */
    int jinterval;  /* second pointer index interval */
    int marked;

    if(headptr==NULL) return headptr; /* return if no elements */
    if(headptr->next==NULL) return headptr; /* return if one element only */

START:
    iptr = headptr;       /* initialize to start at head of list */

    while(iptr!=NULL) /* loop from start to end of list*/
    {
		iinterval = iptr->leftlast - iptr->leftfirst + 1 + iptr->rightlast - iptr->rightfirst + 1;

        jptr=iptr->next;
        prevjptr=iptr;
        while(jptr!=NULL) /* loop until not enough overlap */
        {
            jinterval = jptr->leftlast - jptr->leftfirst + 1 + jptr->rightlast - jptr->rightfirst + 1;
            overlap = IntervalOverlap(iptr, jptr);

            /* if one is  outside of bounds, so will be all the next ones, so don't look at them */
            if(jptr->leftfirst>iptr->leftlast) { 
				break;
			} 

            /* overlap in relation to iinterval */
            if( !(overlap/(double)iinterval<(paramset.redident)) && IsRedundant(iptr,jptr))
            {
                        iptr->redundant = (int)jptr;
                        //printf("%d-%d is red to %d-%d\n",iptr->leftfirst,iptr->leftlast,jptr->leftfirst,jptr->leftlast);
            } else

            /* overlap in relation to jinterval */
            if( !(overlap/(double)jinterval<(paramset.redident)) && IsRedundant(jptr,iptr))
            {

                        jptr->redundant = (int)iptr;
                        //printf("%d-%d is red to %d-%d\n",jptr->leftfirst,jptr->leftlast,iptr->leftfirst,iptr->leftlast);
            }

            /* update and continue to next iteration of inner loop */
            prevjptr = jptr;
            jptr = jptr->next;

        }
        previptr = iptr;
        iptr = iptr->next;
    }

    /* a repeat cannot be redundant to something that is already redundant */
    marked = 0;
    for(iptr = headptr;iptr!=NULL;iptr=iptr->next ) {

		if (iptr->redundant) { 
            jptr = (IL*)iptr->redundant;

            /* only mark for deletion if the repeat we are redundant to, is not, in turn, redundant to anything */
            if (jptr->redundant==0) {
                iptr->markedfordeletion = 1;
                marked = 1;
            } 
		}

	}

	/* now go though the list and remove all the redundant ones */
    for(iptr = headptr;iptr!=NULL; ) {

		if (iptr->markedfordeletion==1) { 

                //printf("\nremoving: %d-%d\n", iptr->leftfirst,iptr->leftlast);

                if(iptr==headptr)
                {
                    headptr=headptr->next;
                    sfree(iptr);
                    iptr=headptr;
                    continue;
                }
                else
                {
                    previptr->next = iptr->next;
                    sfree(iptr);
					iptr = previptr->next;
                    continue;
                }
		}
	
		previptr = iptr;
        iptr->redundant = 0;
		iptr=iptr->next;

	}

    /* now we need to redo the whole process, until we stop getting redundant repeats */
    if (marked)
        goto START;

    return headptr;
}

/* this fuction retuns overlap if two entries overlap */
int     IntervalOverlap(IL* iptr, IL* jptr)
{
    int beg,end, overlap;
    int beg2,end2, overlap2;
	int overlapTotal;

    beg = ( iptr->leftfirst > jptr->leftfirst ) ? iptr->leftfirst : jptr->leftfirst;
    end = ( iptr->leftlast  < jptr->leftlast  ) ? iptr->leftlast  : jptr->leftlast;
    overlap = end - beg + 1;    

    beg2 = ( iptr->rightfirst > jptr->rightfirst ) ? iptr->rightfirst : jptr->rightfirst;
    end2 = ( iptr->rightlast  < jptr->rightlast  ) ? iptr->rightlast  : jptr->rightlast;
    overlap2 = end2 - beg2 + 1;    

	
	overlap=( overlap > 0 ) ? overlap : 0;
	overlap2=( overlap2 > 0 ) ? overlap2 : 0;
	overlapTotal=overlap+overlap2;

    return overlapTotal;
}


/* this fuction retuns 1 if iptr is redundant in reference to jptr */
int     IsRedundant(IL* iptr,IL* jptr)
{
	if ( (iptr->leftlength+iptr->rightlength) <= (jptr->leftlength+jptr->rightlength) &&
        (iptr->score<= 1.1*jptr->score)) return 1; 

    return 0;
}

void   CleanAlignments(IL * headptr, char * alignmentfile)
{
    char string1[MAX_TEXT_FILE_LINE];
    char string2[MAX_TEXT_FILE_LINE];
    char string3[MAX_TEXT_FILE_LINE];
    char tempfile[_MAX_PATH+4];
    FILE *al_fp, *tmp_fp;
    IL*  currptr;
    int moving; /* in loop moving=1 to copy strings to tempfile */

    if (headptr==NULL) return; /* return if empty list */

    /* make name of temporary file */
    strcpy(tempfile,alignmentfile);
    strcat(tempfile,".tmp");

    /* open files */
    al_fp = fopen(alignmentfile, "r");
    tmp_fp= fopen(tempfile, "w");

    moving = 1; /* starts by moving lines to temporary file */
    currptr=headptr;
    while( fgets(string1, MAX_TEXT_FILE_LINE, al_fp) != NULL )
    {
        if( string1[0]=='F' ) /* Ocurrs when "Found at" is encountered */
        {
            fgets(string2, MAX_TEXT_FILE_LINE, al_fp);
            fgets(string3, MAX_TEXT_FILE_LINE, al_fp);
            if(currptr!=NULL)
            {
                if( strstr(string3, currptr->ref)!=NULL )
                {
                    fputs(string1, tmp_fp);
                    fputs(string2, tmp_fp);
                    fputs(string3, tmp_fp);
                    currptr = currptr->next;
                    moving = 1;
                    continue;
                }
            }
            moving = 0;
        }

        if(string1[0]=='D') /* Occurs when "Done." is encountered */
        {
            fgets(string2, MAX_TEXT_FILE_LINE, al_fp);
            fputs(string1, tmp_fp);
            fputs(string2, tmp_fp);
            fputs("\n", tmp_fp);
            break;
        }

        if(moving) fputs(string1, tmp_fp);
    }

    fclose(al_fp);
    fclose(tmp_fp);

    remove(alignmentfile);
    rename(tempfile,alignmentfile);

    return;

}
/***********************  BreakAlignments *********************************************************/
void   BreakAlignments(IL * headptr, char * alignmentfile)
{
    FILE * al_fp;
    FILE * out_fp;
    char outfile[_MAX_PATH+1];
    int nfiles;
    IL* currptr;
    int alignments;
    char headlns[30][MAX_TEXT_FILE_LINE]; /* to hold the heading section of the file */
    char buffer[MAX_TEXT_FILE_LINE]; /* holds one line of alignment data */
    int headcnt;
    char nextchar;
    int i,j;


    /* Find out how many alignments there are and how many files will
       be needed */
    for( alignments=0,currptr=headptr; currptr!=NULL; currptr=currptr->next, alignments++);
    nfiles = alignments / EO_MAX_TBL ;
    if( (alignments % EO_MAX_TBL )> 0 ) nfiles++ ;
    if(nfiles==0) nfiles = 1; /* make sure at least one file is generated */

    /* if only one file will be needed just rename file */
    if(nfiles==1)
    {
        MakeFileName( outfile, alignmentfile,1);

		remove(outfile);

        rename(alignmentfile,outfile);
		
        return;    
    }

    /* get heading lines */
    al_fp = fopen(alignmentfile,"r");
    for (i=0,headcnt=0; i<30; i++)
    {
        /* find if next line is the "Found at" line */
        nextchar = getc(al_fp);
        ungetc(nextchar,al_fp);
        if(nextchar=='F') break;
        fgets( headlns[i], MAX_TEXT_FILE_LINE, al_fp);
        headcnt++;
    }

    /* loop creating files */
    for(i=1;i<=nfiles;i++)
    {
        /* create name of ith file */
        MakeFileName( outfile, alignmentfile,i);
        /* open the file for writing */
        out_fp = fopen(outfile,"w");
        /* output heading */
        for(j=0;j<headcnt;j++) fputs(headlns[j],out_fp);
        /* output n of N identifier */
        fprintf(out_fp,"File %d of %d\n\n\n\n",i,nfiles);

        /*move alignments to the current output file*/
        for(j=0; j<EO_MAX_TBL; j++)
        {
            /* copy the first line of the alignment */
            fgets(buffer, MAX_TEXT_FILE_LINE, al_fp);
            fputs(buffer, out_fp);
            /* copy successive lines checking for the "Found at" and
               the "Done." lines */
            while(1)
            {
                nextchar = getc(al_fp);
                ungetc(nextchar,al_fp);
                if(nextchar=='F'||nextchar=='D') break;
                fgets(buffer, MAX_TEXT_FILE_LINE, al_fp);
                fputs(buffer, out_fp);
            }
            if(nextchar=='F') continue;
            if(nextchar=='D') break;
        }

        /* Output closing Lines */
        fprintf(out_fp,"\nDone.\n</PRE></BODY></HTML>\n");

        /*close file*/
        fclose(out_fp);
    }

    /*close original file*/
    fclose(al_fp);

    /* delete original file */
    remove(alignmentfile);

    return;
}

void    MakeDataFile(IL * headptr,char * datafile,int data)
{
    FILE* fp;
    IL* lpointer;
    int charcount;

    /* if data = 1 then produce new datafile overwriting old one */
    /* otherwise delete old file */
    if(data)
    {
       fp = fopen(datafile,"w");
        if (paramset.ngs != 1) {
          fprintf(fp,"Inverted Repeats Finder Program writen by:\n\nGary Benson\nBioinformatics Program\nBoston University\nVersion %s\n\n\n%s\n\n\n%s\n\n",versionstring, hsequence, hparameters);
        }
        for(lpointer=headptr;lpointer!=NULL;lpointer=lpointer->next)
        {
 
			/* print to the dat file */
  			fprintf(fp,"%d %d %d %d %d %d %d %.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %d %d ",
                    lpointer->leftfirst, lpointer->leftlast, lpointer->leftlength,
                    lpointer->rightfirst, lpointer->rightlast, lpointer->rightlength,
                    lpointer->separation, lpointer->matches, lpointer->indels,
                    lpointer->score, lpointer->atcount, lpointer->gccount,
                    lpointer->atpairs, lpointer->gcpairs, lpointer->gtpairs, lpointer->center, lpointer->avecenter );		


			/* left sequence */
            for(charcount=lpointer->leftfirst; charcount<=lpointer->leftlast;charcount++)
                fprintf(fp,"%c", Sequence[charcount]);

			fprintf(fp," ");
	
			/* right sequence */
            for(charcount=lpointer->rightfirst; charcount<=lpointer->rightlast;charcount++)
                fprintf(fp,"%c", Sequence[charcount]);

			//fprintf(fp," ");
	
			/* loop */
            //for(charcount=lpointer->leftlast+1; charcount<=lpointer->rightfirst-1;charcount++)
            //    fprintf(fp,"%c", Sequence[charcount]);

            fprintf(fp,"\n");
        }


        fclose(fp); 
    }
    else remove(datafile);
    return;
}

void    MakeMaskedFile(IL* headptr,int masked,char*  Sequence,char* maskfile)
{
    int count,printcr;
    int masker;
    FILE *fp;
    IL* lpointer;

    if(masked)
    {
        fp = fopen(maskfile,"w");

        /* Ouput sequence description from global variable to file*/
        fprintf(fp,">%s",&hsequence[10]);

        for(lpointer=headptr;lpointer!=NULL;lpointer=lpointer->next)
        {
            for(masker=lpointer->leftfirst; masker<=lpointer->leftlast; masker++)
				Sequence[masker]='N';
            for(masker=lpointer->rightfirst; masker<=lpointer->rightlast; masker++)
				Sequence[masker]='N';
        }
        printcr=0;
        for(count=1;Sequence[count]!='\0';count++)
        {
            fputc(Sequence[count], fp);
            printcr++;
            if(printcr>=60)
            {
                printcr=0;
                fputc('\n', fp);
            }
        }
        fputc('\n', fp);
        fputc('\n', fp);
        fclose(fp);
    }
    
    return;
}


#endif


