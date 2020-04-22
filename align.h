/*************************************************************************
*
*   align.h :   modified alignment routines done for the purposes of IRF 
*				algorithm. 
*				06/26/2002
*
*
*************************************************************************/

#ifndef ALIGN_H
#define ALIGN_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <malloc.h>
//changed 4.22.20 Gary Benson
#include <limits.h>

#define GLOBAL 0
#define LOCAL 1
#define TRUE 1
#define FALSE 0
#define RECENTERCRITERION 3

/* going in max between this and delta[g] */
#define MINBANDRADIUS 6

/* going out max between this and delta[g] */
//#ifdef IRF_MEGA
#define MINBANDRADIUS_OUTER 100
//#else
//#define MINBANDRADIUS_OUTER 20
//#endif


#define VLNN -1000
#define CUTOFFSCORE 0  /* alignment stops when this is hit */


struct pairalign{
  int length;
  int score;
  char *textprime,*textsecnd;
  int *indexprime,*indexsecnd;
} AlignPair;

        /*
        if(objptr==NULL)\
	    {\
	        debugerror("\n functionname: Out of memory!");\
	        exit(0);\
        }\ */

#define new1Darrayfunc(type,functionname,length)\
	 type *functionname(length)\
{\
        type *objptr=(type *)scalloc((length),sizeof(type));\
        memory_stats_print("\n functionname: requesting ", (length) * sizeof(type) );\
        return(objptr);\
 }
new1Darrayfunc(char,newAlignPairtext,length)

new1Darrayfunc(char,newLine,length)

new1Darrayfunc(int,newAlignPairindex,length)




int Maxlength;
int Reportmin=30;
int Maxrealrow, Maxrow, Maxcol, Largestjump;
int Maxscore;
int Maxrealcol=0;
int Classlength;
int Leftcount;
int Rightcount;
int Separation;
int pwidth=75;

#ifdef IRF_GT_MATCH_YES
#define match(a, b) (SM[256*((a)&0xDF)+(b)])
#else
#define match(a, b) ((((a)&0xDF)==(b))?Alpha:Beta)
#endif

 /* returns match mismatch matrix value */

int *Bandcenter = NULL;
/* int S[MAXINVERTEDREPEATSIZE+1][MAXINVERTEDREPEATSIZE+1];*/
int **S;

/* version 2A adds max3 and max2 */
#define max2(a,b) (((a)>=(b))?(a):(b))
/* returns max of 2 in order a,b */

#define max3(a,b,c) (((a)>=(b))?(((a)>=(c))?(a):(c)):(((b)>=(c))?(b):(c)))
/* returns max of 3 in order a,b,c */

#define max4(a,b,c,d)(((a)>=(b))?(((a)>=(c))?(((a)>=(d))?(a):(d)):(((c)>=(d))?(c):(d))):(((b)>=(c))?(((b)>=(d))?(b):(d)):(((c)>=(d))?(c):(d)) )) 
/* returns max of 4 in order a,b,c,d */

#define min3(a,b,c) (((a)<=(b))?(((a)<=(c))?(a):(c)):(((b)<=(c))?(b):(c)))
/* returns max of 3 in order a,b,c */


//debugmessage("\nc1: %c     c2: %c    l: %d   i: %d   j:%d\n",c1,c2,l,i,j);return;
#define fill_align_pair(c1,c2,l,i,j)\
  AlignPair.textprime[l]=c1;\
  AlignPair.textsecnd[l]=c2;\
  AlignPair.indexprime[l]=i;\
  AlignPair.indexsecnd[l]=j

/* points to the second string (the "pattern") inside the SequenceComplement */
char *EC;

int ECstart;

double Rows=0;




/*******************************************************************/

/***************************** narrowbandwrap()  **************************/

/*******************************************************************/
#define	ALIGNIN  1
#define	ALIGNOUT 2

void narrowbandwrap(int center,int start,int bandradius,int tuplesize, int option, int reportmin, int g, int h, int foundat)
	 
	 /* start is end of pattern in text */
	 /* tuplesize is the size of tuple used for this pattern size */
	 
	 /* started Feb. 7 1997 */
	 /* simplified from the original version to work with IRF on  July. 5 2002 */

{
  register int *pup, *pdiag, *pcurr,*pleft;
  int c,r,realr,end_of_trace,passed_starting_point=0,
  maxscore,maxrow,maxcol;
  int maxrealrow;
  char currchar;
  int w, matches_in_diagonal,matchatmax_col,i,k,maxrowscore,width,bigwidth,largestjump,
  lastmatchatmax_col,match_yes_no;
	

  w=bandradius;
  width=2*w+1;
  bigwidth=4*w+2;
  r=0;
  maxscore=0;
  realr=start+1;
  Bandcenter[r]=0;
  matches_in_diagonal=0;
  matchatmax_col=-2;
  

#ifdef IS_ALIGNMENT_PROFILER_ON
  /* alignment profiler */
	AlignmentProfiler.alignmentTotal[g].interval[h]++;
#endif

  /* to make sure we are not going to overflow our buffers */
  if((MAXBANDWIDTH*2)<bigwidth)
  {
	debugerror("\nIn narrowbandwrap, MAXBANDWIDTH*2: %d exceeded by bigwidth: %d\n",
	   MAXBANDWIDTH*2,bigwidth);
	exit(-1);
  }

		 
  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  /**********************************************************************/
  

  /* initialize top row */
  pcurr=&S[r][0];
  for(i=0;i<bigwidth;i++)
  {
	  if (i<=2*w || i>=3*w)
		*pcurr=VLNN;
      else 
		*pcurr=(i-2*w-1)*DeltaIndel;

	  pcurr++;
  }


  /* compute until end of trace */      
  end_of_trace=FALSE;  
  while ((realr>1) && (r<MAXWRAPLENGTH))
  {
/*	  if (c>MAXINVERTEDREPEATSIZECONSTANT) {

		  debugerror("Repeat too large. Stopping alignment!");
		  break;
		 
	  }
*/
	  r++;
	  realr--;
	  Rows++;
	  maxrowscore=VLNN;
	  lastmatchatmax_col=matchatmax_col;
	  currchar=Sequence[realr];  

	  /* test to recenter the bandcenter */
	  if(matches_in_diagonal>=tuplesize)
	  {
          // remember largest jump, Gelfand for version 3.03
          largestjump = Bandcenter[r] - matchatmax_col;
          if (largestjump<0) largestjump=-largestjump;
          if (largestjump>Largestjump) Largestjump = largestjump;

		  // recenter band 
		  Bandcenter[r]=matchatmax_col+1; // or should it be moved by 2?
		


	  } 
	  else 

	  {
		  // don't recenter 
		  Bandcenter[r]=Bandcenter[r-1]+1;
	  }
	  
		

	  /* init current row */
      pcurr=&S[r][0];
      for(i=0;i<bigwidth;i++)
	  {
	    *pcurr=VLNN;
	    pcurr++;
	  }
	 
 
	  /* change of bandcenter determines which inputs go into which cells */
	  k=Bandcenter[r]-Bandcenter[r-1];
	  pcurr=&S[r][w+1];
	  pleft=pcurr-1;
	  pdiag=&S[r-1][w+k];
	  pup=pdiag+1;
	  c=Bandcenter[r]-w-1;
	  
	  
	  /* test for tow centers meeting when alignment is going inward */
	  if (option==ALIGNIN && realr<ECstart+c) {
		  break;
	  }
	  
	  
	  /* test for reaching the end of the string when alignment is going outward */
	  if (c+ECstart>Length) {	
		  break;
	  }
	  
	  
	  for(i=0;i<=2*w;i++) 
	  {


			  /* because we don't start at the beginning*/
			  if (c<0) 
				  {match_yes_no=Beta; }
			  else {

				  /* because we don't want to compare characters that overlap */
				  if (option==ALIGNIN && c+ECstart>=realr) 
					{match_yes_no=Beta;  }
				  else {

					  /* because we don't want to look outside the sequence boundaries */
					  if (option==ALIGNOUT && c+ECstart>Length) 
						{match_yes_no=Beta;  }
					  else
						match_yes_no=match(currchar,EC[c]);

				  }
		
			  }
			 
			  
			  /* calculate the value of the cell */
			  *pcurr=max3(*pdiag+match_yes_no,*pup+DeltaIndel,*pleft+DeltaIndel);
			  

#ifdef IS_ALIGNMENT_PROFILER_ON
  /* alignment profiler */
	AlignmentProfiler.cellsTotal[g].interval[h]++;
	AlignmentProfiler.cellsTemp++;  // this if the alignment is successfull will this value be used
#endif

	
			  /* test_trace_and_forward_maxscore (used to be a macro) */		  
			  if(*pcurr>maxscore) {
			  //if(*pcurr>=maxscore) {
				  maxscore=*pcurr;
				  maxrealrow=realr;
				  maxcol=c;
				  maxrow=r;
			  }

			 /* test_maxrowscore_with_match (used to be a macro) */
			 if(*pcurr>maxrowscore)
			 {
				 maxrowscore=*pcurr;
				 if((*pcurr==(*pdiag+Alpha))&&(match_yes_no==Alpha))
					 matchatmax_col=c;
                 else if(Gtmatchyesno&&(*pcurr==(*pdiag+GTAlpha))&&(match_yes_no==GTAlpha)) /* added by Gelfand to make sure alignment jumps correcly when a run of GT matches is found */
                     matchatmax_col=c;
				 else matchatmax_col=-2;
			 }
			 
			 /* move on to the next character of the pattern */
			 pcurr++;
			 pleft++;
			 pdiag++;
			 pup++;
			 c++;
	 }  

	/* a condition to stop the alignment*/
    //if (option==ALIGNIN || passed_starting_point) {
	    if (maxrowscore<=CUTOFFSCORE) {
			    end_of_trace=TRUE;
			    break;
	    }
    //} else {

        /* added in version 3.02 */
        /* This sometimes happens when there is a large gap and the OUTWARD alignment cannot breach it  */
        /*                                     c                                 i                      */
        /* ***************************---------------------*******     *******************              */                                   
        /* Let's NOT test for  CUTOFFSCORE going OUT until we pass i                                    */                                                                             
        //passed_starting_point = ( (Bandcenter[r]-w-1+ECstart) >= foundat );

    //}


	/* test bor band shift */
	if((matchatmax_col-lastmatchatmax_col)==1)
	  matches_in_diagonal++;
	else matches_in_diagonal=0;



	} 
	
	/* test for report */
	if (maxscore>=reportmin)
	{
		Maxrealrow=maxrealrow;
		Maxrow=maxrow;
		Maxcol=maxcol;
		Maxscore=maxscore;
		
	}

//if (foundat==6049328&&center==12098505&&option==ALIGNIN) printAlignmentMatrix(start,w,r);
//if (foundat==6049328&&center==12098505&&option==ALIGNOUT) printAlignmentMatrix(start,w,r);

}


/*******************************************************************/

/*********************  printAlignmentMatrix()  ***********/

/*******************************************************************/

int printAlignmentMatrix(int start, int w, int r) {


	int j,z,i;
	int *pcurr, bigwidth,realr;

	bigwidth=4*w+2;
	realr=start+1;

	fprintf(stderr,"\n");
	
	for (j=0;j<=r;j++) {
		pcurr=&S[j][0];
		
		fprintf(stderr,"%d,",realr);
		for (z=0; z<Bandcenter[j]; z++)
			fprintf(stderr,",");
		
		
		for (i=0;i<bigwidth;i++)
		{
			fprintf(stderr,"%d",*pcurr);
			
			if (i<bigwidth-1) fprintf(stderr,",");
			pcurr++;
		}
		realr--;
		fprintf(stderr,"\n");
	}
	
	
	//exit(0);
    return 0;

}

/*******************************************************************/

/*********************  get_narrowband_pair_alignment()  ***********/

/*******************************************************************/


void get_narrowband_pair_alignment(int bandradius)
	 
	 
	 /* for a repeat of EC, do a traceback alignment, */
	 /* ending at row Maxrow and column Maxcol */
	 
	 
	 
	 
#define test_match_mismatch \
	 if (S[r][i]==S[r-1][upi-1]+match(x[realr], y[c]))\
{\
   length++;\
	   fill_align_pair(x[realr],y[c],length,realr,c);c--;\
	   realr++;\
		 r--;\
		   i=upi-1;\
 }\
else

#define test_up \
if (S[r][i]==S[r-1][upi]+DeltaIndel)\
{\
   length++;\
	 fill_align_pair(x[realr],'-',length,realr,c+1);\
	 realr++;\
	   r--;\
		 i=upi;\
 }\
else

#define test_left \
if(i==0)\
{\
   debugerror("\nget_pair_alignment_with_copynumber: error in trace back");\
	 debugerror("\nattempted to compute left branch when i==0");\
	   debugerror("\nS[%d][%d]=%d",r,i,S[r][i]);\
	 break;\
 }\
else if (S[r][i]==S[r][i-1]+DeltaIndel)\
{\
   length++;\
	   fill_align_pair('-',y[c],length,realr+1,c);c--;\
	   i=i-1;\
 }\
else

#define report_error_match_up_left \
{\
   debugerror("\nget_pair_alignment_with_copynumber: error in trace back");\
	 debugerror("\nr=%d i=%d c=%d S: row=%d  column=%d  upi=%d  Sequence: realrow=%d  EC: realcol=%d",\
		r,i,c,r,i,upi,realr,c);\
		  debugerror("\nS=%d  Sleft=%d  Sup=%d  Sdiag=%d  match=%d",\
			 S[r][i],S[r][i-1],S[r-1][upi],\
			 S[r-1][upi-1],match(x[realr], y[c]));\
			   break;\
 }





{
  
  int i,length;
  char *x,*y;
  int realr,r,c,w;
  int upi,k;

  
	w=bandradius;
	x=Sequence;
	y=EC; 
  



	realr=Maxrealrow;
	r=Maxrow;
	c=Maxcol;
	k=Maxcol-Bandcenter[r];
	i=2*w+k+2;			
  
	AlignPair.score=S[r][i]; 
	length=0;

	for (;;)
	{	





	  /* stop at zeros or -1000 for local */
	  if (S[r][i]<=VLNN || c<0 || r<=0)
	  { 
			AlignPair.length=length;
			return;
	  }
	  else
	  {
			k=c-Bandcenter[r-1];
			upi=2*w+k+2;
	
	
/* DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#ifdef IRFDEBUG
  debugmessage("S[%d][%d]=%d   y[%d]=%c x[%d]=%c  upi=%d S[r-1][upi-1]=%d match=%d\n",r,i,S[r][i],c,y[c],realr,x[realr],upi,S[r-1][upi-1],match(x[realr], y[c]));
#endif
/* DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


			test_match_mismatch
				test_up
					test_left
						report_error_match_up_left; 
		

	  }



	}



}

/************************************************************/ 
int d_range(int d) 
{
  return((int)floor(2.3*sqrt(Pindel*d)));
}
   
//#ifdef dsfaasdfasdfasd
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

/*******************************************************************/
/*******************************************************************/
void Init_centerseenlist()
{
  Centerseenlist->next=NULL;
}
/*******************************************************************/
/*******************************************************************/

void free_centerseenlist()
{
  struct centerlistelement *entry, *entrylast;
  entry=Centerseenlist->next;
  Centerseenlist->next=NULL;
  while (entry!=NULL)
  {
	entrylast=entry;
	entry=entry->next;
	sfree(entrylast);
  }
}

/*******************************************************************/
/*******************************************************************/
void add_to_centerseenlist(int center,int rightstart, int rightend,int  score, int g, int h)
	 
/* center:				calculated center of the align pair */
/* rightstart:			beginning of the right pair */
/* rightend:			ending of the right pair */
	 

{
  struct centerlistelement *entry,*entrylast,*temp,*templast;


  entry=(struct centerlistelement *)
	scalloc(1,sizeof(struct centerlistelement));

  memory_stats_print("\n add_to_centerseenlist: requesting", 1 * sizeof(struct centerlistelement) );
  
/*  if(entry==NULL){debugerror("\nAdd_to_centerseenlist: Out of memory!");exit(0);} */


  entry->index=center;
  entry->rightstart=rightstart;
  entry->rightend=rightend;
  entry->delta=d_range(rightend-rightstart+1);
  entry->g = g;
  entry->h = h;
  //entry->searchedThrough=0;

  /* insert at the appropriate place, sorted by the index + delta descending */
  templast=Centerseenlist;
  for (temp=Centerseenlist->next; temp!=NULL; temp=temp->next) {
	    if (temp->index+2*temp->delta<entry->index+2*entry->delta) break;
		templast=temp;
  }
  entry->next=temp;
  templast->next=entry;
#ifdef IS_ALIGNMENT_PROFILER_ON
  Centerseenlist_Size++;
#endif
  /* now go until the end of the list and prune the old ones */
  entry=Centerseenlist->next;
  entrylast=Centerseenlist;
  


  while (entry!=NULL) {
		

		/* remove if too old */
		if (entry->index+2*entry->delta<Runningmincenter3[NTS]) 
			{
				entrylast->next=entry->next;
				temp=entry;
				entry=entrylast;
				sfree(temp);
#ifdef IS_ALIGNMENT_PROFILER_ON
                Centerseenlist_Size--;
#endif
			}

		/* go to next entry */
		entrylast=entry;  
		entry=entry->next;
  }

	

}

/*******************************************************************/
/*******************************************************************/

int search_for_center_match_in_centerseenlist(int center, int i, int g, int h)

/* finding a match means ruling 
out redoing the alignment for that center  */

/* center:		center between two matches */     

{
	struct centerlistelement *entry;
	int deltalow,deltahigh;
	entry=Centerseenlist->next;
	
	
	while (entry!=NULL)
	{	
		
		/* determine lowest and highest possible centers in range */
		/* note we add or subtract delta, but don't count tc */
		/* so width of range is delta+1 */
		deltalow=center-2*entry->delta;    
		deltahigh=center+2*entry->delta;
		

		/* check if we have a match */
		if (deltalow <= entry->index) { 
			if (entry->index <= deltahigh) /* if centers coencise */
			{	
			
				/* if projected size falls insided the known align pair */
				/* this check is to make sure that we don't eliminate the possibility */
				/* of an irepeat having the same center but but having a much larger/smaller radius. */
				if (entry->rightstart <= i && i <= entry->rightend) { 

                    /* added by Gelfand on 11/15/2004 for version 3.03 */
                    /* to allow doing nearby centers found at a higher interval */
                    if (paramset.intcheck==1) {

                       /* check larger intervals */
                       if ( g < entry->g || ( g == entry->g && h <= entry->h ) )
                           return(TRUE);

                    } else if (paramset.intcheck==2) {

                       /* check all intervals */
                       if ( g == entry->g && h == entry->h ) 
                           return(TRUE);

                    } else {

                           return(TRUE);
                    }

				}


                /* MAYBE DO THAT THING DR. BENSON SAID ABOUT REMOVING MATCHES SO WE WON'T CALL IT AGAIN??? */
				/* new */
				//entry->searchedThrough++;

				//if (entry->searchedThrough >= 3 ) {
				
					//recalulate
					//entry->rightend = (int)( i + ( i - entry->rightstart ) *.1 );
					//entry->delta=d_range(entry->rightend-entry->rightstart+1);
					//entry->searchedThrough=0;
				//}
				/* end of new */
			}
		} else return(FALSE); /* if lower, means we passed the point of interest (list is sorted) 
		and should indicate that we found nothing */

		
		/* go to next entry */
		entry=entry->next;
	}
	return(FALSE);
	
}


/*******************************************************************/
/*******************************************************************/

int search_for_center_match_in_centerseenlist_exact(int center, int rightstart, int rightend, int g, int h)

/* finding a match means no need to write to harddrive  */

/* center:		center between two matches */     

{
	struct centerlistelement *entry;
	int deltalow,deltahigh;
	entry=Centerseenlist->next;
	
	
	while (entry!=NULL)
	{	
		
		/* determine lowest and highest possible centers in range */
		/* note we add or subtract delta, but don't count tc */
		/* so width of range is delta+1 */
		deltalow=center-2*entry->delta;    
		deltahigh=center+2*entry->delta;
		

		/* check if we have a match */
		if (deltalow <= entry->index) { 
			if (entry->index <= deltahigh) /* if centers coencise */
			{	
			
				/* if projected size falls insided the known align pair */
				/* this check is to make sure that we don't eliminate the possibility */
				/* of an irepeat having the same center but but having a much larger/smaller radius. */ 
				if (entry->rightstart==rightstart && entry->rightend==rightend) { 

                    /* added by Gelfand on 11/15/2004 for version 3.03 */
                    /* to allow doing nearby centers found at a higher interval */
                    if (paramset.intcheck==1) {

                       /* check larger intervals */
                       if ( g < entry->g || ( g == entry->g && h <= entry->h ) )
                           return(TRUE);

                    } else if (paramset.intcheck==2) {

                       /* check all intervals */
                       if ( g == entry->g && h == entry->h ) 
                           return(TRUE);

                    } else {

                           return(TRUE);
                    }

				}


			}
		} else return(FALSE); /* if lower, means we passed the point of interest (list is sorted) 
		and should indicate that we found nothing */

		
		/* go to next entry */
		entry=entry->next;
	}
	return(FALSE);
	
}


/************************************************************/ 

/*******************************************************************/

/***************************** statistics() **************************/

/*******************************************************************/

void get_3mod_centers(int *r1,int *r2, int *r3) {
	
	int valuelowest=INT_MAX, valuehighest=-1 ,i,center,size;
	int *centersarray,c1,c2,c3,v1,v2,v3;
	
	/* find the losest and highest centers */
	for (i=1; i<AlignPair.length; i++) {
		if ((AlignPair.textprime[i]==AlignPair.textsecnd[i])&&
			(AlignPair.textprime[i]!='-') &&
			(AlignPair.textsecnd[i]!='-')) {

			center=AlignPair.indexprime[i]+AlignPair.indexsecnd[i]+ECstart;

			if (center>valuehighest) valuehighest=center; 

			if (center<valuelowest) valuelowest=center;

		}
	}


	/* allocate buffer */
	size=valuehighest-valuelowest+1;
	centersarray=(int*) scalloc(size+100,sizeof(int));

    memory_stats_print("\n get_3mod_centers: requesting", (size+100) * sizeof(int) );

/*	if (centersarray==NULL ) {
		debugerror("ERROR(get_3mod_centers): cannot allocate %d bytes of memory. Aborting!",(size+100)*(sizeof(int)));
		exit(0);
	}
*/

	/* calculate all the centers */
	for (i=1; i<AlignPair.length; i++) {
		if ((AlignPair.textprime[i]==AlignPair.textsecnd[i])&&
			(AlignPair.textprime[i]!='-') &&
			(AlignPair.textsecnd[i]!='-')) {

			// get center
			center=AlignPair.indexprime[i]+AlignPair.indexsecnd[i]+ECstart;

			// normalize it
			center-=valuelowest;

			// store it
			centersarray[center]++;

		}
	}


	/* get 3 most frequent */
	c1=c2=c3=0;
	v1=v2=v3=0;
	for (i=0; i<size; i++) {
		if (centersarray[i]>c1) { c1=centersarray[i]; v1=i; }
	}
	for (i=0; i<size; i++) {
		if (centersarray[i]>c2&&centersarray[i]!=c1) { c2=centersarray[i]; v2=i; }
	}
	for (i=0; i<size; i++) {
		if (centersarray[i]>c3&&centersarray[i]!=c1&&centersarray[i]!=c2) { c3=centersarray[i]; v3=i; }
	}


	/* unnormalize the centers and copy and return positions*/
	*r1=(c1>0)?v1+valuelowest:0;
	*r2=(c2>0)?v2+valuelowest:0;
	*r3=(c3>0)?v3+valuelowest:0;


	/* free memory */
	sfree(centersarray);
	
}

/**********************       print_alignment_headings     *********************************************/
int print_alignment_headings(void)
{
	
  int c1, c2, c3;

  /* headings */
  fprintf(Fptxt,"\n\n<A NAME=\"%d--%d,%d,%d--%d,%d\">",
		 AlignPair.indexprime[1],AlignPair.indexprime[AlignPair.length],Leftcount,
		 ECstart+AlignPair.indexsecnd[AlignPair.length],ECstart+AlignPair.indexsecnd[1],Rightcount);
  fprintf(Fptxt,"</A>");


    #if defined(WINDOWSGUI)
        fprintf(Fptxt,"<P>See <FONT COLOR=\"#0000FF\">Alignment Explanation</FONT> in Inverted Repeats Finder Help</P><BR>\n");
    #elif defined(WINDOWSCONSOLE)
        fprintf(Fptxt,"<A HREF=\"http://tandem.bu.edu/irf/irf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n");
    #elif defined(UNIXGUI)
        error: Unix GUI code not implemented
    #elif defined(UNIXCONSOLE)
        fprintf(Fptxt,"<A HREF=\"http://tandem.bu.edu/irf/irf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n");
    #endif



  fprintf(Fptxt,"    Indices: %d--%d,%d--%d  Loop: %d  Score: %d\n    Left Length: %d  Right Length: %d  Average length: %.1lf\n",
		 AlignPair.indexprime[1],AlignPair.indexprime[AlignPair.length],
		 ECstart+AlignPair.indexsecnd[AlignPair.length],ECstart+AlignPair.indexsecnd[1],Separation,Maxscore,Leftcount,Rightcount,(Leftcount+Rightcount)/2.0);
		 
  /* print 3 most frequent centers */
  get_3mod_centers(&c1,&c2,&c3);

  fprintf(Fptxt,"    3 most common centers (desc):  %d  %d  %d\n\n",c1,c2,c3);

  return c1;
  
} 


void shift_pattern_indices (int patternsize)
{
	int downshift, upshift, l;

	downshift = AlignPair.indexsecnd[1];
	upshift = patternsize-downshift;
	for(l=1;l<=AlignPair.length;l++)
	{
		if(AlignPair.indexsecnd[l]>=downshift)
			AlignPair.indexsecnd[l]-=downshift;
		else
			AlignPair.indexsecnd[l]+=upshift;
	}
}

/**********************       print_alignment     *********************************************/
void print_alignment(void) {

    int i, j, g, h, m,alignmentslab=75;


	/* show a short left flanking sequence */

        if(AlignPair.indexprime[1]!=1)
        {
            m=AlignPair.indexprime[1]-10;
            if (m<1) m=1;
            j=m;
            fprintf(Fptxt,"  %9d >> (LF) ",j);
            for(i=1;i<=10;i++)
            {
                fputc(Sequence[j],Fptxt);
                j++;
                if(j==AlignPair.indexprime[1]) break;
            }
            fprintf(Fptxt,"\n");
        }


     /* show a short right flanking sequence */

          if(AlignPair.indexsecnd[1]+ECstart!=Length)
        {
            m=ECstart+AlignPair.indexsecnd[1]+10;
            if (m>Length) m=Length;
            j=ECstart+AlignPair.indexsecnd[1]+1;
            fprintf(Fptxt,"  %9d << (RF) ",m);
            for(i=1;i<=10;i++) 
            {
                fputc(Sequence[m],Fptxt); 
                m--;
                if(m<j) break;
            }
            fprintf(Fptxt,"\n\n");
        }


	 /* output the alignment, #alignmentslab nucleotiedes per line */
        g=1;
        for (;;)
        {
            j=g;
            fprintf(Fptxt,"            ");
            h=0;
            fputc('\n',Fptxt);
            fprintf(Fptxt,"  %9d >> ",AlignPair.indexprime[j]); 
           
			
			for (i=0;i<alignmentslab; ++i)
            {
                if (j==AlignPair.length+1) break;

                fputc(AlignPair.textprime[j++],Fptxt);
                h++;
            }  

			fprintf(Fptxt," >> %d\n",AlignPair.indexprime[g+h-1]);
            j=g;
            fprintf(Fptxt,"  %9d << ",ECstart+AlignPair.indexsecnd[j]); 
            h=0;
           
			
			for (i=0;i<alignmentslab; ++i)
            {
                if (j==AlignPair.length+1) break;

                if ((AlignPair.textprime[j]==AlignPair.textsecnd[j])&&
                    (AlignPair.textprime[j]!='-') &&
                    (AlignPair.textsecnd[j]!='-'))
						{fputc('*',Fptxt); j++; }
				else
						{
#ifdef IRF_GT_MATCH_YES
                        if (Gtmatchyesno && ((AlignPair.textprime[j]=='G'&&AlignPair.textsecnd[j]=='A') || (AlignPair.textprime[j]=='T'&&AlignPair.textsecnd[j]=='C')))
                            {fputc('!',Fptxt); j++; }
                        else
                            fputc(Complementascii[AlignPair.textsecnd[j++]],Fptxt);
#else                    
                            fputc(Complementascii[AlignPair.textsecnd[j++]],Fptxt);
#endif
                
                        } 

                h++;
            } 

			fprintf(Fptxt," << %d\n\n",ECstart+AlignPair.indexsecnd[g+h-1]);
            g=j;
            if (j==AlignPair.length+1) break;
        }
		  fprintf(Fptxt,"\n");

}


/************************************************************/ 

/******************** flanking sequence *********************/

void print_flanking_sequence(int flank_length)
{

  int m,n,k,i,j;

  m=AlignPair.indexprime[1]-flank_length;
  if (m<1) m=1;
  n=AlignPair.indexsecnd[1]+ECstart+flank_length;
  if (n>Length) n=Length;

  if (m==AlignPair.indexprime[1])
    {
      fprintf(Fptxt,"\nLeft flanking sequence: None");
    }
  else
    {
      fprintf(Fptxt,"\nLeft flanking sequence: Indices %d -- %d\n",
	      m,AlignPair.indexprime[1]-1);
      k=AlignPair.indexprime[1];
      j=m;
      for(;;)
	{
	  for(i=1;i<=pwidth-10;i++)
	    {
	      fputc(Sequence[j],Fptxt);
	      j++;
	      if(j>=k) break;
	    }
	  fprintf(Fptxt,"\n");
	  if(j>=k) break;
	} 
    }
  if (n==(AlignPair.indexsecnd[1]+ECstart))
    {
      fprintf(Fptxt,"\n\nRight flanking sequence: None");
    }
  else
    {
      fprintf(Fptxt,"\n\nRight flanking sequence: Indices %d -- %d\n",
	      AlignPair.indexsecnd[1]+ECstart+1,n);
      j=AlignPair.indexsecnd[1]+ECstart+1;
      for(;;)
	{
	  for(i=1;i<=pwidth-10;i++)
	    {
	      fputc(Sequence[j],Fptxt);
	      j++;
	      if(j>n) break;
	    }
	  fprintf(Fptxt,"\n");
	  if(j>n) break;
	} 
    }
  fprintf(Fptxt,"\n\n");

  return;

}

/*******************************************************************/

/**************************  reverse()  ****************************/

/*******************************************************************/


void reverse()
	 
	 /* reverses the alignment in AlignPair */
	 
	 
	 
	 
{
  int j, tempi,ml;
  char temp;
  
  ml=AlignPair.length;

  for (j=1;j<=ml;j++)
  {
	  
	temp=AlignPair.textprime[j];
	AlignPair.textprime[j]=Complementascii[AlignPair.textsecnd[ml-j+1]];
	AlignPair.textsecnd[ml-j+1]=Complementascii[temp];
	
	tempi=AlignPair.indexprime[j];
	AlignPair.indexprime[j]=AlignPair.indexsecnd[ml-j+1] + ECstart;
	AlignPair.indexsecnd[ml-j+1]=tempi - ECstart;
	  
  }
}


/**********************       get_statistics     *********************************************/
void get_statistics()
{
 
	int ACGTcount[26],count,i,match,mismatch,indel,x,lp,atpairs,gcpairs,gtpairs,most_common_center,dummy,center;
	
	match=0;
	mismatch=0;
	indel=0;
	x=1;
    dummy = 0;
	atpairs=0;
	gcpairs=0;
    gtpairs=0;
    most_common_center=0;


	/* get the ACGTcount */
	ACGTcount['A'-'A']=0;
	ACGTcount['C'-'A']=0;
	ACGTcount['G'-'A']=0;
	ACGTcount['T'-'A']=0;
	count=0;
	Leftcount=0;
	for(i=1;i<=AlignPair.length;)
	{
		if(AlignPair.textprime[i]!='-')
		{
			ACGTcount[AlignPair.textprime[i]-'A']++;
			count++;
			Leftcount++;
		}
		i++;
	}
	Rightcount=0;
	for(i=1;i<=AlignPair.length;)
	{
		if(AlignPair.textsecnd[i]!='-')
		{
			ACGTcount[AlignPair.textsecnd[i]-'A']++;
			count++;
			Rightcount++;
		}
		i++;
	}	

	/* get the matches mismatches */
	lp=1;
	
	while (lp<=AlignPair.length) {


		if((AlignPair.textprime[lp]!='-')&&(AlignPair.textsecnd[lp]!='-')) {
			
            if(AlignPair.textprime[lp]==AlignPair.textsecnd[lp]
#ifdef IRF_GT_MATCH_YES
                || (Gtmatchyesno && ((AlignPair.textprime[lp]=='G'&&AlignPair.textsecnd[lp]=='A') || (AlignPair.textprime[lp]=='T'&&AlignPair.textsecnd[lp]=='C')))
#endif
            ) /* ver 3.02, for GT match */
			{
				match++;	

#ifdef IRF_GT_MATCH_YES     
                /* more complicated check if GT pair is allowed */
                if ((AlignPair.textprime[lp]=='A'&&AlignPair.textsecnd[lp]=='A') || (AlignPair.textprime[lp]=='T'&&AlignPair.textsecnd[lp]=='T'))
					atpairs++;
				else if ((AlignPair.textprime[lp]=='G'&&AlignPair.textsecnd[lp]=='G') || (AlignPair.textprime[lp]=='C'&&AlignPair.textsecnd[lp]=='C'))
					gcpairs++;
                else
                    gtpairs++;

#else
				if (AlignPair.textprime[lp]=='A' || AlignPair.textprime[lp]=='T')
					atpairs++;
				else
					gcpairs++;
#endif


			} else
			{
				mismatch++;
			}

			
		} else if(AlignPair.textsecnd[lp]!='-') {
			indel++;
		} else if(AlignPair.textprime[lp]!='-') {
			indel++;
		} else {
			debugerror("this should not be happening!!!!!\n");
			exit(-1);
		}

		lp++;
		
	}
	x=match+mismatch+indel;
	//printf("indel: %d   match: %d  sparation: %d\n",indel,match,Separation);

if (!paramset.HTMLoff) {

	/* output the alignment */
	most_common_center = print_alignment_headings();

	
	/* print the alignment itself */
	print_alignment();


	/* add some statistics to aligment */
	fprintf(Fptxt,"\nStatistics");
	fprintf(Fptxt,"\nMatches: %d (%0.2f%%),  Mismatches: %d (%0.2f%%), Indels: %d (%0.2f%%)\n",
          match,    100*(float)match/x,
          mismatch, 100*(float)mismatch/x,
          indel,    100*(float)indel/x);
	fprintf(Fptxt,"\nATpairs: %0.2f%%,  CGpairs: %0.2f%%",100*(float)atpairs/match,100*(float)gcpairs/match);
#ifdef IRF_GT_MATCH_YES
    fprintf(Fptxt,",  GTpairs: %0.2f%%",100*(float)gtpairs/match);
#endif
	fprintf(Fptxt,"\n\n\n");

			  
	
	if(paramset.flankingsequence)
	{
		  //reverse();
		  print_flanking_sequence(paramset.flankinglength);
		  //reverse();
	}     

} else get_3mod_centers(&most_common_center,&dummy,&dummy);

     center = AlignPair.indexprime[AlignPair.length] + (ECstart+AlignPair.indexsecnd[AlignPair.length]);

	/* print to the dat file */
  	 fprintf(Fpdat,"%d--%d,%d,%d--%d,%d %d %d %d %d %d %d %d %.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %d %d\n",
		
		// AlignPair.indexprime[1],AlignPair.indexprime[AlignPair.length]
		//		  best_match_distance,Copynumber,Classlength,( int ) OUTPUTcount,


		 AlignPair.indexprime[1],AlignPair.indexprime[AlignPair.length],Leftcount,
		 ECstart+AlignPair.indexsecnd[AlignPair.length],ECstart+AlignPair.indexsecnd[1],Rightcount,

		 
		 AlignPair.indexprime[1],AlignPair.indexprime[AlignPair.length],Leftcount,
		 ECstart+AlignPair.indexsecnd[AlignPair.length],ECstart+AlignPair.indexsecnd[1],Rightcount,

		
		 Separation,

		 (100*(float)match/x),(100*(float)indel/x),AlignPair.score,
		 (100*((float)ACGTcount['A'-'A']+(float)ACGTcount['T'-'A'])/count),
		 (100*((float)ACGTcount['C'-'A']+(float)ACGTcount['G'-'A'])/count),


		 (100*(float)atpairs/match),(100*(float)gcpairs/match),(100*(float)gtpairs/match),/* Gelfand, Nov 30, 2004, added GTPairs */ 
         
         
         center, most_common_center /* Gelfand, Nov 30, 2004, added center and most common center */ ); 

	 //fflush(Fpdat);

}

 

#endif
