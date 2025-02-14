/***************************************************************/
/*                                                             */
/*       Inverted Repeats Finder                               */
/*       written by Gary Benson                                */
/*                                                             */
/*       started 02-08-01                                      */
/*                                                             */
/*                                                             */
/***************************************************************/

#include <stdio.h>
#include <stddef.h> /* has size_t definition */
#include <stdlib.h> /* has calloc definition */
#include <ctype.h>  /* includes toupper(c) */
#include <string.h> /* includes strncat() */
#include <math.h>   /* for ceil function */

//#define EASY_LIFE_DEBUG

#include "easylife.h"
#include "memorycheck.h"
#include "test.interval.sums.h"  /* added by Gary Benson 12/04/02 for testing */
#include "irf.h"
#include "centerlist.h" /* routines for the centerlists */
#include "krunsums2Anewt.h"
#include "intervals.h"  /* routines for intervals */
#include "centerlist.11.03.h" /* routines for match element linked lists 
                                 Gary Benson 11/05/03 */
#include "centertag.h"  /* routines for tags */
#include "new.centerlists.h"
#include "irfrun.h"
#include "print.interval.sums.h"  /* added by Gary Benson 12/04/02 for testing */
#include "profiler.h"
#include "align.h"


/*******************************************************************/

/***************************** debugmessage function ***************/

/*******************************************************************/

//#define IRFDEBUG
//#define IRFDEBUG2
    
void doCriticalErrorAndQuit(const char *format, ... )	{ 

	va_list				argp;
	

	if (format==NULL) return;

	va_start(argp, format);
	
	vfprintf(stderr,format, argp); 
	fflush(stderr);

	va_end(argp);

	exit(-1);

}

/* end of added on 10/30/03 */

void debugerror(const char *format, ... )	{ 
	
	va_list				argp;
	

	if (format==NULL) return;

	va_start(argp, format);
	
	vfprintf(stderr,format, argp); 
	//vfprintf(Fptxt, format, argp); 
	fflush(stderr);
	//fflush(Fptxt);

	va_end(argp);

}
/*******************************************************************/
/**************************  init_index()  *************************/
/*******************************************************************/
void init_index(void)

{
  int i;

  /* index has 256 entries so that finding the entries for A, C, G and T */
  /* require no calculation */
  Index=(int *)scalloc(256,sizeof(int));
  memory_stats_print("\n init_index: requesting", 256 * sizeof(int) );
/*  if(Index==NULL){debugerror("\nInit_index: Out of memory!");exit(0);} */

  for (i=0; i<256; i++)
	Index[i]=-1;

  Index['A']=0;
  Index['C']=1;
  Index['G']=2;
  Index['T']=3;
}

/******************************************************************************************/
/* Added by Gelfand on 10/14/2004 to search for "similar" tuples that allow GT match      */
/* This function goes through all the codes that are SIMILAR in the scheme of Gs matching 
   Ts for RNA.
   Call it with a buffer containing integer codes of DNA letters and the size of the tuple.
   Call it "pow(2,N)-1" times where N is the number of Gs and Cs in the string, or it will
   go for the next iteration.
*/
int all_GT_tuple_codes(int *dnastring, int *dnastringOriginal, char *intcomp, int tuple) {

  int j,tp,rccode;

  /* find all matching tuples */
  tp = tuple-1;

  /* find G or T */
  for (; tp>=0; tp--)
      if (dnastringOriginal[tp]>1) break;
  dnastring[tp]-=2; // change to A or C


  /* binary addition simmulation */
  for (; tp>=0; tp--) {
    if ( dnastring[tp] <0 ) {

        dnastring[tp]+=4;  // change back to G or T

        /* don't carry the first caracter */
        if (tp==0) break;

        /* find the next G or T again */
        for (j=tp-1; j>=0; j--)
            if (dnastringOriginal[j]>1) break;

        /* if not found, quit */
        if (dnastringOriginal[j]<=1)
            break;

        /* change tp */
        tp=j+1;

        dnastring[tp-1]-=2; // change to A or C
        continue;  
    }
    break;
  }

  /* debug 
  for (j=0; j<tuple; j++)
      fprintf(stderr,"%d ",dnastring[j]);
  fprintf(stderr,"\n");*/

  /* calculate code */
  rccode=0;
  for (j=tuple-1;j>=0;j--) {
    rccode=rccode*4+Complement[intcomp[dnastring[j]]]; /* important to use it like that to make sure MR works ok */
  }

  return rccode;
}

/*******************************************************************/
/**************************  Init_rc_codes()  **********************/
/*******************************************************************/
/* added on 09/28/04 by Gelfand to speed up RC codes ***************/
void Init_rc_codes(void)

{
  int j,tp,tuple,s,code,rccode,k,simcode,ncodes,u,g;
  int p4t;
  char intcomp[4];
  int dnastring[20]; /* must be bigger then the largest tuple */
  int dnastringcopy[20];

  /* important to use it like that to make sure MR works ok */
  intcomp[0]='A';
  intcomp[1]='C';
  intcomp[2]='G';
  intcomp[3]='T';

  /* set RCcodes  */
  {
    RCcodes = (int *)scalloc(four_to_the[Tuplesize[NTS]],sizeof(int));
    memory_stats_print("\ninit_rc_codes: requesting", four_to_the[Tuplesize[NTS]] * sizeof(int) );
    //if(RCcodes==NULL){debugerror("\ninit_rc_codes: Out of memory!");exit(0);}


    /* Added by Gelfand on 10/14/2004 to search for "similar" tuples 
    that allow GT match      */
        for(g=1;g<=NTS;g++) {
	    /* G. Benson 12/27/22 fixed error sizeof(int) --> sizeof(int *) */
            RCcodesSimilar[g] = (int **)scalloc(four_to_the[Tuplesize[g]],sizeof(int *));
            memory_stats_print("\ninit_rc_codes: requesting", four_to_the[Tuplesize[g]] * sizeof(int) );
            //if(RCcodesSimilar[g]==NULL){debugerror("\ninit_rc_codes: Out of memory!");exit(0);}

#ifdef IRF_GT_MATCH_YES
            /* this part ensures the code will work when GTmatch is not selected are not present */
            if (!paramset.gtmatchyesno) {
#endif

                  /* set the arrays to -1 */
                      for (j=0; j<four_to_the[Tuplesize[g]]; j++) {
                            RCcodesSimilar[g][j] = (int*)scalloc(1,sizeof(int));
                            memory_stats_print("\ninit_rc_codes: requesting", 1 * sizeof(int) );
                            //if(RCcodesSimilar[g][j]==NULL){debugerror("\ninit_rc_codes: Out of memory!");exit(0);}
                            RCcodesSimilar[g][j][0]=-1;
                      }
#ifdef IRF_GT_MATCH_YES
            }
#endif
        }

    /* init integer string */
    tuple = Tuplesize[NTS];
    for (j=0; j<tuple; j++) 
        dnastring[j]=0;
    dnastring[j-1]=-1;    

    /* debug 
    fprintf(stderr,"\n"); */

	/* go through all compositions */
    p4t=four_to_the[tuple];
	for (j=0; j<p4t; j++) {

        tp = tuple-1;
        dnastring[tp]++;

	    /* binary addition simmulation */
	    for (; tp>0; tp--) {
		    if ( dnastring[tp] > 3 ) {dnastring[tp]=0; dnastring[tp-1]++; continue;  }
		    break;
	    }

        /* calculate code and rccode */
        code=0;
        for (s=0;s<tuple;s++) {
            code=code*4+dnastring[s]; 
        }
        rccode=0;
        for (s=tuple-1;s>=0;s--) {
            rccode=rccode*4+Complement[intcomp[dnastring[s]]]; /* important to use it like that to make sure MR works ok */
        }
        RCcodes[code]=rccode;


        /* debug 
        for (s=0;s<tuple;s++)
            fprintf(stderr,"%d ",dnastring[s]);
        fprintf(stderr,"   c: %d   rc: %d\n",code,rccode);
        */
    }

    /* debug 
    fprintf(stderr,"\n\n");
    */
		
  }

#ifdef IRF_GT_MATCH_YES

  /* Added by Gelfand on 10/14/2004 to search for "similar" tuples 
  that allow GT match      */
  if (paramset.gtmatchyesno) {

      /* now do the tables for smaller tuples */
      for(g=NTS;g>=1;g--)
      {

        /* init integer string */
        tuple = Tuplesize[g];
        for (j=0; j<tuple; j++) 
            dnastring[j]=0;
        dnastring[j-1]=-1;   


	    /* go through all compositions */
        p4t=four_to_the[tuple];
	    for (j=0; j<p4t; j++) {

            tp = tuple-1;
            dnastring[tp]++;

	        /* binary addition simmulation */
	        for (; tp>0; tp--) {
		        if ( dnastring[tp] > 3 ) {dnastring[tp]=0; dnastring[tp-1]++; continue;  }
		        break;
	        }

            /* calculate code and rccode */
            code=0;
            for (s=0;s<tuple;s++) {
                code=code*4+dnastring[s]; 
            }
            rccode=0;
            for (s=tuple-1;s>=0;s--) {
                rccode=rccode*4+Complement[intcomp[dnastring[s]]]; /* important to use it like that to make sure MR works ok */
            }


            memcpy(dnastringcopy,dnastring,sizeof(int)*tuple);
    
            /* calculate number of codes for this tuple, depends on the # of Gs and Cs */
            ncodes=0;
            for (k=0; k<tuple; k++) {
                if (dnastringcopy[k]>1) ncodes++;
            }
            ncodes = (int)pow(2,ncodes)-1;

            /* debug 
            fprintf(stderr,"\n\n");
            fprintf(stderr,"   code: %d   rcode: %d\n",code,rccode);
              for (u=0; u<tuple; u++)
                fprintf(stderr,"\t%d",dnastring[u]);
            fprintf(stderr,"\n\n");
            fprintf(stderr,"\tsimilar:\n"); */
            /* end debug */

            /* go through the codes and build the "GT similar" table */
            for (k=0;k<ncodes;k++) {

                simcode=all_GT_tuple_codes(dnastringcopy,dnastring,intcomp,tuple);

                { 

                    if (NULL==RCcodesSimilar[g][rccode]) {
                        RCcodesSimilar[g][rccode] = (int*)scalloc(ncodes+1,sizeof(int));
                        memory_stats_print("\ninit_rc_codes: requesting", (ncodes+1) * sizeof(int) );
                        //if(RCcodesSimilar[g][rccode]==NULL){debugerror("\ninit_rc_codes: Out of memory!");exit(0);}
                    }

                    RCcodesSimilar[g][rccode][k]=simcode;
                }   

                /* debug 
                for (u=0; u<tuple; u++)
                    fprintf(stderr,"\t%d ",dnastringcopy[u]);
                fprintf(stderr,"\n");
                fprintf(stderr,"\t\t%d:\t%d (%d)\n\n",k,simcode,RCcodesSimilar[g][rccode][k]); 
                 */
                /* end debug */
        

            }
            if (NULL!=RCcodesSimilar[g][rccode]) 
                RCcodesSimilar[g][rccode][k]=-1;
            else {
                RCcodesSimilar[g][rccode] = (int*)scalloc(1,sizeof(int));
                memory_stats_print("\ninit_rc_codes: requesting", 1 * sizeof(int) );
                //if(RCcodesSimilar[g][rccode]==NULL){debugerror("\ninit_rc_codes: Out of memory!");exit(0);}
                RCcodesSimilar[g][rccode][0]=-1;
            }


        } /* for (j=0; j<p4t; j++) {*/

      } /* end of for(g=NTS;g>=1;g--) */

  } /* end of if (paramset.gtmatchyesno) */

#endif

}


/*******************************************************************/
/**************************    init_sm()   *************************/
/*******************************************************************/
/* Only to be used with IRF_GT_MATCH_YES defined!!!                */
void init_sm(int match,int gtmatch, int mismatch)

{
  int i,j,*currint;

  /* SM has 256*w56 entries to map MATCH-MISMATCH matrix */
  SM=(int *)scalloc(256*256,sizeof(int));
  memory_stats_print("\n init_index: requesting", 256*256 * sizeof(int) );
//  if(SM==NULL){debugerror("\nInit_sm: Out of memory!");exit(0);}


  /* generate 256x256 into matrix */
  for(i=0,currint=SM;i<=255;i++)
  {
    for(j=0;j<=255;j++,currint++)
    {
        if(i==j) *currint=match;
        else *currint=mismatch;
    }
  }

  /* add the GT match  */
  /* Note that the second character is complemented!!! */
  SM['G'*256+'A']=gtmatch;
  SM['T'*256+'C']=gtmatch;

}

/*******************************************************************/
/**************************  init_complement()  ********************/
/***************** added 02-08-01 reverse complement ***************/
/*******************************************************************/

void init_complement(void)

{
  int i;

  /* complement has 256 entries so that finding the entries for A, C, G and T */
  /* require no calculation */
  Complement=(int *)scalloc(256,sizeof(int));
  memory_stats_print("\n init_complement: requesting", 256 * sizeof(int) );
  //if(Complement==NULL){debugerror("\nInit_complement: Out of memory!");exit(0);}

  for (i=0; i<256; i++)
	Complement[i]=-1;

  /* for mirror repeats */
  if (MRyesno) {
    Complement['A']=0;
    Complement['C']=1;
    Complement['G']=2;
    Complement['T']=3;
  } else {
    Complement['A']=3;
    Complement['C']=2;
    Complement['G']=1;
    Complement['T']=0;
  }

}

/*******************************************************************/
/**************************  Clear_Histories()  ********************/
/***************** added 02-08-29 ****************** ***************/
/*******************************************************************/

void Clear_Histories(void) {
	  
	  int g;

	  for(g=1;g<=NTS;g++)
	  {
		  memset( History[g], 0, sizeof(struct historyentry) * Historysize[g] );
		  memset(Tuplehash[g], 0,  four_to_the[Tuplesize[g]] * sizeof(int));
		  Nextfreehistoryindex[g]=1;  /* set all to 1 because 0 indicates */
		  /* Tuplehash points to nothing */
	  }
	  
}

/*******************************************************************/
/*******************  init_complement_ascii()   ********************/
/********** added 06-27-02 reverse complement ascii ****************/
/*******************************************************************/

void init_complement_ascii(void)

{
  /* complement has 256 entries so that finding the entries for A, C, G and T which are complement ascii values */
  /* require no calculation */
  int i;
  Complementascii=(char *)scalloc(256,sizeof(int));
  
  memory_stats_print("\n init_complement_ascii: requesting", 256 * sizeof(int) );
  //if(Complementascii==NULL){debugerror("\ninit_complement_ascii: Out of memory!");exit(0);}
  
  for (i=0; i<256; i++)
	Complementascii[i]=i;

  /* for mirror repeats */
  if (MRyesno) {
      Complementascii['A']='A';
      Complementascii['C']='C';
      Complementascii['G']='G';
      Complementascii['T']='T';
      Complementascii['a']='a';
      Complementascii['c']='c';
      Complementascii['g']='g';
      Complementascii['t']='t';
  } else {
      Complementascii['A']='T';
      Complementascii['C']='G';
      Complementascii['G']='C';
      Complementascii['T']='A';
      Complementascii['a']='t';
      Complementascii['c']='g';
      Complementascii['g']='c';
      Complementascii['t']='a';
  }

}

/*******************************************************************/
/************* InitTuplesizeTuplemaxdistance() *********************/
/*******************************************************************/

void InitTuplesizeTuplemaxdistance()
{

 
  //debugmessage("\ntuple sizes 0,4,5,7"); 
  Tuplemaxdistance[0]=0;
  Tuplemaxdistance[1]=154;
  Tuplemaxdistance[2]=813;
  Tuplemaxdistance[3]=14800;

  
  /* debugmessage("\ntuple distances 0, 29, 159, MAXDISTANCE");  */
}

/******************************************************************/
/************* InitK_RunsSumsAndRandomWalk()  *********************/
/*******************************************************************/

void InitK_RunsSumsAndRandomWalk() 
{
  int g;
  
  for(g=1;g<=NTS;g++)
  {
    /* K runs sums values */
    
    /* Random Walk values */
    Delta[g]=Intervals[g].delta[Intervals[g].num_intervals]; 
  }
}


/******************************************************************/
/**************************  process_sequence() ********************/
/************************02-08-01 reverse complement ***************/
/*******************************************************************/


void process_sequence(FASTASEQUENCE* pseq)
{
  
  char datstring[_MAX_PATH],txtstring[_MAX_PATH],htmlstring[_MAX_PATH],maskstring[_MAX_PATH],paramstring[_MAX_PATH],gtextra[20]="";
 
  int build_entire_code,g,badcharindex,found,progbarpos,percentincrease,onepercent,lastlocation/*, pfd=0*/,scoreIn,simCode,nsc;
  int mintuplesize,maxtuplesize;
  int code,y,i,h,d,j,z,zz,rccode,p,k,LJ;
  int index;
  int startforsecondalignment;
  int cindex;
  int intervalsize;
  int rnkp;
  int jtemp;
  int success_at_larger_interval_size;
	double frequent_center_ratio;
  char ltest;
  double  im3; /*debug, remove laters*/
  
#ifdef IS_ALIGNMENT_PROFILER_ON
  /* alignment profiler */
	ClearAlignmentProfiler();
	AlignmentProfiler.sequenceStarted=time(NULL);
#endif


  /* generate the parameter string to be used in file names */
  sprintf(paramstring,"%d.%d.%d.%d.%d.%d.%d.%d",
		paramset.match,paramset.mismatch,paramset.indel,
		paramset.PM,paramset.PI,paramset.minscore,paramset.maxlength,paramset.maxloop);
  
  
  /* print the names of the files */
  sprintf(htmlstring,"%s.%s.html",paramset.outputprefix,paramstring);
  sprintf(txtstring ,"%s.%s.txt.html",paramset.outputprefix,paramstring);
  sprintf(datstring ,"%s.%s.dat",paramset.outputprefix,paramstring);
  sprintf(maskstring,"%s.%s.mask",paramset.outputprefix,paramstring);
  
 

  /* open files */
  Fpdat=fopen(datstring,"w");
		
 
if (!paramset.HTMLoff) {

  /* start html file */
  Fptxt=fopen(txtstring,"w");

  fprintf(Fptxt,"<HTML>");
  fprintf(Fptxt,"<HEAD>");
  fprintf(Fptxt,"<TITLE>%s</TITLE>",txtstring);
  fprintf(Fptxt,"</HEAD>");
  fprintf(Fptxt,"<BODY bgcolor=\"#FBF8BC\">");
  fprintf(Fptxt,"<PRE>");
  
  fprintf(Fptxt,"\nInverted Repeats Finder Program written by:");
  fprintf(Fptxt,"\n\n   Gary Benson, Yevgeniy Gelfand");
  fprintf(Fptxt,"\n      Bioinformatics Program");
  fprintf(Fptxt,"\n          Boston University");
  fprintf(Fptxt,"\n\nVersion %s", versionstring);
  fprintf(Fptxt,"\n\nPlease cite:\nP. E. Warburton, J. Giordano, F. Cheung, Y. Gelfand and G. Benson. \nInverted Repeat Structure of the Human Genome: The X-Chromosome \nContains a Preponderance of Large, Highly Homologous Inverted \nRepeats That Contain Testes Genes, \nGenome Research, 14:1861-1869, 2004. 10.1101/gr.2542904.\n");

  if (paramset.gtmatchyesno) 
      sprintf(gtextra," -gt %d", paramset.gtmatch);

  /*
  fprintf(Fptxt,"\n\nSequence: %s\n\nParameters: %d %d %d %d %d %d %d %d%s%s%s%s%s\n\n",
		pseq->name,paramset.match,paramset.mismatch,paramset.indel,
        paramset.PM,paramset.PI,paramset.minscore,paramset.maxlength,paramset.maxloop,paramset.datafile?" -d":"",paramset.maskedfile?" -m":"",paramset.flankingsequence?" -f":"",paramset.lowercase?" -l":"",paramset.mryesno?" -mr":"",gtextra);
  */
  fprintf(Fptxt,"\n\nSequence: %s\nParameters: %s\n",pseq->name,paramset.parameters);

  fprintf(Fptxt,"Length: %d",Length);
  fprintf(Fptxt,"\nACGTcount: A:%3.2f, C:%3.2f, G:%3.2f, T:%3.2f, N:%3.2f\n\n",
		(double)pseq->composition['A'-'A']/Length,
		(double)pseq->composition['C'-'A']/Length,
		(double)pseq->composition['G'-'A']/Length,
		(double)pseq->composition['T'-'A']/Length,
		1.0 - (double)(pseq->composition['A'-'A']+pseq->composition['C'-'A']+pseq->composition['G'-'A']+pseq->composition['T'-'A'])/Length);
}  
	
  /* start dat file */
  fprintf(Fpdat,"Inverted Repeats Finder Program written by:\n\n");
  fprintf(Fpdat,"Gary Benson\n");
  fprintf(Fpdat,"Bioinformatics Program\n");
  fprintf(Fpdat,"Boston University\n");
  fprintf(Fpdat,"Version %s\n", versionstring);
  

  /*
  fprintf(Fpdat,"\n\nSequence: %s\n\nParameters: %d %d %d %d %d %d %d %d%s%s%s%s%s\n\n",
		pseq->name,paramset.match,paramset.mismatch,paramset.indel,
		paramset.PM,paramset.PI,paramset.minscore,paramset.maxlength,paramset.maxloop,paramset.datafile?" -d":"",paramset.maskedfile?" -m":"",paramset.flankingsequence?" -f":"",paramset.lowercase?" -l":"",paramset.mryesno?" -mr":"",gtextra);
  */
  fprintf(Fpdat,"\n\nSequence: %s\n\n\nParameters: %s\n\n",pseq->name,paramset.parameters);
  fflush(Fpdat);
	
  /* clean all the datastructures */
  //ClearCenterlists();
  ClearCenterlists3();
  Clear_MaxMincenters3();
  Clear_Histories();
//  Clear_Tags();
  free_centerseenlist();
  mintuplesize=Tuplesize[1];
  maxtuplesize=Tuplesize[NTS];
	
	


  onepercent = Length/100;
  percentincrease = 0;
  progbarpos = 0;

  /* start processing */
  //debugmessage("\n\nProcess sequence");
  build_entire_code=1;
  for(i=0;i<=Length;i++)   /* Sequence starts at 1, not zero */
  {

    /* if percent changed then set indicator */
    percentincrease++;
    if(percentincrease==onepercent)
    {
        percentincrease = 0;
        progbarpos++;
        paramset.percent = progbarpos;
        fprintf(stderr,"."); fflush(stderr);
    }

	/* adjust the running max and min centers */
	for(k=1;k<=NTS;k++)
		Update_Runningmaxmincenters3(i,k);

    if((i==0)  /* before start of sequence or */
      ||(strchr("ACGT",Sequence[i])==NULL))  /* not one of A,C,G,T */
    {
      badcharindex=i;
      build_entire_code=1;
      /* find first good string of mintupsize characters */
      g=0;
      while((g<mintuplesize)&&(i<Length))
      {
        i++;
	 	/* adjust the running max and min centers */
		for(k=1;k<=NTS;k++)
		Update_Runningmaxmincenters3(i,k);
        if(strchr("ACGT",Sequence[i])==NULL)
        {
          badcharindex=i;
          g=0;
        }
        else g++;
      }
      if(g<mintuplesize) break;  /* i=Length and minimum tuple not found */
    }
    if(build_entire_code)
    {
      code=0;
      for(g=badcharindex+1;g<=i;g++)
      {
        code=code*4+Index[Sequence[g]]; 
				
        /* debugmessage("\nSequence[%d]=%c, Index=%d, code=%d",g,Sequence[g],Index[Sequence[g]],code);*/
      }

#ifdef old_stuff_nut_used_sdfasdfasdfasd      
      /* added 02-08-01 reverse complement code */
      rccode=0;
      for(g=i;g>=badcharindex+1;g--)
      {
        rccode=rccode*4+Complement[Sequence[g]];
        /* debugmessage("\nSequence[%d]=%c, Complement=%d, rccode=%d",g,Sequence[g],Complement[Sequence[g]],rccode); */
      }
      rccode=rccode*four_to_the[maxtuplesize-(i-badcharindex)];  /* must increase length of rccode when not as long as
      maxtuplesize */
#endif

      rccode=RCcodes[code]; /* added by Gelfand on 09/28/04 */


      if(i-badcharindex>=maxtuplesize){build_entire_code=0;}

    }
		
    else
    {
      code=(code%four_to_the[Tuplesize[NTS]-1])*4+Index[Sequence[i]];
      
#ifdef old_stuff_nut_used_sdfasdfasdfasd 
      /* added 02-08-01 reverse complement code */
      rccode=(rccode/4)+(four_to_the[Tuplesize[NTS]-1]*Complement[Sequence[i]]);
#endif

      rccode=RCcodes[code]; /* added by Gelfand on 09/28/04 */


    }
    Tuplecode[NTS]=code;
    /* debugmessage("\nTuplecode[%d] %d\n",NTS,Tuplecode[NTS]);*/
    /* added 02-08-01 reverse complement code */
    Tuplerccode[NTS]=rccode;
    /* debugmessage("\nTuplerccode[%d] %d\n",NTS,Tuplerccode[NTS]); */
    
    for(h=NTS-1;h>=1;h--)
    {
      Tuplecode[h]=code%four_to_the[Tuplesize[h]];
      /* debugmessage("\nTuplecode[%d] %d\n",h,Tuplecode[h]); */
      
      /* added 02-08-01 reverse complement code */
      Tuplerccode[h]=rccode/four_to_the[Tuplesize[NTS]-Tuplesize[h]];
      /* debugmessage("\nTuplerccode[%d] %d\n",h,Tuplerccode[h]);*/
    }
    
    /* test that tuplecodes are being processed correctly */
    /*
	debugmessage("\n\n i=  %d",i);
		//if (i==2000) exit(0);
		
    debugmessage("\n   tupsize index           code                 rccode");
    for(h=1;h<=NTS;h++)
    {
      // debugmessage("\nTuplecode[%d] %d\n",h,Tuplecode[h]);
      debugmessage("\n        %d          ",h);
      for(g=1;g<=Tuplesize[NTS]-Tuplesize[h];g++){ debugmessage("  ");}
      prntcode=Tuplecode[h];
      for(g=1;g<=Tuplesize[h];g++)
      {
        debugmessage("%d ",prntcode/four_to_the[Tuplesize[h]-g]);
        prntcode=prntcode%four_to_the[Tuplesize[h]-g];
      }
      
      // debugmessage("\nTuplerccode[%d] %d\n",h,Tuplerccode[h]); 
      debugmessage("          ");
      for(g=1;g<=Tuplesize[NTS]-Tuplesize[h];g++){ debugmessage("  ");}
      prntcode=Tuplerccode[h];
      for(g=1;g<=Tuplesize[h];g++)
      {
        debugmessage("%d ",prntcode/four_to_the[Tuplesize[h]-g]);
        prntcode=prntcode%four_to_the[Tuplesize[h]-g];
      }
    }
    debugmessage("\n");
    */
    /* process index i using all the tuplesizes */
//    fprintf(stderr,"\ni=%d",i);
    


    g=1;
    ltest = Sequence[i+1]; /* for the "look ahead character test", Gelfand */
    while((g<=NTS)&&(i-badcharindex>=Tuplesize[g])) /* check for valid tuple */
    {
			
      
      h=Nextfreehistoryindex[g];      /* next free index in history list */
      j=h+1;                          /* advance next free index         */
      if(j>=Historysize[g]) j=1;      /* we use a circular history list  */
      
			/*  if((History[g][j].location!=0)   if the next entry has already been used  */
      /*  &&(j==Tuplehash[g][History[g][j].code]))  check Tuplehash. If it still */
			/* Tuplehash[g][History[g][j].code]=0;   points here, zero it out.     */
      
			
			
			
			
      if(History[g][j].location!=0)                /* if the next entry has already been used  */
			{
				if(j==Tuplehash[g][History[g][j].code])    /* check Tuplehash. If it still */
					Tuplehash[g][History[g][j].code]=0;         /* points here, zero it out.     */
				else                                        /* otherwise, scan through list and find */
				{                                                                                                                                                             
					/* entry that points to j and make it point */
					jtemp=Tuplehash[g][History[g][j].code];   /*to zero */
					while((jtemp!=j)&&(jtemp!=0))             /* exits loop either after pointer to entry j */
					{                                                                                                                                                           
						/* is cleared or finds zero pointer */
						if(History[g][jtemp].previous==j)
						{
							History[g][jtemp].previous=0;
							jtemp=0;
						}
						else jtemp=History[g][jtemp].previous;
					}
					
				}  
				
			}	  
			
  	  Nextfreehistoryindex[g]=j;      /* store next free index */

																			
			y=Tuplehash[g][Tuplecode[g]]; /* index in history list of last
      occurrence of code */
			
         


      /* move BEFORE the z extraction. Gelfand, 11/10/04 */
      Tuplehash[g][Tuplecode[g]]=h;   /* store index of current tuple   */
      History[g][h].location=i;       /* store info about current tuple */
      History[g][h].previous=y;
      History[g][h].code=Tuplecode[g];
            
      /* added 02-08-01                                     */
      /* for inverted repeats, z is index in history        */
      /* of last occurrence of reverse complement of code   */
      //z=Tuplehash[g][Tuplerccode[g]]; /* index in history list of last occurrence of reverse complement code */
	  

   /* added on 10/14/04 by Gelfand to for GTmatch tuple checking */
   simCode = Tuplerccode[g];
   for (nsc=-1; simCode!=-1; simCode=RCcodesSimilar[g][Tuplerccode[g]][nsc]) {

       //debug
     /*     if (nsc==-1) fprintf(stderr,"\nRegular! %d ",simCode);
          else fprintf(stderr,"\nSimcode! %d",simCode);
              fflush(stderr);
      */

//       {
//           fprintf(stderr,"g=%d, Tuplesize[g]=%d, simCode=%d, four_to_the[Tuplesize[g]=%d\n",g, Tuplesize[g], simCode, four_to_the[Tuplesize[g]]);
//           exit(1);
//       }

      z=Tuplehash[g][simCode];




      /* following modified 02-08-01        */
      /* includes while loop                */
      /* reverse complement                 */
      zz=0;        /* zz holds entry which points to z */
      /* if Tuplehash points to z, zz=0   */
      /* debugmessage("\nz:%3d   zz:%3d",z,zz); */      
      while(z!=0)
      {


        d=i-History[g][z].location; /* d=distance between matching tuples */
                                    /** note that d must be at least as large as one tuple so that the tuples are
                                    not overlapping
        **/
        lastlocation = i - d - Tuplesize[g]; /* need for LookAhead and GTMatch detection, Gelfand 11/18/2004 */
     	
		//fprintf(stderr,"i: %d     y: %d     z: %d     d: %d   g:  %d    Tuplecode[g]: %d    Tuplerccode[g]:  %d    History[g][j].location: %d   History[g][j].code:  %d\n",i,y,z,d,g, Tuplecode[g],Tuplerccode[g],History[g][j].location,History[g][j].code);
			
        /*fprintf(stderr,"\n    Matching Tuples:  tuplesize index:%d  last reverse complement index:%d  distance:%d"
          ,g,History[g][z].location,d); */

        if(d>Tuplemaxdistance[g])   /* if d exceeds Tuplemaxdistance, then  */
        {                           /* make the previous location 0.        */
          /* We are no longer interested in the z */
          //debugmessage("\n     !distance between tuples exceeds Tuplemaxdistance[%d]: %d!",g,Tuplemaxdistance[g]);
          if(zz==0) {
              
                /* Gelfand, 11/10/04 */
                /* NOTE: it appears code and rccode can be the same. But it does not seem to be a problem */

                Tuplehash[g][simCode]=0; /* if no previous z, reset tuplehash         */

          } else History[g][zz].previous=0;           /* else zero history pointer at zz.          */
          /* Entry at z will be zeroed out when reused */ 
          z=0;                      /* set z to zero to stop while loop */
        } 
        else if (d<Tuplesize[g])  /* if d is too small, tuple and reverse complement overlap */
        {
          //debugmessage("\n     !tuple and reverse complement overlap!");
          zz=z;
          z=History[g][z].previous; /* get next matching tuple */
        }
        else /* process tuple */
        {


		  

          zz=z;
          z=History[g][z].previous; /* get next matching tuple */   /*******stopped here 02-08-01 *******/
					
					//if (g==1) fprintf(stderr,"\nzz=%d  z=%d  d=%d  locationOfZ: %d \n",zz,z,d,History[g][z].location);
          /* process */
          

          /* calculate correct Center index */
          index=2*i-d-Tuplesize[g]+1; /* index is true center times 2 to avoid 1/2 values */
          
          /* calculate entry in circular list for center index */
          cindex=index%Centerlistsize[g];
          
          //debugmessage("\n     Center[%d][%d].end_of_list:%d",g,cindex,Center[g][cindex].end_of_list);
          /* add tuple to Center[g][cindex] */
          /*debugmessage("\n      add tuple match i:%d  true centerindex:%d  cindex:%d  Tuplesize:%d  tuplesize index:%d",
            i, index, cindex, Tuplesize[g], g);
			*/		
		              
           
              //if ((Sequence[i]=='G' && ((simCode/four_to_the[Tuplesize[g]-1])==0)) || (Sequence[i]=='T' && ((simCode/four_to_the[Tuplesize[g]-1])==1)))							
              //if ((Sequence[i]=='G' && Sequence[rcloc=(History[g][zz].location-Tuplesize[g]+1)]=='T') || (Sequence[i]=='T' && Sequence[rcloc]=='G')) {							
              //    
              //    add_tuple_match_to_Center3(i,Tuplesize[g],index,g,cindex, rcloc, GTDetectMatch ); 
              //} else */



//          add_tuple_match_to_Center(i,Tuplesize[g],index,g,cindex,z,y);
          //add_tuple_match_to_Center2(i,Tuplesize[g],index,g,cindex);

          add_tuple_match_to_Center3(i,Tuplesize[g],index,g,cindex, Gtmatchyesno ? (lastlocation + 1) : 0 );


					
									
          /*debugmessage("\n         after add_tuple_match Center[%d][%d].end_of_list: %d",
            g,cindex,Center[g][cindex].end_of_list);        */
					
          /* print for testing purposes */
          //PrintCenterMatchlist(g,cindex);


          /* Check if the next character also matches, if yes, we will pick it up on the next try, leave for now */
          /* Changed in v3.02 */
          /* Two changes are made to this code, Gelfand, July 14, 2004 */
          /* a. Forgot that complemented sequence was completely uppercased, therefore I cannot use it to predict 
                value of the next tuple. Using first sequence now.
             b. Now checking if Tuplemaxdistance[g] is not exceeded, otherwise the tuple will nevel be tested.
          */
          if (paramset.la && lastlocation>=1 && (ltest == Complementascii[Sequence[lastlocation]]) && Index[ltest]!=-1 && (Tuplemaxdistance[g] > (d+Tuplesize[g]))) {
              //Center3[g][cindex].waitingToAlign = 1;
              continue;
          }

									
          /* using intervals to test the criteria after a match is found */
          h=Intervals[g].num_intervals;
          success_at_larger_interval_size=0;  /* prevents redundant alignment if IR already found with larger interval */
          while((h>=1)&&(!success_at_larger_interval_size))
          {
              
              intervalsize=Intervals[g].size[h];
              
     
              /* debugmessage("\n           Tuplesize index: %d  Interval index: %d  Interval size: %d",
              g,h,intervalsize);*/
              
              rnkp=Intervals[g].RNKP[h];
              frequent_center_ratio=Intervals[g].frequent_center_ratio[h];
              /* restrict testing to those matches which are far enough from */
              /* the center to have enough matches to meet the RNKP criteria */
              /* suggestion, at least 75% of the RNKP criteria               */
              /* this assumes that the RNKP criteria is satisfied by one long */
              /* run of matches                                               */
              /* if(i-(index/2)>=rnkp*.75)*/ /* index is center*2 *//* must redo this
              to get interval correct */
              
              /* see if interval size is too large, that is, the position in the sequence 
              is too close to the center to meet the rnkp criterion for this interval size */
              
              if((2*i-index>=rnkp*2)
                  &&((im3=intervalmatches3(i,g,cindex,h))>=
                  frequent_center_ratio*rnkp)) /* test changed by Gary Benson 8/2003 
                                               to reflect lowest ratio of most frequent center tuple matches to RNKP value 
                                               in 95% of simulation trials */
                                               /*.3*rnkp))*/ /* test added by Gary Benson 12/05/02 */
                                               /* requires a minimum number of matches in center */ 
            { 
//							if(TestIntervalCriteria3(g,index,intervalsize,h,i)) {
//							if(TestIntervalCriteria4(g,index,h,i)) {


// test to find something not being passed by TestIntervalCriteria34
/*
if ((i==120354) && (g==2) && (h==1) && (index==240197)) {
        
                                                // we must run TestIntervalCriteria34(g,index,h,i) anyway
                                                // to fill the values correctly and move the lo and hi interval pointers 
                                                TestIntervalCriteria34(g,index,h,i);

                                                // to view the center info, Gelfand, 07/30.2004 
                                                fprintf(stderr,"\n:i = %d", i);
                                                fprintf(stderr,"\n:index = %d", index);
                                                fprintf(stderr,"\n tupleindex: %d",TestMaxIntervalSum.tupleindex);
                                                fprintf(stderr,"\n intervalindex: %d",TestMaxIntervalSum.intervalindex);
                                                fprintf(stderr,"\n matches_in_interval: %d",TestMaxIntervalSum.matches_in_interval);
                                                fprintf(stderr,"\n lowcenter: %d",TestMaxIntervalSum.lowcenter);
                                                fprintf(stderr,"\n highcenter: %d",TestMaxIntervalSum.highcenter);
                                                fprintf(stderr,"\n testcenter: %d",TestMaxIntervalSum.testcenter);
                                                fprintf(stderr,"\n im3: %d",im3);
                                                fprintf(stderr,"\n frequent_center_ratio*rnkp: %.4lf\n",frequent_center_ratio*rnkp);
                                                
}
*/  

							if(TestIntervalCriteria34(g,index,h,i)) {
								


                                
                            

								
								Classlength=i-index/2;
								
                                /*
                                if (Center3[g][cindex].waitingToAlign) {
                                    
                                    Center3[g][cindex].waitingToAlign = 0;
                                    found = 0;
                                } else
                                */

									/* check if this center has already been processed */
									found=search_for_center_match_in_centerseenlist(index,i,g,h);
								

								if (!found) {
									


#ifdef IS_ALIGNMENT_PROFILER_ON
	/* alignment profiler */
	AlignmentProfiler.cellsTemp=0;
#endif

									/* try to align in starting from outside */
									Maxscore=0;
									ECstart= index - i;
									EC=SequenceComplement+ECstart;
									narrowbandwrap( index, i, max(MINBANDRADIUS,Delta[g]), RECENTERCRITERION, ALIGNIN, Alpha*5,g,h,i);
									
									
									//debugmessage("\nindex: %d  i: %d  Classlength: %d\n\n", index, i, Classlength);
									//exit(0);
									
									/* This test is redundant because the narrowbandwrap function above should not have even
									been called if the were so little matches. But let's just make sure we have at least 
									5 matches, else output an error. This should never happen.
									*/
									//if (Maxscore>=Alpha*5) {
									if (Maxscore>=Reportmin) {
										
										
										
										/* trace back and get the align pair */
										get_narrowband_pair_alignment(max(MINBANDRADIUS,Delta[g]));
										
										/*{
											// print the pair again 
											int charcounter;
											debugmessage("\nlength: %d\n\n", AlignPair.length);
											
											debugmessage("left sequence: %s\n", &AlignPair.textprime[1]);
											debugmessage("left, first index: %d\n", AlignPair.indexprime[1]);
											debugmessage("left, second index: %d\n\n", AlignPair.indexprime[AlignPair.length]);
											
											debugmessage("right sequence: ");
											for (charcounter=1; charcounter<=AlignPair.length; charcounter++)
												debugmessage("%c",Complementascii[AlignPair.textsecnd[AlignPair.length-charcounter+1]]);
											
											debugmessage("\nright, first index: %d\n", ECstart+AlignPair.indexsecnd[AlignPair.length]);
											debugmessage("right, second index: %d\n\n", ECstart+AlignPair.indexsecnd[1]);
											//exit(0); 
										}*/
										
										/* try to align out starting from the inside */
										Maxscore=0;
                                        Largestjump = 0;
                                        scoreIn=AlignPair.score;
										// startforsecondalignment=ECstart+AlignPair.length-1;
										startforsecondalignment=ECstart+AlignPair.indexsecnd[1];
										Classlength=max(startforsecondalignment,Length-AlignPair.indexprime[1]+1);
										
										/*debugmessage("\n ???: %d ECstart(old): %d AlignPair.length: %d  startforsecondalignment: %d  Classlength: %d\n\n", ECstart+AlignPair.indexsecnd[1],ECstart, AlignPair.length,startforsecondalignment, Classlength);
										{
											int charcounter;
											debugmessage("right sequence: ");
											for (charcounter=1; charcounter<=AlignPair.length; charcounter++)
												debugmessage("%c",Complementascii[AlignPair.textsecnd[AlignPair.length-charcounter+1]]);
										}
										debugmessage("Length: %d  AlignPair.indexprime[1]: %d  Classlength: %d\n",Length,AlignPair.indexprime[1],Classlength); 
										*/
										//exit(0);
										
										ECstart=AlignPair.indexprime[1];
										EC=SequenceComplement+ECstart;
										narrowbandwrap(index, startforsecondalignment, max(MINBANDRADIUS_OUTER,Delta[g]),RECENTERCRITERION, ALIGNOUT, Reportmin,g,h,i);
										


										/* check if the score is high enought to report the repeat */

										if (Maxscore>=Reportmin) {
											

											/* trace back and get the align pair */
											get_narrowband_pair_alignment(max(MINBANDRADIUS_OUTER,Delta[g]));



                                            /* print the ones that have not reached the starting point, should not happen after the fix
                                               to the alignment routing in v 3.02 */
                                            /*if ((ECstart+AlignPair.indexsecnd[1])<i) {
                                                int _op, c1, c2, c3;
                                                get_3mod_centers(&c1,&c2,&c3);
                                                fprintf(stderr,"\nPremature finish detected at i=%d! (diff: %d) (intervalsize: %d) (oldscore: %d   newscore: %d) (center found: %d   most common centers: %d,%d,%d)",i,i-(ECstart+AlignPair.indexsecnd[1]),intervalsize,scoreIn,AlignPair.score, index, c1, c2, c3); fflush(stderr);
                                                pfd++;
                                                fprintf(stderr,"leftarm: ");
                                                for (_op=1; _op<=AlignPair.length; _op++) {
                                                    fprintf(stderr,"%c",AlignPair.textprime[_op]);
                                                }
                                                fprintf(stderr," rightflank: ");
                                                for (_op=1; _op<=100; _op++) {
                                                    fprintf(stderr,"%c",Sequence[AlignPair.indexprime[AlignPair.length]+_op]);
                                                }
                                            }*/
                                            
#ifdef IRF_THIRD_ALIGNMENT                     
                                        if (paramset.a3) {                                               
                                            /* Third alignment. Because some repeats seem to have incorrect (or unfinished) endpoints near the loop */
                                            /* Changed in v3.03 */
                                            /* Gelfand,  Nov 23, 2004 */  
                                            Maxscore=0;
                                            Largestjump=0;
                                            startforsecondalignment = ECstart+AlignPair.indexsecnd[1]; //ECstart+AlignPair.indexsecnd[1];
                                           	ECstart= AlignPair.indexprime[1];
									        EC=SequenceComplement+ECstart; 
                                            if (paramset.a3==1&&Largestjump!=0)
                                                LJ=min(MINBANDRADIUS_OUTER,Largestjump+5); //This speeds up the program by about 10-20 %
                                            else 
                                                LJ=MINBANDRADIUS_OUTER;
                                            narrowbandwrap(index, startforsecondalignment, LJ,RECENTERCRITERION, ALIGNIN, Reportmin,g,h,i);
                                            get_narrowband_pair_alignment(LJ);
                                            //fix coordinates after the third pass to look like we only did 3 passes
                                            reverse();
                                            ECstart = AlignPair.indexprime[1];
                                        }
#endif



                                            /* check if in the list, do not output if it is */
									        found=search_for_center_match_in_centerseenlist_exact(index,
                                                                                                  ECstart+AlignPair.indexsecnd[AlignPair.length], /* this used to be wrong before, Gelfand Feb 06, 2006*/
                                                                                                  ECstart+AlignPair.indexsecnd[1], g, h);
                                            if (!found) {




#ifdef IS_ALIGNMENT_PROFILER_ON
	/* alignment profiler */
	AlignmentProfiler.repeatsBeforeRedundancy[g].interval[h]+=1;
#endif

											    /* change the letters of the second align string to uppercase */
											    for (p=1; p<=AlignPair.length; p++)
												    if (AlignPair.textsecnd[p]>='a'&&AlignPair.textsecnd[p]<='z') AlignPair.textsecnd[p]+=('A'-'a');

											    for (p=1; p<=AlignPair.length; p++)
												    if (AlignPair.textprime[p]>='a'&&AlignPair.textprime[p]<='z') AlignPair.textprime[p]+=('A'-'a');


											    /* output repeat */
											    AlignPair.textprime[AlignPair.length+1]='\0';
											    AlignPair.textsecnd[AlignPair.length+1]='\0';
											    Separation=ECstart+AlignPair.indexsecnd[AlignPair.length]-AlignPair.indexprime[AlignPair.length]-1;
                                                
                                                if (!paramset.HTMLoff) 
											        fprintf(Fptxt,"\n\n\nFound at i: %d    center: %d    interval size: %d    g: %d",i,index,intervalsize,g);

                                                /* to view the center info, Gelfand, 07/30.2004 */
                                                /*fprintf(Fptxt," tupleindex: %d",TestMaxIntervalSum.tupleindex);
                                                fprintf(Fptxt," intervalindex: %d",TestMaxIntervalSum.intervalindex);
                                                fprintf(Fptxt," matches_in_interval: %d",TestMaxIntervalSum.matches_in_interval);
                                                fprintf(Fptxt," lowcenter: %d",TestMaxIntervalSum.lowcenter);
                                                fprintf(Fptxt," highcenter: %d",TestMaxIntervalSum.highcenter);
                                                fprintf(Fptxt," testcenter: %d",TestMaxIntervalSum.testcenter); */



											    /* next line added by Gary Benson, 12/04/02 for testing */
											    
											   

											    get_statistics();
											    
                                            if (!paramset.HTMLoff) {
                                                // to view the center info, Gelfand, 07/30.2004 
                                                fprintf(Fptxt,"\n:i = %d", i);
                                                fprintf(Fptxt,"\n tupleindex: %d",TestMaxIntervalSum.tupleindex);
                                                fprintf(Fptxt,"\n intervalindex: %d",TestMaxIntervalSum.intervalindex);
                                                fprintf(Fptxt,"\n matches_in_interval: %.4lf",TestMaxIntervalSum.matches_in_interval);
                                                fprintf(Fptxt,"\n lowcenter: %d",TestMaxIntervalSum.lowcenter);
                                                fprintf(Fptxt,"\n highcenter: %d",TestMaxIntervalSum.highcenter);
                                                fprintf(Fptxt,"\n testcenter: %d",TestMaxIntervalSum.testcenter);
                                                fprintf(Fptxt,"\n im3: %.4lf",im3);
                                                fprintf(Fptxt,"\n frequent_center_ratio*rnkp: %.4lf\n",frequent_center_ratio*rnkp);
                                            }
            
                                                

											    //PrintMaxIntervalSum(i);

											    
											    /* print the pair again */
											    /*debugmessage("\nlength: %d\n\n", AlignPair.length);
											    
											    debugmessage("left sequence: %s\n", &AlignPair.textprime[1]);
											    debugmessage("left, first index: %d\n", AlignPair.indexprime[1]);
											    debugmessage("left, second index: %d\n\n", AlignPair.indexprime[AlignPair.length]);
											    
											    debugmessage("right sequence: ");
											    for (charcounter=1; charcounter<=AlignPair.length; charcounter++)
												    debugmessage("%c",Complementascii[AlignPair.textsecnd[AlignPair.length-charcounter+1]]);
											    
											    debugmessage("\nright, first index: %d\n", ECstart+AlignPair.indexsecnd[AlignPair.length]);
											    debugmessage("right, second index: %d\n\n", ECstart+AlignPair.indexsecnd[1]);
											    //exit(0);
											    */

											    /* add to centerseen */
											    add_to_centerseenlist(index,ECstart+AlignPair.indexsecnd[AlignPair.length],ECstart+AlignPair.indexsecnd[1],AlignPair.score,g,h);
										    

                                            }

#ifdef IS_ALIGNMENT_PROFILER_ON
	/* alignment profiler */
	AlignmentProfiler.alignmentSuccessfull[g].interval[h]+=2;
	AlignmentProfiler.cellsSuccessfull[g].interval[h]+=AlignmentProfiler.cellsTemp;
#endif

											/* note that you found something at this interval size */
											success_at_larger_interval_size=1;
										}
										
										
										
									} 
									/*
									else {
										debugmessage("ERROR: could not match the minimum number of matches. TestCriteria function is probobly incorrect.");
									} */
									
									
									
					} /* end of found */
					
					
					//exit(0);
					
					
					
					
					
				} /* end of TestIntervalCriteria3 */
            }


           h--; 
          }  /* end of loop for different interval sizes */
					
			
          
        } /* end of process tuple */

      } /* end of z loop */

     nsc++;
     } /* end of the loop for "similar GT" tuple matches */



      g++;
    } /* end of g loop */
  }
	

 /*  {
  struct centerlistelement *temp;

  for (temp=Centerseenlist; temp!=NULL; temp=temp->next)
	  fprintf(stderr,"\nindex: %d      rightstart: %d     rightend: %d",temp->index,temp->rightstart,temp->rightend);

	} */
 

#ifdef IS_ALIGNMENT_PROFILER_ON
  /* print the alignment profiler */
  {
	char profilerstring[_MAX_PATH];
	FILE * profile;
    sprintf(profilerstring ,"%s.%s.profiler.xls",paramset.outputprefix,paramstring);
	profile=fopen(profilerstring,"w");
	if (profile==NULL) {
		fprintf(stderr,"Cannot open the profile file for writing. Aborting!");
		exit(-1);
	}
	AlignmentProfiler.sequenceEnded=time(NULL);
	fprintf(profile,"%s\n\n",pseq->name);
	PrintAlignmentProfiler(profile);
	fclose(profile);
  }
#endif 

  /* debug */
  //printf("\ntotal pfd: %d",pfd);

if (!paramset.HTMLoff) {
  /* close the dat file */
  fprintf(Fptxt,"Done.\n");
  fprintf(Fptxt,"</PRE>"); 
  fprintf(Fptxt,"</BODY>");
  fprintf(Fptxt,"</HTML>");
  fclose(Fptxt);
}
  fclose(Fpdat); 

  
  /* this function define on trfclean.h */
  fprintf(stderr,"\nResolving output...");
  fflush(stderr);
  IRFClean( datstring, htmlstring, txtstring,paramset.datafile,paramset.maxlength,paramset.maxloop,paramset.maskedfile,maskstring);
	
}

/***********************************************
*   This routine can act on a multiple-sequence
*   file and calls process_sequence() routine as 
*	many times as it needs to. First, it allocates
*   all program structures and arrays.
************************************************/

int IRFControlRoutine(void) {

  FASTASEQUENCE seq;
  char  line[1000];
  int *stemp;
  FILE *srcfp,*outmfp,*outdfp,*destmfp,*destdfp;
  char  source[_MAX_PATH],prefix[_MAX_PATH],destm[_MAX_PATH],destd[_MAX_PATH],
        desth[_MAX_PATH],paramstring[_MAX_PATH],outh[_MAX_PATH],input[_MAX_PATH],outm[_MAX_PATH],outd[_MAX_PATH];
  int a,i,rc,g,foundsome=0,bufferSize=0;
  FILE  *desthfp;



  /* set algorithm's parameters */
  Alpha = paramset.match;  
  Beta  = -paramset.mismatch;
  DeltaIndel = -paramset.indel;
  if (Beta>0) Beta=-Beta;
  if (DeltaIndel>0) DeltaIndel=-DeltaIndel;
  if (paramset.gtmatchyesno) {
      Gtmatchyesno = 1;
      GTAlpha = paramset.gtmatch;
  } else {
      Gtmatchyesno = 0;
      GTAlpha = Beta;
  }
  if (paramset.mryesno) {
      MRyesno = 1;
  } else {
      MRyesno = 0;
  }
  PM    = paramset.PM;
  PI    = paramset.PI;
  Pindel=(float)PI/100;
  Pmatch=(float)PM/100;
  Reportmin = paramset.minscore;
  GTDetectMatch = 1.0;
#ifdef IRF_GT_MATCH_YES
  if (Gtmatchyesno) {
    //GTDetectMatch = ((double)GTAlpha) / ((double)Alpha);
    GTDetectMatch = ((double)(GTAlpha - Beta))  /  ((double)(Alpha - Beta));
  }
  //GTDetectMatch = .001;
#endif

  /* init some stuff */
  init_complement_ascii();
  init_index();
  init_complement();


  /***********************************************************/
  /******************* load the sequence  ********************/
  /***********************************************************/
  
  /* save names locally so they can be replaced later */
  strcpy(source,paramset.inputfilename);
  strcpy(prefix,paramset.outputprefix);
  
  
  /* open input file for reading */
  srcfp = fopen(source,"rb");
  if(srcfp==NULL)
  {
          sprintf(line,"%s not found",source);
          fprintf(stderr, "%s", line);
          paramset.endstatus = CTRL_BADFNAME;
          paramset.running = 0;
          return CTRL_BADFNAME;
  }


  /* get file size */
  fseek( srcfp, 0, SEEK_END );
  bufferSize = ftell( srcfp );
  fseek( srcfp, 0, SEEK_SET );


  /* to save memory, let's make sure our maximum stem loop is smaller or equal to the fize size (assuming maximum will occur when only one sequence in a file, the whole sequence is one stemploop split in the middle and are lot of insertions on both sides) */
  MAXWRAPLENGTH = min( MAXWRAPLENGTH, bufferSize );


  fprintf(stderr,"\nLoading sequence...");
  fflush(stderr);
  rc=LoadSequenceFromFile(&seq,srcfp);
  if(rc==-1)
  {
    fprintf(stderr,"\nThe file is not in FASTA format!\n\n");
    paramset.endstatus = CTRL_BADFORMAT; /* ok for now */
    exit(CTRL_BADFORMAT);
  }
  if(rc==-2)
  {
    fprintf(stderr,"\nCould not allocate memory for sequence!\n\n");
    paramset.endstatus = CTRL_BADFORMAT; /* ok for now */
    exit(CTRL_BADFORMAT);
  }
  
  /* set the sequence pointer. more global vars! */ 
  Sequence = seq.sequence-1;					/* start one character before */
  SequenceComplement=seq.sequencecomplement-1;	/* start one character before */
  Length=    seq.length;


  //debug
  /*Sequence[Length+1]=SequenceComplement[Length+1]=0;
  fprintf(stderr,"buffer1: %s\n",Sequence+1);
  fprintf(stderr,"buffer2: %s\n",SequenceComplement+1);
  exit(0);
	*/

  fprintf(stderr,"\nAllocating Memory...");
  fflush(stderr);

  /***********************************************************/
  /***** allocate memory for the scoring matrix  *************/
  /***********************************************************/

  S = (int **) scalloc(MAXWRAPLENGTH+1,sizeof(int *));
  memory_stats_print("\n IRFControlRoutine(S): requesting", (MAXWRAPLENGTH+1) * sizeof(int *) );
/*  if(S==NULL)
  {
      debugerror("Unable to allocate memory for S array");
      exit(0);
  }*/
  stemp = (int *) scalloc((MAXWRAPLENGTH+1)*(2*MAXBANDWIDTH+2),sizeof(int));
  memory_stats_print("\n IRFControlRoutine(stemp): requesting", (MAXWRAPLENGTH+1)*(2*MAXBANDWIDTH+2)*sizeof(int) );
  /*if(stemp==NULL)
  {
      debugerror("Unable to allocate memory for stemp array");
      exit(0);
  }*/
  for(i=0;i<=MAXWRAPLENGTH;i++)
  {
      S[i]= stemp;
      stemp+=2*MAXBANDWIDTH+2;
  }
  S[0][0]=1;


  /***********************************************************/
  /** allocate memory for bandcenter (v3.14 change, Gelfand) */
  /***********************************************************/
  Bandcenter = (int *) scalloc(MAXWRAPLENGTH+1,sizeof(int));
  memory_stats_print("\n IRFControlRoutine(stemp): requesting", (MAXWRAPLENGTH+1)*sizeof(int) );

  /* AlignPair holds the characters and alignments of the current */
  /* primary and secondary sequences  */
  AlignPair.textprime=newAlignPairtext(2*MAXWRAPLENGTH);
  AlignPair.textsecnd=newAlignPairtext(2*MAXWRAPLENGTH);
  AlignPair.indexprime=newAlignPairindex(2*MAXWRAPLENGTH);
  AlignPair.indexsecnd=newAlignPairindex(2*MAXWRAPLENGTH);

  /* new in 3.02 */
#ifdef IRF_GT_MATCH_YES
  init_sm(Alpha,GTAlpha,Beta); 
#endif

  fprintf(stderr,"\nInitializing data structures...");
  fflush(stderr);


  /* v3.02 08/18/04, I had to bring this outside, because I had to move the 
  InitIntervals() function before the InitTuplesizeTuplemaxdistance() function */

  /* Waiting time calculations */
  /* if(PM==80) */
  Tuplesize[0]=0;
  Tuplesize[1]=4;
  Tuplesize[2]=5;
  Tuplesize[3]=7;

   /***********************************************************/
  /*************** initialize intervals info  ****************/
  /***********************************************************/
  
  InitIntervals(); 

  /***********************************************************/
  /*** initialize tuplesize and tuplemaxdistance info ********/
  /* tuplemaxdistance is how far back a match at tuple size  */
  /* can be tested                                           */
  /***********************************************************/
  
  InitTuplesizeTuplemaxdistance();


  /* changed on 01/28/04 to find repeats with separation that is close to loopsize */
  
  /* overwrite with advanced parameters ,v3.02 08/18/04 */
  if (paramset.t4!=-1) Tuplemaxdistance[1]=paramset.t4+2*Intervals[1].size[Intervals[1].num_intervals];
  if (paramset.t5!=-1) Tuplemaxdistance[2]=paramset.t5+2*Intervals[2].size[Intervals[2].num_intervals];
  if (paramset.t7!=-1) Tuplemaxdistance[3] = paramset.t7 + 2*MAXINTERVALSIZE; 

  /* debug */
  fprintf(stderr,"\n\nTuples"
    "\n  tupsize index   tuplesize   tuplemaxdistance");
  for(g=1;g<=NTS;g++)
  {
    fprintf(stderr,"\n     %d              %d               %d",g,Tuplesize[g],Tuplemaxdistance[g]);
  }

  /***********************************************************/
  /********** initialize tuplehash and history ***************/
  /* tuplehash holds possible tuple codes                    */
  /* history holds pointers to previous occurrence of same   */
  /* code                                                    */
  /***********************************************************/
  
  
  for(g=1;g<=NTS;g++)
  {
    Tuplehash[g]=(int *)scalloc(four_to_the[Tuplesize[g]],sizeof(int));
    Historysize[g]=2*(Tuplemaxdistance[g]+1)+2; /* The idea here is that no */
    /* previous history pointer points back */
    /* more than Tuplemaxdistance.  Then, when */
    /* History entry is reused, following */
    /* links from the current will exceed the */
    /* maxdistance before reaching the reused */
    /* entry. */
    History[g]=(struct historyentry *)
      scalloc(Historysize[g], sizeof(struct historyentry));
    Nextfreehistoryindex[g]=1;  /* set all to 1 because 0 indicates */
    /* Tuplehash points to nothing */

    memory_stats_print("\n IRFControlRoutine(Tuplehash[g]): requesting", four_to_the[Tuplesize[g]]* sizeof(int) );
    //if( NULL == Tuplehash[g]) {debugerror("\nUnable to allocate memory to tuplehash[%d]. Aborting!", g ); exit(0); }

    memory_stats_print("\n IRFControlRoutine(History[g]): requesting", Historysize[g] * sizeof(struct historyentry) );
    //if( NULL == History[g]) {debugerror("\nUnable to allocate memory to History[%d]. Aborting!", g ); exit(0); }

  }


  
  /***********************************************************/
  /*************** initialize centerlist info ****************/
  /***********************************************************/
  
  Init_MaxMincenters3();
  //InitCenterlists();
  InitCenterlists3();
  /* Init_links(); */ 
  //Init_Tags();
  Init_centerseenlist();  
  Init_rc_codes();


  /***********************************************************/
  /*** initialize k run sums and random walk values   ********/
  /* k run sums is the number of matches that indicate       */
  /* a repeat.                                               */
  /* random walk indicates drift of center due to indels     */
  /***********************************************************/
  InitK_RunsSumsAndRandomWalk();



  //debugmessage("\nSequence: %s\n", seq.sequence); exit(0);
  //debugmessage("SequenceComplement: %s\n", SequenceComplement);
  //debugmessage("\nLength: %d\n",Length);
  
  
  
  
  
  /* based on number of sequences in file use different approach */
  if(rc==0) /* only one sequence in file */
  {
	  paramset.multisequencefile = 0;
	  paramset.sequenceordinal = 1;
	  /* call trf and return */
	  
	  /* progress */
	  if(paramset.multisequencefile)
	  {
		  sprintf(line,"\nScanning Sequence %d...",
			  paramset.sequenceordinal);
		  fprintf(stderr,"%s", line);
	  }
	  else
	  {
		  fprintf(stderr,"\nScanning...");
	  }
	  fflush(stderr);
	  
	  process_sequence(&seq);
	  sfree(seq.sequence);
      sfree(seq.sequencecomplement);
	  fclose(srcfp);
	  paramset.running = 0;
      i=1;


      /* print the file to stdout */
      if(paramset.ngs && paramset.datafile)
        {

	       fprintf(stdout,"@%s\n",seq.name);

               /* generate the parameter string to be used in file names */
               sprintf(paramstring,"%d.%d.%d.%d.%d.%d.%d.%d",
                  paramset.match,paramset.mismatch,paramset.indel,
                  paramset.PM,paramset.PI,paramset.minscore,paramset.maxlength,paramset.maxloop);

                sprintf(outd,"%s.%s.dat",prefix,paramstring);
                outdfp = fopen(outd,"r");
                if (NULL==outdfp) {fprintf(stderr,"\nUnable to open data file %s. Aborting!",outd); exit(1); }
                while(1)
                 {
                         a = getc(outdfp);
                         if(a==EOF) break;
                         putc(a,stdout);
                }
                fclose(outdfp);
                
                /* remove intermediary file */
                remove(outd);
         }


	paramset.endstatus = CTRL_SUCCESS;
  }
  else if (rc>0)/* multiple sequences in file */
  {


	  paramset.multisequencefile = 1;
	  paramset.sequenceordinal = 1;
	  
	  
	  /*********************************************************
	  *   if there are more files need to produce sumary-style
	  *   output.
	  **********************************************************/
	  
	  /* generate the parameter string to be used in file names */
	  sprintf(paramstring,"%d.%d.%d.%d.%d.%d.%d.%d",
		  paramset.match,paramset.mismatch,paramset.indel,
		  paramset.PM,paramset.PI,paramset.minscore,paramset.maxlength,paramset.maxloop);
	  
	  /* open sumary table file */
if (!paramset.HTMLoff) {
	  sprintf(desth,"%s.%s.summary.html",prefix,paramstring);
	  desthfp = fopen(desth,"w");
             if(desthfp==NULL)
                 {
                          fprintf(stderr,"Unable to open html file for writing in IRFControlRoutine routine!");
                          exit(-1);
                 }

}	  
	  /* open masked file if requested */
	  if(paramset.maskedfile)
	  {
		  sprintf(destm,"%s.%s.mask",prefix,paramstring);
		  destmfp = fopen(destm,"w");
                  if(destmfp==NULL)
                      {
                              fprintf(stderr,"Unable to open mask file for writing in IRFControlRoutine routine!");
                              exit(-1);
                      }
	  }
	  /* open datafile if requested */
	  if(paramset.datafile)
	  {

                   if (paramset.ngs) {
                           destdfp = stdout;
                   } else {
                           sprintf(destd,"%s.%s.dat",prefix,paramstring);
                           destdfp = fopen(destd,"w");
                           if(destdfp==NULL)
                           {
                                   fprintf(stderr,"Unable to open data file for writing in IRFControlRoutine routine!");
                                   exit(-1);
                           }
                   }

	  }

if (!paramset.HTMLoff)  {
	  /* start output of sumary file */
	  fprintf(desthfp,"<HTML>");
	  fprintf(desthfp,"<HEAD>");
	  fprintf(desthfp,"<TITLE>Output Summary</TITLE>");
	  /* fprintf(desthfp,"<BASE TARGET=\"Table\">"); */
	  fprintf(desthfp,"</HEAD>");
	  fprintf(desthfp,"<BODY bgcolor=\"#FBF8BC\">");
	  fprintf(desthfp,"<PRE>");
	  
	  fprintf(desthfp,"\nInverted Repeats Finder Program written by:<CENTER>");
	  fprintf(desthfp,"\nGary Benson");
	  fprintf(desthfp,"\nBioinformatics Program");
	  fprintf(desthfp,"\nBoston University");
	  fprintf(desthfp,"\nVersion %s</CENTER>\n", versionstring);
  	  fprintf(desthfp,"\n\nPlease cite:\nP. E. Warburton, J. Giordano, F. Cheung, Y. Gelfand and G. Benson. \n\"Inverted Repeat Structure of the Human Genome: The X-Chromosome \nContains a Preponderance of Large, Highly Homologous Inverted \nRepeats That Contain Testes Genes\", \nGenome Research, 14:1861-1869, 2004. 10.1101/gr.2542904.\n");
	  
	  
	  fprintf(desthfp,"\n\n<B>Multiple Sequence Summary</B>\n\n");
	  fprintf(desthfp,"Only sequences containing repeats are shown!\n\n");
	  fprintf(desthfp,"Click on sequence description to view repeat table.\n\n");
	  
	  /* print beginning of table */
	  fprintf(desthfp,"<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n");
	  fprintf(desthfp,"<TR><TD WIDTH=80><CENTER>Sequence\nIndex</CENTER></TD>"
		  "<TD WIDTH=400><CENTER>Sequence\nDescription</CENTER></TD>"
		  "<TD WIDTH=80><CENTER>Number of\nRepeats</CENTER></TD>"
		  "</TR>\n");
}	  
	  
		  /******************************************
		  *   process every sequence in file
	  *******************************************/
	  i=1;
	  for(;;)
	  {

	  //debugmessage("\nSequence: %s\n", seq.sequence); 
	  //debugmessage("\nLength: %d\n",Length);
	  //getch();

		  /* set the prefix to be used for naming of output */
		  sprintf(input,"%s.s%d",prefix,i);
		  strcpy(paramset.inputfilename,input);
		  strcpy(paramset.outputprefix,input);
		  
		  /* progress */
		  if(paramset.multisequencefile)
		  {
			sprintf(line,"\nScanning Sequence %d...",
					paramset.sequenceordinal);
			fprintf(stderr, "%s", line);
		  }
		   else
		  {
			fprintf(stderr,"\nScanning...");
		  }
		 fflush(stderr);

		  /* call the inverted repeats finder routine */
		  process_sequence(&seq);
		  
		  /* append new output to destination files */
		  if(paramset.maskedfile)
		  {
			  /* recreate the name of the masked sequence file */
			  sprintf(outm,"%s.s%d.%s.mask",prefix,i,paramstring);
			  outmfp = fopen(outm,"r");
	                  if (NULL==outmfp) {fprintf(stderr,"\nUnable to open masked file %s. Aborting!",outd); exit(1); }

			  /* copy until end of file */
			  while(1)
			  {
				  a = getc(outmfp);
				  if(a==EOF) break;
				  putc(a,destmfp);
			  }
			  fclose(outmfp);
			  
			  /* remove intermediary file */
			  remove(outm);
		  }
		  if(paramset.datafile)
		  {

			  if (paramset.ngs) { fprintf(destdfp,"@%s\n",seq.name); }

			  sprintf(outd,"%s.s%d.%s.dat",prefix,i,paramstring);
			  outdfp = fopen(outd,"r");
	                  if (NULL==outdfp) {fprintf(stderr,"\nUnable to open data file %s. Aborting!",outd); exit(1); }
        
        // Change by G. Benson 2.14.25
        // added &&(!paramset.ngs)
        // to fix bug with -ngs output for multiple sequences in the same input file. The table entries for the first
        // sequence were reported correctly, but for subsequent sequences, the first 6 entries (or all if there were
        // or fewer) were excluded.  
			  if((i!=1)&&(!paramset.ngs)) /* discard header if not on first sequence, but only if also not ngs */
			  {
				  for(a=0;a<6;a++) fgets(line,256,outdfp);
			  }
			  while(1)
			  {
				  a = getc(outdfp);
				  if(a==EOF) break;
				  putc(a,destdfp);
			  }
			  fclose(outdfp);
			  
			  /* remove intermediary file */
			  remove(outd);
		  }
		  
		  /* print table rows based on repeat count */
		  sprintf(outh,"%s.%s.1.html",input,paramstring);

		  if(paramset.outputcount>0)
		  {
if (!paramset.HTMLoff) {
			  /* print a table raw to the summary table */
			  fprintf(desthfp,"<TR><TD><CENTER>%d</CENTER></TD>"
				  "<TD><CENTER><A TARGET=\"%s\" HREF=\"%s\">%s</A>"
				  "</CENTER></TD><TD><CENTER>%d</CENTER></TD></TR>",
				  i,outh,outh,seq.name,paramset.outputcount);
}
			  foundsome=1;
		  }
		  else
		  {
			  /* remove html files if no output in it */
			  remove(outh);
			  sprintf(line,"%s.%s.1.txt.html",input,paramstring);
			  remove(line);
		  }
		  
		  /* free the data associated with the sequence */
		  sfree(seq.sequence);
		  sfree(seq.sequencecomplement);
		  
		  /* if more sequences load and repeat */
		  if(rc>0)
		  {
			  rc = LoadSequenceFromFile(&seq,srcfp);
			  paramset.sequenceordinal++;
			  i++;
			    
			  /* set the sequence pointer. more global vars! */ 
			  Sequence = seq.sequence-1;					/* start one character before */
			  SequenceComplement=seq.sequencecomplement-1;	/* start one character before */
			  Length=    seq.length;

		  }
		  else
		  {

			  if(rc==-1)
			  {
				  fprintf(stderr,"\nThe file is not in FASTA format!\n\n");
				  exit(0);
			  }
			  if(rc==-2)
			  {
				  fprintf(stderr,"\nCould not allocate memory for sequence!\n\n");
				  exit(0);
			  }

			  break;
		  }
    }

    /* close table and html body */
if (!paramset.HTMLoff) {
    fprintf(desthfp,"\n</TABLE>\n");
    if(!foundsome)
    {
        fprintf(desthfp, "\nNo Repeats Found!<BR>");
    }
    fprintf(desthfp,"\n</BODY></HTML>\n");
}

    /* close files */
	fclose(srcfp);
    if(paramset.maskedfile) fclose(destmfp);
    if(paramset.datafile && !paramset.ngs) fclose(destdfp);

if (!paramset.HTMLoff)
    fclose(desthfp);

    /* set output file name to the summary table */
    strcpy(paramset.outputfilename,desth);

  }
  
  
  
  /* cleanup */
  fprintf(stderr,"\nFreeing Memory...");
  fflush(stderr);

  // There is some memory corruption somewhere in the program that crashes the program during memory release 
  // so we comment for now
  /*
  sfree(S[0]); 
  sfree(S);
  sfree(Bandcenter);

  sfree(AlignPair.textprime);
  sfree(AlignPair.textsecnd);
  sfree(AlignPair.indexprime);
  sfree(AlignPair.indexsecnd);
  free_centerseenlist();
  sfree(Maxcenter3);
  sfree(Mincenter3);
  sfree(Runningmaxcenter3);
  sfree(Runningmincenter3);
  sfree(Oldrunningmaxcenter3);
  sfree(Oldrunningmincenter3);
  ClearCenterlists3();
  for(g=1;g<=NTS;g++) { 
    sfree(Center3[g]);
    sfree(Tuplehash[g]);
    sfree(History[g]);
    sfree(RCcodes);
  }
  sfree(Index);
  sfree(SM);
  sfree(Complement);
  sfree(Complementascii);
*/
  /* finish message */
  fprintf(stderr,"\nDone\n\n");

  /* debug */
#ifdef EASY_LIFE_DEBUG
  sMemorySummary(1);
#endif

  paramset.endstatus = CTRL_SUCCESS;
  return i;

}
/************************************	PrintBanner	*****************************************/
void PrintBanner(void)
{
    fprintf(stderr,"\nInverted Repeats Finder, Version %s", versionstring);
    fprintf(stderr,"\nCopyright (C) Dr. Gary Benson 2002-2023. All rights reserved.\n");

    return;
}

/************************************	GetNamePartAddress	*********************************/
char* GetNamePartAddress(char* name)
{
    int i;
    char *pname;
    #if defined(UNIXCONSOLE)||defined(UNIXGUI)
    char dirsymbol = '/';
    #elif defined(WINDOWSCONSOLE)||defined(WINDOWSGUI)
    char dirsymbol = '\\';
    #endif


    i = strlen(name)-1;
    pname = &name[i];
    while(i>0&&*pname!=dirsymbol)
    {
        pname--;
        i--;
    }
    if(*pname==dirsymbol) pname++;
    return pname;
}

#define COMMAND_PARSE(a) \
        if(!stricmp(av[a],"-t4")) {\
            if (ac>(a+1)) {\
              paramset.t4 = atoi(av[a+1]);\
              if (paramset.t4==0&&av[a+1][0]!='0') ac=100; {/* to force error message */\
                if (paramset.t4<0) paramset.t4=0;\
                /*if (paramset.t4>10000) paramset.t4=10000;*/\
              }\
            } else ac=100; /* to force error message */\
        }\
        else\
        if(!stricmp(av[a],"-t5")) {\
            if (ac>(a+1)) {\
              paramset.t5 = atoi(av[a+1]);\
              if (paramset.t5==0&&av[a+1][0]!='0') ac=100; {/* to force error message */\
                if (paramset.t5<0) paramset.t5=0;\
                /*if (paramset.t5>10000) paramset.t5=10000;*/\
              }\
            } else ac=100; /* to force error message */\
        }\
        else\
        if(!stricmp(av[a],"-t7")) {\
            if (ac>(a+1)) {\
              paramset.t7 = atoi(av[a+1]);\
              if (paramset.t7==0&&av[a+1][0]!='0') ac=100; {/* to force error message */\
                if (paramset.t7<0) paramset.t7=0;\
                /*if (paramset.t7>500000) paramset.t7=500000;*/\
              }\
            } else ac=100; /* to force error message */\
        }\
        else\
        if(!stricmp(av[a],"-i1")) paramset.intcheck=1;\
        else\
        if(!stricmp(av[a],"-i2")) paramset.intcheck=2;\
        else\
        if(!stricmp(av[a],"-a3")) paramset.a3=1;\
        else\
        if(!stricmp(av[a],"-a4")) paramset.a3=2;\
        else\
        if(!stricmp(av[a],"-r0")) paramset.redalg=0;\
        else\
        if(!stricmp(av[a],"-r2")) paramset.redalg=2;\
        else\
        if(!stricmp(av[a],"-h")) paramset.HTMLoff=1;\
        else\
        if(!stricmp(av[a],"-m")) paramset.maskedfile=1;\
        else\
        if(!stricmp(av[a],"-f")) paramset.flankingsequence=1;\
        else\
        if(!stricmp(av[a],"-d")) paramset.datafile=1;\
        else\
        if(!stricmp(av[a],"-l")) paramset.lowercase=1;\
        else\
        if(!stricmp(av[a],"-la")) paramset.la=1;\
        else\
        if(!stricmp(av[a],"-mr")) paramset.mryesno=1;\
        else\
	if(!strcmp(av[a],"-ngs") || !strcmp(av[a],"-Ngs") || !strcmp(av[a],"-NGS") || !strcmp(av[a],"-N") || !strcmp(av[a],"-n")) paramset.ngs=1;\
        else\
        if(!stricmp(av[a],"-gt")) {\
            paramset.gtmatchyesno=1;\
            if (ac>(a+1)) {\
              paramset.gtmatch = atoi(av[a+1]);\
              if (paramset.gtmatch==0&&av[a+1][0]!='0') ac=100; /* to force error message */\
            } else ac=100; /* to force error message */\
        }\
        else\
        if(!stricmp(av[a],"-r")) {\
            if (ac>(a+1)) {\
              paramset.redident = (double)atoi(av[a+1])/100.0;\
                if (paramset.redident<.6) paramset.redident=.6;\
                if (paramset.redident>1.0) paramset.redident=1.0;\
            } else ac=100; /* to force error message */\
        }\
        else \
        if (atoi(av[a])!=0||av[a][0]=='0') /* OK, number field*/\
            {} /* OK */\
        else\
            ac=100; /* to force error message */
        

/************************************	main	*********************************************/
int main(int ac, char** av)

{
    int i;
    char *pname;

    PrintBanner();

    /* NOTE: before release, make sure all user input is validated!!! That is not the case now! */

    paramset.datafile = 0;
    paramset.maskedfile = 0;
    paramset.flankingsequence = 0;
	paramset.lowercase = 0;
    paramset.gtmatchyesno = 0;
    paramset.mryesno = 0;
    paramset.flankinglength = 500;
    paramset.t4=-1;
    paramset.t5=-1;
    paramset.t7=-1;
    paramset.HTMLoff = 0;
    paramset.redident = .9;
    paramset.redalg = 1;
    paramset.a3 = 0;
    paramset.la = 0;
    paramset.intcheck = 0;
    paramset.parameters[0]='\0';
    paramset.ngs=0;

    if (ac>=11 && ac <=99)
        for (i=10; i<ac; i++) {
            COMMAND_PARSE(i);
            /* forced from inside the macro */
            if (ac==100) break;
        }


#if (defined(UNIXGUI)+defined(UNIXCONSOLE))<1
        paramset.ngs = 0; /* this is for unix systems only */
#endif


    if  (paramset.ngs == 1) {
                paramset.datafile=1;
    }

    if (ac<10||ac>99)
    {
        fprintf(stderr,"\n\nPlease use: %s File Match Mismatch Delta PM PI Minscore Maxlength MaxLoop [options]\n", av[0]);
        fprintf(stderr,"\nWhere: (all weights, penalties, and scores are positive)");
        fprintf(stderr,"\n  File = sequences input file");
        fprintf(stderr,"\n  Match  = matching weight"); 
        fprintf(stderr,"\n  Mismatch  = mismatching penalty"); 
        fprintf(stderr,"\n  Delta = indel penalty");
        fprintf(stderr,"\n  PM = match probability (whole number)");
        fprintf(stderr,"\n  PI = indel probability (whole number)");
        fprintf(stderr,"\n  Minscore = minimum alignment score to report");
        fprintf(stderr,"\n  MaxLength = maximum stem length to report (10,000 minimum and no upper limit, but system will run out memory if this is too large)");
        fprintf(stderr,"\n  MaxLoop = filters results to have loop less than this value (will not give you more results unless you increase -t4,-t4,-t7 as well)");
        fprintf(stderr,"\n  [options] = one or more of the following :");
        fprintf(stderr,"\n               -m    masked sequence file");
        fprintf(stderr,"\n               -f    flanking sequence");
        fprintf(stderr,"\n               -d    data file");
        fprintf(stderr,"\n               -h    suppress HTML output\n");

	fprintf(stderr,"\n               -l    lowercase letters do not participate in a k-tuple match, but can be part of an alignment");
	fprintf(stderr,"\n               -gt   allow the GT match (gt matching weight must follow immediately after the switch)");
        fprintf(stderr,"\n               -mr   target is mirror repeats");
        fprintf(stderr,"\n               -r    set the identity value of the redundancy algorithm (value 60 to 100 must follow immediately after the switch)\n");

        fprintf(stderr,"\n               -la   lookahead test enabled. Results are slightly different as a repeat might be found at a different interval. Faster.");
        fprintf(stderr,"\n               -a3   perform a third alignment going inward. Produces longer or better alignments. Slower.");
        fprintf(stderr,"\n               -a4   same as a3 but alignment is of maximum narrowband width. Slightly better results than a3. Much slower.");    
        fprintf(stderr,"\n               -i1   Do not stop once a repeat is found at a certain interval and try larger intervals at nearby centers. Better(?) results. Slower.");
        fprintf(stderr,"\n               -i2   Do not stop once a repeat is found at a certain interval and try all intervals at same and nearby centers. Better(?) results. Much slower.");
        fprintf(stderr,"\n               -r0   do not eliminate redundancy from the output");
        fprintf(stderr,"\n               -r2   modified redundancy algorithm, does not remove stuff which is redundant to redundant. Slower and not good for TA repeat regions, would not leave the largest, but a whole bunch.\n");

        fprintf(stderr,"\n               -t4   set the maximum loop separation for tuple of length4 (default 154, separation <=1,000 must follow)");
        fprintf(stderr,"\n               -t5   set the maximum loop separation for tuple of length5 (default 813, separation <=10,000 must follow)");
        fprintf(stderr,"\n               -t7   set the maximum loop separation for tuple of length7 (default 14800, limited by your system's memory, make sure you increase maxloop to the same value)");
#if (defined(UNIXGUI)+defined(UNIXCONSOLE))>=1
        fprintf(stderr,"\n               -ngs  more compact .dat output on multisequence files, returns 0 on success. ");
#endif
        fprintf(stderr,"\n");
        fprintf(stderr,"\nNote the sequence file should be in FASTA format:");
        fprintf(stderr,"\n");
        fprintf(stderr,"\n>Name of sequence");
        fprintf(stderr,"\n   aggaaacctg ccatggcctc ctggtgagct gtcctcatcc actgctcgct gcctctccag");
        fprintf(stderr,"\n   atactctgac ccatggatcc cctgggtgca gccaagccac aatggccatg gcgccgctgt");
        fprintf(stderr,"\n   actcccaccc gccccaccct cctgatcctg ctatggacat ggcctttcca catccctgtg");
        fprintf(stderr,"\n");

        exit(0);
    }


    /* store the parameter string */
    for (i=2; i<ac; i++) {
        strcat(paramset.parameters,av[i]); 
        strcat(paramset.parameters," ");
    }
    if (strlen(paramset.parameters)>_MAX_PATH) {
        fprintf(stderr,"\n\nInput string too large. This should not happen.");
        exit(0);
    }


    /* using these two options together is: 1. confusing  2. does not display correct output */
    /* disallow for now */
    if (paramset.gtmatchyesno && paramset.mryesno) {
        fprintf(stderr,"\n\nSorry, you cannot use -mr and -gt options together.");
        exit(0);
    }

    /* get other input parameters */
    strcpy(paramset.inputfilename,av[1]);



    pname = GetNamePartAddress(av[1]);
    strcpy(paramset.outputprefix,pname);
    paramset.match    = atoi(av[2]);
    paramset.mismatch = atoi(av[3]);
    paramset.indel    = atoi(av[4]);
    paramset.PM = atoi(av[5]);
    paramset.PI = atoi(av[6]);
    paramset.minscore = atoi(av[7]);
    paramset.maxlength = atoi(av[8]);
    paramset.maxloop = atoi(av[9]);
    paramset.guihandle=0;


    /* set the matrix length*/
    if (paramset.maxlength<10000) paramset.maxlength=10000;
    MAXWRAPLENGTH = paramset.maxlength;
    /* Note1: we do another check inside the control routine to make sure MAXWRAPLENGTH 
       is not lorger then the sequence length, might save us memory */
    /* Note2: there might be multiple sequence so we'll use the file size instead */

    /* loop must be larger than 10,000 */
    //paramset.maxloop = max(10000,paramset.maxloop);
    //paramset.maxloop = min(100000,paramset.maxloop);

    /* call the fuction that controls execution */
    if (paramset.ngs) {

            int rc = IRFControlRoutine();

            if (rc>=1)
                    return 0;
            if (rc==0)
                    return CTRL_NOTHINGPROCESSED;
            else
                    return rc;
    } else {
            return IRFControlRoutine();
    }


}


