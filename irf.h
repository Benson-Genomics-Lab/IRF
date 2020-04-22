#ifndef _MAX_PATH
#define _MAX_PATH 260
#endif

/* uncomment only one platform target identifier */
//#define WINDOWSGUI
//#define WINDOWSCONSOLE
//#define UNIXGUI
#define UNIXCONSOLE

#ifndef WINDOWSCONSOLE
int stricmp( char *str1, char* str2) {
    int result;
    char *sdup1,*sdup2,*src;

    sdup1 = strdup(str1);
    sdup2 = strdup(str2);

    for (src=sdup1; *src!='\0'; src++) {
        *src = toupper(*src);
    }
    for (src=sdup2; *src!='\0'; src++) {
        *src = toupper(*src);
    }

    result = strcmp(sdup1,sdup2);

    free(sdup1);
    free(sdup2);

    return result;
}

#define max(a,b) (((a)>=(b))?(a):(b))
#define min(a,b) (((a)<=(b))?(a):(b))

#endif

/* uncomment this one to force MEGA IRF (larger and wider ) alignment matrix and larger MINBANDRADIUS_OUTER  */
#define IRF_MEGA

/* uncomment the line below to add the profiler .xls file in the output */
//#define IS_ALIGNMENT_PROFILER_ON

/* uncomment the line below to allow GT matching, u still need to call the switch */
#define IRF_GT_MATCH_YES

/* uncomment the line below to perform a third alignment (may get better endpoints), u still need to call the switch */
#define IRF_THIRD_ALIGNMENT

#define MAXTUPLESIZES 10
#define MAXNUMINTERVALS 10
#define MAXBANDWIDTH 220    /* must be equal to at least 2*MINBANDRADIUS_OUTER */ 

int MAXWRAPLENGTH = 500000;
//#ifdef IRF_MEGA         
//#define MAXWRAPLENGTH 500000
//#else
//#define MAXWRAPLENGTH 100000
//#endif




/* version */
#define versionstring "3.07"

/* DNA sequence */
char *Sequence; 

/* a COMPLEMENTED copy of Sequence */
char *SequenceComplement; 

int Length;
FILE *Fptxt,*Fpin;
FILE *Fpdat;

#define CTRL_SUCCESS    0
#define CTRL_BADFNAME   -1
#define CTRL_BADFORMAT  -2
#define CTRL_NOTHINGPROCESSED  -3

//temp debug global variable
int THE_CENTER=0;

int four_to_the[]={1,4,16,64,256,1024,4096,16384,65536,262144,1048576};

/* #define Number_tuple_sizes  4 */
int NTS=3;                   /* number of different tuple sizes to use;
                                   preset for all distances , (assigned here since ver 3.01 )*/
/* int Tuplesize[NTS+1]={0,2,3,5,7};*/ /* what the different sizes are */
/* int Tuplesize[NTS+1]={0,4,5,6,7};*/ /* what the different sizes are */
/* int Tuplesize[NTS+1]={0,3,4,5,7};*/

int Tuplesize[MAXTUPLESIZES+1];
int Tuplemaxdistance[MAXTUPLESIZES+1];
/* int Tuplemaxdistance[MAXTUPLESIZES+1]={0,30,80,200,MAXDISTANCE};*/ /* upper distance for each tuplesize */
/* int Tuplemaxdistance[MAXTUPLESIZES+1]={0,29,83,159,MAXDISTANCE};*/
/* int Tuplemaxdistance[MAXTUPLESIZES+1]={0,29,159,MAXDISTANCE};*/


int Tuplecode[MAXTUPLESIZES+1];       /* this is where the actual tuple codes,
                                              from different size tuples,
                                              encountered at a sequence location
                                              are stored. */
int Tuplerccode[MAXTUPLESIZES+1];     /* this is where the actual reverse complement 
                                              tuple codes, from different size tuples,
                                              encountered at a sequence location
                                              are stored. */
int *Tuplehash[MAXTUPLESIZES+1];      /* points to last location of code
                                              in history list */
int Historysize[MAXTUPLESIZES+1];     /* size of history lists */

int Nextfreehistoryindex[MAXTUPLESIZES+1]; /*next free location in history index*/

int *RCcodes = NULL;                   /* points to lookup table for RC codes (padded for largest tuple), added on 09/28/04 by Gelfand */
int **RCcodesSimilar[MAXTUPLESIZES+1]; /* points to lookup tables of similar codes for each code for each tuplesize, added on 10/14/04 by Gelfand */

struct historyentry{
  int location,
  previous,
  code;
} *History[MAXTUPLESIZES+1];

int *Maxcenter, *Mincenter, *Runningmaxcenter, *Runningmincenter;
int *Oldrunningmincenter, *Oldrunningmaxcenter;
int *Maxcenter3, *Mincenter3, *Runningmaxcenter3, *Runningmincenter3;
int *Oldrunningmincenter3, *Oldrunningmaxcenter3;



int *ACGT_History[4];
int Nextfree_ACGT_History[4];

int *Index, *Complement, *SM = NULL;
char *Complementascii;

int Delta[MAXTUPLESIZES+1];  /* stores delta values for matching centers to test */




int Alpha;							/* match bonus */
int Gtmatchyesno;                   /* GT match */
int GTAlpha;					    /* GT match bonus (v3.02, for GT match, 0 no match)*/
int Beta;							/* mismatch penalty */
int DeltaIndel;						/* indel penalty */  
int MRyesno;                        /* mirror repeat on */
double  GTDetectMatch;              /* how much to give to GT tuple during detection */

/* note Delta global var is already used for something else */

double	PI;
double  PM;
double  Pindel;						/* expected probability of a single character indel in 
									the worst case inverted repeat. Pindel should be tied
									to the indel cost parameter */
double  Pmatch;




char Line[300];

typedef struct
{
    int match;
    int gtmatch;
    int mismatch;
    int indel;
    int minscore;
    int maxlength;
    int maxloop;
    int PM;
    int PI;
    int datafile;
    int maskedfile;
    int flankingsequence;
    int flankinglength;
    int lowercase;
    int gtmatchyesno;
    int mryesno;
    int redalg;
    int t4;
    int t5;
    int t7;
    int HTMLoff;
    int ngs;
    int endstatus;
    int intcheck;
    int a3;
    int la;
    double redident;
    char inputfilename[_MAX_PATH+1]; /* constant defined in stdlib */
    char outputprefix[_MAX_PATH+1];
    char outputdirectory[_MAX_PATH+1];
    char outputfilename[_MAX_PATH+1];
    char parameters[_MAX_PATH+1];
    int  multisequencefile; /* flags if file has more than one sequence */
    int  sequenceordinal; /* holds seq. index starting on 1 */
    int  outputcount; /* repeats found */
    int  guihandle; /* this variable is only used in the GUI version */
    int  running;
	int  percent; /* for progress report */
} IRFPARAMSET;


/* the following structure is used to pass a sequence to the algorithm */
#define MAXSEQNAMELEN 200

typedef struct
{
    int length;
    int composition[26];
    int nucleotides;
    char name[MAXSEQNAMELEN];
    char* sequence;
	char* sequencecomplement;

} FASTASEQUENCE;

/* to remember the centers which were already seen */
struct centerlistelement{
  int index;  // index is the average sum of all the matching letters (again, to prevent 1/2 it is twice the actual center)
  int rightstart;
  int rightend;
  int delta;
  char g; // declared as char to save space, Gelfand 11/15/04
  char h; // declared as char to save space, Gelfand 11/15/04
  //int searchedThrough;
  struct centerlistelement *next;
}Centerseenlist[1];

IRFPARAMSET paramset; /* this global controls the algorithm */

void debugerror(const char *format, ... );
//void debugmessage(const char *format, ... );
