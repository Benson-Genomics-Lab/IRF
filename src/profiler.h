#ifndef _PROFILER_H
#define _PROFILER_H

#ifdef IS_ALIGNMENT_PROFILER_ON

#ifdef IS_ALIGNMENT_PROFILER_ON
__int64 Centerseenlist_Size = 0;
#endif

#include "time.h"


typedef struct apinterval{
  __int64 interval[MAXNUMINTERVALS];                  
} APINTERVAL;

struct alignmentProfiler
{
    APINTERVAL alignmentTotal[MAXTUPLESIZES];
	APINTERVAL alignmentSuccessfull[MAXTUPLESIZES];
	APINTERVAL cellsTotal[MAXTUPLESIZES];
	APINTERVAL cellsSuccessfull[MAXTUPLESIZES];
	APINTERVAL repeatsBeforeRedundancy[MAXTUPLESIZES];
	__int64 cellsTemp;
	time_t sequenceStarted;
	time_t sequenceEnded;

} AlignmentProfiler;


/*******************************************************/
/* ClearAlignmentProfiler clearts the alignment prof   */
/*******************************************************/

void ClearAlignmentProfiler(void)
{
  int g,h;

  for(g=1;g<=NTS;g++)
  {

    for(h=1;h<=Intervals[g].num_intervals;h++)
    {

		AlignmentProfiler.alignmentTotal[g].interval[h]=0;
		AlignmentProfiler.alignmentSuccessfull[g].interval[h]=0;
		AlignmentProfiler.cellsTotal[g].interval[h]=0;
		AlignmentProfiler.cellsSuccessfull[g].interval[h]=0;
		AlignmentProfiler.repeatsBeforeRedundancy[g].interval[h]=0;
    }
  }

	
  AlignmentProfiler.cellsTemp=0;


}

/*******************************************************/
/* PrintAlignmentProfiler prints the alignment prof    */
/*******************************************************/

void PrintAlignmentProfiler(FILE* out)
{
	double  aPercent, cPercent;
	__int64 alignmentTotalSum=0,alignmentSuccessfullSum=0,cellsTotalSum=0,cellsSuccessfullSum=0,repeatsBeforeRedundancy=0;
	int g,h,gstart;


	fprintf(out,"\n\n/********************************************/\n");
	fprintf(out,"/***********  ALIGNMENT PROFILER  ***********/\n");
	fprintf(out,"********************************************/");

	fprintf(out,"\n\nTime taken:   %d seconds",AlignmentProfiler.sequenceEnded-AlignmentProfiler.sequenceStarted);
	fprintf(out,"\n\n\t\t\t\t\t\tALIGNMENTS\t\t\t\t\t\t\tCELLS\n");
	fprintf(out,"Tuple Size\t\tInterval\t\tCalled\t\tSuccess\t\t%%\t\tCalled\t\tSuccess\t\t%%\t\trepeatsBeforeRedundancy\n\n");

	for(g=1;g<=NTS;g++)
	{
		gstart=0;
		for(h=1;h<=Intervals[g].num_intervals;h++)
		{
			aPercent=AlignmentProfiler.alignmentTotal[g].interval[h]?((double)AlignmentProfiler.alignmentSuccessfull[g].interval[h]/(double)AlignmentProfiler.alignmentTotal[g].interval[h]*100):0;
			cPercent=AlignmentProfiler.cellsTotal[g].interval[h]?((double)AlignmentProfiler.cellsSuccessfull[g].interval[h]/(double)AlignmentProfiler.cellsTotal[g].interval[h]*100):0;
			
			if (gstart==0)
				fprintf(out,"%d",Tuplesize[g]);
			else
				fprintf(out,"");

			fprintf(out,"\t\t%d\t\t%I64d\t\t%I64d\t\t%.7lf%%\t\t%I64d\t\t%I64d\t\t%.7lf%%\t\t%I64d\n",
				Intervals[g].size[h],
				AlignmentProfiler.alignmentTotal[g].interval[h],
				AlignmentProfiler.alignmentSuccessfull[g].interval[h],
				aPercent,
				AlignmentProfiler.cellsTotal[g].interval[h],
				AlignmentProfiler.cellsSuccessfull[g].interval[h],
				cPercent,
                AlignmentProfiler.repeatsBeforeRedundancy[g].interval[h]);

			alignmentTotalSum+=AlignmentProfiler.alignmentTotal[g].interval[h];
			alignmentSuccessfullSum+=AlignmentProfiler.alignmentSuccessfull[g].interval[h];
			cellsTotalSum+=AlignmentProfiler.cellsTotal[g].interval[h];
			cellsSuccessfullSum+=AlignmentProfiler.cellsSuccessfull[g].interval[h];
            repeatsBeforeRedundancy+=AlignmentProfiler.repeatsBeforeRedundancy[g].interval[h];

			gstart++;
		}
	}

	fprintf(out,"\n\nTOTALS:\t\t\t\t%I64d\t\t%I64d\t\t\t\t%I64d\t\t%I64d\n",alignmentTotalSum,alignmentSuccessfullSum,cellsTotalSum,cellsSuccessfullSum);

    fprintf(out,"\nrepeatsBeforeRedundancy:\t\t%I64d",repeatsBeforeRedundancy);
    fprintf(out,"\nCenterseenlist_Size:\t\t%I64d",Centerseenlist_Size);

}


#endif //IS_ALIGNMENT_PROFILER_ON
#endif //_PROFILER_H