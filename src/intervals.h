#ifndef _INTERVALS_H
#define _INTERVALS_H

/* Intervals stores the number and sizes of intervals for each tuple size */
/* an interval is the sliding window which totals matches with a common */
/* center for testing with the criteria */

/*

Interval testing

1. Each tuple size has a {\bf Tuplemaxdistance} which is the largest
tested separation between matching tuples (in the case of inverted
repeats, matching reverse complement tuples).  {\bf Tuplemaxdistance} is
essentially the largest allowed separation between parts of an inverted
repeat and suggests that the largest inverted repeat allowed would be of
size Tuplemaxdistance/2.

2. Since we don't know the size of the central hairpin in the stem loop
inverted repeat structure, I want to permit rather large loops, say up
to 1000 characters.

3. As an example, let Tuplemaxdistance = 5000 and the center loop can be
up to 1000.  So there are possibly up to 2000 characters on each side
of the central loop that are inverted repeats.

4. The length of the inverted repeat is also unknown.  So here we use
the idea of intervals.  An interval is a possible size for an inverted
repeat. Each tuple size has several interval sizes.  In our example,
interval sizes might be 100, 200, 400, 800, 1600, 3200.  In general, I
like the idea of using powers of 2 spacing between interval sizes.

5. Given a {\bf Center} that we are testing, we can check each interval
size in turn for an inverted repeat.  Start with the largest size
interval and work down.  An inteverted repeat of size 600 might not be
found in the interval of size 800, but will be found in the interval of
size 400.

6. Since the largest inverted repeat allowed (Tuplemaxdistance/2) is
larger than the interval size, the interval has to slide along the list
of matches found for the given center, sliding from the center out.    
Along the way, the number of matches in the interval is tested.  This
includes matches in the given center and matches at nearby centers, also
within the interval.

7.  All intervals for a center slide ahead when a new match for the center is found.
Two pointers in the center are maintained for this.  They are called
{\bf interval_lo} and {\bf interval_hi}.  When a new match is found,
interval_hi is set to point to the entry containing this new match. 
Then interval_lo is dragged ahead as far as it needs to go.

8.  When a nearby center is tested, it also has to consider the
interval.  In this case, the interval might be shifted left or right. 
This is done with the same pointers.  So there is some back and forth
shifting of these pointers.

*/

struct intervallist{
	int num_intervals;         /* number of intervals for the specific tuple size */
	int size[MAXNUMINTERVALS]; /* size of the intervals for the specific tuple size */
	int RNKP[MAXNUMINTERVALS]; /* RNKP criteria for the intervals for specific tuple size */ 
  int delta[MAXNUMINTERVALS];                 /* delta is +/- centers for the interval size  */
	/* next added by Gary Benson 8/2003 */
	double frequent_center_ratio[MAXNUMINTERVALS]; /* fraction of RNKP criteria.  
															 This is taken from simulation of the ratio of matches found in 
															 k-tuples classified as belonging to the most frequent center, 
															 i.e., the center with the most k-tuple matches, over the RNKP value 
															 for the interval size.  This fraction is the 5% cutoff in those
															 simulations, meaning 95% of the time a higher ratio is found. */
} Intervals[MAXTUPLESIZES];

/* added in v3.01, this number is how much further the distance would look */
#define MAXINTERVALSIZE 2400

/*******************************************************************/
/*********************** InitIntervals() ***************************/
/************************ 05-3-01 started *************************/
/*******************************************************************/

void InitIntervals(void)
{

  struct distribution_parameters distparam; /* distparam is the object containing the average
                                               (expected) number of matches and the variance for
                                               this number */

  int h;


  Intervals[1].num_intervals=2;
  Intervals[2].num_intervals=2;
  Intervals[3].num_intervals=4;



  /* g=1 */
  for(h=1;h<=Intervals[1].num_intervals;h++)
  {
		Intervals[1].size[h]=20*pow(2,h-1);

		distparam=distribution_k_run_sums_version4(Pmatch, Intervals[1].size[h], Tuplesize[1]); /* call to procedure for expectation and variance */
		Intervals[1].RNKP[h]=distparam.exp-1.65*sqrt(distparam.var); /* formula for calculating 5% */

		Intervals[1].delta[h]=(int)floor(2.3*sqrt(Pindel* Intervals[1].size[h]));
		//printf("\nPmatch: %lf  n: %d  k: %d  exp: %.8lf var: %.8lf\n",Pmatch,Intervals[1].size[h],Tuplesize[1],distparam.exp,distparam.var);
  }
  Intervals[1].frequent_center_ratio[1]=0.80;
	Intervals[1].frequent_center_ratio[2]=0.51;

  /* g=2 */
  for(h=1;h<=Intervals[2].num_intervals;h++)
  {
		Intervals[2].size[h]=80*pow(2,h-1);

		distparam=distribution_k_run_sums_version4(Pmatch, Intervals[2].size[h], Tuplesize[2]); /* call to procedure for expectation and variance */
		Intervals[2].RNKP[h]=distparam.exp-1.65*sqrt(distparam.var); /* formula for calculating 5% */

		Intervals[2].delta[h]=(int)floor(2.3*sqrt(Pindel* Intervals[2].size[h]));
        //printf("\nPmatch: %lf  n: %d  k: %d  exp: %.8lf var: %.8lf\n",Pmatch,Intervals[2].size[h],Tuplesize[2],distparam.exp,distparam.var);
  }
  Intervals[2].frequent_center_ratio[1]=0.44;
	Intervals[2].frequent_center_ratio[2]=0.31;

  /* g=3 */
  for(h=1;h<=Intervals[3].num_intervals;h++)
  {
		Intervals[3].size[h]=300*pow(2,h-1);

		distparam=distribution_k_run_sums_version4(Pmatch, Intervals[3].size[h], Tuplesize[3]); /* call to procedure for expectation and variance */
		Intervals[3].RNKP[h]=distparam.exp-1.65*sqrt(distparam.var); /* formula for calculating 5% */

		Intervals[3].delta[h]=(int)floor(2.3*sqrt(Pindel* Intervals[3].size[h]));
        //printf("\nPmatch: %lf  n: %d  k: %d  exp: %.8lf var: %.8lf\n",Pmatch,Intervals[3].size[h],Tuplesize[3],distparam.exp,distparam.var);
  }
  Intervals[3].frequent_center_ratio[1]=0.27;
	Intervals[3].frequent_center_ratio[2]=0.19;
  Intervals[3].frequent_center_ratio[3]=0.13;
	Intervals[3].frequent_center_ratio[4]=0.09;


  /* debug */
 /* debugmessage("\n\nIntervals"
    "\n  tupsize index   tuplesize   intervalsize index   intervalsize    delta       RNKP Pmatch=%lf",Pmatch);

  for(g=1;g<=NTS;g++)
  {
    debugmessage("\n");
    for(h=1;h<=Intervals[g].num_intervals;h++)
    {
      debugmessage("\n    %d                %d                 %d               %4d         %3d         %3d",
        g, Tuplesize[g],h,Intervals[g].size[h],Intervals[g].delta[h],Intervals[g].RNKP[h]);
    }
  }


  getch();
  exit(0);*/

}



/*******************************************************************/
/*********************** TestIntervalCriteria4() ********************/
/************************ 12-05-02 started *************************/
/*******************************************************************/



int TestIntervalCriteria4(int g, int t, int h, int i)
/* g is tuplesize index */
/* t is true center  */
/* h is the index of the interval */
/* i is the location of last matched character in sequence */
{
	
	int s,c,delta,t_matches,fl,flc,fh,fhc,sum,m,fc;
	struct centerlist *centergc, *centergfc;
	struct matchentry *firstmatch, *currentmatch;
	int index,maxindex,intervalend;
	struct intervallist *intervalsg;
	int old_h_index;
	int sum1; /* to correct error using sum for two variables Gary Benson 11/18/03 */
	
	s=Centerlistsize[g];
	c=t%s;
	centergc=&Center[g][c];
	intervalsg=&Intervals[g];
	delta=intervalsg->delta[h];

	//printf("\n\nTestIntervalCriteria4(%d,%d,%d,%d)",g,t,h,i);
	
	/* get matches for t */
	t_matches=centergc->intervalmatches[h];
	//printf("\n  t_matches=%d",t_matches);

	/* get matches in delta range below t */
	fl=max(t-2*delta,Mincenter[g]);
	flc=fl%s;
	
	sum=0;
	fc=flc;
	while(fc!=c)
	{
		centergfc=&Center[g][fc];
		if(centergfc->intervalmatches[h]!=0) /* only process if there are matches recorded */
		{
			/* test here that the matches recorded in intervalmatches[h] 
			fall into the current interval using the interval boundary 
			variable interval_lo[h].  If not, then recalculate the 
			interval matches for this value of g */
			
			if(centergfc->interval_lo[h]<i-intervalsg->size[h]+1) 
			{
				/* recalculate intervalmatches[h] and interval_lo[h] for 1..h */
				firstmatch=centergfc->match;
				currentmatch=firstmatch+centergfc->end_of_list;
				old_h_index=currentmatch->h_index; /* keeps track of the high end of the lowest match in an interval */
				
				sum1=0;		
				index=1;
				maxindex=intervalsg->num_intervals;
				intervalend=i-intervalsg->size[1]+1;
				
				while(currentmatch>=firstmatch)
				{
				while(currentmatch->h_index<intervalend) /* we are counting any match run that is partially or wholly
					within the interval */
				{
					centergfc->intervalmatches[index]=sum1;
					centergfc->interval_lo[index]=old_h_index; /* store h_index for lo end match in this interval for this count of matches */
					index++;
					if(index>maxindex) break; /* break out of while loop */
					intervalend=i-intervalsg->size[index]+1;
				}
				if(index>maxindex) break;
				sum1+=currentmatch->num_matches;
				old_h_index=currentmatch->h_index;
				currentmatch--;
				}
				while(index<=maxindex)
				{
					centergfc->intervalmatches[index]=sum1;
					centergfc->interval_lo[index]=old_h_index;
					index++;
				}
				
			}
			m=centergfc->intervalmatches[h];
		  //printf("\n  m=%d, (%d,%d,%d,%d)",m,i,g,fc,h);


			if(m>t_matches) return(0);  /* t is not best center */
			sum+=m;
		}
		//fc=(++fc%s);
		fc++;
		fc=fc%s;
	}
	
	sum+=t_matches;
	
	/* get matches in delta range above t */
	fh=min(t+2*delta,Maxcenter[g]);
	fhc=fh%s;
	fc=(c+1)%s;
	while(fc!=fhc)
	{
		centergfc=&Center[g][fc];
		if(centergfc->intervalmatches[h]!=0) /* only process if there are matches recorded */
		{
			/* test here that the matches recorded in intervalmatches[h] 
			fall into the current interval using the interval boundary 
			variable interval_lo[h].  If not, then recalculate the 
			interval matches for this value of g */
			
			if(centergfc->interval_lo[h]<i-intervalsg->size[h]+1) 
			{
				/* recalculate intervalmatches[h] and interval_lo[h] for 1..h */
				firstmatch=centergfc->match;
				currentmatch=firstmatch+centergfc->end_of_list;
				old_h_index=currentmatch->h_index; /* keeps track of the high end of the lowest match in an interval */
				
				sum1=0;		
				index=1;
				maxindex=intervalsg->num_intervals;
				intervalend=i-intervalsg->size[1]+1;
				
				while(currentmatch>=firstmatch)
				{
				while(currentmatch->h_index<intervalend) /* we are counting any match run that is partially or wholly
					within the interval */
				{
					centergfc->intervalmatches[index]=sum1;
					centergfc->interval_lo[index]=old_h_index; /* store h_index for lo end match in this interval for this count of matches */
					index++;
					if(index>maxindex) break; /* break out of while loop */
					intervalend=i-intervalsg->size[index]+1;
				}
				if(index>maxindex) break;
				sum1+=currentmatch->num_matches;
				old_h_index=currentmatch->h_index;
				currentmatch--;
				}
				while(index<=maxindex)
				{
					centergfc->intervalmatches[index]=sum1;
					centergfc->interval_lo[index]=old_h_index;
					index++;
				}
				
			}
			m=centergfc->intervalmatches[h];
		  //printf("\n  m=%d, (%d,%d,%d,%d)",m,i,g,fc,h);
			if(m>t_matches) return(0);  /* t is not best center */
			sum+=m;
		}
		//fc=(++fc%s);
		fc++;
		fc=fc%s;
	}
	
	//printf("\ntotal=%d, test=%d",sum,intervalsg->RNKP[h]);
	/* if sum is greater than rnkp criteria, return true */
	if(sum>=intervalsg->RNKP[h]) 
	{
		/* following 7 lines added by Gary Benson, 12/04/02 for testing */
		TestMaxIntervalSum.tupleindex=g;
		TestMaxIntervalSum.intervalindex=h;
		TestMaxIntervalSum.matches_in_interval=sum;
		TestMaxIntervalSum.lowcenter=fl;
		TestMaxIntervalSum.highcenter=fh;
		TestMaxIntervalSum.testcenter=t;
		return(1);
	}
	else return(0);
}
 
#endif
