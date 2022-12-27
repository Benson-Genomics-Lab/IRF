#ifdef removethis_stuff


void PrintMaxIntervalSum(int i)
/* i is the location of the tuple match in the sequence */
{

  int f,r,t,g,h,s,matches,interval,interval_lo,c,cc,k;
  
  f=TestMaxIntervalSum.lowcenter;
  r=TestMaxIntervalSum.highcenter;
  t=TestMaxIntervalSum.testcenter;

  g=TestMaxIntervalSum.tupleindex;
  h=TestMaxIntervalSum.intervalindex;

 	s=Centerlistsize[g];
  interval=Intervals[g].size[h];
	
	interval_lo=i-interval;
	matches=TestMaxIntervalSum.matches_in_interval;

  
  fprintf(Fptxt,"\nInterval_Lo=%d  Interval_Hi=%d  Matches_in_Interval=%d",
    interval_lo,i,matches);
    
  
  fprintf(Fptxt,"\nCenter        Matchcount");
  for(c=f;c<=r;c++)
	{
		cc=c%s;
		matches=Center[g][cc].intervalmatches[h];
		if(matches!=0)
		{
			if(c==t)
				fprintf(Fptxt,"\n*%12d    %5d",c,matches);
			else
  			fprintf(Fptxt,"\n %12d    %5d",c,matches);
			for(k=0;k<=Center[g][cc].end_of_list;k++)
				if(Center[g][cc].match[k].h_index>interval_lo)
  					fprintf(Fptxt,"\n                     %d  %d %d",Center[g][cc].match[k].num_matches,
						Center[g][cc].match[k].l_index,Center[g][cc].match[k].h_index);
		}
  }

  fprintf(Fptxt,"\n");
}

#endif

