
void Clear_MaxMincenters(void) {

   int g;

  for(g=1;g<=NTS;g++)
  {
    Maxcenter[g]=2*Length-2*Tuplesize[g]+1;
    Mincenter[g]=2*Tuplesize[g]+1;
    Runningmaxcenter[g]=0;
    Runningmincenter[g]=0;
    Oldrunningmaxcenter[g]=0;
    Oldrunningmincenter[g]=0;
	/*debugmessage("\n\nInit_MaxMincenters\n  g  Maxcenter  Mincenter  Runningmax  Runningmin  Oldrunmax  Oldrunmin");
    debugmessage("\n%3d %6d  %9d  %9d   %9d   %9d",g,Maxcenter[g],Mincenter[g],
      Runningmaxcenter[g],Runningmincenter[g],Oldrunningmincenter[g]);
	*/
  }

}



void Init_MaxMincenters(void)
{

  Maxcenter=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Mincenter=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Runningmaxcenter=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Runningmincenter=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Oldrunningmaxcenter=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  Oldrunningmincenter=(int *)scalloc(MAXTUPLESIZES+1,sizeof(int));
  
  memory_stats_print("\n Init_MaxMincenters(Maxcenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters(Mincenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters(Runningmaxcenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters(Runningmincenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters(Oldrunningmaxcenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );
  memory_stats_print("\n Init_MaxMincenters(Oldrunningmincenter): requesting", (MAXTUPLESIZES+1) * sizeof(int) );

  /* added by Gelfand on 10/30/03 
  if ( NULL == Maxcenter || NULL == Mincenter || NULL == Runningmaxcenter || NULL == Runningmincenter || NULL == Oldrunningmaxcenter || NULL == Oldrunningmincenter ) 
    { debugerror("\nMemory error in Init_MaxMincenters. Aborting!"); }
*/

  Clear_MaxMincenters();

}


void Update_Runningmaxmincenters2(int i, int g)
/* moves the maxcenter and mincenter as the index of the sequence moves */
/* i is index in sequence */
/* g is in 1..NTS */
{
	int h,ormc,cormc;
	
	/* define the new Runningmaxcenter and Runningmincenter */
	
	Runningmaxcenter[g]=2*i-2*Tuplesize[g]+1;  
	Runningmincenter[g]=2*i-Tuplemaxdistance[g]-Tuplesize[g]+1;  
	//printf("\n    g=%3d  Runningmaxcenter=%3d  Runningmincenter=%3d  Oldrunningmaxcenter=%3d  Oldrunningmincenter=%3d",
	//	g,Runningmaxcenter[g],Runningmincenter[g],Oldrunningmaxcenter[g],Oldrunningmincenter[g]);
	
	/* clear centers Oldruningmincenter to Runningmincenter-1 */
	if(Runningmincenter[g]>=Mincenter[g]) /* Mincenter is the smallest possible center */
	{
		/* ormc is a true center */
		ormc=Oldrunningmincenter[g];
		while(ormc!=Runningmincenter[g])
		{
			/* cormc is an index in the centerlist */
			cormc=ormc%Centerlistsize[g];
			
			/* remove information stored in center */
			
			Center[g][cormc].index=-2; /* -2 means no index */
			
			Center[g][cormc].numentries=0;
			Center[g][cormc].nummatches=0;
			Center[g][cormc].end_of_list=-1;
			for(h=0;h<MAXNUMINTERVALS;h++)
				Center[g][cormc].intervalmatches[h]=0;

			ormc++;
		}
		Oldrunningmincenter[g]=ormc;
	}
	
	/* no clearing necessary for runningmaxcenter */

	Oldrunningmaxcenter[g]=Runningmaxcenter[g];
}
