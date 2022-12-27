/****************************************************/
/* add_tuple_match_to_Center2                        */
/****************************************************/

void add_tuple_match_to_Center2(int location, int size, int d, int g, int c) 

/* location is index of tuple        */
/* size is size of tuple             */
/* d is the true center              */
/* g is in [1..NTS]                  */
/* c is index in circular centerlist */

{
  struct centerlist *centergc;
  struct matchentry *lastmatch, *currentmatch, *firstmatch;
  int sum, maxindex, index, intervalend, old_h_index;
	
  centergc=&Center[g][c];
	
  /* test if this list is currently used for entries at center d    */
  /* if not, initialize the information, updaterunningminmaxcenters */
	/* has already cleared the entry                                  */
	
  if(centergc->index!=d) centergc->index=d;
	
 	lastmatch=centergc->match+centergc->end_of_list;

  if((centergc->numentries!=0)
    &&(lastmatch->h_index==location-1)) 
		/* there are entries and this is an adjacent tuple, */
		/* just add on to last entry */
  {
		//debugmessage("\n        add to existing entry");
		lastmatch->h_index++;
		lastmatch->tot_matches++;
		lastmatch->num_matches++;
		centergc->nummatches++;
  }
  else          /* need a new entry here */
  {
		
		
		//debugmessage("\n        new entry");
		centergc->numentries++;
		centergc->end_of_list++;
		
		lastmatch++;
		lastmatch->h_index=location;
		lastmatch->l_index=location-size+1;
		lastmatch->num_matches=size;
		
		centergc->nummatches+=size;
		lastmatch->tot_matches=centergc->nummatches;
		
  }

	
	/* get totals for different interval sizes */
	
  firstmatch=centergc->match;
  currentmatch=lastmatch;
  old_h_index=currentmatch->h_index; /* keeps track of the high end of the lowest match in an interval */
  
  
  sum=0;
  
  index=1;
  maxindex=Intervals[g].num_intervals;
  intervalend=location-Intervals[g].size[1]+1;
  
  while(currentmatch>=firstmatch)
  {
	while(currentmatch->h_index<intervalend) /* we are counting any match run that is partially or wholly
												within the interval */
	{
	  centergc->intervalmatches[index]=sum;
	  centergc->interval_lo[index]=old_h_index; /* store h_index for lo end match in this interval 
	  for this count of matches */
	  index++;
	  if(index>maxindex) break; /* break out of while loop */
	  intervalend=location-Intervals[g].size[index]+1;		
	}
	if(index>maxindex) break;
	sum+=currentmatch->num_matches;
	old_h_index=currentmatch->h_index;
	currentmatch--;
  }
  
  while(index<=maxindex)
  {
	  centergc->intervalmatches[index]=sum;
	  centergc->interval_lo[index]=old_h_index;
	  index++;
  }
	
}


