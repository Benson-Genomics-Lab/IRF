#ifndef _MYMEMLIB
#define _MYMEMLIB

#define HIDING_BUFFER_SIZE	200000
#define HIDING_PLACE		100000

char	HidingBuffer[HIDING_BUFFER_SIZE];



/* added by Gelfand on 10/30/03 */
//#define MEMORY_TRACKER
void __memory_stats_print(char * message, int size )	{ 
    
    static int total = 0;

    total+=size;

    fprintf(stderr,"%s %d bytes (total = %d MB )", message, size, total / (1024*1024) );
    fflush(stderr);

}


#ifdef MEMORY_TRACKER
    #define memory_stats_print(_x, _y) __memory_stats_print(_x, _y)
#else
    #define memory_stats_print(_x, _y) ;   
#endif

    

/* The idea is that we take a copy of our memory and periodically check if it has not changed
The problem is that our copy is also succeptible to "unathorized foreign changes." To counter
that, let's hide our copy inside a large preallocated memory buffer. This way, the chance of
a random "writethrough" is minimized */ 


/********************************************************************************************/
void MYMEMLIB_createMemorySnapshot(void * mptr, int size) {
/* mptr  -  pointer to memory that we need to remember */
/* size  -  size of that memory */

	int		*sizePtr;			// remember the size
	void	**locationPtr; 		// remember the location of the source
	char	*memPtr;			// remember the data

	// resolve our variable space
	sizePtr		=	(int*)		(HidingBuffer+HIDING_PLACE);
	locationPtr	=	(void**)	(HidingBuffer+HIDING_PLACE+sizeof(int));
	memPtr		=	(char*)		(HidingBuffer+HIDING_PLACE+sizeof(int)+sizeof(void*));


	// take a snapshot
	*sizePtr=size;
	*locationPtr=mptr;
	memcpy(memPtr,mptr,size);
}
/********************************************************************************************/
int MYMEMLIB_checkMemoryIntegrity(void) {
	
	int		*sizePtr;			// remember the size
	void	**locationPtr; 		// remember the location of the source
	char	*memPtr;			// remember the data

	// resolve our variable space
	sizePtr		=	(int*)		(HidingBuffer+HIDING_PLACE);
	locationPtr	=	(void**)	(HidingBuffer+HIDING_PLACE+sizeof(int));
	memPtr		=	(char*)		(HidingBuffer+HIDING_PLACE+sizeof(int)+sizeof(void*));

	// check memory integrity
	return memcmp(*locationPtr,memPtr,*sizePtr);
}

/* create a mem snapshot, for debugging, remove this later!!! 
  if (Centerseenlist->next!=NULL) MYMEMLIB_createMemorySnapshot((void*)Centerseenlist->next,sizeof(struct centerlistelement));

  //Centerseenlist->rightend=7;


  if (Centerseenlist->next!=NULL)
	  if (MYMEMLIB_checkMemoryIntegrity()!=0) {
		printf("ERROR: (MYMEMLIB_checkMemoryIntegrity) memory was corrupted. Aborting!!!\n");
		exit(0);
	  }

*/
#endif
