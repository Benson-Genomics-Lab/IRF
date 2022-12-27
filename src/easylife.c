#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include <stdarg.h> 


#include "easylife.h"


/******************************************************************************************

EASY_LIFE - data structures, memory checks, useful code snippets

*   Written by Yevgeniy Gelfand, April 3, 2003
*   Last Revision Date  :   January 18, 2005

*   This library is free software; you can redistribute it and/or
*   modify it under the terms of the GNU Lesser General Public
*   License as published by the Free Software Foundation; either
*   version 2.1 of the License, or (at your option) any later version.
*
*   This library is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*   Lesser General Public License for more details.
*
*   You should have received a copy of the GNU Lesser General Public
*   License along with this library; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
*   Please direct all questions related to this program or bug reports or updates to
*   Yevgeniy Gelfand (ygelfand@bu.edu)


DESCIPTION:

 *  This module provides a lot of usefull structures for the mannipulation of:
      arrays, hashes, lists, stacks, queues, memory checks, etc...

 *  This module is thread safe.

 
 *  These people provided help or ideas with this module:
        - Patrick Coskren
        - Alfredo Rodriguez
        - Douglas Ehlenberger

 *  qsort and insertion sort (which I modified to fit my data structures) come from "Mastering Algorithms with C"  
    by Kyle Loudon, published by O'Reilly & Associates. Code used is surrounded by special comments.
    That is a great book, I highly recommend it.


SIDE EFFECTS, LIMITATIONS: 

 *  This library was tested on a Windows platform only, but with minor modifications
    it should work on any platform.

 *  Note: doCriticalErrorAndQuit must be defined in order to use
    this module. See buttom of this page for prototype. That
    function is called on memory allocation failures,
    illegal lookups and incorrect input sent into functions.
    All structures here grow and shrink dynamically (except hash, 
    didn't get around to that yet).

 *  It would be nice to make the hash dynamic. Also, integer hash 
    is based on a string hash now, which is slow, would be nice to rewrite it.
    A balanced tree is something that I might add to this module. 
    I was asked to write a srealloc function as well.

VERSIONS:

 *  Version  :   1.02.
        - added library to check memory allocations (smalloc,scalloc,sfree)

 *  Version  :   1.01.
        - found a bug in the list copy function: 
          delete function wasn't set correctly
        - added a list append function
 
 *  Version  :   1.00
        - first release version

*******************************************************************************************/


/******************************** internal *******************************\
*
* The functions below are used internally by Easy Life
*
\*************************************************************************/

int EasyLifeGetPrime(int number);
unsigned int EasyLifeElfUp(char *text);



/******************************* streamlibs ******************************\
*
* The functions below provide C++ like character stream capabilites
*
\*************************************************************************/


/* --------------------------------------------------------------------- */
/* works similarly to fgets, see fgets help for more info				 */
BSTREAM* bopen( char *data ) {

	BSTREAM *bp = (BSTREAM*)malloc( sizeof(BSTREAM) );

	if ( NULL == bp ) return NULL; 

	bp->lastTBFDataPos = data;

	return bp;
}

/* --------------------------------------------------------------------- */
void bclose( BSTREAM *bp ) {

	free(bp);
}

/* --------------------------------------------------------------------- */
char* bgets( char *buffer, int lineLength, BSTREAM *bp ) {


	char* src;
	int count = 0;

	/* exit condition */
	if ( NULL == bp || NULL == bp->lastTBFDataPos) return NULL;
	src = bp->lastTBFDataPos;

	/* Read one line at a time until the end of the data. */
    while(*src!=13&&*src!=10&&*src!='\0'&&count<=lineLength)  {  buffer[count] = *src; src++; count++; }
	buffer[ ( count > lineLength ) ? lineLength : count ]=0;
	if ( *src == 0 || count > lineLength ) {
		bp->lastTBFDataPos = NULL;
		return buffer;
	}

	/* get to the next printable character */
    if ( *src == 13 && *(src+1) == 10 )
        src+=2;
	else if( *src == 13 || *src == 10 ) 
        src++;

	if ( *src == 0 )  {
		bp->lastTBFDataPos = NULL;
		return NULL;
	}

	/* set the stoppage at next good character */
	bp->lastTBFDataPos = src;
	return buffer;

}



/******************************* EasyArray *******************************\
*
* The functions below let you maintain a dynamic array (doubles each time)
*
* Note: the last item of the array points to NULL
*
* Note: the size of an array is the number of items not including
* the item pointing to NULL
*
* Note: it is unsafe to use CopyEasyArray functions on an array
* that contains variable length strings. Make sure each element
* is of the same size if you intend to use use the CopyEasyArray 
* function. Otherwise you might recieve an "illegal memory access" 
* error. Or better yet, use CopyEasyArrayOfCStrings.
*
\*************************************************************************/

#define		EASY_ARRAY_INITIAL_SIZE		500 

/* --------------------------------------------------------------------- */
EASY_ARRAY* EasyArrayCreate( int initialSize, void* (*copy)(const void *data), void (*destroy)(void *item) ) {

	EASY_ARRAY *newEasyArray;

	newEasyArray = (EASY_ARRAY*) calloc( 1, sizeof(EASY_ARRAY) );
	if ( NULL == newEasyArray ) {
		doCriticalErrorAndQuit("\nCreateEasyArray: Memory error. Aborting!");
		exit(1);
	}

	newEasyArray->size = 0;
	newEasyArray->copy = copy;
	newEasyArray->destroy = destroy;
    newEasyArray->reserved = initialSize > 0 ? initialSize : EASY_ARRAY_INITIAL_SIZE;
	newEasyArray->array = (void**) calloc( newEasyArray->reserved + 1 , sizeof(void*) );
	if ( NULL == newEasyArray->array) {
		doCriticalErrorAndQuit("\nCreateEasyArray: Memory error. Aborting!");
		exit(1);
	}

	return newEasyArray;
}

/* --------------------------------------------------------------------- */
EASY_ARRAY* EasyArrayCopy( EASY_ARRAY* easyArray, int copyData ) {

  int			i, size, reserved; 
  EASY_ARRAY	*newEasyArray;
  void			*newData;


    /* check if not NULL (note, empty means at least one element exists) */
    if ( NULL == easyArray ) {
        return NULL;
    }


    /* allocate new structure */
    newEasyArray = ( EASY_ARRAY* ) calloc( 1, sizeof( EASY_ARRAY ) );
    if ( NULL == newEasyArray ) {
	    doCriticalErrorAndQuit("\nEasyArrayCopy: Memory error1. Aborting!");
	    exit(1);
    }

    reserved = newEasyArray->reserved = easyArray->size * 2;
    size = newEasyArray->size = easyArray->size;
    newEasyArray->array = ( void** ) calloc( reserved + 1 , sizeof( void* ) );
    if ( NULL == newEasyArray->array ) {
        doCriticalErrorAndQuit("\nEasyArrayCopy: Memory error2. Aborting!");
        exit(1);
    }


    /* if copyData, create copies of items and place into the new array */
    if (copyData) {
     

        /* assign copy and destroy functions */
        newEasyArray->copy = easyArray->copy;
        newEasyArray->destroy = easyArray->destroy;       


        /* check if the copy function is provided */
        if ( NULL == newEasyArray->copy) {
	        doCriticalErrorAndQuit("\nEasyArrayCopy: No copy function provided. Aborting!");
	        exit(1);
        }


        /* create copies of items and place into the new array */
	    for ( i = 0; i < size; i++ ) {

            newData = newEasyArray->copy( easyArray->array[i] );

            /* check if memory alloc did not fail */
	        if ( NULL == newData ) {
		        doCriticalErrorAndQuit("\nEasyArrayCopy: Memory error. Aborting!");
		        exit(1);
	        }

		    newEasyArray->array[i] = newData;
	    }


    /* else, copy list item pointers into the new array */
    } else {


        /* assign copy and destroy functions */
        newEasyArray->copy = NULL;
        newEasyArray->destroy = NULL;   
        
        
        /* copy list item pointers */
	    for ( i = 0; i < size; i++ ) {
		    newEasyArray->array[i] = easyArray->array[i];
	    }

    }


    /* a check */
    if ( i != newEasyArray->size ) { /* one extra for NULL */
	   doCriticalErrorAndQuit("\nEasyArrayCopy: number of elements does not match the size field of the incoming list. Aborting!");
       exit(1);
    }


    /* set last element to NULL */
    newEasyArray->array[i] = NULL;

    
    return newEasyArray;

}

/* --------------------------------------------------------------------- */

/**************************************************************************
* This function creates an Easy Array from an Easy List. If the copyData
* flag is set to 0, data is not copied, only pointers to data (basically
* array shares data with the list. Also COPY and DESTROY functions of the
* new array are set to NULL. If copyData is 1, data is copied and COPY
* and DESTROY functions are also inherrited from the list. 
**************************************************************************/
EASY_ARRAY* EasyArrayCreateFromList ( EASY_LIST* easyList, int copyData ) {

    EASY_NODE *eListIterator;
	EASY_ARRAY *newEasyArray;
    int i;
    void **newData;


	/* check if not NULL (note, empty means at least one element exists) */
	if ( NULL == easyList ) {
		return NULL;
	}


    /* create array structure */
	newEasyArray = (EASY_ARRAY*) calloc( 1, sizeof(EASY_ARRAY) );
	if ( NULL == newEasyArray ) {
		doCriticalErrorAndQuit("\nEasyArrayCreateFromList: Memory error. Aborting!");
		exit(1);
	}


    /* create array elements (the same number as items in the easy list passed in ) */
	newEasyArray->size = easyList->size;
	newEasyArray->reserved = easyList->size  * 2;
	newEasyArray->array = (void**) calloc( newEasyArray->reserved + 1 , sizeof(void*) );
	if ( NULL == newEasyArray->array) {
		doCriticalErrorAndQuit("\nEasyArrayCreateFromList: Memory error. Aborting!");
		exit(1);
	}


    /* if copyData, create copies of items and place into the new array */
    if (copyData) {
     

        /* assign copy and destroy functions */
        newEasyArray->copy = easyList->copy;
        newEasyArray->destroy = easyList->destroy;       


        /* check if the copy function is provided */
        if ( NULL == newEasyArray->copy) {
	        doCriticalErrorAndQuit("\nEasyArrayCreateFromList: No copy function provided. Aborting!");
	        exit(1);
        }


        /* create copies of items and place into the new array */
	    for ( i = 0, eListIterator = easyList->head; eListIterator != NULL;  eListIterator = eListIterator->next, i++ ) {

            newData = (void**)newEasyArray->copy( eListIterator->item );

            /* check if memory alloc did not fail */
	        if ( NULL == newData ) {
		        doCriticalErrorAndQuit("\nEasyArrayCreateFromList: Memory error. Aborting!");
		        exit(1);
	        }

		    newEasyArray->array[i] = newData;
	    }


    /* else, copy list item pointers into the new array */
    } else {


        /* assign copy and destroy functions */
        newEasyArray->copy = NULL;
        newEasyArray->destroy = NULL;   
        
        
        /* copy list item pointers */
	    for ( i = 0, eListIterator = easyList->head; eListIterator != NULL;  eListIterator = eListIterator->next, i++ ) {
		    newEasyArray->array[i] = eListIterator->item;
	    }

    }


    /* a check */
    if ( i != newEasyArray->size ) { /* one extra for NULL */
		doCriticalErrorAndQuit("\nEasyArrayCreateFromList: number of elements does not match the size field of the incoming list. Aborting!");
		exit(1);
    }


    /* set last element to NULL */
    newEasyArray->array[i] = NULL;
    

    /* return the new structure */
    return newEasyArray;

}

/* --------------------------------------------------------------------- */
void EasyArrayInsert( EASY_ARRAY* easyArray, void *item ) {

	int	 arraySize=0, newReserved, i;


	/* check if not NULL (note, empty means at least one element exists) */
	if ( NULL == easyArray ) {
		doCriticalErrorAndQuit("\nAddItemToEasyArray: NULL value received. Aborting!");
		exit(1);
	}


	/* get size of the array + 1 */
	arraySize = easyArray->size;


	/* check if not too large for our data, else double the size of our data */
	if ( arraySize >= easyArray->reserved ) {

		newReserved = 2 * easyArray->reserved;
		easyArray->array = (void**)realloc( easyArray->array, ( newReserved + 1 ) * sizeof( void* ) );
		if ( NULL == easyArray->array ) {
			doCriticalErrorAndQuit("\nAddItemToEasyArray: Memory error. Aborting!");
			exit(1);		
		}
		easyArray->reserved = newReserved;


		/* zero out the rest of the array */
		for ( i=arraySize; i<=newReserved; i++ )
			easyArray->array[i] = NULL;
	}


	/* assign the pointer */
	easyArray->array[arraySize] = item;
	easyArray->size = arraySize + 1;


}

/* --------------------------------------------------------------------- */
void EasyArrayDestroy(EASY_ARRAY* easyArray) {

  int	i = 0;
  void	*eachItem;


  /* check if not NULL (note, empty means at least one element exists) */
  if ( NULL == easyArray ) {
	return;
  }


  /* Iterate through the elemets until a NULL value is encountered.*/
  if ( NULL != easyArray->destroy ) {

      while ( NULL != ( eachItem = easyArray->array[i] )) {
    
	      easyArray->destroy( eachItem );

	      i++;

      }

  }

  /* delete the pointers to items */
  free( easyArray->array );


  /* No operations are allowed now, but clear the structure as a precaution. */
  memset( easyArray, 0, sizeof(easyArray) );


  /* delete the structure */
  free( easyArray );
  
}

/* --------------------------------------------------------------------- */
/* Note: you can simply access inside the array instead of calling this
function. The only advantage of this function is that it checks bounds */
void* EasyArrayItem( EASY_ARRAY* easyArray, int position ) {

  void **array;

  if ( position < 0 || position >= easyArray->reserved || NULL == ( array = (void**)easyArray->array[position] ) ) {
		doCriticalErrorAndQuit("\nGetItemFromEasyArray: Illegal lookup OR Value Not Set (reserved range:0-%d, asked for %d. Aborting!",easyArray->reserved-1,position);
		exit(1);
  }
	
  return array;

}

/* --------------------------------------------------------------------- */
void EasyArrayRemove( EASY_ARRAY* easyArray, int position ) {

  int j = 0;
  void **array;


  /* check if not NULL (note, empty means at least one element exists) */
  if ( NULL == easyArray || NULL == ( array = easyArray->array ) || position < 0 || position >= easyArray->size || NULL == array[position] ) {
		doCriticalErrorAndQuit( "\nRemoveItemFromEasyArray: Illegal lookup. Aborting!" );
		exit( 1 );
  }


  /* free it */
  if ( NULL != easyArray->destroy )
	  easyArray->destroy( array[position] ); 


  /* move all elements back */
  for ( j = position+1; array[j] != NULL; j++ ) {
	  array[j - 1] = array[j];
  }


  /* the last one is now free */
  array[j - 1] = NULL;
  easyArray->size--;


} 



/******************************* EasyList ********************************\
*
* The functions below let you maintain a dynamic doubly linked list
*
* Note: copy and destroy functions in the constructor are optional.
* If destroy is omitted, data is not freed.
* If copy is omitted, EasyArrayCopy operation will fail.
*
\*************************************************************************/

/* ---------------------------------------------------------------------- */
EASY_LIST* EasyListCreate( void* (*copy)(const void *data), void (*destroy)(void *item) ) {

	EASY_LIST *newEasyList;

	newEasyList = (EASY_LIST*) calloc( 1, sizeof(EASY_LIST) );
	if ( NULL == newEasyList ) {
		doCriticalErrorAndQuit("\nEasyListCreate: Memory error. Aborting!");
		exit(1);
	}

	newEasyList->size = 0;
	newEasyList->head = NULL;
	newEasyList->tail = NULL;
	newEasyList->copy = copy;
	newEasyList->destroy = destroy;

	return newEasyList;
}

/* ---------------------------------------------------------------------- */
void EasyListAppend( EASY_LIST* easyListTo, EASY_LIST* easyListFrom, int copyData ) {

	EASY_LIST *newEasyList;

    /* check if both exist */
	if (NULL == easyListTo || NULL == easyListFrom) return;

    /* create new from last */
    newEasyList = EasyListCopy( easyListFrom, copyData );

    /* add */
    (easyListTo->size)+=(newEasyList->size);
    if (easyListTo->tail) {
        
        easyListTo->tail->next = newEasyList->head;
        
        if (newEasyList->head) 
            newEasyList->head->prev = easyListTo->tail;

        if (newEasyList->tail)
            easyListTo->tail = newEasyList->tail;

    } else {

        easyListTo->head = newEasyList->head;
        easyListTo->tail = newEasyList->tail;
    }

    /* free structure, but not contents */
    free(newEasyList);
}

/* ---------------------------------------------------------------------- */
EASY_LIST* EasyListCopy	( EASY_LIST* easyList, int copyData ) {

	EASY_LIST *newEasyList;
	EASY_NODE *eListIterator;
    void **newData;

    /* create a new empty list */
	newEasyList = EasyListCreate ( easyList->copy, easyList->destroy );


    /* if copyData, create copies of items and place into the new array */
    if (copyData) {
     

        /* assign copy and destroy functions */
        newEasyList->copy = easyList->copy;
        newEasyList->destroy = easyList->destroy;       


        /* check if the copy function is provided */
        if ( NULL == newEasyList->copy) {
	        doCriticalErrorAndQuit("\nEasyListCopy: No copy function provided. Aborting!");
	        exit(1);
        }


        /* create copies of items and place into the new array */
	    for ( eListIterator = easyList->tail; eListIterator != NULL;  eListIterator = eListIterator->prev ) {

            newData =  (void**)newEasyList->copy( eListIterator->item );

            /* check if memory alloc did not fail */
	        if ( NULL == newData ) {
		        doCriticalErrorAndQuit("\nEasyListCopy: Memory error. Aborting!");
		        exit(1);
	        }

	    	EasyListInsertHead( newEasyList, newData );
	    }


    /* else, copy list item pointers into the new array */
    } else {


        /* assign copy and destroy functions */
        newEasyList->copy = NULL;
        newEasyList->destroy = NULL;   
        
        
        /* copy list item pointers */
        for ( eListIterator = easyList->tail; eListIterator != NULL;  eListIterator = eListIterator->prev ) {
        
            EasyListInsertHead( newEasyList, eListIterator->item );   
        
        }

    }


	return newEasyList;
}

/* ---------------------------------------------------------------------- */
void EasyListDestroy( EASY_LIST* easyList ) {

  /* Check if not NULL. */
  if ( NULL == easyList ) {
	return;
  }


  /* Iterate through the elemets until a NULL value is encountered.*/
  while ( EasyListSize( easyList )) {
	EasyListRemoveHead( easyList );
  }


  /* No operations are allowed now, but clear the structure as a precaution. */
  memset( easyList, 0, sizeof(easyList) );


  /* delete the structure */
  free( easyList );

}

/* ---------------------------------------------------------------------- */
void EasyListInsertAfter ( EASY_LIST* easyList, EASY_NODE* easyNode, void *item ) {

	EASY_NODE	*newEasyNode;

	/* Allocate storage for the new element. */
	newEasyNode = (EASY_NODE*) calloc( 1, sizeof(EASY_NODE) );
	if ( NULL == newEasyNode ) {
		doCriticalErrorAndQuit("\nEasyListInsertAfter: Memory error. Aborting!");
		exit(1);
	}

	/* Insert data. */
	newEasyNode->item = item;

	if ( easyNode == NULL ) {
	
		/* Handle insertion at the head of the list. */
		if ( easyList->size == 0) 
			easyList->tail = newEasyNode;
		else
			easyList->head->prev = newEasyNode;

		newEasyNode->next = easyList->head;
		newEasyNode->prev = NULL;
		easyList->head = newEasyNode;

	} else {

		/* Handle insertion somewhere other than at the head. */
		if ( easyNode->next == NULL )
			easyList->tail = newEasyNode;

		newEasyNode->next = easyNode->next;
		newEasyNode->prev = easyNode;

		easyNode->next = newEasyNode;
		if ( newEasyNode->next != NULL ) 
			newEasyNode->next->prev = newEasyNode;

	}

	/* Adjust the size. */
	easyList->size++;


}


/* ---------------------------------------------------------------------- */
void EasyListRemoveNode( EASY_LIST* easyList, EASY_NODE* easyNode ) {

	EASY_NODE	*oldEasyNode;

	/* Do not allow removal from an empty list.  */
	if ( EasyListSize( easyList ) == 0 ) {
		doCriticalErrorAndQuit("\nEasyListRemoveNode: Illegal operation. Aborting!");
		exit(1);
	}

	/* Remove the element from the list. */				
	if ( easyNode == NULL ) {

		/* Handle removal from the head of the list. */
		oldEasyNode = easyList->head;
		easyList->size--;


		if ( easyList->size == 0) {

			easyList->tail = NULL;
			easyList->head = NULL;

		} else {

			easyList->head = easyList->head->next;
			
			if ( easyList->head!= NULL ) 
				easyList->head->prev = NULL;
		}

	} else {

		oldEasyNode = easyNode;
		easyList->size--;

		if ( oldEasyNode->prev!= NULL ) 
			oldEasyNode->prev->next = oldEasyNode->next;
		else
			easyList->head = oldEasyNode->next;	

		if ( easyNode->next!= NULL ) 
			oldEasyNode->next->prev = oldEasyNode->prev;
		else
			easyList->tail = oldEasyNode->prev;	

	}


	/* Free the element: call user defined function if provided. */
	if ( NULL != easyList->destroy )
	  easyList->destroy( oldEasyNode->item );

	free( oldEasyNode );

}



/******************************** sorting ********************************\
*                                                                            
* Various sorting algorithms                                       
*                                                                            
**************************************************************************/


/*
 *  This insertion sort and qsort algortihms come from the "Mastering Algorithms with C"  by    
 *  Kyle Loudon, published by O'Reilly & Associates 
 *
 */

void EasyArrayInsertionSort( EASY_ARRAY* easyArray, int (*compare)( const void *item1, const void *item2 ) )  {

	void               **a = easyArray->array;
	void               *key;
	int                i, j;


	/* Repeatedly insert a key element among the sorted elements. */
	for ( j = 1; j < easyArray->size; j++ ) {

		key = a[j];
		i = j - 1;

		/* Determine the position at which to insert the key element. */
		while ( i >= 0 && compare( a[i], key ) > 0) {
		    a[i + 1] = a[i];
			i--;
		}

		a[i + 1] = key;

	}

}

/* ---------------------------------------------------------------------- */
void EasyListInsertionSort( EASY_LIST* easyList, int (*compare)( const void *item1, const void *item2 ) )  {

	EASY_NODE          *j, *i, *prev;
    void               *key;


    if ( easyList->size < 2 ) return;

	/* Repeatedly insert a key element among the sorted elements. */
	for ( j = easyList->head->next; j != NULL; j =  j->next ) {

        key = j->item;
		i = j->prev;
        prev = j;

		/* Determine the position at which to insert the key element. */
		while ( i != NULL && compare( i->item, key ) > 0 ) {
		    i->next->item = i->item;
            prev = i;
			i = i->prev;
		}

		prev->item = key;

	}

}

/* ---------------------------------------------------------------------- */
static int partition( EASY_ARRAY* easyArray, int i, int k, int (*compare)( const void *item1, const void *item2 ) ) {

    void               **a = easyArray->array;
    void               *pval, *temp;
    int                r[3],tempInt;


    /*  Use the median-of-three method to find the partition value. */
    r[0] = (rand() % (k - i + 1)) + i;
    r[1] = (rand() % (k - i + 1)) + i;
    r[2] = (rand() % (k - i + 1)) + i;
    
    if ( r[0] > r[1] ) { tempInt = r[1]; r[1] = r[0]; r[0] = tempInt; }
    if ( r[1] > r[2] ) { tempInt = r[2]; r[2] = r[1]; r[1] = tempInt; }

    pval = a[ r[1] ];

    /* qsort way (select middle value and swap it to the beginning of the array ) */
    //tempInt= ( k ) / 2;
    //temp = a[ tempInt ]; a[ tempInt ] = a[i]; a[i] = temp;
    //pval = a[ tempInt ];

    /* Create two partitions around the partition value. */
    i--;
    k++;

    while (1) {


       /* Move left until an element is found in the wrong partition. */
       do {

          k--;

       } while ( compare( a[k], pval ) > 0 ) ;


       /* Move right until an element is found in the wrong partition. */
       do {

          i++;

       } while ( compare( a[i], pval ) < 0 );

       if ( i >= k ) {

            /* Stop partitioning when the left and right counters cross. */
            break;

       } else {

            /* Swap the elements now under the left and right counters. */
            temp = a[i];
            a[i] = a[k];
            a[k] = temp;

       }

    }

    /* Return the position dividing the two partitions. */
    return k;

}

/* ---------------------------------------------------------------------- */
static void EasyArrayQuickSortHelper( EASY_ARRAY* easyArray, int i, int k, int (*compare)( const void *item1, const void *item2 ) ) {

    int j;

    /* Stop the recursion when it is not possible to partition further. */
    while (i < k) {

       /* Determine where to partition the elements. */
       j = partition( easyArray, i, k, compare );


       /* Recursively sort the left partition.*/
       if ( ( j - i ) >= 9 )
            EasyArrayQuickSortHelper( easyArray, i, j, compare );
       else {
        
           EASY_ARRAY tempArray;

           /* let insertion sort finish the job */
           tempArray.array = easyArray->array + i;
           tempArray.size = j - i + 1;
           EasyArrayInsertionSort( &tempArray, compare );

       }

       /* Iterate and sort the right partition.*/
       i = j + 1;

    }


}


/*
 *  End of code from the "Mastering Algorithms with C"  by Kyle Loudon,  
 *  published by O'Reilly & Associates 
 *
 */


/* ---------------------------------------------------------------------- */
void EasyArrayQuickSort( EASY_ARRAY* easyArray, int (*compare)( const void *item1, const void *item2 ) ) {

       
    EasyArrayQuickSortHelper( easyArray, 0, easyArray->size - 1, compare );
    //qsort( easyArray->array, easyArray->size, sizeof(void**), compare );

}

/* ---------------------------------------------------------------------- */
void EasyListQuickSort( EASY_LIST* easyList, int (*compare)( const void *item1, const void *item2 ) ) {


    EASY_NODE *eListIterator;
    int i;


    /* create an array out of list (without copying the data, just pointers) */
    EASY_ARRAY* easyArray = EasyArrayCreateFromList( easyList, 0 );


    /* sort pointers with quick sort */
    EasyArrayQuickSortHelper( easyArray, 0, easyArray->size - 1, compare );


    /* copy the pointers back to the list */
	for ( i = 0, eListIterator = easyList->head; eListIterator != NULL;  eListIterator = eListIterator->next, i++ ) {
		    eListIterator->item = EasyArrayItem( easyArray, i );
	}

    
    /* delete the array */
    EasyArrayDestroy( easyArray );
 
}

/* ---------------------------------------------------------------------- */
int EasyArraySearch( EASY_ARRAY* easyArray, int isSorted, const void *target, int (*compare)( const void *item1, const void *item2 ) ) {

    int                left, middle, right;


    left  = 0;
    right = easyArray->size - 1;


    /* Use binary search if the array is sorted. */
    if (isSorted) {


        /* Continue searching until the left and right indices cross. */
        while ( left <= right ) {

            middle = (left + right) / 2;

            switch ( compare( EasyArrayItem( easyArray, middle ), target ) ) {

                case -1:

                    /* Prepare to search to the right of the middle index. */
                    left = middle + 1;
                    break;


                case 1:

                    /* Prepare to search to the left of the middle index. */
                    right = middle - 1;
                    break;


                case 0:

                    /* Return the exact index where the data has been found. */
                    return middle;


                default:

                    /* This will happen if the compare function returned an illegal value. */
                    doCriticalErrorAndQuit("\nEasyArrayBiSearch: Illegal operation. Aborting!");
		            exit(1);
	       

            } /* end of switch */
    
        } /* end of search loop */


    /* If not sorted, use linear search. */
    } else {

        /* Continue searching until the left and right indices cross. */
        while ( left <= right ) {

            if ( 0 == compare( EasyArrayItem( easyArray, left ), target ) )
                return left;

            left++;
        } /* end of search loop */

    }


    /* Return that the data was not found. */
    return -1;

}

/* ---------------------------------------------------------------------- */
void* EasyArrayComputeMin( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) ) {

    int  size;
    void *smallest;
    void **src;


    /* If empty retun NULL */
    if ( 0 == (size = EasyArraySize(easyArray)) )
        return NULL;


    /* set the smallest to the first element */
    smallest = EasyArrayItem( easyArray, 0 );


    /* Return first value if the array is sorted or only 1 value is in the array. */
    if ( isSorted || 1 == size ) {
        
        return smallest;


    /* If not sorted, use linear search. */
    } else {


        /* search for min value */
        for ( src = (easyArray->array + 1); *src != NULL; src++ ) {

            if ( -1 == compare( *src, smallest ) )
                smallest = *src;

        } /* end of search loop */

    }


    /* return smallest */
    return smallest;

}

/* ---------------------------------------------------------------------- */
void* EasyArrayComputeMax( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) ) {

    int  size;
    void *largest;
    void **src;


    /* If empty retun NULL */
    if ( 0 == (size = EasyArraySize(easyArray)) )
        return NULL;


    /* set the largest to the last element */
    largest = EasyArrayItem( easyArray, size - 1 );


    /* Return first value if the array is sorted or only 1 value is in the array. */
    if ( isSorted || 1 == size ) {
        
        return largest;


    /* If not sorted, use linear search. */
    } else {


        /* search for min value */
        for ( src = (easyArray->array + 1); *src != NULL; src++ ) {

            if ( 1 == compare( *src, largest ) )
                largest = *src;

        } /* end of search loop */

    }


    /* return largest */
    return largest;

}

/* ---------------------------------------------------------------------- */
int EasyArrayComputeFrequency( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ), void *value ) {

    int count = 0, size;
    void **src;


    /* If empty return 0 */
    if ( 0 == (size = EasyArraySize(easyArray)) )
        return 0;


    /* search for min value */
    for ( src = easyArray->array; *src != NULL; src++ ) {

            if ( 0 == compare( *src, value ) )
                count++;

    } /* end of search loop */


    /* return the count of items equal to value */
    return count;
}

/* ---------------------------------------------------------------------- */
void* EasyArrayComputeMedian( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) ) {
    
    int size;
    void **result;


    /* If empty return 0 */
    if ( 0 == (size = EasyArraySize(easyArray)) )
        return NULL;


    /* If sorted, find the median value immediately. */
    if ( isSorted || 1 == size ) {

        result = (void**)EasyArrayItem( easyArray, size / 2 );

    /* If not sorted, create a copy of the array and sort it. */
    } else {

        EASY_ARRAY tempArray, *newEasyArray;

        /* create a copy of the array (but make sure it does not create 
         * new data and do not delete old data ) */
        tempArray.array = easyArray->array;
        tempArray.size = easyArray->size;
        newEasyArray = EasyArrayCopy( &tempArray, 0 );
        

        /* sort it */
        EasyArrayQuickSort( newEasyArray, compare );


        /* find the median value */
        result = (void**)EasyArrayItem( newEasyArray, size / 2 );


        /* delete the copied array */
        EasyArrayDestroy( newEasyArray );
    }

    
    /* return the middle value */
    return result;

}

/* ---------------------------------------------------------------------- */
void* EasyArrayComputeMode( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) ) {
    
    int max = 0, largestMax = 0, size;
    void **src;
    void *last, *mode;


    /* If empty return 0 */
    if ( 0 == (size = EasyArraySize(easyArray)) )
        return NULL;


    /* Return first value if only 1 value is in the array. */
    if ( 1 == size ) {
        return EasyArrayItem( easyArray, 0 );
    }


    /* If sorted, find the mode value immediately. */
    if ( isSorted ) {


        /* search for mode value */
        mode = last = EasyArrayItem( easyArray, 0 );
        for ( src = easyArray->array+1; *src != NULL; src++ ) {

                if ( 0 == compare( *src, last ) )
                    max++;
                else {
                    if ( max > largestMax ) {
                        largestMax = max;
                        mode = last;
                        
                        max = 0;
                    }
                }

                last = *src;

        } /* end of search loop */
        
        if ( max > largestMax ) mode = last;

    /* If not sorted, create a copy of the array and sort it. */
    } else {

        EASY_ARRAY tempArray, *newEasyArray;

        /* create a copy of the array (but make sure it does not create 
         * new data and do not delete old data ) */
        tempArray.array = easyArray->array;
        tempArray.size = easyArray->size;
        newEasyArray = EasyArrayCopy( &tempArray, 0 );
        

        /* sort it */
        EasyArrayQuickSort( newEasyArray, compare );


        /* search for mode value */
        mode = last = EasyArrayItem( newEasyArray, 0 );
        for ( src = newEasyArray->array+1; *src != NULL; src++ ) {

                if ( 0 == compare( *src, last ) )
                    max++;
                else {
                    if ( max > largestMax ) {
                        largestMax = max;
                        mode = last;
                        last = *src;
                        max = 0;
                    }
                }

        } /* end of search loop */

        if ( max > largestMax ) mode = last;


        /* delete the copied array */
        EasyArrayDestroy( newEasyArray );
    }

    
    /* return the mode */
    return mode;

}



/******************************* EasyStringHash *****************\
*
* Functions below let you maintain a hash that accepts a string as a 
* key. This string is used to generate a unique key based on some bit
* shifting function. Collisions are resolved via a linked list.
*
* Note: copy and destroy functions in the constructor are optional.
* If destroy is omitted, the data will not be freed.
* If copy is omitted, EasyHashCopy operation will fail.
*
* The initialSize specifies the initial size of the hash. A closest
* prime number is actually used instead. The hash dynamically grows
* and shrinks based on the fill factor.
*
* Note: If the set function is called with a NULL value for an item,
* the entry corresponding to the key will be removed, rather than
* being set to NULL. This way we don't need a separate function to
* cleanup empty slots later. Thefore, if the get function returns
* NULL, it means the entry does not exist in all cases.
*
\*****************************************************************/

#define EASY_LIFE_INIT_STRING_HASH_SIZE 14591

/* ---------------------------------------------------------------------- */
EASY_STRING_HASH* EasyStringHashCreate( int initialSize, void* (*copy)(const void *data), void (*destroy)(void *item) ) {

	EASY_STRING_HASH *newEasyStringHash;


    /* allocate the hash item */
	newEasyStringHash = (EASY_STRING_HASH*) calloc( 1, sizeof(EASY_STRING_HASH) );
	if ( NULL == newEasyStringHash ) {
		doCriticalErrorAndQuit("\nEasyStringHashCreate: Memory error. Aborting!");
		exit(1);
	}


    /* if zero size use default size */
    if( initialSize <= 0 ) 
        initialSize = EASY_LIFE_INIT_STRING_HASH_SIZE;
    else
        initialSize = EasyLifeGetPrime( initialSize );


    /* allocate the rack */
    newEasyStringHash->rack = (SHITEM**) calloc( initialSize, sizeof(SHITEM*) );
    if( NULL == newEasyStringHash->rack )
    {
		doCriticalErrorAndQuit("\nEasyStringHashCreate: Memory error. Aborting!");
		exit(1);
    }


    /* set other wars */
	newEasyStringHash->size = initialSize;
	newEasyStringHash->copy = copy;
	newEasyStringHash->destroy = destroy;

	return newEasyStringHash;
}

/* ---------------------------------------------------------------------- */
void removeSHITEM( EASY_STRING_HASH* easyStringHash, SHITEM* hitem ) {

    /* only remove item if the destroy function is provided */
    if ( NULL != easyStringHash->destroy )
        easyStringHash->destroy( hitem->item );

    /* free the key and the slot */
    free( hitem->key );
    free( hitem );
}

/* ---------------------------------------------------------------------- */
void EasyStringHashDestroy( EASY_STRING_HASH* easyStringHash )
{
    int i;
    SHITEM *curr;

    /* step trough and free any items in rack */
    for(  i = 0; i < easyStringHash->size; i++ )
    {
        /* free all items in that slot of the rack and slots themselves */
        while( easyStringHash->rack[i] )
        {
            curr = easyStringHash->rack[i];
            easyStringHash->rack[i] = curr->next;
            removeSHITEM( easyStringHash, curr );
        }
    }

    /* free rack and structure */
    free(easyStringHash->rack);
    free(easyStringHash);

    return;
}

/* ---------------------------------------------------------------------- */
void EasyStringHashSet( EASY_STRING_HASH* easyStringHash, char *key, void *item ) {
   
    SHITEM *hitem, *prev;
    int hkey;
    unsigned int ukey;


    /* generate an integer key from the string */
    ukey = EasyLifeElfUp(key);


    /* compute the hash key */
    hkey = ukey % easyStringHash->size;


    /* find the value */
    prev = NULL;
    for( hitem = easyStringHash->rack[hkey]; hitem != NULL;  ) {
        if( 0 == strcmp( key, hitem->key ))  break; 
        prev = hitem;
        hitem = hitem->next;
    }


    /* if incoming value is null, remove the hitem */
    if ( NULL == item ) {

        if ( NULL != hitem ) {

            if( NULL == prev )
            {
                easyStringHash->rack[hkey] = hitem->next;
            }
            else
            {
                prev->next = hitem->next;
            }
                
            removeSHITEM( easyStringHash, hitem );
        }

        return;
    }


    /* if hitem is not in list, create a new THITEM */
    if ( NULL == hitem ) {
        hitem = (SHITEM*) malloc(sizeof(SHITEM));
        if ( NULL == hitem ) {
		    doCriticalErrorAndQuit("\nEasyStringHashSet: Memory error. Aborting!");
		    exit(1);
        }
        hitem->next = easyStringHash->rack[hkey];
        easyStringHash->rack[hkey] = hitem;
        hitem->key = strdup(key);
    }


    /* set the value of the hitem */
    hitem->item = item; 

}

/* ---------------------------------------------------------------------- */
void* EasyStringHashGet( EASY_STRING_HASH* easyStringHash, char *key ) {

    SHITEM* hitem;
    int hkey;
    unsigned int ukey;


    /* generate an integer key from the string */
    ukey = EasyLifeElfUp(key);    


    /* compute the hash key */
    hkey = ukey % easyStringHash->size;


    /* get the value */
    for( hitem = easyStringHash->rack[hkey]; hitem != NULL; hitem = hitem->next )
    {
        if( 0 == strcmp( key, hitem->key )) break;
    }


    /* return it */
    return hitem ? hitem->item : NULL;

}

/* ---------------------------------------------------------------------- */
void EasyStringHashPrint( EASY_STRING_HASH* easyStringHash ) {

    SHITEM *hitem;
    int i;

    /* compute the hash key */
    for( i = 0; i < easyStringHash->size; i++ )
    {
        /* print rack id */
        fprintf(stderr,"\n%d :", i );

        for( hitem = easyStringHash->rack[i]; hitem != NULL; hitem = hitem->next )
        {
            fprintf(stderr,"  %s (%d)", hitem->key, hitem->item ); /* prints the address of the item, not item itself */
        }
    }

    return;

}

/******************************* EasyIntsHash ***************\
*
* Functions below let you maintain a hash that accepts one or more 
* integer keys. The keys are then written to a character buffer, 
* which is used as a key to the EasyStringHash. This is an example
* of how easy a hash based on multiple values can be created using
* the EASY_STRING_HASH.
*  
\*****************************************************************/

#define easy_ints_hash_convert_buf_size 1024

/* ---------------------------------------------------------------------- */
void EasyIntsHashSet( EASY_INTS_HASH* easyIntsHash, void *item, int key1, ... ) {

    /* increase this based on the the maximum allowed number of parameters */
    char text[easy_ints_hash_convert_buf_size] = ""; 

    char temp[35];
	va_list marker;
    int i, count=0;


    /* Initialize variable arguments. */
    i = key1;
    va_start( marker, key1 );           
   

    /* get and store all the keys */
    while( i != -1 )
    {
      sprintf( temp, "%d-", i );
      count+=strlen(temp);
      if ( count >= easy_ints_hash_convert_buf_size ) {
            doCriticalErrorAndQuit("\nEasyIntsHashSet: Buffer overflow. Aborting!");
		    exit(1);
      }
      strcat( text, temp );
      i = va_arg( marker, int );
    }

    /* Reset variable arguments. */
    va_end( marker );            


    /* set the item */
    EasyStringHashSet( easyIntsHash, text, item );

}

/* ---------------------------------------------------------------------- */
void* EasyIntsHashGet( EASY_INTS_HASH* easyIntsHash, int key1, ... ) {

    /* increase this based on the the maximum allowed number of parameters */
    char text[easy_ints_hash_convert_buf_size] = ""; 

    char temp[35];
	va_list marker;
    int i,count=0;


    /* Initialize variable arguments. */
    i = key1;
    va_start( marker, key1 );           
   

    /* get and store all the keys */
    while( i != -1 )
    {
      //printf("\n%d", i );
      sprintf( temp, "%d-", i );
      count+=strlen(temp);
      if ( count >= easy_ints_hash_convert_buf_size ) {
            doCriticalErrorAndQuit("\nEasyIntsHashSet: Buffer overflow. Aborting!");
		    exit(1);
      }
      strcat( text, temp );
      i = va_arg( marker, int );
    }
    
    /* Reset variable arguments. */
    va_end( marker );   

    
    /* set the item */
    return EasyStringHashGet( easyIntsHash, text );

}



/******************************* Internal ************************\
*
* These are various internal functions that are used by Easy Life: 
*
\*****************************************************************/




/******************************************************************
 * This function produces an integer value by shifting the bits
 * of the passed in string. Used for hashing.
 *
 *******************************************************************/
unsigned int EasyLifeElfUp(char *text)
{

    unsigned int h=0,g;

    while(*text)
    {
        h = (h<<4)+ *text;
        if(g = h & 0xF0000000)
            h ^= g >>24;
        h &= ~g;
        text++;
    }

    return h;
}

/******************************************************************
* the following list of primes is used by the getprime routine. It
* starts at 1103 and continues with primes that are about 20%
* greater than the previous one up to roughly 10 million.
* Read the description of EasyLifeGetPrime for more information.
*******************************************************************/
int EasyLifeSomePrimes[] =
{
    1103,
    1327,
    1597,
    1931,
    2333,
    2801,
    3371,
    4049,
    4861,
    5839,
    7013,
    8419,
    10103,
    12143,
    14591,
    17519,
    21023,
    25229,
    30293,
    36353,
    43627,
    52361,
    62851,
    75431,
    90523,
    108631,
    130363,
    156437,
    187751,
    225307,
    270371,
    324449,
    389357,
    467237,
    560689,
    672827,
    807403,
    968897,
    1162687,
    1395263,
    1674319,
    2009069,
    2410897,
    2893081,
    3471697,
    4166047,
    4999273,
    5999129,
    7198963,
    8638807,
    10366597,
    0       /* last must be zero to exit loop */
};

/******************************************************
* Seach for prime greater than or equal to the size
* requested. If size is too large returns the largest
* prime in the list
*******************************************************/
int EasyLifeGetPrime(int number)
{
    int *p;

    if(number<EasyLifeSomePrimes[0]) return EasyLifeSomePrimes[0];

    for(p = EasyLifeSomePrimes;*p!=0;p++)
    {
        if(*p>=number) return *p;
    }
    return *(p-1);
}



/******************************* SimpleArray ******************************\
*
* This is just a wrapper around regular ANSI funcs to maintain a regular
* array. In debug mode (EASY_LIFE_DEBUG defined), they map to special 
* bound checking functions, other time they map to C arrays with type
* checking  and allocation checking  but no bound checking.
*
\**************************************************************************/

static EASY_INTS_HASH* __EasyLifeSimpleArrayAddressHash = NULL;
static long long __EasyLifeSimpleArrayAllocations=0;
static long long __EasyLifeSimpleArrayBytes=0;

typedef struct {
    int elements;
    int esize;
    int line;
    char *file;
} __EL_SIMPLE_ARRAY_DEBUG_INFO;


void EsadiFree(void *mem) {

    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi=(__EL_SIMPLE_ARRAY_DEBUG_INFO*)mem;;

    if (NULL!=mem&&esadi->file) {
        free(esadi->file);
        esadi->file=NULL;
    }

    free(mem);
}

/* ---------------------------------------------------------------------- */
void* simpleArrayCreateDebug( int ELEMENTS, int ESIZE, char *file, int line ) {

	void *array;
    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi;


    /* allocate a void* array */
	array = (void*) calloc( ELEMENTS, ESIZE );
	if ( NULL == array ) {
		doCriticalErrorAndQuit("\nsimpleArrayCreateDebug: Memory error. Aborting!");
		exit(1);
	}


    /* check if address hash exists */
    if ( !__EasyLifeSimpleArrayAddressHash ) {

        __EasyLifeSimpleArrayAddressHash = EasyIntsHashCreate( 10000000, NULL, EsadiFree );
    }


    /* remember array size and esize */
    esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)malloc( sizeof(__EL_SIMPLE_ARRAY_DEBUG_INFO) );
	if ( NULL == esadi ) {
		doCriticalErrorAndQuit("\nsimpleArrayCreateDebug: Memory error on line %d in file %s. Aborting!", line, file);
		exit(1);
	}
    esadi->elements = ELEMENTS;
    esadi->esize = ESIZE;
    esadi->line = line;
    esadi->file = strdup(file);
	if ( NULL == (esadi->file) ) {
		doCriticalErrorAndQuit("\nsimpleArrayCreateDebug: Memory error on line %d in file %s. Aborting!", line, file);
		exit(1);
	}


    /* check if this address is not already used, if used, that means
    the user has not delocated it. Bad user!!! Bum out. */
    if (EasyIntsHashGet( __EasyLifeSimpleArrayAddressHash, (int)array, -1 )!=NULL) {
		doCriticalErrorAndQuit("\nsimpleArrayCreateDebug: Address %d already allocated. Line %d in file %s. You have a block somewhere in your code not delocated properly (make sure to use sfree for scalloc and saDestroy for SaCreate, etc. Aborting!", line, file);
		exit(1);
    }


    /* use the address of the array as the hash key */
    EasyIntsHashSet( __EasyLifeSimpleArrayAddressHash, esadi, (int)array, -1 ); 


    /* remember number of allocations and bytes */
    __EasyLifeSimpleArrayAllocations+=1;
    __EasyLifeSimpleArrayBytes+=(((long long)ELEMENTS)*((long long)ESIZE));

    return array;
}

/* ---------------------------------------------------------------------- */
void* simpleArrayCreateRelease( int ELEMENTS, int ESIZE ) {

	void *array;

    /* allocate a void* array */
	array = (void*) calloc( ELEMENTS, ESIZE );
	if ( NULL == array ) {
		doCriticalErrorAndQuit("\nsimpleArrayCreateRelease: Memory error. Aborting!");
		exit(1);
	}

    return array;
}

/* ---------------------------------------------------------------------- */
char* simpleArrayAtDebug( void* inArray , int POS ) {

    int elements, esize;
    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi;
    char* ARRAY = (char*)inArray;

    //printf("got to simpleArrayAtDebug!!!array=%d\n", (int)inArray );

    /* get the array debug info */
    if ( NULL == __EasyLifeSimpleArrayAddressHash || NULL == ( esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)EasyIntsHashGet( __EasyLifeSimpleArrayAddressHash, (int)ARRAY, -1 ) ) ) {
        
        //if (__EasyLifeSimpleArrayAddressHash)
        //    EasyIntsHashPrint( __EasyLifeSimpleArrayAddressHash );

        doCriticalErrorAndQuit("\nsimpleArrayAtDebug: Could not lookup debug info for address %d! Aborting!", ARRAY );        
        exit(1);
    }
    elements = esadi->elements;
    esize = esadi->esize;


    /* check bounds */
	if ( POS >= elements || POS < 0  ) {
		doCriticalErrorAndQuit("\nsimpleArrayAtDebug: Illegal memory access! Aborting!\n"
            "\tblock start: %d\n "
            "\tblock end: %d\n "
            "\tesize: %d\n "
            "\tarraysize: %d\n"
            "\tposition accessed: %d\n "
            "\taddress range accessed: %d-%d\n ", ARRAY, ARRAY + elements * esize - 1, esize, elements, POS, ARRAY + POS * esize, ARRAY + ( POS + 1 ) * esize - 1);
		exit(1);
	}


    /* return the position asked for */
    return ARRAY + POS * esize;

}

/* ---------------------------------------------------------------------- */
void simpleArrayDestroyDebug( void* ARRAY, char *file, int line ) {

    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi;

    /* check if debug info hash exists */
    if ( NULL == __EasyLifeSimpleArrayAddressHash || NULL == (esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)EasyIntsHashGet( __EasyLifeSimpleArrayAddressHash, (int)ARRAY, -1 )) ) {
        
        //if (__EasyLifeSimpleArrayAddressHash)
        //    EasyIntsHashPrint( __EasyLifeSimpleArrayAddressHash );

        doCriticalErrorAndQuit("\nsimpleArrayDestroyDebug: Address %d is not present in the address hash table! Line %d in file %s. Aborting!", ARRAY, line, file );
        exit(1);
    }

    /* clear the memory to prevent zombie objects from using it (hopefully) */
    memset( ARRAY, 0, (esadi->elements) * (esadi->esize) );

    /* clear the address */
    EasyIntsHashSet( __EasyLifeSimpleArrayAddressHash, NULL, (int)ARRAY, -1 ); 

    /* free the memory */
    free( ARRAY );

}

/* ---------------------------------------------------------------------- */
void simpleMemorySummaryDebug( int printUnfreed ) {

    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi;
    SHITEM *hitem;
    int i; 
    long long bytesLeft = 0, allocationsLeft = 0;

    /* check if debug info hash exists */
    if ( NULL == __EasyLifeSimpleArrayAddressHash) {
        
        doCriticalErrorAndQuit("\nsimpleMemorySummaryDebug: No hash table present! Aborting!" );
        exit(1);
    }

    /* compute the hash key */
    for( i = 0; i < __EasyLifeSimpleArrayAddressHash->size; i++ )
    {

        for( hitem = __EasyLifeSimpleArrayAddressHash->rack[i]; hitem != NULL; hitem = hitem->next )
        {
            esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)hitem->item; 
            bytesLeft+=( ((long long)(esadi->elements))*((long long)(esadi->esize)) );
            allocationsLeft+=1;
        }
    }

    /* print */
    fprintf(stderr,""
        "total allocations: \t\t\t%I64d\n"
        "allocations not deleted: \t\t%I64d\n"
        "total bytes allocated: \t\t\t%I64d\n"
        "bytes unallocated: \t\t\t%I64d\n",
        __EasyLifeSimpleArrayAllocations, allocationsLeft, __EasyLifeSimpleArrayBytes, bytesLeft );


    /* print unfreed allocations */
    if (printUnfreed){

        fprintf(stderr,"\n\nUNFREED ALLOCATIONS:\n");
        for( i = 0; i < __EasyLifeSimpleArrayAddressHash->size; i++ )
        {

            for( hitem = __EasyLifeSimpleArrayAddressHash->rack[i]; hitem != NULL; hitem = hitem->next )
            {
                esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)hitem->item; 
                fprintf(stderr,"Line %d in file \"%s\"\n",esadi->line,esadi->file);
            }
        }        

    } else fprintf(stderr,"Pass 1 to this function to view the unfreed allocations.\n"); 
}

/* ---------------------------------------------------------------------- */
void simpleMemorySummaryRelease( void ) {

    doCriticalErrorAndQuit("\nsimpleMemorySummaryRelease: This operation is not supported in the release mode. Aborting!");
    exit(1);
}

/******************************* SimpleMatrix ****************************\
*
* This is just a wrapper around regular ANSI funcs to maintain a regular
* matrix. In debug mode (EASY_LIFE_DEBUG defined), they map to special 
* bound checking functions, other time they map to C arrays with type
* checking and allocation checking but no bound checking.
*
\*************************************************************************/

/* ---------------------------------------------------------------------- */
void* simpleMatrixCreateDebug( int ROWS, int COLS, int ESIZE ) {

    int i;
	char *matrix;
    void **rows;
    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi;


    /* allocate a void* array */
	matrix = (char*)calloc( ROWS * COLS, ESIZE );
	rows = (void**)calloc( ROWS, sizeof(void*) );
	if ( NULL == matrix || NULL == rows ) {
		doCriticalErrorAndQuit("\nsmCreate: Memory error. Aborting!");
		exit(1);
	}


    /* check if address hash exists */
    if ( !__EasyLifeSimpleArrayAddressHash ) {

        __EasyLifeSimpleArrayAddressHash = EasyIntsHashCreate( 50000, NULL, free );
    }


    /* remember address of the rows array */
	if ( NULL == (esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)malloc( sizeof(__EL_SIMPLE_ARRAY_DEBUG_INFO))) ) {
		doCriticalErrorAndQuit("\nsimpleArrayCreateDebug: Memory error. Aborting!");
		exit(1);
	}
    esadi->elements = ROWS;
    esadi->esize = sizeof(void*);
    EasyIntsHashSet( __EasyLifeSimpleArrayAddressHash, esadi, (int)rows, -1 ); 
    //printf("rows: %d\n", (int)rows );


    /* assign rows */
    for ( i = 0; i < ROWS; i++ ) {

            /* link the rows to the matrix */
            rows[i] = matrix + i * COLS * ESIZE;

            /* remember array size and esize and store in the address hash */
	        if ( NULL == (esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)malloc( sizeof(__EL_SIMPLE_ARRAY_DEBUG_INFO))) ) {
		        doCriticalErrorAndQuit("\nsimpleArrayCreateDebug: Memory error. Aborting!");
		        exit(1);
	        }
            esadi->elements = COLS;
            esadi->esize = ESIZE;
            EasyIntsHashSet( __EasyLifeSimpleArrayAddressHash, esadi, (int)(matrix + i * COLS * ESIZE ), -1 ); 
            //printf("matrix: %d\n", (int)(matrix + i * COLS * ESIZE ) );
    }


    return rows;
}

/* ---------------------------------------------------------------------- */
void* simpleMatrixCreateRelease( int ROWS, int COLS, int ESIZE ) {

    int i;
	char *matrix;
    void **rows;


    /* allocate a void* array */
	matrix = (char*)calloc( ROWS * COLS, ESIZE );
	rows = (void**)calloc( ROWS, sizeof(void*) );
	if ( NULL == matrix || NULL == rows ) {
		doCriticalErrorAndQuit("\nsmCreate: Memory error. Aborting!");
		exit(1);
	}


    /* assign rows */
    for ( i = 0; i < ROWS; i++ ) {

            /* link the rows to the matrix */
            rows[i] = matrix + i * COLS * ESIZE;
    }


    return rows;
}

/* ---------------------------------------------------------------------- */
void simpleMatrixDestroyDebug( void* MATRIX ) {

    void **rows = (void**)MATRIX;
    __EL_SIMPLE_ARRAY_DEBUG_INFO *esadi;
    int i, ROWS;


    /* check if debug info hash exists */
    if ( NULL == __EasyLifeSimpleArrayAddressHash || NULL == ( esadi = (__EL_SIMPLE_ARRAY_DEBUG_INFO*)EasyIntsHashGet( __EasyLifeSimpleArrayAddressHash, (int)MATRIX, -1 ) ) ) {

        doCriticalErrorAndQuit("\nsimpleMatrixDestroyDebug: Address %d is not present in the address hash table! Aborting!", MATRIX );
        exit(1);
    }
    ROWS = esadi->elements;


    /* Remove all entries from the address hash for this matrix     */
    /* Note: We are not removing the first entry because it will    */
    /* be removed by simpleArrayDestroyDebug( rows[0] ) which is    */
    /* called below.                                                */
    for ( i = 1; i < ROWS; i++ ) {

            EasyIntsHashSet( __EasyLifeSimpleArrayAddressHash, NULL, (int)rows[i], -1 ); 

    }


    /* delete matrix */
    simpleArrayDestroyDebug( rows[0],  __FILE__, __LINE__ );


    /* delete rows */
    simpleArrayDestroyDebug( rows,  __FILE__, __LINE__ );

}

/* ---------------------------------------------------------------------- */
void simpleMatrixDestroyRelease( void* MATRIX ) {

    void **rows = (void**)MATRIX;


    /* delete matrix */
    free( rows[0] );


    /* delete rows */
    free( rows );

}
