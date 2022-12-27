#ifndef _EASY_LIFE_H
#define _EASY_LIFE_H

//undefine this to execute in the debug mode!
//#define EASY_LIFE_DEBUG

typedef struct {
    char            *lastTBFDataPos;
} BSTREAM;


typedef struct {
    void            **array;
    int             size;
    void*           ( *copy )( const void *item );
    void            ( *destroy )( void *item );
    int             reserved; /* not to be used by the user */
} EASY_ARRAY;


typedef struct tagEASY_NODE {
    void                    *item;
    struct tagEASY_NODE     *next;
    struct tagEASY_NODE     *prev;
} EASY_NODE;


typedef struct {
    int             size;
    EASY_NODE       *head;
    EASY_NODE       *tail;
    void*           ( *copy )( const void *item );
    void            ( *destroy )( void *item );
} EASY_LIST;


typedef struct tagSHITEM {
    char            *key;
    void            *item;
    struct tagSHITEM *next;
} SHITEM;


typedef struct
{
    int             size;
    void*           ( *copy )( const void *item );
    void            ( *destroy )( void *item );
    SHITEM**        rack;
} EASY_STRING_HASH;

/******************************* StreamLibs ******************************\
*
* The functions below provide C++ like character stream capabilites
*
* Note: this input source for bopen is a '\0' terminated string 
*
\*************************************************************************/


BSTREAM*             bopen                                      ( char *data );

void                 bclose                                     ( BSTREAM *bp );

char*                bgets                                      ( char *buffer, int lineLength, BSTREAM *bp );



/******************************* EasyArray *******************************\
*
* Functions below let you maintain a dynamic array (doubles each time)
*
* Note: the last item of the array points to NULL
*
* Note: the size of an array is the number of items not including
* the item pointing to NULL
*
* Note: copy and destroy functions in the constructor are optional.
* If destroy is omitted, the data will not be freed.
* If copy is omitted, EasyArrayCopy operation will fail.
*
\*************************************************************************/


EASY_ARRAY*             EasyArrayCreate                         ( int initialSize, void* (*copy)(const void *data), void (*destroy)(void *item) );

EASY_ARRAY*             EasyArrayCreateFromList                 ( EASY_LIST* easyList, int copyData ); 

EASY_ARRAY*             EasyArrayCopy                           ( EASY_ARRAY* easyArray, int copyData );

void                    EasyArrayDestroy                        ( EASY_ARRAY* easyArray );

void                    EasyArrayInsert                         ( EASY_ARRAY* easyArray, void *item );

void                    EasyArrayRemove                         ( EASY_ARRAY* easyArray, int position );

void*                   EasyArrayItem                           ( EASY_ARRAY* easyArray, int position );

#define                 EasyArraySize( easyArray  )             ( (easyArray)->size )



/************************ Sorting/Searching ******************************\
*                                                                            
* Various array sorting and searching algorithms 
*
* Some conventions:
*
*    when algorithm returns an index, -1 means the search value not found
*                                                                            
**************************************************************************/


void					EasyArrayInsertionSort					( EASY_ARRAY* easyArray, int (*compare)( const void *item1, const void *item2 ) );

void					EasyArrayQuickSort						( EASY_ARRAY* easyArray, int (*compare)( const void *item1, const void *item2 ) );

int						EasyArraySearch					    	( EASY_ARRAY* easyArray, int isSorted, const void *target, int (*compare)( const void *item1, const void *item2 ) ); 

int                     EasyArrayComputeFrequency               ( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ), void *value );	

void*                   EasyArrayComputeMin                     ( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) );	

void*                   EasyArrayComputeMax                     ( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) );	

void*                   EasyArrayComputeMode                    ( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) );

void*                   EasyArrayComputeMedian                  ( EASY_ARRAY* easyArray, int isSorted, int (*compare)( const void *item1, const void *item2 ) );                       



/******************************* EasyList ********************************\
*
* The functions below let you maintain a dynamic doubly linked list
*
* Note: copy and destroy functions in the constructor are optional.
* If destroy is omitted, the data will not be freed.
* If copy is omitted, EasyListCopy operation will fail.
*
\*************************************************************************/


EASY_LIST*              EasyListCreate                          ( void* (*copy)(const void *data), void (*destroy)(void *item) );

EASY_LIST*              EasyListCopy                            ( EASY_LIST* easyList, int copyData );

void                    EasyListAppend                          ( EASY_LIST* easyListTo, EASY_LIST* easyListFrom, int copyData );

void                    EasyListDestroy                         ( EASY_LIST* easyList );

#define                 EasyListSize( easyList )                ( (easyList)->size )

#define                 EasyListIsEmpty( easyList    )          ( (easyList)->size == 0 )

#define                 EasyListHead( easyList )                ( (easyList)->head )

#define                 EasyListTail( easyList )                ( (easyList)->tail )

#define                 EasyListIsHead( easyList , easyNode )   ( (easyNode) == (easyList)->head ? 1 : 0 )

#define                 EasyListIsTail( easyNode )              ( (easyNode)->next == NULL ? 1 : 0 )

#define                 EasyListItem( easyNode )                ( (easyNode)->item )

#define                 EasyListNext( easyNode )                ( (easyNode)->next )    

#define                 EasyListPrevious( easyNode )            ( (easyNode)->prev )


void                    EasyListInsertAfter                     ( EASY_LIST* easyList, EASY_NODE* easyNode, void *item );

#define                 EasyListInsertBefore                    ( EasyListInsertAfter( easyList, NULL != (easyNode) ? ( (easyNode)->prev ) : NULL , item ) )

#define                 EasyListInsertHead( easyList, item )    ( EasyListInsertAfter( easyList, NULL, item ) )

#define                 EasyListInsertTail( easyList, item )    ( EasyListInsertAfter( easyList, ( (easyList)->tail ) , item ) )    


void                    EasyListRemoveNode                      ( EASY_LIST* easyList, EASY_NODE* easyNode );

#define                 EasyListRemoveHead( easyList )          ( EasyListRemoveNode( easyList, NULL ) )    

#define                 EasyListRemoveTail( easyList )          ( EasyListRemoveNode( easyList, ( (easyList)->tail ) ) )            



/************************ Sorting/Searching ******************************\
*                                                                            
* Various link list sorting and searching algorithms                                       
*                                                                            
**************************************************************************/


void EasyListInsertionSort( EASY_LIST* easyList, int (*compare)( const void *item1, const void *item2 ) );

void EasyListQuickSort( EASY_LIST* easyList, int (*compare)( const void *item1, const void *item2 ) );


            
/******************************* EasyStack *******************************\
*                                                                            
*  Implement stacks as linked lists.                                         
*                                                                            
**************************************************************************/


typedef EASY_LIST EASY_STACK;


#define                 EasyStackCreate                         ( EasyListCreate )

#define                 EasyStackDestroy                        ( EasyListDestroy )

#define                 EasyStackPush( easyStack, item )        ( EasyListInsertAfter( easyStack, NULL, item ) )

#define                 EasyStackPop( easyStack )               ( EasyListRemoveNode( easyStack, NULL ) )

#define                 EasyStackPeek( easyStack )              ( (easyStack)->head == NULL ? NULL : (easyStack)->head->item )

#define                 EasyStackSize( easyStack )              ( (easyStack)->size )


            
/******************************* EasyQueue *******************************\
*                                                                            
*  Implement queues as linked lists.                                         
*                                                                            
**************************************************************************/


typedef EASY_LIST EASY_QUEUE;


#define                 EasyQueueCreate                         ( EasyListCreate )

#define                 EasyQueueDestroy                        ( EasyListDestroy )

#define                 EasyQueueInsert( easyQueue, item )      ( EasyListInsertAfter( easyQueue, NULL, item ) )

#define                 EasyQueueRemove( easyQueue )            ( EasyListRemoveNode( easyQueue, ( (easyQueue)->tail ) ) )

#define                 EasyQueuePeek( easyQueue )              ( (easyQueue)->tail == NULL ? NULL : (easyQueue)->tail->item )

#define                 EasyQueueSize( easyQueue )              ( (easyQueue)->size )



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
* and shrinks based on the fill factor (well not yet :)
*
* Note: If the set function is called with a NULL value for an item,
* the entry corresponding to the key will be removed, rather than
* being set to NULL. This way we don't need a separate function to
* cleanup empty slots later. Thefore, if the get function returns
* NULL, it means the entry does not exist in all cases.
*
\*****************************************************************/


EASY_STRING_HASH*       EasyStringHashCreate                    ( int initialSize, void* (*copy)(const void *data), void (*destroy)(void *item) );

void                    EasyStringHashDestroy                   ( EASY_STRING_HASH* easyStringHash );

void                    EasyStringHashSet                       ( EASY_STRING_HASH* easyStringHash, char *key, void *item );

void*                   EasyStringHashGet                       ( EASY_STRING_HASH* easyStringHash, char *key );

void                    EasyStringHashPrint                     ( EASY_STRING_HASH* easyStringHash );

#define                 EasyStringHashSize( easyStringHash  )   ( (easyStringHash)->size )



/******************************* EasyIntsHash ***************\
*
* Functions below let you maintain a hash that accepts one or more 
* integer keys. The keys are then written to a character buffer, 
* which is used as a key to the EasyStringHash. This is an example
* of how easy a hash based on multiple values can be created using
* the EASY_STRING_HASH.
*  
\*****************************************************************/


typedef EASY_STRING_HASH EASY_INTS_HASH;


#define                 EasyIntsHashCreate                      ( EasyStringHashCreate )

#define                 EasyIntsHashDestroy                     ( EasyStringHashDestroy )

/* USE -1 as terminator, or else!!! */
void                    EasyIntsHashSet                         ( EASY_INTS_HASH* easyIntsHash, void *item, int key1, ... );

/* USE -1 as terminator, or else!!! */
void*                   EasyIntsHashGet                         ( EASY_INTS_HASH* easyStringHash, int key1, ... );

#define                 EasyIntsHashPrint                       ( EasyStringHashPrint )

#define                 EasyIntsHashSize( easyIntsHash  )       ( (easyIntsHash)->size )



/******************************* SimpleArray *****************************\
*
* This is just a wrapper around regular ANSI funcs to maintain a regular
* array. In debug mode (EASY_LIFE_DEBUG defined), they map to special 
* bound checking functions, other time they map to C arrays with type
* checking and allocation checking but no bound checking.
*
\*************************************************************************/


void*                   simpleArrayCreateDebug                  ( int ELEMENTS, int ESIZE, char *file, int line );

void*                   simpleArrayCreateRelease                ( int ELEMENTS, int ESIZE );

char*                   simpleArrayAtDebug                      ( void* ARRAY, int POS );

void                    simpleArrayDestroyDebug                 ( void* ARRAY, char *file, int line );


#ifdef EASY_LIFE_DEBUG

    #define             saCreate( ELEMENTS, ESIZE )             ( simpleArrayCreateDebug( ELEMENTS, ESIZE, __FILE__, __LINE__ ) )

    #define             saAt( ARRAY , POS, TYPE )               ( *(TYPE*)simpleArrayAtDebug( ARRAY , POS ) )

    #define             saAtPtr( ARRAY , POS, TYPE )            ( (TYPE*)simpleArrayAtDebug( ARRAY , POS ) )

    #define             saDestroy( ARRAY )                      ( simpleArrayDestroyDebug( ARRAY , __FILE__, __LINE__ ) )

#else

    #define             saCreate( ELEMENTS, ESIZE )             ( simpleArrayCreateRelease( ELEMENTS, ESIZE ) )

    #define             saAt( ARRAY , POS, TYPE )               ( ARRAY[POS] )

    #define             saAtPtr( ARRAY , POS, TYPE )            ( ARRAY + POS )

    #define             saDestroy( ARRAY )                      ( free( ARRAY ) )

#endif



/******************************* SimpleMemory *****************************\
*
* Keeps track of your memory allocations. Does not let you Free wrong memory.
* Allows you to view if you have unallocated stuff at the end and other 
* statistics. In release mode maps to basic C functions but allocation
* checking is still performed.
*
\*************************************************************************/


void                    simpleMemorySummaryDebug                ( int printUnfreed );

void                    simpleMemorySummaryRelease              ( void );


#ifdef EASY_LIFE_DEBUG

    #define             scalloc( ELEMENTS, ESIZE )              ( simpleArrayCreateDebug( ELEMENTS, ESIZE,  __FILE__, __LINE__ ) )

    #define             smalloc( ESIZE )                        ( simpleArrayCreateDebug( 1, ESIZE,  __FILE__, __LINE__ ) )

    #define             sfree( MEMORY )                         ( simpleArrayDestroyDebug( MEMORY,  __FILE__, __LINE__ ) )

    #define             sMemorySummary( PRINT_ALLOCATIONS )     ( simpleMemorySummaryDebug(PRINT_ALLOCATIONS) )

#else

    #define             scalloc( ELEMENTS, ESIZE )              ( simpleArrayCreateRelease( ELEMENTS, ESIZE ) )
//    #define             scalloc( ELEMENTS, ESIZE )              ( calloc( ELEMENTS, ESIZE ) )

    #define             smalloc( ESIZE )                        ( simpleArrayCreateRelease( 1, ESIZE ) )
//    #define             smalloc( ESIZE )                        ( malloc( ESIZE ) )

    #define             sfree( MEMORY )                         ( free( MEMORY ) )

    #define             sMemorySummary( PRINT_ALLOCATIONS )     ( simpleMemorySummaryRelease() )

#endif



/******************************* SimpleMatrix ****************************\
*
* This is just a wrapper around regular ANSI funcs to maintain a regular
* matrix. In debug mode (EASY_LIFE_DEBUG defined), they map to special 
* bound checking functions, other time they map to C arrays with type
* checking and allocation checking  but no bound checking.
*
\*************************************************************************/


void*                   simpleMatrixCreateDebug                 ( int ROWS, int COLS, int ESIZE );

void*                   simpleMatrixCreateRelease               ( int ROWS, int COLS, int ESIZE );

void                    simpleMatrixDestroyDebug                ( void* MATRIX );

void                    simpleMatrixDestroyRelease              ( void* MATRIX );


#ifdef EASY_LIFE_DEBUG

    #define             smCreate( ROWS, COLS, ESIZE )           ( simpleMatrixCreateDebug( ROWS, COLS, ESIZE ) )

    #define             smAt( MATRIX , ROW , COL, TYPE )        ( *(TYPE*)simpleArrayAtDebug( *(TYPE**)simpleArrayAtDebug( MATRIX , ROW ) , COL ) )

    #define             smAtPtr( MATRIX , ROW , COL, TYPE )     ( (TYPE*)simpleArrayAtDebug( *(TYPE**)simpleArrayAtDebug( MATRIX , ROW ) , COL ) )

    #define             smDestroy( MATRIX )                     ( simpleMatrixDestroyDebug( MATRIX ) )

#else

    #define             smCreate( ROWS, COLS, ESIZE )           ( simpleMatrixCreateRelease( ROWS, COLS, ESIZE ) )

    #define             smAt( MATRIX , ROW , COL, TYPE )        ( MATRIX[ ROW ][ COL ] )
    
    #define             smAtPtr( MATRIX , ROW , COL, TYPE )     ( MATRIX[ ROW ] + COL )

    #define             smDestroy( MATRIX )                     ( simpleMatrixDestroyRelease( MATRIX ) )

#endif



/******************************* doCriticalErrorAndQuit ******************\
*
* This function must be defined by your application to provide correct
* program flow.
*
\*************************************************************************/
extern void doCriticalErrorAndQuit(const char *format, ... );

#endif