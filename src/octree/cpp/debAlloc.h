#ifdef DEBUG_MEMORY

#define MEM_VERBOSE_NO	0
#define MEM_VERBOSE_LOG	1
#define MEM_VERBOSE_ALL	2


/* FABIO: Debugging helpers... */
extern "C" {
void *MyMalloc(size_t, char *, int);
void *MyCalloc(size_t, size_t, char *, int);
void *MyRealloc(void *, size_t, char *, int);
void MyFree(void *, char *, int);

/* Helper funcs */
void MyAllocStat(const char *optFileOut);
void MyCheckAllocInUse(const char *optFileOut , long minIndex );
void MyAllocSet(int verboseAlloc, int verboseFree);
void MyAllocNoCheckAtExit();
}

#define malloc(a)	MyMalloc(a,  __FILE__,__LINE__)
#define calloc(a,b)	MyCalloc(a,b,__FILE__,__LINE__)
#define free(ptr)	MyFree  (ptr,__FILE__,__LINE__)
#define realloc(ptr,sz)	MyRealloc(ptr, sz, __FILE__,__LINE__)

/* Shortways for typed allocs */
#define MYNEWC(type)			( (type *) malloc(sizeof(type)) )
#define MYNEWCV(type, count)	( (type *) malloc(sizeof(type)*(count)) )

#endif	// DEBUG_MEMORY
