#include <stdio.h>
#include <stdlib.h>

/* FABIO: Debugging helpers... */
#include <assert.h>

#define LISTSIZE	1024*1024
static long byCount=0;
static long verboseMalloc = 0;	/* If not set, no debugging output on malloc */
static long verboseFree  = 0;	/* If not set, no debugging output on free   */
static long listCount =0;
static void *ToBeCleared[LISTSIZE];

#define headStr "CHECK123456HEAD"
#define HEADSZ  sizeof(headStr)
#define tailStr "CHECK123456TAIL"
#define TAILSZ  sizeof(tailStr)

typedef struct {
	char *fileName;
	int fileLineNum, posInList, nameLen;
	double pad;
	size_t origAmt;		/* More restrictive type for alignment */
	char headChk[HEADSZ];
}	myStructPre;


// Seriale me==0
// In MPI, me == rank
int me=0;
void *MyMalloc(size_t amt, char *file, int line)
{
	char *ptr, *p = NULL;
	size_t newA, nameLen, preSz;
	myStructPre *theAlloc;


	byCount += amt;


	nameLen = strlen(file)+1;
	if (nameLen % 8 != 0)
		nameLen = ((nameLen / 8) + 1)*8;

	preSz = sizeof(myStructPre) + nameLen;
	newA = preSz + amt + TAILSZ;

	ptr = (char *) malloc(newA);
	if (ptr != NULL){
		strcpy(ptr, file);						/* File Name */
		p = ptr + nameLen;


		theAlloc = (myStructPre *) p;

		theAlloc->fileName = ptr;
		theAlloc->fileLineNum = line;			/* File Line num */
		theAlloc->origAmt = amt;				/* Amt req'ed */
		theAlloc->nameLen = nameLen;	/* Header Size */
		strcpy(theAlloc->headChk, headStr);		/* Head Check zone */

		p += sizeof(myStructPre);

		if (listCount < LISTSIZE) {				/* Pos in in-use-list */
			ToBeCleared[listCount] = p;
			theAlloc->posInList = listCount++;
		} else {
			theAlloc->posInList = LISTSIZE + 1;
		}


		strcpy(p + amt, tailStr);				/* Tail Check zone */

		if (verboseMalloc)
		if (me<10)
			fprintf(stderr, "MyMalloc) me:%d %s:%d - Got a ptr(=%020ld) of %d bs.\n", me ,file, line, listCount-1, amt);

	}


	return (void *)p;
}
void *MyCalloc(size_t a, size_t b, char *file, int line)
{
	char *p;
	size_t newA = a*b;


	p = (void *) MyMalloc(newA, file, line);
	memset(p, 0, newA);

	return (void *)p;
}
static char *MyDecodePtr(void *ptr, char*msg, char*fromFile, int fromLine,
					int *res, char **fileAlloc,
					int *lineAlloc, long *listPos, size_t *amt)
{
	char *p, *tailPos;
	myStructPre *theAlloc;
	int res1;


	p = (char*)ptr - sizeof(myStructPre);
	theAlloc = (myStructPre *) p;

	res1 = 10 * (strcmp(theAlloc->headChk, headStr) != 0);

	*fileAlloc = theAlloc->fileName;			/* File Name */
	*lineAlloc = theAlloc->fileLineNum;			/* File Line num */
	*amt = theAlloc->origAmt;					/* Amt req'ed */

	*listPos   = theAlloc->posInList;			/* Pos in in-use-list */
	if (*listPos == LISTSIZE + 1)
		*listPos = -1;


	tailPos = (char*)ptr + *amt;
	res1 += strcmp(tailPos , tailStr) != 0;
	*res = res1;

	return (p - theAlloc->nameLen);
}
const char *ExpandCheckCode(int res)
{
		static const char *statusCode;

		switch (res) {
			case 1:
				statusCode = "After";
				break;
			case 10:
				statusCode = "Before";
				break;
			case 11:
				statusCode = "Before+After";
				break;
			default:
				statusCode = "OK";
		}
		return statusCode;
}
void MyFree(void *ptr, char *file, int line)
{
	size_t amt;
	long	listPos;
	char *p, *fileAlloc;
	int lineAlloc;
	int res;


	if (ptr == NULL){
		if (1 || verboseFree)
			fprintf(stderr, "MyFree) me:%d %s:%d - Got NULL.\n", me, file, line);
		return;
	}

	p = MyDecodePtr(ptr, "MyFree", file, line, &res, &fileAlloc, &lineAlloc, &listPos, &amt);
	if (res != 0) {
		fprintf(stderr, "MyFree) me:%d Status:%s %s:%d - Bad ptr.\n", me, ExpandCheckCode(res), file, line);
		abort();
	}

	byCount -= amt;
	if (listPos < LISTSIZE) {
			ToBeCleared[listPos] = NULL;
	}


	if (verboseFree)
		if (me<10)
			fprintf(stderr, "MyFree) me:%d %s:%d - Freed ptr(=#%020ld) of %d bs from %s:%d.\n", me, file, line, listPos, amt, fileAlloc, lineAlloc);

	if (p)
		free(p);
}
void *MyRealloc(void *ptr, size_t size, char *file, int line)
{
	if (ptr == NULL)
		return MyMalloc(size, file, line);

	if (size == 0){
		MyFree(ptr, file, line);
		return NULL;
	}

	/* Standard case: shrinking/expanding */
	{
		void *newP;
		size_t oldsize,commSz;
		long	listPos;
		char *p, *fileAlloc;
		int lineAlloc;
		int res;


		p = MyDecodePtr(ptr, "MyFree", file, line, &res, &fileAlloc, &lineAlloc, &listPos, &oldsize);


		commSz = (size<oldsize)? size : oldsize;
		newP = MyMalloc(size,file,line);
		if (newP) {
			memcpy(newP, ptr, commSz);
			MyFree(ptr, file, line);
			return newP;
		}

		return NULL;
	}
}

/* Helper routines  */
void MyCheckAllocInUse(const char* fileOut, long imin)
{
	size_t amt;
	long	listPos, i;
	int lineAlloc, res;
	void *ptr;
	char *p, *fileAlloc;
	FILE *fpOut;

	if (fileOut == NULL)
		fpOut = stderr;
	else {
		fpOut = fopen(fileOut, "a");
		if (!fpOut)
				fpOut = stderr;
	}

	for (i=(imin!=-1 ? imin : 0) ; i < listCount; i++){
		ptr = ToBeCleared[i];
		if (ptr) {
			p = MyDecodePtr(ptr, "MyCheckAllocInUse", "", 0, &res, &fileAlloc, &lineAlloc, &listPos, &amt);

			fprintf(fpOut, "MyCheckInUse) me:%d ptr(=#%020ld) status:%s dims:%d bs from %s:%d.\n", me,i, ExpandCheckCode(res), amt, fileAlloc, lineAlloc);
		}
	}

	if (fpOut != stderr)
		fclose(fpOut);
}
void MyAllocSet(int newVM, int newVF)
{
		verboseMalloc = newVM;
		verboseFree = newVF;

		fprintf(stderr, "cMyAllocSet) me:%d verboseMalloc=%d verboseFree=%d.\n", me, newVM, newVF);
}
static long MyAllocLive()
{
		long i, c=0;
		for (i=0; i<listCount; i++)
				if (ToBeCleared[i])
						c++;

		return c;
}
void MyAllocStat(const char* fileOut)
{
	FILE *fpOut;
	long live = MyAllocLive();

	if (fileOut == NULL)
			fpOut = stderr;
	else {
			fpOut = fopen(fileOut, "a");
			if (!fpOut)
					fpOut = stderr;
	}
	fprintf(fpOut, "MyAllocStat) me:%d Allocs:%ld (live:%ld) Count=%12ld(=%8ldKbs = %4ldMbs).\n",
			me, listCount, live, byCount, byCount>>10,  byCount>>20);

	if (fpOut != stderr)
			fclose(fpOut);
}
