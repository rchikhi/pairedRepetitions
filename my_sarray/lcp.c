#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sarray.h"

/* 
   int *lcp(const int *a, const char *s, int n)
   Precondition: a is suffix array for string s of
      length n *including* the terminating '\0'. 
   Return value: longest-common-prefix array; 0 on error.
   Reference: T. Kasai, G. Lee, H. Arimura, S.Arikawa
   and K. Park, "Linear-time longest-common-prefix
   computation in suffix arrays and its applications",
   Proc 12th Annual Conference on Combinatorial Pattern
   Matching, Springer, LNCS 2089 (2001) 181-192.

   lcp[x] is the length of the longest common prefix of
   suffixes s[a[x-1]..] and s[a[x]..].

   The algorithm determines the elements of lcp in the
   order that the suffixes occur in s.  It uses this fact:
   If the lcp for suffix s[i..] has length h, where h>0,
   then the lcp for suffix s[i+1..] is at least h-1.

   Proof.
   Let the immediate lexicographic predecessor of suffix
      s[i..] be s[j..], i.e. lex[i]=lex[j]+1.
   If s[i..] and s[j..] have a common prefix of length h,
      where h>0, then s[i+1..] and s[j+1..] have a common
      prefix of length h-1.
   Since s[i+1..] and s[j+1..] differ from s[i..] and s[j..]
      respectively only by the deletion of a common first
      letter, the two pairs must be similarly ordered.
      Hence s[j+1..] lexicographically precedes s[i+1..].
   Since s[i+1..] shares a common prefix of length h-1 with
      some lexicographic predecessor, namely s[j+1..], it
      must share a common prefix of length at least h-1 with
      its immediate predecessor.  Otherwise the suffix array 
      would be out of order.

   Running time is O(n).

   Proof.
   h is bounded by n; and h is decreased by 1 at most
   n times. Hence h is increased at most 2n times.
   This bounds the number of executions of the inner loop.
*/

/* 
   inv is the inverse of a: if inv[i]=x then a[x]=i.
   In other words, inv[i] is the index x of the
   pointer (in array a) to suffix s[i..].
*/

my_lc*
lcp(const my_int *a, const char *s, my_int n) 
{
	my_lc *lcp = (my_lc*)malloc(n*sizeof(my_lc));

	if(lcp == 0)
		return 0;
	if(lcpa(a, s, lcp, n) == 0) {
		free(lcp);
		return 0;
	}
	return lcp;
}

/* lcpa is used by the java native method */

int
lcpa(const my_int *a, const char *s0, my_lc *lcp, my_int n)
{
	my_int i;
	my_lc h;
	uchar *s = (uchar*)s0;
	printf(".");fflush(stdout);
	my_int *inv = (my_int*)malloc(n*sizeof(my_int));
	my_int x,j;

	if(inv == 0)
		return 0;
	for(i=0; i<n; i++)
		inv[a[i]] = i;
	printf(".");fflush(stdout);

	h = 0;			/* visit in string order */
	for(i=0; i<n-1; i++) {  /* omit last, least suff */
		x = inv[i]; /* i,j,x,h as in intro */
		j = a[x-1];
		/*uchar *p1 = s + i + h;
		uchar *p0 = s + j + h;
		while(*p1++ == *p0++)
			h++;*/
		while (*((uchar *)(s+i+h))==*((uchar *)(s+j+h))) 
			h++;
		lcp[x] = h;
		if(h > 0)
			h--;
	}
	printf(".");fflush(stdout);
	lcp[0] = 0;	/* least suffix has no predecessor */
	free(inv);
	return 1;
}
