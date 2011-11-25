/****************************************************************
Copyright (C) Lucent Technologies 1997
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name Lucent Technologies or any of
its entities not be used in advertising or publicity pertaining
to distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/
/*
	Suffix sort
	Peter M. McIlroy
	M. Douglas McIlroy

Prototype
	int ssarray(int a[]);

Purpose
	Return in a[] a suffix array for the original
	contents of a[].  (The original values in a[]
	are typically serial numbers of distinct tokens
	in some list.)

Precondition
	Array a[] holds n values, with n>=1.  Exactly k 
	distinct values, in the range 0..k-1, are present.
	Value 0, an endmark, appears exactly once, at a[n-1].

Postcondition
	Array a[] is a copy of the internal array p[]
	that records the sorting permutation: if i<j
	then the original suffix a[p[i]..n-1] is
	lexicographically less than a[p[j]..n-1].

Return value
	-1 on error.
	Otherwise index i such that a[i]==0, i.e. the
        index of the whole-string suffix, used in
	Burrows-Wheeler data compression.

Cost
	Sorting takes time O(n log m), where m is the longest
	duplicate (possibly overlapped) substring, and uses n words
	of temporary memory.  There are at most log n passes,
	each involving essentially 5 serial and 6 random
	traversals of size-n arrays (a[] and p[]).

	A tradeoff desribed after step (3) below would save
	about 20% in time at the cost of another n words.

Terminology
	"The h-gram at a[i]" is the string a[i], a[i+1], ...,
	a[i+h-1], with subscripts taken mod n, except that
	when h=0 the "0-gram" at a[i] is a[i].  (More succinctly,
	the h-gram at a[i] is a[i], ..., a[i+max(h,1)-1].)
	"The h-successor of the h-gram at a[i]" is the h-gram
	at a[i+h].  
	
Method
	Order 2h-grams by doing steps (1)-(5) below
	for h = 0, 1, 2, 4, ...  Because of the unique endmark,
	ordering n-grams (or longer grams) is exactly the same as 
	ordering the suffixes of a[].

	This is a fillip on radix sort.  The i-th stage of a
	traditional radix sort distributes on the i-th "character"
	counting from the right of the keys, thus arranging the keys
	so that their i-character suffixes are in lexicographic order.
	Here, instead, we distribute on h-gram prefixes of 2h-grams,
	where h = 2^i, thus arranging 2h-grams in order.  Each 
	h-gram is coded into one element of a[].
	This is possible within a constant word width because
	there are only n h-grams for any h.  With 2^i-grams being
	placed in order at stage i, and each stage taking time O(n),
	sorting is complete after log n stages for a total
	time of O(n log n).

History
	The general radix-sorting notion is due to Manber and
	Myers, SODA '91, who used a big-endian radix sort.  They
	also made the observation that by cleverly multiplexing
	uses of the two arrays, one can get by with O(1) extra
	space beyond arrays a[] and p[].  The circular-array idea,
	which greatly simplifies the program, is due to P. McIlroy.
	Steps (1) and (2) below are based on M&M's method for
	the first cycle (h=0); they do something entirely different
	for the h-gram doubling cycles and develop in a[] the
	inverse of permutation p[].  

	The state of the art has improved since this code was
	written.  Variants developed independently by (at
	least) Larssen at Lund, Sadakane at Tokyo and Quinlan at
	Bell Labs is several times as fast; see Larssen and Sadakane 
	LU-CS-TR:99-214, Dept. of CS, Lund Univ. Sweden.

Data used in the sort
	a[]   array of h-gram codes.
	al[]  linked-list area, overlaid on a[].
	p[]   permutation.
	pl[]  list heads, overlaid on p[].
	ORIG  mark at the end of lists, originally on p[].
	BUCK  mark at start of "buckets", stretches of p[]
	      that point to identical codes in a[].

	Were the overlaid arrays pulled apart, the algorithm would
	simplify only slightly, mainly by dropping the list-chasing
	inner loop from steps (1) and (2).

Description of steps

(0)	(0a) Initialize p[] to contain the identity permutation.
	(0b) Place mark ORIG on every element of p[].

State before step 1
	Each element a[i] encodes the h-gram at a[i].
	The h-grams are coded 0,..,k-1 in increasing
	lexicographic order.  When h>0, permutation p[] 
	lists h-grams in lexicographic order: if i<j then
	a[p[i]] <= a[p[j]].  When h=0, p[] contains
	the input order 0,..,n-1.

(1)	Construct in al[] linked lists of like-valued
	codes in a[] ordered in reverse of the order
	of their h-successors as given by p[].
	(Done by iterating over elements of a[] in p[]-order
	and operating on their h-predecessors.)

	Place in pl[j] the head of the list of elements with
	value j.  Thus, if a[7]=a[3]=a[9]=2, and these
	elements were visited in this order according to p[],
	then pl[2]=9, al[9]=3, al[3]=7.  The list is developed
	by pushing previous contents of pl[j] onto the list.
	In particular the original content of pl[j] (which
	was p[j]) gets pushed to the end of the list (al[7]
	in the example), still bearing the mark ORIG.  The mark
	serves both to note the end of the list and to flag the
	fact that this element contains a quantity of a different
	kind.

	During step (1), a[] contains three kinds of data:
	unprocessed elements of a[], list links of al[]
	and list ends, marked ORIG, that contain original
	data from p[].  p[] contains two kinds of data:
	elements of p[], marked ORIG, and elements of pl[].

	The elements of p[] that have been displaced by pl[] 
	are retrieved by chasing the lists; this is done by
	a for loop.  Each list will be chased at most once in
	step (1) and once in step (2) and the total size of
	all lists will never exceed n, so the overhead of list
	chasing is at most linear in n.

State before step (2)
	Array a[] has been wiped out and replaced by
	the lists of al[], which together with pl[] give
	exactly the same information.  The first k elements
	of p[] are occupied by pl[] and the remaining elements
	are occupied by ORIG-marked elements of p[].

(2)	Make p[] point to the (former) elements of a[]
	in order of increasing value.  This is done in reverse:
	biggest values first, and for each value running out the
	list created in (1), which is already in reverse order.
	New values fill in the size-n array p[] from the top,
	while the size-k array pl[] of list heads shrinks
	toward the bottom.  Because no list is empty, the 
	two uses of p[] cannot collide.

	Place mark BUCK at each bucket start, namely the
	element of p[] last transferred from each list.

State before step 3
	Permutation p[] lists h-grams in lexicographic order.
	It is bucketed by h-grams.  Within each bucket h-grams
	are ordered by their h-successors.

(3)	Reconstruct codes in a[] from the buckets.  The
	element of a[] pointed to from the first bucket is 0;
	the elements pointed to from the second bucket are 1;
	and so on.

	This step would be unnecessary if al[] did not share
	storage with a[].  Then steps (2) and (4) could
	be combined.

State before step 4
	Permutation p[] is as before step 3, and code array
	a[] is as before step 1.

(4)	Refine buckets by values of h-successors in a[] by
	placing a BUCK mark on each p[i] where the h-successors
	of a[p[i]] and a[p[i-1]] differ.

State before step 5
	Array p[] holds an ordered, bucketed list of 2h-grams.

(5)	Recode a[] according to the new buckets, exactly as
	a[] was reconstructed in step (3).  Count the
	buckets to determine a new value of k.  If k==n, every
	bucket is a singleton and sorting is complete.


Bad input
	Tests that jump to "out" are solely defensive.  They
	may be dropped if you like to live dangerously.
*/

#include <stdlib.h>
#include "sarray.h"

enum {
	ORIG = ~(~0u>>1),			/* sign bit */
	BUCK = ~(~0u>>1)
};

#define pred(i, h) ((t=(i)-(h))<0?  t+n: t)
#define succ(i, h) ((t=(i)+(h))>=n? t-n: t)

int
ssarray(int a[])
{
	int h, i, j, l, n, t;
	int k = 0;				/* initialized for lint */
	int *p = 0;
	int result = -1;
#	define al a
#	define pl p


	for(j=n=0; a[n]>0; n++)			/* find n */
		if(a[n] > j)
			j = a[n];		/* and max element */
	if(a[n++]<0 || j>=n)
		goto out;
	p = malloc(n*sizeof(int));
	if(p == 0)
		goto out;

	for(i=0; i<n; i++)			/* (0) initialize */
		p[i] = i | ORIG;

	for(h = 0; ; h = h==0? 1: 2*h) {
		for(i=0; i<n; i++) {		/* (1) link */
			for(j=p[i]; !(j&ORIG); j=al[j]);
			j = pred(j&~ORIG, h);
			l = a[j];
			al[j] = pl[l];
			pl[l] = j;
		}

		if(h == 0) {			/* find k */
			for(k=0; k<n; k++)
				if(pl[k]&ORIG)
					break;

			for(i=k; i<n; i++)	/* check input */
				if(!(pl[i]&ORIG))
					goto out;
		}

		for(i=n; --k>=0; ) {		/* (2) order */
			j = pl[k];
			do
				p[--i] = j;
			while(!((j=al[j]) & ORIG));
			p[i] |= BUCK;
		}

		for(i=0; i<n; i++) {		/* (3) reconstruct */
			if(p[i] & BUCK)
				k++;
			a[p[i]&~BUCK] = k;
		}

		for(i=0, j=-1; i<n; i++, j=l) {	/* (4) refine */
			l = a[succ(p[i]&~BUCK, h)];
			if(l != j)
				p[i] |= BUCK;

		}

		for(i=0, k=-1; i<n; i++) {	/* (5) recode */
			if(p[i] & BUCK)
				k++;
			a[p[i]&~BUCK] = k;
			p[i] |= ORIG;		/* (0b) reinitialize */
		}
		if(++k >= n)
			break;
	}

	for(i=0; i<n; i++) {
		a[i] = p[i] & ~ORIG;
		if(a[i] == 0)
			result = i;
	}
	
out:
	free(p);
	return result;
}
