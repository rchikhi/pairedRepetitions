/*
Suffix arrays
A suffix array is a simple data structure that enables lookup of any substring of a text and identification of repeated substrings. It is more compact than a suffix tree and is amenable to storage in secondary memory.

A suffix array is (notionally) a sorted list of the suffixes of a given string. (A suffix is a substring that extends to the end of the given string.) The sorted list is presented as an array of integers that identify the suffixes in order. The suffix array for the 11-letter text ABRACADABRA, with an endmark # that collates low, is given below, together with the corresponding suffixes.

11  #
10  A#
 7  ABRA#
 0  ABRACADABRA#
 3  ACADABRA#
 5  ADABRA#
 8  BRA#
 1  BRACADABRA#
 4  CADABRA#
 6  DABRA#
 9  RA#
 2  RACADABRA#

Together, the suffix array and string enable binary search for any substring, e.g. CAD.

A useful auxiliary data structure is an `LCP array', an array of lengths of the longest common prefix between each substring and its predecessor in the suffix array. The second column below gives the LCP array for the previous example. The third column shows that a suffix array may also be interpreted as an ordered list of circular shifts.

11  0  #ABRACADABRA
10  0  A#ABRACADABR
 7  1  ABRA#ABRACAD
 0  4  ABRACADABRA#
 3  1  ACADABRA#ABR
 5  1  ADABRA#ABRAC
 8  0  BRA#ABRACADA
 1  3  BRACADABRA#A
 4  0  CADABRA#ABRA
 6  0  DABRA#ABRACA
 9  0  RA#ABRACADAB
 2  2  RACADABRA#AB

The last column of letters, ARD#RCAAAABB, is the Burrows-Wheeler transform, which is central to B-W data compression. 

*/

/* $Id$  */

#include <stdio.h>
#include <string.h>
#include "sarray.h"

int main()
{
    char* s = "ABRACADABRA#";
    my_int n = strlen(s)-1;
    my_int* a = (my_int*)malloc((n+1)*sizeof(my_int));
    my_lc* b = (my_lc*)malloc((n+1)*sizeof(my_lc));
    my_int i;
    
    bsarray(s,a,n);
    lcpa(a,s,b,n+1);
    
    printf("string:%s\n",s);

    printf("suffix array:\n");
    for(i=0;i<n+1;i++)
	printf("%d:\t%s\n",a[i],s+a[i]);
    printf("lcp array:\n");
    for(i=0;i<n+1;i++)
	printf("%d\t%d\n",a[i],b[i]);

    free(a);
    free(b);
    return 0;
}
