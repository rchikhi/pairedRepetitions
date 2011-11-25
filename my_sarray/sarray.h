#include <stdint.h>

#ifndef SARRAY_H
#define SARRAY_H
typedef unsigned char uchar;

#ifndef MY_INT
#define MY_INT
typedef long int my_int;
typedef uint32_t my_lc;
typedef uint32_t my_r;
#endif

int ssarray(int *a);
my_int sarray(my_int *a, my_int n);
my_int bsarray(const uchar *b, my_int *a, my_int n);
my_lc *lcp(const my_int *a, const char *s, my_int n);
int lcpa(const my_int *a, const char *s, my_lc *b, my_int n);
my_int *scode(const char *s);
uchar *codetab(const uchar *s);
uchar *inverse(const uchar *t);
#endif
