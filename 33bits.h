#include "my_sarray/sarray.h"
#include <stdint.h>
#include <Judy.h>

typedef struct
{ 
  uint32_t* array;
  unsigned char *boolv;
}
my_33bits;

my_33bits* alloc_33bits(my_int len_sa);
void free_33bits();
uint64_t get_33bits(my_33bits* array, my_int idx);
int set_33bits(my_33bits* array, my_int idx,uint64_t value);
my_int size_33bits(my_33bits* array, my_int len_sa);

