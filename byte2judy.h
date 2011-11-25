#include "my_sarray/sarray.h"
#include <stdint.h>
#include <Judy.h>

typedef struct
{ 
  uint8_t* array;
  Pvoid_t judy;
}
my_bj;

my_bj* alloc_bj(my_int len_sa);
void free_bj();
Word_t get_bj(my_bj* array, my_int idx);
int set_bj(my_bj* array, my_int idx,Word_t value);
my_int size_bj(my_bj* array, my_int len_sa);

