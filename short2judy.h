#include "my_sarray/sarray.h"
#include <stdint.h>
#include <Judy.h>

typedef struct
{ 
  uint16_t* array;
  Pvoid_t judy;
}
my_sj;

my_sj* alloc_sj(my_int len_sa);
void free_sj();
Word_t get_sj(my_sj* array, my_int idx);
int set_sj(my_sj* array, my_int idx,Word_t value);
my_int size_sj(my_sj* array, my_int len_sa);

