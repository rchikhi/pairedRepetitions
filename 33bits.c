
#include "33bits.h"
#include "my_sarray/sarray.h"

uint64_t retval=0;

my_33bits* alloc_33bits(my_int len_33bits)
{
	my_33bits *new=(my_33bits*)malloc(sizeof(my_33bits));
	new->array=(uint32_t*)malloc((len_33bits)*sizeof(uint32_t));
	if (new->array==NULL) error("can't alloc memory for 33-bits structure");
        new->boolv = (unsigned char *) malloc((len_33bits+8) / 8);
	if (new->boolv==NULL) error("can't alloc bool array for 33-bits structure");
	memset(new->boolv, 0, (len_33bits +8)/ 8);
	return new;
}

my_int size_33bits(my_33bits* array, my_int len_33bits)
{
	return (len_33bits)*sizeof(uint32_t)+((len_33bits+8)/8);
}

void free_33bits(my_33bits* array)
{
	free(array->array);
	free(array->boolv);
	free(array);
}

uint64_t get_33bits(my_33bits* array, my_int idx)
{
	if (((array->boolv[idx/8]>>(idx%8))&1)==0)
		retval=(array->array)[idx];
	else 
		retval=(uint64_t)((uint64_t)((uint64_t)(array->array)[idx])&0xFFFFFFFF)+(uint64_t)((uint64_t)1<<32);
	return retval;


}

int set_33bits(my_33bits *array, my_int idx,uint64_t value)
{
	if (((value>>32)&1) == 1)
	{
		array->boolv[idx / 8] |= ( 1) << (idx % 8);
	}
	else 
	{ 
		if ((((array->boolv[idx/8])>>(idx%8))&1)==1)
			array->boolv[idx / 8] ^= ( 1) << (idx % 8);
	}
	array->array[idx]=value&(0xFFFFFFFF);
	
	return 1;
}


