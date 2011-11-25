
#include "byte2judy.h"
#include "my_sarray/sarray.h"
#include <Judy.h>

my_bj* alloc_bj(my_int len_sa)
{
	my_bj *new=(my_bj*)malloc(sizeof(my_bj));
	uint8_t *array;
	array=(uint8_t*)malloc((len_sa+1)*sizeof(uint8_t));
	if (array==NULL) error("can't alloc memory for byte judy structure");
	new->array=array;
	new->judy=NULL;
	return new;
}

my_int size_bj(my_bj* array, my_int len_sa)
{
	int Rc_word;
	JLMU(Rc_word,array->judy);
	return (len_sa+1)*sizeof(uint8_t)+Rc_word;
}

void free_bj(my_bj* array)
{
	Word_t Bytes;
	JLFA(Bytes,array->judy);
	free(array->array);
	free(array);
}

Word_t get_bj(my_bj* array, my_int idx)
{
	PWord_t pval;
	if ((array->array)[idx]==UINT8_MAX)
	{
		JLG(pval,array->judy,idx);
		return *pval;
	}
	return (array->array)[idx];


}

int set_bj(my_bj *array, my_int idx,Word_t value)
{
	PWord_t pval;
	int Rcint;
	if (value >= UINT8_MAX)
	{
		JLI(pval,array->judy,idx);
		*pval=value;
		(array->array)[idx]=UINT8_MAX;
	}
	else 
	{
		//JLD(Rcint,(array->judy),idx);
		(array->array)[idx]=(uint8_t)value;
	}
	return 0;
}


