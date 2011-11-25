
#include "short2judy.h"
#include "my_sarray/sarray.h"
#include <Judy.h>

my_sj* alloc_sj(my_int len_sa)
{
	my_sj *new=(my_sj*)malloc(sizeof(my_sj));
	uint16_t *array;
	array=(uint16_t*)malloc((len_sa+1)*sizeof(uint16_t));
	if (array==NULL) error("can't alloc memory for short judy structure");
	new->array=array;
	new->judy=NULL;
	return new;
}

my_int size_sj(my_sj* array, my_int len_sa)
{
	int Rc_word;
	JLMU(Rc_word,array->judy);
	return (len_sa+1)*sizeof(uint16_t)+Rc_word;
}

void free_sj(my_sj* array)
{
	Word_t Bytes;
	JLFA(Bytes,array->judy);
	free(array->array);
	free(array);
}

Word_t get_sj(my_sj* array, my_int idx)
{
	PWord_t pval;
	if ((array->array)[idx]==0)
	{
		JLG(pval,array->judy,idx);
		return *pval;
	}
	return (array->array)[idx];


}

int set_sj(my_sj *array, my_int idx,Word_t value)
{
	PWord_t pval;
	int Rcint;
	if (value >= UINT16_MAX || value==0)
	{
		JLI(pval,array->judy,idx);
		*pval=value;
		(array->array)[idx]=0;
	}
	else 
	{
		//JLD(Rcint,(array->judy),idx);
		(array->array)[idx]=(uint16_t)value;
	}
	return 0;
}


