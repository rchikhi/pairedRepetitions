
#include <Judy.h>
#include "my_sarray/sarray.h"

Pvoid_t judy;
Pvoid_t judy2;
PWord_t pval;
Word_t idx;


void init_p(my_int len_sa)
{
}

void free_p()
{
}

void init_p_entries(my_int len_sa)
{
	judy = (PWord_t) NULL;
}

void free_p_entries(my_int len_sa)
{
	my_int i;
	Word_t    Bytes;
	idx=0;
	JLF(pval,judy,idx);		
	while (1)
	{
		if (pval==NULL) break;	
		J1FA(Bytes,*pval);
		JLN(pval,judy,idx);
	}
	JLFA(Bytes,judy);
}

my_int p_size(my_int len_sa)
{
	my_int i;
	Word_t    Bytes;
	Word_t total=0;
	idx=0;
	JLF(pval,judy,idx);		
	while (1)
	{
		if (pval==NULL) break;	
		J1MU(Bytes,*pval);
		total+=Bytes;
		JLN(pval,judy,idx);
	}
	JLMU(Bytes,judy);
	total+=Bytes;
	return total;
}



void set_p(my_int j,my_int k,int v)
{
// using encoding: 10 -> 2, 01 -> 1, 00 -> 2
	int Rc_int;
	JLI(pval,judy,j);
	if (v&1==1){
		J1S(Rc_int,*pval,2*k);}
	else	J1U(Rc_int,*pval,2*k);
	int v2=v>>1;
	if (v2&1==1){
		J1S(Rc_int,*pval,2*k+1);}
	else	J1U(Rc_int,*pval,2*k+1);

}

int return_p(my_int j,my_int k)
{
	int Rc_int;
	JLG(pval,judy,j);
	if (pval==NULL)
		return 0;
	J1T(Rc_int,*pval,2*k+1);
	if (Rc_int==1)
		return 2;
	J1T(Rc_int,*pval,2*k);
	if (Rc_int==1)
		return 1;
   	return 0;
}

