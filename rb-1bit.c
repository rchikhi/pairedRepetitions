
#include <Judy.h>
#include "my_sarray/sarray.h"

Pvoid_t judyrb;
Word_t idxrb;
int Rc_intrb;

void init_rb()
{
	judyrb = (PWord_t) NULL;
}
void rb_refresh()
{
	J1FA(idxrb,judyrb);
}

void set_rb(my_int k)
{
	J1S(Rc_intrb,judyrb,k);
}

int return_rb(my_int k)
{
	J1T(Rc_intrb,judyrb,k);
	if (Rc_intrb==1)
		return 1;
   	return 0;
}

