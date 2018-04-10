#include "tools.h"
#include "param.h"
#include "arrays.h"

#include "derivatives.h"

void ah_geometry(void)
{
	// Extra mixed derivatives necessary for bicubic interpolation.
	if (ahfind_bicubic)
	{
		diff2rz(Drz_A_a, A_a, 1, 1);
		diff2rz(Drz_A_b, A_b, 1, 1);
		diff2rz(Drz_A_h, A_h, 1, 1);
		diff2rz(Drz_A_c, A_c, 1, -1);
		diff2rz(Drrz_a, Dr_a, -1, 1);
		diff2rz(Drrz_b, Dr_b, -1, 1);
		diff2rz(Drrz_h, Dr_h, -1, 1);
		diff2rz(Drrz_c, Dr_c, -1, -1);
		diff2rz(Drrz_phi, Dr_phi, -1, 1);
		diff2rz(Drzz_a, Dz_a, 1, -1);
		diff2rz(Drzz_b, Dz_b, 1, -1);
		diff2rz(Drzz_h, Dz_h, 1, -1);
		diff2rz(Drzz_c, Dz_c, 1, 1);
		diff2rz(Drzz_phi, Dz_phi, 1, -1);

		if (!(strcmp(slicing, "maximal") == 0))
		{
			diff2rz(Drz_K, K, 1, 1);
		}
	}

	return;
}