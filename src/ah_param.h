#ifdef MAIN_FILE
double *array_th;
double *array_H;
double *array_dH;
double *array_a;
double *array_b;
double *array_h;
double *array_c;
double *array_phi;
int *ah_mask_array;

int NthTotal;

double drr;
double dth;

int ah_found;
double start;

double ah_max;
double ah_min;

int ah_count;

FILE *f_H;
#else
extern double *array_th;
extern double *array_H;
extern double *array_dH;
extern double *array_a;
extern double *array_b;
extern double *array_h;
extern double *array_c;
extern double *array_phi;
extern int *ah_mask_array;

extern int NthTotal;

extern double drr;
extern double dth;

extern int ah_found;
extern double start;

extern double ah_max;
extern double ah_min;

extern int ah_count;

extern FILE *f_H;
#endif
