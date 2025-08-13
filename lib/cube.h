
/* FILE dll.h */
#ifdef BUILD_DLL
/* DLL export */
#define EXPORT __declspec(dllexport)



/* Declare our Add function using the above definitions. */
EXPORT double twoe(int a, int b, int c,int d,double *alpha_a,double *alpha_b,double *alpha_c, double *alpha_d, double *Da,double *Db,double *Dc, double *Dd, int *la, int *lb, int *lc,int *ld, double *xa, double *xb, double *xc, double *xd, int lenLa, int lenLb, int lenLc, int lenLd);

#endif    // DLL_H
  

