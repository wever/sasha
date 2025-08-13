#ifndef DLL_H
#define DLL_H

/* FILE dll.h */
#ifdef BUILD_DLL
/* DLL export */
#define EXPORT __declspec(dllexport)
#else
/* EXE import */
#define EXPORT __declspec(dllimport)
#endif



/* Declare our Add function using the above definitions. */
EXPORT void overlap(double *D, double *A, double *S, double *Center, int *CartAng, int N, int *lenL, int *locL);
EXPORT void dipole(double *D, double *A, double *Mu, double *Center, int *CartAng, int N, int *lenL, int *locL);
EXPORT void kei(double *D, double *A, double *KEI, double *Center, int *CartAng, int N, int *lenL, int *locL);
EXPORT void nai(double *D, double *A, double *NAI, double *Center, int *CartAng, double *R, int *Z, int N, int nnucleus, int *lenL, int *locL);
EXPORT void loopSCF(double *G, double *ERI, double *P, int N);




#endif    // DLL_H
  

