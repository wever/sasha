#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <locale.h>
#include "cube.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif


	void computecube(double *x_plot, double *y_plot, double *z_plot, double *A, double *D, double *C, double *Center, int *CartAng, int *lenL,int *listlocL, int N, int n, int nptsx, int nptsy, int nptsz, double xmin, double ymin, double zmin, double voxelx, double voxely, double voxelz,int nos_initial, int nos_final, int *Z, double *R, int plot_type, char *inputName) {

		setlocale(LC_ALL, "en_US.utf8");
		FILE *fp;
		fp = fopen(inputName,"w");

		int x,y,z,cgf,mi,ni,state,Nprimitive,Nprimitive2,n_plot,locL1,locL2;
		char *program_name;
		program_name = "SASHA PROGRAM";
		double cube,L1,L2;
		printf("%d\n",plot_type);

		double bfs(int mini, int loc, int Np,  double xp,double yp,double zp) {
			return pow((xp - Center[3*mini]),(CartAng[3*mini]))*pow((yp - Center[3*mini+1]),(CartAng[3*mini+1]))*pow((zp - Center[3*mini+2]),(CartAng[3*mini+2]))*exp(-A[loc+Np]*(pow((xp - Center[3*mini]),2.) + pow((yp - Center[3*mini+1]),2.) + pow((zp - Center[3*mini+2]),2.)));
		}



		if (plot_type == 0) {
			n_plot = -1*n;
		}
		else {
			n_plot = n;
		}


		void header() {
			int nn,nstate;

			fprintf(fp,"%s %s %s\n","Created by",program_name,"package");
			fprintf(fp,"%s\n","-----------------------------------");
			fprintf(fp,"%5d %11.6f %11.6f %11.6f\n", n_plot,xmin,ymin,zmin);
			fprintf(fp,"%5d %11.6f %11.6f %11.6f\n", nptsx,voxelx,0.,0.);
			fprintf(fp,"%5d %11.6f %11.6f %11.6f\n", nptsy,0.,voxely,0.);
			fprintf(fp,"%5d %11.6f %11.6f %11.6f\n", nptsz,0.,0.,voxelz);
			for (nn = 0; nn < n; nn++) {
				fprintf(fp,"%5d %11.6f %11.6f %11.6f %11.6f\n", Z[nn],1.*Z[nn],R[3*nn + 0],R[3*nn + 1],R[3*nn + 2]);
			}

			if (plot_type == 0) {
				fprintf(fp,"%5d", (nos_final - nos_initial));
				for (nstate = nos_initial; nstate < nos_final; nstate++) {
					fprintf(fp,"%5d",nstate);
				}
				fprintf(fp,"\n");
			}      
			return;
		}

		if (plot_type == 0) {
			nos_initial++;
			nos_final++;
			nos_final++;
			header();
			int count_line = 0;
			int count_tot = 0;
			for (x = 0; x < nptsx; x++) {
				for (y = 0; y < nptsy; y++) {
					for (z = 0; z < nptsz; z++) {
						for (state = nos_initial; state < nos_final; state++) {
							cube = 0;
							for (cgf = 0; cgf < N; cgf++) {
								L1 = lenL[cgf];
								locL1 = listlocL[cgf];
								for (Nprimitive = 0; Nprimitive < L1; Nprimitive++) {
									cube += C[cgf + N*state]*D[locL1+Nprimitive]*bfs(cgf,locL1,Nprimitive,x_plot[x],y_plot[y],z_plot[z]);

								}   
							}   
							if (fabs(cube) < 1e-12) {
								cube = 0;
							} 
							fprintf(fp,"%13.5e", cube);
							count_line ++; 
							count_tot ++; 
							if ((count_line % 6 == 0) || (count_tot % ((nos_final-nos_initial)*nptsz)) == 0) {
								fprintf(fp,"\n");
								count_line = 0;
							}   
						}   
					}   
				}   
			}
		}


		if (plot_type == 1) {

			/* ADJUST TO PRINT FROM 1 (NOT 0) TO N, AND AVOIDING NINIT-NFINAL = 0 */
			nos_initial++;
			nos_final++;
			nos_final++;
			header();
			int count_line = 0;
			int count_tot = 0;
			for (x = 0; x < nptsx; x++) {
				for (y = 0; y < nptsy; y++) {
					for (z = 0; z < nptsz; z++) {
						cube = 0.;
						for (mi = 0; mi < N; mi++) {
							L1 = lenL[mi];
							locL1 = listlocL[mi];
							for (ni = 0; ni < N; ni++) {
								L2 = lenL[ni];
								locL2 = listlocL[ni];
								for (Nprimitive = 0; Nprimitive < L1; Nprimitive++) {
									for (Nprimitive2 = 0; Nprimitive2 < L2; Nprimitive2++) {

										cube += C[mi+ni*N]*D[locL1+Nprimitive]*D[locL2+Nprimitive2]*bfs(mi,locL1,Nprimitive,x_plot[x],y_plot[y],z_plot[z])*bfs(ni,locL2,Nprimitive2,x_plot[x],y_plot[y],z_plot[z]);


									}
								}
							}
						}
						if (fabs(cube) < 1.e-12) {
							cube = 0.;
						}
						fprintf(fp,"%14.5e", cube);
						count_line ++; 
						count_tot ++; 
						if ((count_line % 6 == 0) || (count_tot % nptsz) == 0) {
							fprintf(fp,"\n");
							count_line = 0;
						}  
					}
				}
			}
		}

		/* Fechar o ficheiro */
		fclose(fp);
	}


#ifdef __cplusplus
} // __cplusplus defined.
#endif
