
/*
    Written by Weverson R. Gomes, 2015.
    Copyright (c) 2015, Weverson R. Gomes.

    This file is part of Sasha Software.

    Sasha Software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sasha Software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sasha Software.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "HPints.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

int a,b,c,d,ab,cd,abtest,cdtest,p,q,r,s,L;
int e,f,g,h;
int la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd;
double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz;
double Dp,Dq,Na,Nb,Nc,Nd;
double Px,Py,Pz,Qx,Qy,Qz,Prc;
double AB,CD,PQ,U;
double alpha_ab, alpha_cd, alpha_p, alpha_q, T, alpha_pq, Up, Uq, onecenterRR;
double m, T;
double Kp,Kq,omega;
double AB2, CD2, PQ2, Wx, Wy, Wz;
int i,j,k,l;
int mtot, im;
double eri;
double ERI;

  double absolute(double value) {
    if (value < 0) {
      return -value;
    }
    else {
      return value;
    }
  }

  int Index(int i,int j,int k,int l) {
    int tmp,ij,kl;
    if (i<j) {
      return Index(j,i,k,l);}
    if (k<l) {
      return Index(i,j,l,k);}
        ij = i*(i+1)/2+j;
        kl = k*(k+1)/2+l;
	if (ij<kl) {
	  tmp = ij;
	  ij = kl;
	  kl = tmp;
	}
    return ij*(ij+1)/2+kl;
  }

  long int Factorial2(int n) {
    if (n > 1) {
      long int f = 1;
      while (n > 1) {
        f *= n;
        n = n - 2;
      }
      return f;
    }
    else {
      return 1;
    }
  }
 
  long int factorial(int n) {
    int c;
    long int fact = 1;
    for (c = 1; c <= n; c++)
      fact *= c;
  return fact;
  }

  double fgama(int ni, long double T) {
  double c,F1,F2,fgam;
  int k;
  if (T<1e-6) {
    fgam =  1./(2.*ni+1.);
  }
  else {
    c = Factorial2(2*ni-1)/pow(2.*T,ni);
    F1 = sqrt(M_PI/(4.*T))*erf(sqrt(T));
    F2 = 0;
    for (k = 1; k < (ni+1); k++) {
      F2 += pow(2.*T,k-1)/Factorial2(2*k-1);
      //printf("%d %8.5Lf %8.5f %8.5f\n",ni,T,c,F2);
    }
    fgam = c*(F1-exp(-T)*F2);
  }
  //printf("%d %8.5Lf %8.5f\n",ni,T,fgam);
  return fgam;
  }

  int binom(int n, int p) {
  return factorial(n)/(factorial(p)*factorial(n-p));
  }
  
  double fk(int k,int la,int lb,double r1,double r2) {

    double f = 0.;
    int i,j;
      for (i = 0; i < (la+1); i++) {
        for (j = 0; j < (lb+1); j++) {
          if (i+j==k) {
            f += binom(la,i)*binom(lb,j)*pow(r1,la-i)*pow(r2,lb-j);
	  }
        }
      }
  return f;
  }

   double radius_pair(double x1,double y1,double z1,double x2,double y2,double z2) {
   return sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
  }
  
   double norm(double a, int l, int m, int n) {
        return pow(pow(2.,2.*(l+m+n)+1.5)*pow(a,l+m+n+1.5)/(Factorial2(2*l-1)*Factorial2(2*m-1)*Factorial2(2*n-1)*pow(M_PI,1.5)),0.5);  
    }  
  

  double twoe(int a, int b, int c,int d,double *alpha_a,double *alpha_b,double *alpha_c, double *alpha_d, double *Da,double *Db,double *Dc, double *Dd, int *la, int *lb, int *lc,int *ld, double *xa, double *xb, double *xc, double *xd, int lenLa, int lenLb, int lenLc, int lenLd) {


    double Fgamma(int ni, long double T) {
      double c,F1,F2,fgam;
      int k;
      if (T<1e-6) {
	fgam =  1./(2.*ni+1.);
      }
      else {
	c = Factorial2(2*ni-1)/pow(2.*T,ni);
	F1 = sqrt(M_PI/(4.*T))*erf(sqrt(T));
	F2 = 0;
	for (k = 1; k < (ni+1); k++) {
	  F2 += pow(2.*T,k-1)/Factorial2(2*k-1);
	  //printf("%d %8.5Lf %8.5f %8.5f\n",ni,T,c,F2);
	}
	fgam = c*(F1-exp(-T)*F2);
      }
    return fgam;
    }


    double vrr(double Dah, double Axh, double Ayh, double Azh, double alpha_ah, int axh, int ayh, int azh, double Dbh, double Bxh, double Byh, double Bzh, double alpha_bh, int bxh, int byh, int bzh, double Dch, double Cxh, double Cyh, double Czh, double alpha_ch, int cxh, int cyh, int czh, double Ddh, double Dxh, double Dyh, double Dzh, double alpha_dh, int dxh, int dyh, int dzh, int M) {


    AB = radius_pair(Ax,Ay,Az,Bx,By,Bz);
    CD = radius_pair(Cx,Cy,Cz,Dx,Dy,Dz);
    AB2 = pow(AB,2);
    CD2 = pow(CD,2);

    alpha_ab = alpha_ah + alpha_bh;
    alpha_cd = alpha_ch + alpha_dh;

    Px = (alpha_ah*Ax+alpha_bh*Bx)/alpha_ab;
    Py = (alpha_ah*Ay+alpha_bh*By)/alpha_ab;
    Pz = (alpha_ah*Az+alpha_bh*Bz)/alpha_ab;

    Qx = (alpha_ch*Cx+alpha_dh*Dx)/alpha_cd;
    Qy = (alpha_ch*Cy+alpha_dh*Dy)/alpha_cd;
    Qz = (alpha_ch*Cz+alpha_dh*Dz)/alpha_cd;

    alpha_pq = alpha_ab + alpha_cd;

    Wx = (alpha_ab*Px+alpha_cd*Qx)/alpha_pq;
    Wy = (alpha_ab*Py+alpha_cd*Qy)/alpha_pq;
    Wz = (alpha_ab*Pz+alpha_cd*Qz)/alpha_pq;

    PQ = radius_pair(Px,Py,Pz,Qx,Qy,Qz);
    PQ2 = pow(PQ,2);

    T = alpha_ab*alpha_cd/(alpha_ab+alpha_cd)*PQ2;

    Kp = pow(2.,0.5)*pow(M_PI,1.25)/alpha_ab*exp(-alpha_ah*alpha_bh/alpha_ab*AB2);
    Kq = pow(2.,0.5)*pow(M_PI,1.25)/alpha_cd*exp(-alpha_ch*alpha_dh/alpha_cd*CD2);

    mtot = axh+ayh+azh+cxh+cyh+czh+M;

    Dp = Da[p]*Db[q];
    Dq = Dc[r]*Dd[s];

    double vrr_terms[axh+1][ayh+1][azh+1][cxh+1][cyh+1][czh+1][mtot+1];
    

    for (im = 0; im < mtot+1; im++) {
        //vrr_terms[0][0][0][0][0][0][im] = Na*Nb*Nc*Nd*Kp*Kq/sqrt(alpha_ab+alpha_cd)*Fgamma(im,T);
        vrr_terms[0][0][0][0][0][0][im] = Kp*Kq/sqrt(alpha_ab+alpha_cd)*Fgamma(im,T);
        //vrr_terms[0][0][0][0][0][0][im] = Kp*Kq/sqrt(alpha_ab+alpha_cd)*Fgamma(im,T);
	//printf("%d %f\n",im,vrr_terms[0][0][0][0][0][0][im]);
	//printf("%f %f %f %f\n",alpha_ah,alpha_bh,alpha_ch,alpha_dh);
    }

    for (i = 0; i < axh; i++) {
	for (im = 0; im < mtot-i; im++) {
	    vrr_terms[i+1][0][0][0][0][0][im] = (Px-Ax)*vrr_terms[i][0][0][0][0][0][im] + (Wx-Px)*vrr_terms[i][0][0][0][0][0][im+1];
            if (i) {
		vrr_terms[i+1][0][0][0][0][0][im] += (i/(2.*alpha_ab)*(vrr_terms[i-1][0][0][0][0][0][im] - alpha_cd/alpha_pq*vrr_terms[i-1][0][0][0][0][0][im+1]));
	    }
	}
    }

	//printf("%d %f\n",axh,vrr_terms[axh][0][0][0][0][0][0]);
    for (j = 0; j < ayh; j++) {
	for (i = 0; i < axh+1; i++) {
	    for (im = 0; im < mtot-i-j; im++) {
		vrr_terms[i][j+1][0][0][0][0][im] = (Py-Ay)*vrr_terms[i][j][0][0][0][0][im] + (Wy-Py)*vrr_terms[i][j][0][0][0][0][im+1];
		if (j) {
		    vrr_terms[i][j+1][0][0][0][0][im] += (j/(2*alpha_ab)*(vrr_terms[i][j-1][0][0][0][0][im] - alpha_cd/alpha_pq*vrr_terms[i][j-1][0][0][0][0][im+1]));
		}
	    }
	}
    }

    
	//printf("%d %d %f\n",axh,ayh,vrr_terms[axh][ayh][0][0][0][0][0]);
    for (k = 0; k < azh; k++) {
	for (j = 0; j < ayh+1; j++) {
	    for (i = 0; i < axh+1; i++) {
		for (im = 0; im < mtot-i-j-k; im++) {
		    vrr_terms[i][j][k+1][0][0][0][im] = (Pz-Az)*vrr_terms[i][j][k][0][0][0][im] + (Wz-Pz)*vrr_terms[i][j][k][0][0][0][im+1];
		    if (k) {
			vrr_terms[i][j][k+1][0][0][0][im] += k/2./alpha_ab*(vrr_terms[i][j][k-1][0][0][0][im] - alpha_cd/alpha_pq*vrr_terms[i][j][k-1][0][0][0][im+1]);
		    }
		}
	    }
	}
    }
	//printf("%d %d %d %f\n",axh,ayh,azh,vrr_terms[axh][ayh][azh][0][0][0][0]);

    for (q = 0; q < cxh; q++) {
	for (k = 0; k < azh+1; k++) {
	    for (j = 0; j < ayh+1; j++) {
		for (i = 0; i < axh+1; i++) {
		    for (im = 0; im < mtot-i-j-k-q; im++) {
			vrr_terms[i][j][k][q+1][0][0][im] = (Qx-Cx)*vrr_terms[i][j][k][q][0][0][im] + (Wx-Qx)*vrr_terms[i][j][k][q][0][0][im+1];
			if (q) {
			    vrr_terms[i][j][k][q+1][0][0][im] += (q/2./alpha_cd*(vrr_terms[i][j][k][q-1][0][0][im] - alpha_ab/alpha_pq*vrr_terms[i][j][k][q-1][0][0][im+1]));
			}
			if (i) {
			    vrr_terms[i][j][k][q+1][0][0][im] += (i/2./alpha_pq*(vrr_terms[i-1][j][k][q][0][0][im+1]));
			}
		    }
		}
	    }
	}
    }

    //printf("%f %f %f\n",Wx,Wy,Wz);

	//printf("%d %d %d %d %f\n",axh,ayh,azh,cxh,vrr_terms[axh][ayh][azh][cxh][0][0][0]);
    for (r = 0; r < cyh; r++) {
	for (q = 0; q < cxh+1; q++) {
	    for (k = 0; k < azh+1; k++) {
		for (j = 0; j < ayh+1; j++) {
		    for (i = 0; i < axh+1; i++) {
			for (im = 0; im < mtot-i-j-k-q-r; im++) {
			    vrr_terms[i][j][k][q][r+1][0][im] = (Qy-Cy)*vrr_terms[i][j][k][q][r][0][im] + (Wy-Qy)*vrr_terms[i][j][k][q][r][0][im+1];
			    if (r) {
				vrr_terms[i][j][k][q][r+1][0][im] += (r/2./alpha_cd*(vrr_terms[i][j][k][q][r-1][0][im] - alpha_ab/alpha_pq*vrr_terms[i][j][k][q][r-1][0][im+1]));
			    }
			    if (j) {
				vrr_terms[i][j][k][q][r+1][0][im] += (j/2./alpha_pq*(vrr_terms[i][j-1][k][q][r][0][im+1]));
			    }
			}
		    }
		}
	    }
	}
    }

    //printf(" %f %f %f %f %f %f\n",Qy,Cy,Wy,Qy,alpha_cd,alpha_ab);
	//printf("%d %d %d %d %d %f\n",axh,ayh,azh,cxh,cyh,vrr_terms[axh][ayh][azh][cxh][cyh][0][0]);

    for (s = 0; s < czh; s++) {
	for (r = 0; r < cyh+1; r++) {
	    for (q = 0; q < cxh+1; q++) {
		for (k = 0; k < azh+1; k++) {
		    for (j = 0; j < ayh+1; j++) {
			for (i = 0; i < axh+1; i++) {
			    for (im = 0; im < mtot-i-j-k-q-r-s; im++) {
				vrr_terms[i][j][k][q][r][s+1][im] = (Qz-Cz)*vrr_terms[i][j][k][q][r][s][im] + (Wz-Qz)*vrr_terms[i][j][k][q][r][s][im+1];
				if (s) {
				    vrr_terms[i][j][k][q][r][s+1][im] += (s/2./alpha_cd*(vrr_terms[i][j][k][q][r][s-1][im] - alpha_ab/alpha_pq*vrr_terms[i][j][k][q][r][s-1][im+1]));
				}
				if (k) {
				    vrr_terms[i][j][k][q][r][s+1][im] += (k/2./alpha_pq*(vrr_terms[i][j][k-1][q][r][s][im+1]));
				}
			    }
			}
		    }
		}
	    }
	}
    }
    //printf("%f\n",vrr_terms[axh][ayh][azh][cxh][cyh][czh][0]);
   //printf("%d %d %d %d %d %d %.8f\n" ,axh,ayh,azh,cxh,cyh,czh,vrr_terms[axh][ayh][azh][cxh][cyh][czh][0]);

//printf("%f\n",vrr_terms[i][j][k][q][r][s][M]);
//printf("%d %d %d %d %d %d %f\n",axh,ayh,azh,cxh,cyh,czh,vrr_terms[axh][ayh][azh][cxh][cyh][czh][M]);
    return vrr_terms[axh][ayh][azh][cxh][cyh][czh][M];//vrr_terms[i][j][k][q][r][s][M];
    }

    double comp_contract(double *Dah, double Axh, double Ayh, double Azh, double *alpha_ah, int axh, int ayh, int azh, double *Dbh, double Bxh, double Byh, double Bzh, double *alpha_bh, int bxh, int byh, int bzh, double *Dch, double Cxh, double Cyh, double Czh, double *alpha_ch, int cxh, int cyh, int czh, double *Ddh, double Dxh, double Dyh, double Dzh, double *alpha_dh, int dxh, int dyh, int dzh) {

	eri = 0;
	for (e = 0; e < lenLa; e++) {
	    for (f = 0; f < lenLb; f++) {
		for (g = 0; g < lenLc; g++) {
		    for (h = 0; h < lenLd; h++) {
		        Na = norm(alpha_ah[e],ax,ay,az);
		        Nb = norm(alpha_bh[f],bx,by,bz);
		        Nc = norm(alpha_ch[g],cx,cy,cz);
		        Nd = norm(alpha_dh[h],dx,dy,dz);


				eri += Na*Nb*Nc*Nd*Dah[e]*Dbh[f]*Dch[g]*Ddh[h]*vrr(Dah[e],Axh,Ayh,Azh,alpha_ah[e],axh,ayh,azh,Dbh[f],Bxh,Byh,Bzh,alpha_bh[f],bxh,byh,bzh,Dch[g],Cxh,Cyh,Czh,alpha_ch[g],cxh,cyh,czh,Ddh[h],Dxh,Dyh,Dzh,alpha_dh[h],dxh,dyh,dzh,0);
		    }
		}
	    }
	}
	//printf("%s %f\n", "FINAL: ", eri);
	//printf("%d %d %d %d %f\n", a,b,c,d, eri);

    return eri;
    }
    
    double hrr(double *Dah, double Axh, double Ayh, double Azh, double *alpha_ah, int axh, int ayh, int azh, double *Dbh, double Bxh, double Byh, double Bzh, double *alpha_bh, int bxh, int byh, int bzh, double *Dch, double Cxh, double Cyh, double Czh, double *alpha_ch, int cxh, int cyh, int czh, double *Ddh, double Dxh, double Dyh, double Dzh, double *alpha_dh, int dxh, int dyh, int dzh) {

    //printf("%d %d %d %d %d %d %d %d %d %d %d %d\n", axh, ayh, azh, bxh, byh, bzh, cxh, cyh, czh, dxh, dyh, dzh);
    if (bxh > 0) {
        return (hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh+1, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh-1, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh) + (Axh-Bxh)*hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh-1, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh));
    }
    
    else if (byh > 0) {
        return (hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh+1, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh-1, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh) + (Ayh-Byh)*hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh-1, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh)); 
    }
    
    else if (bzh > 0) {
        return (hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh+1, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh-1, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh) + (Azh-Bzh)*hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh-1, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh)); 
    }
        
    else if (dxh > 0) {
        return (hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh+1, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh-1, dyh, dzh) + (Cxh-Dxh)*hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh-1, dyh, dzh)); 
    }
        
    else if (dyh > 0) {
        return (hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh+1, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh-1, dzh) + (Cyh-Dyh)*hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh-1, dzh));
    }
    
    else if (dzh > 0) {
        return (hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh+1, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh-1) + (Czh-Dzh)*hrr(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh-1));
    }
       
    return comp_contract(Dah, Axh, Ayh, Azh, alpha_ah, axh, ayh, azh, Dbh, Bxh, Byh, Bzh, alpha_bh, bxh, byh, bzh, Dch, Cxh, Cyh, Czh, alpha_ch, cxh, cyh, czh, Ddh, Dxh, Dyh, Dzh, alpha_dh, dxh, dyh, dzh);
    }


    

    /* CARTESIAN AXIS OF A, B, C AND D */
    Ax = xa[0];
    Ay = xa[1];
    Az = xa[2];
    Bx = xb[0];
    By = xb[1];
    Bz = xb[2];
    Cx = xc[0];
    Cy = xc[1];
    Cz = xc[2];
    Dx = xd[0];
    Dy = xd[1];
    Dz = xd[2];

    /* CARTESIAN ANGULAR MOMENTUM OF A, B, C AND D */
    ax = la[0];
    ay = la[1];
    az = la[2];
    bx = lb[0];
    by = lb[1];
    bz = lb[2];
    cx = lc[0];
    cy = lc[1];
    cz = lc[2];
    dx = ld[0];
    dy = ld[1];
    dz = ld[2];
 
    //Na = pow(2.*alpha_ah/M_PI,0.75)*pow(4.*alpha_ah,((axh+ayh+azh)/2.))/pow(Factorial2(2*axh-1)*Factorial2(2*ayh-1)*Factorial2(2*azh-1),0.5);
    //Nb = pow(2.*alpha_bh/M_PI,0.75)*pow(4.*alpha_bh,((bxh+byh+bzh)/2.))/pow(Factorial2(2*bxh-1)*Factorial2(2*byh-1)*Factorial2(2*bzh-1),0.5);
    //Nc = pow(2.*alpha_ch/M_PI,0.75)*pow(4.*alpha_ch,((cxh+cyh+czh)/2.))/pow(Factorial2(2*cxh-1)*Factorial2(2*cyh-1)*Factorial2(2*czh-1),0.5);
    //Nd = pow(2.*alpha_dh/M_PI,0.75)*pow(4.*alpha_dh,((dxh+dyh+dzh)/2.))/pow(Factorial2(2*dxh-1)*Factorial2(2*dyh-1)*Factorial2(2*dzh-1),0.5);
    //printf("%d %d %d %d %f\n", a,b,c,d,ERI);
    ERI = 0.;
    ERI = hrr(Da, Ax, Ay, Az, alpha_a, ax, ay, az,Db, Bx, By, Bz, alpha_b, bx, by, bz, Dc, Cx, Cy, Cz, alpha_c, cx, cy, cz, Dd, Dx, Dy, Dz, alpha_d, dx, dy, dz);                   
	//printf("%s %f\n","FINAL: ",ERI);

	return ERI;

}

#ifdef __cplusplus
} // __cplusplus defined.
#endif




