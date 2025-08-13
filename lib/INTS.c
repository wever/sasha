#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ints.h"

#ifdef __cplusplus
extern "C"
{
#endif
  
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

int a,b,c,d,ab,cd,p,q,r,s,i,j,k,l,k1,k2;
int l1,m1,n1,l2,m2,n2,l3,m3,n3,l4,m4,n4;
int nn;
double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz;
double Na,Nb,Nc,Nd;
double Px,Py,Pz,Qx,Qy,Qz,Prc;
double sum_alpha,sum_beta,e1,e2,E12,rate_alpha,rate_beta,AB,CD,PQ,U;

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
  
  
  
  void overlap(double *D, double *A, double *S, double *Center, int *CartAng, int N, int *lenL, int *locL) {

  double Oxyz(int la, int lb, double xa, double xb, double xd) {

  double O = 0;
  for (k1 = 0; k1 <= la; k1++) {
    for (k2 = 0; k2 <= lb; k2++) {
      if ((k1+k2)%2 == 0) {
        O += binom(la,k1)*binom(lb,k2)*pow(xd-xa,la-k1)*pow(xd-xb,lb-k2)*Factorial2(k1+k2-1)/(pow(2*sum_alpha,(k1+k2)/2)*sqrt(M_PI/sum_alpha));
      }
    }
  }
  return sqrt(M_PI/sum_alpha)*O;
  }

  for (a = 0; a < N; a++) {
    l1 = CartAng[3*a];
    m1 = CartAng[3*a+1];
    n1 = CartAng[3*a+2]; 
    Ax = Center[3*a];
    Ay = Center[3*a+1];
    Az = Center[3*a+2];
    for (p = 0; p < lenL[a]; p++) {
      Na = pow((pow(2,(2*(l1+m1+n1)+1.5))*pow(A[locL[a]+p],(l1+m1+n1+1.5))/(Factorial2(2*l1-1)*Factorial2(2*m1-1)*Factorial2(2*n1-1)*pow(M_PI,1.5))),0.5);
      for (b = a+1; b < N; b++) {
        l2 = CartAng[3*b];
        m2 = CartAng[3*b+1];
        n2 = CartAng[3*b+2]; 
	Bx = Center[3*b];
	By = Center[3*b+1];
	Bz = Center[3*b+2];
        for (q = 0; q < lenL[b]; q++) {
          Px = (A[locL[a]+p]*Ax+A[locL[b]+q]*Bx)/(A[locL[a]+p]+A[locL[b]+q]);
          Py = (A[locL[a]+p]*Ay+A[locL[b]+q]*By)/(A[locL[a]+p]+A[locL[b]+q]);
          Pz = (A[locL[a]+p]*Az+A[locL[b]+q]*Bz)/(A[locL[a]+p]+A[locL[b]+q]);
          Nb = pow((pow(2,(2*(l2+m2+n2)+1.5))*pow(A[locL[b]+q],(l2+m2+n2+1.5))/(Factorial2(2*l2-1)*Factorial2(2*m2-1)*Factorial2(2*n2-1)*pow(M_PI,1.5))),0.5);
          sum_alpha = (A[locL[a]+p] + A[locL[b]+q]);
          rate_alpha = A[locL[a]+p]*A[locL[b]+q]/sum_alpha;          
          AB = radius_pair(Ax,Ay,Az,Bx,By,Bz);
          S[N*a + b] += Na*Nb*D[locL[a]+p]*D[locL[b]+q]*pow(M_PI/sum_alpha,1.5)*exp(-rate_alpha*pow(AB,2))*Oxyz(l1,l2,Ax,Bx,Px)*Oxyz(m1,m2,Ay,By,Py)*Oxyz(n1,n2,Az,Bz,Pz);
      

        }
	if (absolute(S[N*a + b]) < 1e-6) {
	  S[N*a + b] = 0.;
	}
      }
    }
  }
}
  
  void dipole(double *D, double *A, double *Mu, double *Center, int *CartAng, int N, int *lenL, int *locL) {

  double Oxyz(int la, int lb, double xa, double xb, double xd) {

  double O = 0;
  for (k1 = 0; k1 <= la; k1++) {
    for (k2 = 0; k2 <= lb; k2++) {
      if ((k1+k2)%2 == 0) {
        O += binom(la,k1)*binom(lb,k2)*pow(xd-xa,la-k1)*pow(xd-xb,lb-k2)*Factorial2(k1+k2-1)/(pow(2*sum_alpha,(k1+k2)/2)*sqrt(M_PI/sum_alpha));
      }
    }
  }
  return sqrt(M_PI/sum_alpha)*O;
  }

  int inc[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  int axis;
  for (axis = 0; axis < 3; axis++) {
  for (a = 0; a < N; a++) {
    l1 = CartAng[3*a];
    m1 = CartAng[3*a+1];
    n1 = CartAng[3*a+2]; 
    //Ax = Center[3*a]*inc[0][axis];
    //Ay = Center[3*a+1]*inc[1][axis];
    //Az = Center[3*a+2]*inc[2][axis];
    Ax = Center[3*a];
    Ay = Center[3*a+1];
    Az = Center[3*a+2];
    for (p = 0; p < lenL[a]; p++) {
      Na = pow((pow(2,(2*(l1+m1+n1)+1.5))*pow(A[locL[a]+p],(l1+m1+n1+1.5))/(Factorial2(2*l1-1)*Factorial2(2*m1-1)*Factorial2(2*n1-1)*pow(M_PI,1.5))),0.5);
      for (b = 0; b < N; b++) {
        l2 = CartAng[3*b];
        m2 = CartAng[3*b+1];
        n2 = CartAng[3*b+2]; 
	//Bx = Center[3*b]*inc[0][axis];
	//By = Center[3*b+1]*inc[1][axis];
	//Bz = Center[3*b+2]*inc[2][axis];
	Bx = Center[3*b];
	By = Center[3*b+1];
	Bz = Center[3*b+2];
        for (q = 0; q < lenL[b]; q++) {
          Px = (A[locL[a]+p]*Ax+A[locL[b]+q]*Bx)/(A[locL[a]+p]+A[locL[b]+q]);
          Py = (A[locL[a]+p]*Ay+A[locL[b]+q]*By)/(A[locL[a]+p]+A[locL[b]+q]);
          Pz = (A[locL[a]+p]*Az+A[locL[b]+q]*Bz)/(A[locL[a]+p]+A[locL[b]+q]);
          Nb = pow((pow(2,(2*(l2+m2+n2)+1.5))*pow(A[locL[b]+q],(l2+m2+n2+1.5))/(Factorial2(2*l2-1)*Factorial2(2*m2-1)*Factorial2(2*n2-1)*pow(M_PI,1.5))),0.5);
          sum_alpha = (A[locL[a]+p] + A[locL[b]+q]);
          rate_alpha = A[locL[a]+p]*A[locL[b]+q]/sum_alpha;          
          AB = radius_pair(Ax,Ay,Az,Bx,By,Bz);
          Mu[N*N*axis+N*a+b] += Na*Nb*D[locL[a]+p]*D[locL[b]+q]*pow(M_PI/sum_alpha,1.5)*exp(-rate_alpha*pow(AB,2))*Oxyz(l1,l2+inc[0][axis],Ax,Bx,Px)*Oxyz(m1,m2+inc[1][axis],Ay,By,Py)*Oxyz(n1,n2+inc[2][axis],Az,Bz,Pz);
      

        }
	if (absolute(Mu[N*N*axis+N*a+b]) < 1e-6) {
	  Mu[N*N*axis+N*a+b] = 0.;
	}
      }
    }
  }
}
}
  
  void kei(double *D, double *A, double *KEI, double *Center, int *CartAng, int N, int *lenL, int *locL) {

  double SKxyz(int la, int lb, double Axyz, double Bxyz, double Pxyz) {
  
  double SK = 0.;
  int ii;
  int labint = (la+lb)/2;
  for (ii = 0; ii <= labint; ii++) {
    k = 2*ii;
    SK += fk(k,la,lb,(Pxyz - Axyz),(Pxyz - Bxyz))*Factorial2(2*ii-1)/pow(2.*sum_alpha,ii);
  }
  //printf("%8.5f\n", SK);
  return SK;
  }	  

  double Ixyz(int la, int lb, double Axyz, double Bxyz, double Pxyz) {
 
  double I;
  if ((la-1) < 0 && (lb-1) >= 0) {
    I = 0.5*(4*A[locL[a]+p]*A[locL[b]+q]*SKxyz(la+1,lb+1,Axyz,Bxyz,Pxyz)-2*A[locL[a]+p]*lb*SKxyz(la+1,lb-1,Axyz,Bxyz,Pxyz));
  }
  else if ((lb-1) < 0 && (la-1) >= 0) {
    I = 0.5*(4*A[locL[a]+p]*A[locL[b]+q]*SKxyz(la+1,lb+1,Axyz,Bxyz,Pxyz)-2*A[locL[b]+q]*la*SKxyz(la-1,lb+1,Axyz,Bxyz,Pxyz));
  }
  else if ((lb-1) < 0 && (la-1) < 0) {
    I = 0.5*(4*A[locL[a]+p]*A[locL[b]+q]*SKxyz(la+1,lb+1,Axyz,Bxyz,Pxyz));
  }
  else {
    I = 0.5*(la*lb*SKxyz(la-1,lb-1,Axyz,Bxyz,Pxyz)+4*A[locL[a]+p]*A[locL[b]+q]*SKxyz(la+1,lb+1,Axyz,Bxyz,Pxyz)-2*A[locL[a]+p]*lb*SKxyz(la+1,lb-1,Axyz,Bxyz,Pxyz)-2*A[locL[b]+q]*la*SKxyz(la-1,lb+1,Axyz,Bxyz,Pxyz));
  }
  return I;
  }

  for (a = 0; a < N; a++) {
    l1 = CartAng[3*a];
    m1 = CartAng[3*a+1];
    n1 = CartAng[3*a+2]; 
    Ax = Center[3*a];
    Ay = Center[3*a+1];
    Az = Center[3*a+2];
    for (p = 0; p < lenL[a]; p++) {
      Na = pow((pow(2,(2*(l1+m1+n1)+1.5))*pow(A[locL[a]+p],(l1+m1+n1+1.5))/(Factorial2(2*l1-1)*Factorial2(2*m1-1)*Factorial2(2*n1-1)*pow(M_PI,1.5))),0.5);
      for (b = a; b < N; b++) {
        l2 = CartAng[3*b];
        m2 = CartAng[3*b+1];
        n2 = CartAng[3*b+2]; 
	Bx = Center[3*b];
	By = Center[3*b+1];
	Bz = Center[3*b+2];
        for (q = 0; q < lenL[b]; q++) {
          Px = (A[locL[a]+p]*Ax+A[locL[b]+q]*Bx)/(A[locL[a]+p]+A[locL[b]+q]);
          Py = (A[locL[a]+p]*Ay+A[locL[b]+q]*By)/(A[locL[a]+p]+A[locL[b]+q]);
          Pz = (A[locL[a]+p]*Az+A[locL[b]+q]*Bz)/(A[locL[a]+p]+A[locL[b]+q]);
          Nb = pow((pow(2,(2*(l2+m2+n2)+1.5))*pow(A[locL[b]+q],(l2+m2+n2+1.5))/(Factorial2(2*l2-1)*Factorial2(2*m2-1)*Factorial2(2*n2-1)*pow(M_PI,1.5))),0.5);
          sum_alpha = (A[locL[a]+p] + A[locL[b]+q]);
          rate_alpha = A[locL[a]+p]*A[locL[b]+q]/sum_alpha;          
          AB = radius_pair(Ax,Ay,Az,Bx,By,Bz);
          KEI[N*a + b] += Na*Nb*D[locL[a]+p]*D[locL[b]+q]*pow(M_PI/sum_alpha,1.5)*exp(-rate_alpha*pow(AB,2))*(Ixyz(l1,l2,Ax,Bx,Px)*SKxyz(m1,m2,Ay,By,Py)*SKxyz(n1,n2,Az,Bz,Pz)+Ixyz(m1,m2,Ay,By,Py)*SKxyz(l1,l2,Ax,Bx,Px)*SKxyz(n1,n2,Az,Bz,Pz)+Ixyz(n1,n2,Az,Bz,Pz)*SKxyz(l1,l2,Ax,Bx,Px)*SKxyz(m1,m2,Ay,By,Py));
      

        }
      }
    }
  }
}
  
  void nai(double *D, double *A, double *NAI, double *Center, int *CartAng, double *R, int *Z, int N, int nnucleus, int *lenL, int *locL) {


  double Axyz() {

    double Xijk(int coord, int la, int lb, double Ai, double Bi, double Pi, double ri) {
    int l, q, Lhalf, t, lhalf;
    double Xijk = 0.;
    double first;
    for (l = 0; l < (la+lb+1); l++) {
      first = fk(l,la,lb,(Pi-Ai),(Pi-Bi))*factorial(l);
      lhalf = l/2.;
      for (q = 0; q < (lhalf+1); q++) {
        Lhalf = (l-2*q)/2.;
        for (t = 0; t < (Lhalf+1); t++) {
          if ((l-2*q-t)==coord) {
            Xijk += first/(factorial(q)*factorial(t)*factorial(l-2*q-2*t))*pow(-1,t)*pow(1./(4.*sum_alpha),q+t)*pow(ri-Pi,l-2*q-2*t);

	  }
        }  
      }
    }
    return Xijk;
  }
    
    double XXX = 0.;
    double Xx = 0.;
    double Xy = 0.;
    for (i = 0; i < (l1+l2+1); i++) {
      Xx = Xijk(i,l1,l2,Ax,Bx,Px,R[3*nn]);
      for (j = 0; j < (m1+m2+1); j++) {
        Xy = Xijk(j,m1,m2,Ay,By,Py,R[3*nn+1]);
        for (k = 0; k < (n1+n2+1); k++) {
          XXX += Xx*Xy*Xijk(k,n1,n2,Az,Bz,Pz,R[3*nn+2])*fgama(i+j+k,U);
          //printf("%8.5f\n", fgama(i+j+k,U));
        }
      }
    }
    return XXX; 
  }

  for (nn = 0; nn < nnucleus; nn++) {
    for (a = 0; a < N; a++) {
      l1 = CartAng[3*a];
      m1 = CartAng[3*a+1];
      n1 = CartAng[3*a+2]; 
      Ax = Center[3*a];
      Ay = Center[3*a+1];
      Az = Center[3*a+2];
    for (p = 0; p < lenL[a]; p++) {
      Na = pow((pow(2,(2*(l1+m1+n1)+1.5))*pow(A[locL[a]+p],(l1+m1+n1+1.5))/(Factorial2(2*l1-1)*Factorial2(2*m1-1)*Factorial2(2*n1-1)*pow(M_PI,1.5))),0.5);
      for (b = a; b < N; b++) {
        l2 = CartAng[3*b];
        m2 = CartAng[3*b+1];
        n2 = CartAng[3*b+2]; 
	Bx = Center[3*b];
	By = Center[3*b+1];
	Bz = Center[3*b+2];
        for (q = 0; q < lenL[b]; q++) {
          Px = (A[locL[a]+p]*Ax+A[locL[b]+q]*Bx)/(A[locL[a]+p]+A[locL[b]+q]);
          Py = (A[locL[a]+p]*Ay+A[locL[b]+q]*By)/(A[locL[a]+p]+A[locL[b]+q]);
          Pz = (A[locL[a]+p]*Az+A[locL[b]+q]*Bz)/(A[locL[a]+p]+A[locL[b]+q]);
          Nb = pow((pow(2,(2*(l2+m2+n2)+1.5))*pow(A[locL[b]+q],(l2+m2+n2+1.5))/(Factorial2(2*l2-1)*Factorial2(2*m2-1)*Factorial2(2*n2-1)*pow(M_PI,1.5))),0.5);
          sum_alpha = (A[locL[a]+p] + A[locL[b]+q]);
          rate_alpha = A[locL[a]+p]*A[locL[b]+q]/sum_alpha;          
          AB = radius_pair(Ax,Ay,Az,Bx,By,Bz);
	  Prc = radius_pair(Px,Py,Pz,R[3*nn],R[3*nn+1],R[3*nn+2]);
	  U = sum_alpha*pow(Prc,2);
          NAI[N*a + b] += -Z[nn]*Na*Nb*D[locL[a]+p]*D[locL[b]+q]*(2*M_PI/sum_alpha)*exp(-rate_alpha*pow(AB,2))*Axyz();
      

        }
      }
    }
  }
}

  }


void loopSCF(double *G, double *ERI, double *P, int N) {
  int count1 = 0;
  for (a = 0; a < N; a++) {
    for (b = 0; b < N; b++) {
	  int count2 = 0;
	  G[N*a+b] = 0.;
      for (c = 0; c < N; c++) {
        for (d = 0; d < N; d++) {
         G[N*a+b] += P[N*c+d]*(ERI[Index(a,b,d,c)] - 0.5*ERI[Index(a,c,d,b)]);
	 count2++;
        }
      }
      count1++;
    }
  }
}

#ifdef __cplusplus
} // __cplusplus defined.
#endif
