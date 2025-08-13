# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
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
"""

'''

AUTHOR: WEVERSON RODRIGUES GOMES
        GOMESWR@GMAIL.COM

************************************************************************************************************* *
****************************************************************************************************************
                                      HARTREE-FOCK SCF METHOD                                                 **
                                                                                                              **
IT IS A GENERIC ALPHA PROGRAM TO CALCULATE HF SCF FROM N ATOMS (JUST TESTED FOR H2O,HeO+ AND O2) USING STO-NG **
(N=1 TO 6). THE PROGRAM USES A ANALYTICAL SOLUTION FOR HIGHER ANGULAR MOMENTUM GAUSSIAN INTEGRALS.            **
                                                                                                              **
THERE IS NO PURPOSE IN ALGORITHMIC EFFICIENCY, BUT ONLY DIDACTIC.                                             **
                                                                                                              **
****************************************************************************************************************
***************************************************************************************************************

INFORMATIONS:
1 - THIS PROGRAM MUST BE ONLY USED FOR EDUCATIONAL PURPOSE
2 - TO RUN THE PROGRAM IS NECESSARY INSTALL PYTHON3 AND THE FOLLOW PACKAGES: NUMPY AND SCIPY
3 - ANY QUESTIONS OR ERROR SEND AN E-MAIL: gomeswr@gmail.com
4 - FOR VERY HIGHER ANGULAR MOMENTUM, IT IS NECESSARY CHANGE  BOYS FUNCTION ALGORITHM
5 - the program is extremely (EXTREMELY!!) inefficient from a computational point of view, but is quite didactic



REFERENCES:

*INTEGRALS:
[1] Clementi, E. & Davis, D. . Electronic structure of large molecular systems. Journal of Computational Physics 1, 223–244 (1966).  # Many printing erros
[2] Saunders, V. R. in Computational Techniques in Quantum Chemistry and Molecular Physics (eds. Diercksen, G. H. F., Sutcliffe, B. T. & Veillard, A.) 347–424 (Springer Netherlands, 1975)
[3] Petersson, T. & Hellsing, B. A detailed derivation of Gaussian orbital-based matrix elements in electron structure calculations. Eur. J. Phys. 31, 37 (2010).
[4] Hô, M. & Hernández-Pérez, J.-M. Evaluation of Gaussian Molecular Integrals. The Mathematica Journal 14, (2012).
[5] Hô, M. & Hernández-Pérez, J.-M. Evaluation of Gaussian Molecular Integrals II. The Mathematica Journal 15, (2013).
[6] Helgaker, T., Jorgensen, P. & Olsen, J. Molecular Electronic-Structure Theory. (Wiley, 2013).
[7] Hô, M. & Hernández-Pérez, J.-M. Evaluation of Gaussian Molecular Integrals III. The Mathematica Journal 16, (2014).

*BOYS FUNCTION:
[8] Matsuoka, O. Field and field gradient integrals based on Gaussian type orbitals. Computer Physics Communications 3, 130–135 (1972). # There is an error in equation 10.
[9] McMurchie, L. E. & Davidson, E. R. One- and two-electron integrals over cartesian gaussian functions. Journal of Computational Physics 26, 218–231 (1978).
[10] Gill, P. M. W., Johnson, B. G. & Pople, J. A. Two-electron repulsion integrals over Gaussian s functions. Int. J. Quantum Chem. 40, 745–752 (1991).
[11] Helgaker, T., Jorgensen, P. & Olsen, J. Molecular Electronic-Structure Theory. (Wiley, 2013).
[12] Ref 2

*SCF:
[13] Szabo, A. & Ostlund, N. S. Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory. (Dover Publications, 1996).
[14] Ramachandran, K. I., Deepa, G. & Namboori, K. Computational Chemistry and Molecular Modeling: Principles and Applications. (Springer, 2008).
[15] Levine, I. N. Quantum Chemistry. (Prentice Hall, 2013).


INTRODUCTION:

The problem is to solve integrals of the general form:

INT(XaOXb) from -inf to inf

Xa(r,alpha,A,a) = (x-Ax)**l*(y-Ay)**m*(z-Az)**n*exp(-alpha|r-A|**2)

Xa = unnormalized Cartesian Gaussian function at nucleus A={Ax,Ay,Az}, alpha is the orbital exponent and the polynomial represents the angular part. The sum of the Cartesian angular momenta ax+ay+az = 0,1,2,... corresponds to functions of type s,p,d,f,....

O = 1 --> Overlap/Density integral (Sab).
O = -1/2Laplaciano**2 --> Energy operator for kinetic energy (Tab).
O = |r-R|**-1 --> Electron-nuclear attraction (Vab).
O = |ra-rb|**-1 --> Electron-electron repulsion (which wold involve double integrals) (ab|cd).
'''

def hartree_fock(N,R,Z,A,D,CartAng,Center,n_electrons,n,lenL,locL,lenAtomN,atoms,mass,nos_initial,nos_final,method,basisset,charge,units,inputName,fnCube,checkprint,plot_type,logger):

    #from scipy.special import binom # THE BINOMIAL COEFFICIENT = n!/(k!(n-k)!)
    #from scipy.misc import factorial2  # DOUBLE FACTORIAL e.g: -1!! = 1
    from numpy import array,zeros,empty,exp,sqrt,pi,set_printoptions,arange,copy,dot,transpose,diag,fill_diagonal,shape,linspace,meshgrid,double,intc,asarray,prod,trace
    from numpy.linalg import eigh
    #from pylab import show,contour,cm,clabel,title,figure
    #import numpy as nm
    from math import factorial
    from time import clock
    import numpy.ctypeslib as npc
    import ctypes as ct
    from itertools import chain
    import sys
    import os
    #from get_basis import inp

    # determine if application is a script file or frozen exe
    if getattr(sys, 'frozen', False):
        exedir= os.path.dirname(sys.executable)
    elif __file__:
        exedir = os.path.dirname(os.path.abspath(__file__))
    #exedir = os.path.dirname(os.path.realpath(__file__))

    if sys.platform.startswith('win32'):
        libPath = exedir+"\lib\_INTS.dll"
        libHPPath = exedir+"\lib\_HPINTS.dll"
        libCubePath = exedir+"\lib\_CUBE.dll"
    elif sys.platform.startswith('linux'):
        libPath = exedir+"/lib/_INTS.so"
        libHPPath = exedir+"/lib/_HPINTS.so"
        libCubePath = exedir+"/lib/_CUBE.so"
    elif sys.platform.startswith('darwin'):
        libPath = exedir+"/lib/_INTS.dylib"
        libHPPath = exedir+"/lib/_HPINTS.dylib"
        libCubePath = exedir+"/lib/_CUBE.dylib"

    #N,R,Z,A,D,CartAng,Center,n_electrons,n,lenL,locL = inp(basisset,units,inputName)

    n_electrons = n_electrons - charge

    if n_electrons < 0:
        logger.error('Wrong Charge Value')
        print("NUMBER OF ELECTRONS: ", n_electrons)
        raise RuntimeError('WRONG CHARGE VALUE')

    if units == "Bohr":
        Rinput = R
    elif units =="Angstroms":
        Rinput = R/1.889726

    print('{0:1s} {1:>55s}'.format("","***********************************************"))
    printformatted = "%s" % ("METHOD:")
    logger.info(printformatted)
    printformatted = "%s %s %s %d\n\n\n" % ("Hartree-Fock/",basisset," charge:",charge)
    logger.info(printformatted)

    print("Input(",units,")")
    for numele in range(n):
        print('{0:10s} {1:10.5f} {2:10.5f} {3:10.5f}'.format(atoms[numele],Rinput[numele,0],Rinput[numele,1],Rinput[numele,2]))

    print('{0:1s} {1:>50s}'.format("\n","***********************************"))
    print("\n\n")



    # Cartesian Gaussian function = (x-Ax)**l*(y-Ay)**m*(z-Az)**n*exp(-ai|r-A|**2)

    def fk(k,la,lb,r1,r2):
        f = 0
        for i in range(la+1):
            for j in range(lb+1):
                if i+j==k:
                    f += binom(la,i)*binom(lb,j)*r1**(la-i)*r2**(lb-j)
                else:
                    pass
        return f

    ############################################### OVERLAPE INTEGRAL #################################################
    def S():
        ts0 = clock()

        printformatted = "%s" % ("OVERLAP INTEGRALS STARTED...")

        logger.info(printformatted)


        S = zeros([N,N],float)

        _overlap =   npc.load_library(libPath,'.')   # . load path
        _overlap.overlap.argtypes = [npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=2),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = intc,ndim=2),ct.c_int,npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1)]
        _overlap.overlap.restype = ct.c_void_p
        CartAng2 = CartAng.astype(intc)
        DD = array(list(chain.from_iterable(D)),float)
        AA = array(list(chain.from_iterable(A)),float)
        lenLL = asarray(lenL,dtype=intc)
        locLL = asarray(locL,dtype=intc)

        _overlap.overlap(DD,AA,S,Center,CartAng2,N,lenLL,locLL)


        S = S + S.T
        fill_diagonal(S,1.)
        printformatted = "%s\n" % ("OVERLAP INTEGRALS DONE")

        logger.info(printformatted)

        global ts
        ts = clock() - ts0
        return S


    ############################################### DIPOLE INTEGRAL #################################################
    def Muint():
        ts0 = clock()

        printformatted = "%s" % ("DIPOLE INTEGRALS STARTED...")

        logger.info(printformatted)


        Mu = zeros([3,N,N],float)

        _dipole = npc.load_library(libPath,'.')   # . load path
        _dipole.dipole.argtypes = [npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=3),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = intc,ndim=2),ct.c_int,npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1)]
        _dipole.dipole.restype = ct.c_void_p
        CartAng2 = CartAng.astype(intc)
        DD = array(list(chain.from_iterable(D)),float)
        AA = array(list(chain.from_iterable(A)),float)
        lenLL = asarray(lenL,dtype=intc)
        locLL = asarray(locL,dtype=intc)

        _dipole.dipole(DD,AA,Mu,Center,CartAng2,N,lenLL,locLL)


        #Mu = Mu + Mu.T
        printformatted = "%s\n" % ("DIPOLE INTEGRALS DONE")

        logger.info(printformatted)

        global tm
        tm = clock() - ts0
        return Mu

############################################### KINETIC-ENERGY INTEGRALS #############################################

    def KEI():
        ts0 = clock()

        printformatted = "%s" % ("KINETIC-ENERGY INTEGRALS STARTED...")

        logger.info(printformatted)


        KEI = zeros([N,N],float)

        _kei =   npc.load_library(libPath,'.')   # . load path
        _kei.kei.argtypes = [npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=2),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = intc,ndim=2),ct.c_int,npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1)]
        _kei.kei.restype = ct.c_void_p
        CartAng2 = CartAng.astype(intc)
        DD = array(list(chain.from_iterable(D)),float)
        AA = array(list(chain.from_iterable(A)),float)
        lenLL = asarray(lenL,dtype=intc)
        locLL = asarray(locL,dtype=intc)

        _kei.kei(DD,AA,KEI,Center,CartAng2,N,lenLL,locLL)
        KEI = KEI + KEI.T - diag(KEI.diagonal())

        printformatted = "%s\n" % ("KINETIC-ENERGY INTEGRALS DONE")

        logger.info(printformatted)

        global tk
        tk = clock() - ts0
        return KEI


############################################### NUCLEAR-ATTRACTION INTEGRAL #############################################

    def NAI():
        ts0 = clock()

        printformatted = "%s" % ("THREE-CENTER NUCLEAR ATTRACTION INTEGRALS STARTED...")

        logger.info(printformatted)

        NAI = zeros([N,N],float)

        _nai =   npc.load_library(libPath,'.')   # . load path
        _nai.nai.argtypes = [npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=2),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = intc,ndim=2),npc.ndpointer(dtype=double,ndim=2),npc.ndpointer(dtype = intc,ndim=1),ct.c_int,ct.c_int,npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1)]
        _nai.nai.restype = ct.c_void_p
        CartAng2 = CartAng.astype(intc)
        Z2 = Z.astype(intc)
        DD = array(list(chain.from_iterable(D)),float)
        AA = array(list(chain.from_iterable(A)),float)
        #RR = array(list(chain.from_iterable(R)),float)
        lenLL = asarray(lenL,dtype=intc)
        locLL = asarray(locL,dtype=intc)

        _nai.nai(DD,AA,NAI,Center,CartAng2,R,Z2,N,n,lenLL,locLL)
        NAI = NAI + NAI.T - diag(NAI.diagonal())

        printformatted = "%s\n" % ("THREE-CENTER NUCLEAR ATTRACTION INTEGRALS DONE")

        logger.info(printformatted)

        global tn
        tn = clock() - ts0
        return NAI

    def ERI2():
        ERI = zeros([N**4],float)

        _twoe = npc.load_library('./_INTS.so','.')   # . load path
        _twoe.twoe.argtypes = [npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=1),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = intc,ndim=2),ct.c_int,npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1)]
        _twoe.twoe.restype = ct.c_void_p
      #c = A.astype(double)
        CartAng2 = CartAng.astype(intc)
        DD = array(list(chain.from_iterable(D)),float)
        AA = array(list(chain.from_iterable(A)),float)
        lenLL = asarray(lenL,dtype=intc)
        locLL = asarray(locL,dtype=intc)

        _twoe.twoe(DD,AA,ERI,Center,CartAng2,N,lenLL,locLL)
       # print(ERI)
        #print('TWO-ELECTRON REPULSION INTEGRALS DONE')
        return ERI


    ###############################################  ELECTRON REPULTION INTEGRAL #############################################

    def ERI():
        ts0 = clock()

        printformatted = "%s" % ("TWO-ELECTRON REPULSION INTEGRALS STARTED...")

        logger.info(printformatted)

        ERI = zeros([N**4],float)


        _twoe = npc.load_library(libHPPath,'.')   # . load path
        _twoe.twoe.argtypes = [ct.c_int,ct.c_int,ct.c_int,ct.c_int,npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),ct.c_int,ct.c_int,ct.c_int,ct.c_int]
        _twoe.twoe.restype = ct.c_double


        #CartAng2 = CartAng.astype(intc)
        lenLL = asarray(lenL,dtype=intc)
        locLL = asarray(locL,dtype=intc)

        for a in range(N):
            for b in range(a+1):
                ab = a*(a+1)/2+b
                for c in range(N):
                    for d in range(c+1):
                        cd = c*(c+1)/2+d;
                        if ab <= cd:
                            ERI[Index(a,b,c,d)] = _twoe.twoe(a,b,c,d,asarray(A[a]),asarray(A[b]),asarray(A[c]),asarray(A[d]),asarray(D[a]),asarray(D[b]),asarray(D[c]),asarray(D[d]),CartAng[a].astype(intc),CartAng[b].astype(intc),CartAng[c].astype(intc),CartAng[d].astype(intc),Center[a],Center[b],Center[c],Center[d],lenL[a],lenL[b],lenL[c],lenL[d])
            #print(a)

        printformatted = "%s\n" % ("TWO-ELECTRON REPULSION INTEGRALS DONE")

        logger.info(printformatted)

        global te
        te = clock() - ts0
        return ERI

    ## BOYS FUNCTION Fm(T) = int(u**(2*m)*exp(-Tu**2)), m = 0,1,2,....
    # FOR BEST ALGORITHM USING CHEBYSHEV INTERPOLATION, SEE REF 10

    # TABLE OF Fni(T_LINE) FOR A SERIES OF VALUES OF THE ARGUMENT, T_line (AT INTERVALS OF 0.1). THE VALUE
    # OF Fni(T) MAY THEN BE OBTAINED (BY FUNCTION fgama(ni,T)) FROM SUCH A TABLE BY INTERPOLATION.
    def fgama_tab():
        mmax = 22
        def fgama_mmax():
            thres = 1e-8
            m = mmax
            Fmmax_T = zeros([mmax*10+1])
            for T_line in arange(0.1,mmax+0.1,0.1):
                i = 0
                summ = 1
                while summ > thres:
                    term = 1
                    for j in arange(2*m+1,2*m+2*i+1+1,2):
                        term *= float(j)
                    summ = (2*T_line)**i/term
                    Fmmax_T[(T_line)*10] += (2*T_line)**i/term
                    i +=1
                Fmmax_T[T_line*10] = exp(-T_line)*Fmmax_T[T_line*10]
            return Fmmax_T
        Ftab = zeros([mmax+1,mmax*10+1])
        Ftab[mmax] = fgama_mmax()
        Ttab = array(arange(0,mmax+0.1,0.1))
        for m in range(mmax-1,-1,-1):
            Ftab[m] = (2*Ttab*Ftab[m+1]+exp(-Ttab))/(2*m+1)
        return Ftab

    # BASED ON REF 9
    def fgama(ni,T):
        if ni >= 16:
          print('ni = ',ni,'HIGH VALUE FOR NI')
        Fni = zeros([ni+1])
        if T<1e-6:
            Fni[ni] = 1/(2*ni+1)
        elif T>=1e-6 and T<=12:
            mmax = 22
            Ttab = array(arange(0,mmax+0.1,0.1))
            value_near = abs(T-Ttab).argmin()
            T_line = Ttab[value_near]
            for k in range(7):
                Fni[ni] += Ftab[ni+k,T_line*10]*(T_line-T)**k/factorial(k)
        elif T>12 and T<30:
            if T<=15:
                g = 0.4999489092 - 0.2473631686/T + 0.321180909/T**2 - 0.3811559346/T**3
            elif T>15 and T<18:
                g = 0.4998436875 - 0.24249438/T + 0.24642845/T**2
            elif T>=18 and T<24:
                g = 0.499093162 - 0.2152832/T
            else:
                g = 0.490
            Fni[0] = 0.5*sqrt(pi/T) - exp(-T)*g/T
            for m in range(1,ni+1):
                Fni[m] = ((2*m-1)*Fni[m-1]-exp(-T))/(2*T)
        else:
            Fni[0] = 0.5*sqrt(pi/T)
            for m in range(1,ni+1):
                Fni[m] = ((2*m-1)*Fni[m-1]-exp(-T))/(2*T)
        return Fni[ni]

    # SAME AS FGAMA(NI,T), BASED ON REF 2 (IT IS SLOWLY THAN fgama(ni,T)).
    def fgama2(ni,T):
        if T<1e-6:
            return  1/(2*ni+1)
        else:
            c = factorial2(2*ni-1)/(2*T)**ni
            F1 = sqrt(pi/(4*T))*erf(sqrt(T))
            F2 = 0
            for k in range(1,ni+1):
                F2 += (2*T)**(k-1)/factorial2(2*k-1)
            return c*(F1-exp(-T)*F2)


    def radius_pair(Ra,Rb):
        return sqrt(sum((Ra-Rb)**2))


    def position_atom(A):
        return sqrt(sum((A)**2))

    def Index(i,j,k,l):
        if (i<j): return Index(j,i,k,l);
        if (k<l): return Index(i,j,l,k);
        ij = i*(i+1)/2+j
        kl = k*(k+1)/2+l
        if (ij<kl):
            tmp = ij
            ij = kl
            kl = tmp
        return ij*(ij+1)/2+kl


    def print_Matrix(x,inputfile):
        fileintegrals = open(inputfile,'w')
        j = N//5+1
        i = o = e = 0
        while o<j:
            if i<N:
                if i+4<N:
                    print('%7i %13i %13i %13i %13i' % (e+1,e+2,e+3,e+4,e+5),file=fileintegrals)
                else:
                    if N%5==1:
                        print('%7i' % (e+1),file=fileintegrals)
                    elif N%5==2:
                        print('%7i %13i' % (e+1,e+2),file=fileintegrals)
                    elif N%5==3:
                        print('%7i %13i %13i' % (e+1,e+2,e+3),file=fileintegrals)
                    elif N%5==4:
                        print('%7i %13i %13i %13i' % (e+1,e+2,e+3,e+4),file=fileintegrals)
                    elif N%5==0:
                        print('%7i %13i %13i %13i %13i' % (e+1,e+2,e+3,e+4,e+5),file=fileintegrals)
            if i<N:
                print('%3i %13.6e' % (e,x[i,e]),file=fileintegrals)
                i+=1
            if i<N:
                print('%13.6e %13.6e' % (x[i,e],x[i,e+1]),file=fileintegrals)
                i+=1
            if i<N:
                print('%13.6e %13.6e %13.6e' % (x[i,e],x[i,e+1],x[i,e+2]),file=fileintegrals)
                i+=1
            if i<N:
                print('%13.6e %13.6e %13.6e %13.6e' % (x[i,e],x[i,e+1],x[i,e+2],x[i,e+3]),file=fileintegrals)
                i+=1
                k = i
            if i<N:
                while k<N:
                    print('%13.6e %13.6e %13.6e %13.6e %13.6e' % (x[k,e],x[k,e+1],x[k,e+2],x[k,e+3],x[k,e+4]),file=fileintegrals)
                    k+=1
            o +=1
            e += 5
            i += 1
        fileintegrals.close()
        return '\n'

    def print_Matrix2(x,inputfile):
        fileintegrals = open(inputfile,'w')
        e = 0
        count = N
        for i in range(N//5):
            print('%13i %15i %15i %15i %15i' % (e+1,e+2,e+3,e+4,e+5),file=fileintegrals)
            for k in range(N):
                print('%3i %15.6e %15.6e %15.6e %15.6e %15.6e' % (k+1,x[k,e],x[k,e+1],x[k,e+2],x[k,e+3],x[k,e+4]),file=fileintegrals)
            e += 5
        for l in range(1):
            if N%5==1:
                print('%13i' % (e+1),file=fileintegrals)
                for k in range(N):
                    print('%3i %15.6e' % (k+1,x[k,e]),file=fileintegrals)
            elif N%5==2:
                print('%13i %15i' % (e+1,e+2))
                for k in range(N):
                    print('%3i %15.6e %15.6e' % (k+1,x[k,e],x[k,e+1]),file=fileintegrals)
            elif N%5==3:
                print('%13i %15i %15i' % (e+1,e+2,e+3),file=fileintegrals)
                for k in range(N):
                    print('%3i %15.6e %15.6e %15.6e' % (k+1,x[k,e],x[k,e+1],x[k,e+2]),file=fileintegrals)
            elif N%5==4:
                print('%13i %15i %15i %15i' % (e+1,e+2,e+3,e+4),file=fileintegrals)
                for k in range(N):
                    print('%3i %15.6e %15.6e %15.6e %15.6e' % (k+1,x[k,e],x[k,e+1],x[k,e+2],x[k,e+3]),file=fileintegrals)
            else:
                pass

        fileintegrals.close()
        return '\n'

    def printtime():
        printformatted = "%s %f %s" % (" OVERLAP INTEGRALS CPU time = ",ts,"seconds")

        logger.info(printformatted)

        printformatted = "%s %f %s" % (" DIPOLE MOMENT  INTEGRALS CPU time = ",tm,"seconds")

        logger.info(printformatted)

        printformatted = "%s %f %s" % (" KINETIC-ENERGY INTEGRALS CPU time = ",tk,"seconds")

        logger.info(printformatted)

        printformatted = "%s %f %s" % (" THREE-CENTER NUCLEAR ATTRACTION INTEGRALS CPU time = ",tn,"seconds")

        logger.info(printformatted)

        printformatted = "%s %f %s" % (" TWO-ELECTRON REPULSION INTEGRALS CPU time = ",te,"seconds")

        logger.info(printformatted)

        printformatted = "%s %f %s" % (" SCF CPU time = ",clock() - t1,"seconds")

        logger.info(printformatted)

        printformatted = "%s %f %s" % (" TOTAL CPU time = ",clock() - t0,"seconds")

        logger.info(printformatted)



    t0 = clock()

    import parameters

    # CONVERGENCE CRITERION FOR_DENSITY MATRIX
    convergence_criterion = parameters.convergence()


    if type(n_electrons) != int or n_electrons%2 != 0:
        logger.error('THE NUMBER OF ELECTRONS MUST BE AN EVEN INTEGER NUMBER')
        #raise RuntimeError('THE NUMBER OF ELECTRONS MUST BE AN EVEN INTEGER NUMBER')
        return False,False,False,False,False,False,True

    # MAXIMUM NUMBER OF CYCLES FOR SCF
    maxcycles = parameters.scfmaxcycles()

    set_printoptions(suppress=True,precision=6)

    #Ftab = fgama_tab()

    # PRINTING TWO ELECTRON INTEGRALS

    S = S()
    if checkprint[0]:
        print_Matrix(S,inputName.split('.')[0]+'_'+basisset+'_OverlapMatrix.txt')
    Mu = Muint()
    T = KEI()
    if checkprint[1]:
        print_Matrix(T,inputName.split('.')[0]+'_'+basisset+'_KineticMatrix.txt')
    V = NAI()
    if checkprint[2]:
        print_Matrix(V,inputName.split('.')[0]+'_'+basisset+'_NuclearAttractionMatrix.txt')
    twoe_matrix = ERI()
    if checkprint[3]:
        fileintegrals = open(inputName.split('.')[0]+'_'+basisset+'_Two_eMatrix.txt','w')
        for a in range(N):
             for b in range(a+1):
                 ab = a*(a+1)/2+b
                 for c in range(N):
                     for d in range(c+1):
                         cd = c*(c+1)/2+d;
                         if ab <= cd:
                             print('%3i %5i %5i %5i, %5.5f' % (a,b,c,d,twoe_matrix[Index(a,b,c,d)]), file=fileintegrals)
                             #print('%3i %5i %5i %5i, %5.5f' % (a,b,c,d,twoe_matrix[Index(a,b,c,d)]))
        fileintegrals.close()
    # COMPUTING CORE-HAMILTONIAL MATRIX --> H_CORE = T + V1 + V2
    H_core =  T + V

    # ORTHOGONAL MATRIX TRANSFORMATION
    val, vec = eigh(S)
    s_transf = diag(val**-0.5)

    # EQUATION 3.167
    # SYMMETRIC ORTHOGONALIZATION
    #X = dot(vec,dot(s_transf,transpose(vec)))

    # EQUATION 3.169
    # CANONICAL ORTHOGONALIZATION
    X = dot(vec,s_transf)

    G = zeros([N,N],float)
    F = zeros([N,N],float)

    # INITIAL GUESS --> DENSITY MATRIX = 0
    # VERY POOR GUESS
    P = zeros([N,N],float)

    count = 1
    delta_error = 1

    # FROM HERE, EACH EQUATION WILL BE REPORT ACCORDING TO SZABO/OSTLUND BOOK

    # EQUATION 3.154 - TWO ELECTRON PART OF THE FOCK MATRIX
    # USE CORE-HAMILTONIAN FOR INITIAL GUESS AT F, I.E. (P=0)
    _loopSCF = npc.load_library(libPath,'.')   # . load path
    _loopSCF.loopSCF.argtypes = [npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=2),ct.c_int]
    _loopSCF.loopSCF.restype = ct.c_void_p

    print("\n\n")
    print("%37s" % "--------------------")
    print("%31s" % "SCF STEPS")
    print("%37s" % "--------------------")

    printformatted = '{0:1s} {1:20s} {2:20s} {3:20s}'.format("","CYCLE","ENERGY","DELTA-E")


    #logger.error(printformatted)
    print(printformatted)

    t1 = clock()
    countlist = []
    ENlist = []



    F = H_core + G
    F_line = dot(transpose(X),dot(F,X))
    e,C_line = eigh(F_line,UPLO='U')
    C = dot(X,C_line)
    # EQUATION 3.145 - DENSITY MATRIX OR CHARGE-DENSITY BOND-ORDER MATRIX
    for ad in range(N):
        for bd in range(N):
            for aa in range(n_electrons//2):
                P[ad,bd] += C[ad,aa]*C[bd,aa]
    P = 2*P


    while delta_error > convergence_criterion:


        _loopSCF.loopSCF(G,twoe_matrix,P,N)
        #print(twoe_matrix)


        # EQUATION 3.154 - FINAL EXPRESSION FOR THE FOCK MATRIX - ONE ELECTRON PART (H_core), WHICH IS FIXED GIVEN THE BASIS SET, AND TWO ELECTRON PART G, WHICH DEPENDS ON THE DENSITY MATRIX P
        F = H_core + G

        # EQUATION 3.274 - VARIATIONAL VALUE OF ELECTRONIC ENERGY
        EN = 0
        for y in range(N):
            for z in range(N):
                EN += P[z,y]*(H_core[y,z] + F[y,z])
        EN = 0.5*EN


        # EQUATION 3.185 - NUCLEAR-NUCLEAR REPULSION PLUS ELECTRONIC ENERGY (E0) EQUAL TO TOTAL ENERGY (ET)
        # NUCLEAR REPULTION
        #NR = 0
        #for AA in range(n):
        #    for BB in range(AA+1,n):
        #        NR += Z[AA]*Z[BB]/radius_pair(R[AA],R[BB])

        # TOTAL ENERGY
        #ENT = EN + NR



        #print(ENT)
        # EQUATION 3.177
        F_line = dot(transpose(X),dot(F,X))

        # EQUATION 3.174 - 'e' IS THE EINGENVALUE AND C_LINE IS THE EIGENVECTOR
        # MATRIX E IS THE ORBITAL ENERGIES FOR OCCUPIED AND VIRTUAL
        # THE FUNCTION EIGH CALCULATES THE EIGENVALUES AND EIGENVECTOR OF A REAL SYMMETRIC OR HERMITIAN MATRIX USING QR ALGORITHM
        # FOR MORE INFORMATION ABOUT QR DECOMPOSITION SEE Newman, M. Computational Physics. (CreateSpace Independent Publishing Platform, 2012)
        e,C_line = eigh(F_line,UPLO='U')

        # EQUATION 3.174
        C = dot(X,C_line)

        oldP = copy(P)
        P = zeros([N,N],float)

        # EQUATION 3.145 - DENSITY MATRIX OR CHARGE-DENSITY BOND-ORDER MATRIX
        for ad in range(N):
            for bd in range(N):
                for aa in range(n_electrons//2):
                    P[ad,bd] += C[ad,aa]*C[bd,aa]
        P = 2*P

        delta_error = 0
        for i in range(N):
            for j in range(N):
                delta_error += (P[i,j] - oldP[i,j])**2

        delta_error = sqrt(delta_error/4)

        if count>maxcycles:

            printtime()
            logger.warning('Maximum number of SCF cycles exceed')
            "TOTAL CPU time = ",clock() - t0,"seconds"
            #raise RuntimeError('Maximum number of SCF cycles exceed')
            return False,False,False,False,False,False,True

##########################################################################
        #printformatted = "%4d %30.6f %30.6f" % (count, ENT, delta_error)
        printformatted = '{0:5d} {1:24f} {2:20f}'.format(count, EN, delta_error)

        #logger.info(printformatted)
        print(printformatted)
###########################################################################

        countlist.append(count)
        ENlist.append(EN)
        count += 1

    print("\n\n\n")
    print("%37s" % "--------------------")
    print("%35s" % "ORBITAL ENERGIES")
    print("%37s" % "--------------------")

    esort = sorted(e)
    print('{0:12s} {1:3s} {2:18s} {3:20s}'.format("","No","","ENERGY/a.u"))
    for i in range(len(esort)):
        if i< n_electrons/2:
            occvirt = 'occ.'
        else:
            occvirt = 'virt.'
        print("%15d %10s %18.5f" % (i,occvirt,esort[i]))


    print("\n\n\n")
    print("%37s" % "--------------------")
    print("%37s" % "MULLIKEN POPULATION")
    print("%31s" % "ANALYSIS")
    print("%37s" % "--------------------")
    mulliken = dot(transpose(P),S)
    diagmulliken = diag(mulliken)


    # EQUATION 3.195
    at = 0
    sumatomcharge = 0.
    for i,atom in enumerate(atoms):
        Matomcharge = 0.
        for l in range(lenAtomN[i]):
            Matomcharge += diagmulliken[at]
            at += 1
        atomcharge = Z[i] - Matomcharge
        sumatomcharge += atomcharge
        print("%18s %18.5f" % (atom,atomcharge))
    print("\n %32s %1.5f" % ("SUM OF ATOMIC CHARGES: ",sumatomcharge))

    print("\n\n\n")
    print("%37s" % "--------------------")
    print("%33s" % "DIPOLE MOMENT")
    print("%37s" % "--------------------")

    # CENTER OF MASS OF MOLECULE
    upfor = 0.
    downfor = 0.
    upfor  = zeros([3],float)
    for axisnuc in range(3):
        for i in range(n):
            upfor[axisnuc] += mass[i]*R[i,axisnuc]

    COM = upfor/sum(mass)

    '''
    COM  = zeros([3],float)
    for axisnuc in range(3):
        for nuci in range(n):
            COM[axisnuc] += Z[nuci]*R[nuci,axisnuc]
    COM = COM/sum(Z)
    print('COM',COM)
    '''

    # CONTRIBUCTION (CLASSICAL) OF THE NUCLEI
    # EQUATION 3.190
    nucdip = zeros([3],float)
    for axisnuc in range(3):
        for i,prot in enumerate(Z):
            nucdip[axisnuc] += prot*(R[i,axisnuc]-COM[axisnuc])
    print("\n")

    # EQUATION 3.192
    Mxyz = zeros([3,N,N],float)
    Dipxyz = zeros([3],float)
    for axisdip in range(3):
        for mi in range(N):
            for ni in range(N):
                Mxyz[axisdip,mi,ni] += Mu[axisdip,mi,ni] + (Center[ni,axisdip] - COM[axisdip])*S[mi,ni]

    # EQUATION 3.190
    for axisdip in range(3):
        for mi in range(N):
            for ni in range(N):
                Dipxyz[axisdip] -= P[mi,ni]*Mxyz[axisdip,ni,mi]
                if abs(Dipxyz[axisdip]) < 1e-6:
                    Dipxyz[axisdip] = 0.

    total_dipxyz = Dipxyz+nucdip
    total_dip = sqrt(sum((Dipxyz+nucdip)**2))

    print('{0:>33} {1:>12} {2:>12}'.format("X","Y","Z"))
    print('{0:1} {1:10.5} {2:12.5} {3:12.5}'.format("Electronic contribution: ",Dipxyz[0],Dipxyz[1],Dipxyz[2]))
    print('{0:1} {1:10.5} {2:12.5} {3:12.5}'.format("Nuclear contribution   : ",nucdip[0],nucdip[1],nucdip[2]))
    print("%s\n" % (                                "                        --------------------------------------"))
    print('{0:1} {1:10.5} {2:12.5} {3:12.5}'.format("Total Dipole Moment    : ",total_dipxyz[0],total_dipxyz[1],total_dipxyz[2]))
    print("%s\n" % (                                "                        --------------------------------------"))
    print('{0:1} {1:10.5}'.format("Magnitude (a.u.)       : ",total_dip))
    print('{0:1} {1:10.5}'.format("Magnitude (Debye)      : ",total_dip/0.393456))

    # EQUATION 3.274 - VARIATIONAL VALUE OF ELECTRONIC ENERGY
    EN = 0
    for bb in range(N):
        for zz in range(N):
            EN += P[bb,zz]*(H_core[bb,zz] + F[bb,zz])
    EN = 0.5*EN

    # EQUATION 3.185 - NUCLEAR-NUCLEAR REPULSION PLUS ELECTRONIC ENERGY (E0) EQUAL TO TOTAL ENERGY (ET)
    # NUCLEAR REPULTION
    NR = 0
    for AA in range(n):
        for BB in range(AA+1,n):
            NR += Z[AA]*Z[BB]/radius_pair(R[AA],R[BB])

    # TOTAL ENERGY
    ENT = EN + NR

    print("\n\n\n")

    printformatted = "%s %f %s %d %s" % ("Electronic Energy = ",EN,"a.u after",count-1,"cycles")

    logger.info(printformatted)

    printformatted = "%s %f %s %s %s %d %s\n\n\n" % ("Total Energy = ",ENT,"a.u"," (Electronic + Nuclear Energy) ","after",count-1,"cycles")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" OVERLAP INTEGRALS CPU time = ",ts,"seconds")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" DIPOLE MOMENT  INTEGRALS CPU time = ",tm,"seconds")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" KINETIC-ENERGY INTEGRALS CPU time = ",tk,"seconds")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" THREE-CENTER NUCLEAR ATTRACTION INTEGRALS CPU time = ",tn,"seconds")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" TWO-ELECTRON REPULSION INTEGRALS CPU time = ",te,"seconds")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" SCF CPU time = ",clock() - t1,"seconds")

    logger.info(printformatted)

    printformatted = "%s %f %s" % (" TOTAL CPU time = ",clock() - t0,"seconds")

    logger.info(printformatted)

    global runstopped
    runstopped = True

    # CUBE COMPUTATION
    if checkprint[5]:

        from cubeCall import cubeCall
        
        print("\n\n\n")

        printformatted = "%s %s" % (plot_type.upper(),"STARTED...")

        logger.info(printformatted)

        cubeCall(R,C,P,A,D,Center,CartAng,lenL,locL,Z,fnCube,N,n,nos_initial,nos_final,plot_type,libCubePath,logger)

        printformatted = "%s %s" % (plot_type.upper(),"DONE")

        logger.info(printformatted)


    return e, count, countlist, EN, ENT, ENlist, P, C, False

        
