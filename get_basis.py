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

from numpy import empty,size,ones,zeros,asarray,size,sum,sqrt,diff
import sys
import os
# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    exedir= os.path.dirname(sys.executable)
elif __file__:
    exedir = os.path.dirname(__file__)

#config_path = os.path.join(application_path, config_name)
#exedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,exedir)
sys.path.insert(0,exedir+'/Basis')
sys.path.insert(0,'Basis')

def basislist():
    #basis = ['sto-3g','sto-6g','3-21g','6-31g','6-31g**','6-31g(d,p)','6-31g**++','6-311g**','6-311g++(2d,2p)','6-311g++(3d,3p)','6-311g++(3df,3pd)','dzvp','cc-pvdz','cc-pvtz','cc-pvqz','cc-pv5z','cc-pv6z','aug-cc-pvdz','aug-cc-pvtz','aug-cc-pvqz','aug-cc-pv5z','aug-cc-pv6z']
    basis = ['sto-3g','sto-6g','3-21g','6-31g','6-31g**','6-31g(d,p)','6-31g**++','6-311g**','dzvp','cc-pvdz','cc-pvtz']
    return basis

def Z_tablelist():
    Z_table = {'H' : 1, 'HE' : 2, 'Li' : 3, 'Be' : 4, 'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8, 'F' : 9, 'NE' : 10, 'NA' : 11, 'MG' : 12, 'AL' : 13, 'SI' : 14, 'P' : 15, 'S' : 16, 'CL' : 17, 'AR' : 18, 'K' : 19, 'CA' : 20, 'SC' : 21, 'TI' : 22, 'V' : 23, 'CR' : 24, 'MN' : 25, 'FE' : 26, 'CO' : 27, 'NI' : 28, 'CU' : 29, 'ZN' : 30, 'GA' : 31, 'GE' : 32, 'AS' : 33, 'SE' : 34, 'BR' : 35, 'KR' : 36, 'RB' : 37, 'SR' : 38, 'Y' : 39, 'ZR' : 40, 'NB' : 41, 'MO' : 42, 'TC' : 43, 'RU' : 44, 'RH' : 45, 'PD' : 46, 'Ag' : 47, 'CD' : 48, 'IN' : 49, 'SN' : 50, 'SB' : 51, 'TE' : 52, 'I' : 53, 'XE' : 54}
    return Z_table

def mass_tablelist():
    mass = [
    0.00,
    1.0008, 4.0026,
    6.941,9.0122,
    10.811,12.011,14.007,15.999,18.998,20.179,
    22.990,24.305,
    26.982,28.086,30.974,32.066,35.453,39.948,
    39.098, 40.078,
    44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845,
    58.9332, 58.6934, 63.546,65.39,
    69.723, 72.61, 74.9216, 78.96, 79.904, 83.80,
    85.4678, 87.62,
    88.90686, 91.224, 92.90638, 95.94, 98, 101.07,
    102.90550, 106.42, 107.8682, 112.411,
    114.818, 118.710, 121.760, 127.60, 126.90447, 131.29,
    132.90545, 137.327, 138.9055, 140.11, 140.90765, 144.24,
    145.0, 150.36, 151.964,
    157.25, 158.92534, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967,
    178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59]
    return mass

def inp(basisset,units,inputName):

    Z_table = Z_tablelist()

    sym2powerlist = {'S' : [0,0,0],'P' : ([1,0,0],[0,1,0],[0,0,1]),'D' : [[2,0,0],[1,1,0],[1,0,1],[0,2,0],[0,1,1],[0,0,2]],'F' : [[3,0,0],[2,1,0],[2,0,1],[1,2,0],[1,1,1],[1,0,2],[0,3,0],[0,2,1],[0,1,2], [0,0,3]]}

    basis_map = {
        '6-31g':'p631',
        '6-31g**':'p631ss',
        '6-31g(d,p)':'p631ss',
        '6-31g**++':'p631ppss',
        '6-31g++**':'p631ppss',
        '6-311g**':'p6311ss',
        '6-311g++(2d,2p)':'p6311pp_2d_2p',
        '6-311g++(3d,3p)':'p6311pp_3d_3p',
        '6-311g++(3df,3pd)':'p6311pp_3df_3pd',
        '3-21g':'p321',
        'sto3g':'sto3g',
        'sto-3g':'sto3g',
        'sto-6g':'sto6g',
        'lacvp':'lacvp',
        
        'ccpvdz':'ccpvdz',
        'cc-pvdz':'ccpvdz',
        'ccpvtz':'ccpvtz',
        'cc-pvtz':'ccpvtz',
        'ccpvqz':'ccpvqz',
        'cc-pvqz':'ccpvqz',
        'ccpv5z':'ccpv5z',
        'cc-pv5z':'ccpv5z',
        'ccpv6z':'ccpv6z',
        'cc-pv6z':'ccpv6z',

        'augccpvdz':'augccpvdz',
        'aug-cc-pvdz':'augccpvdz',
        'augccpvtz':'augccpvtz',
        'aug-cc-pvtz':'augccpvtz',
        'augccpvqz':'augccpvqz',
        'aug-cc-pvqz':'augccpvqz',
        'augccpv5z':'augccpv5z',
        'aug-cc-pv5z':'augccpv5z',
        'augccpv6z':'augccpv6z',
        'aug-cc-pv6z':'augccpv6z',
        
        'dzvp':'dzvp',
        }

    atoms_line = []

    basis = basis_map[basisset]
    basis_data = __import__(basis).basis_data


    inputread = open(inputName,'r')
    inputlines = inputread.readlines()

    #method = []
    l = []
    atoms = []
    count_lines = 0
    for line in  inputlines:
        lines = line.strip()
        if lines:
            atoms.append(lines.split()[0])
            l.append(line.split()[1:4])

            count_lines += 1

    inputread.close()

    R = l

    atoms = list(map(str.upper,atoms))

    Z = empty([len(atoms)],int)
    CartAng = []
    Center = []
    lenAtomN = [0]
    #A = empty([N,L],float)
    #D = empty([N,L],float)
    A = []
    D = []

    N_count = 0
    N2_count = 0
    N3_count = 0
    for atom in atoms:
        Z[N3_count] = Z_table[str(atom)]
        N_count = 0
        for orb in range(len(basis_data[Z_table[atom]])):
            
            #CartAng[N_count:N_count+size(sym2powerlist[str(basis_data[Z_table[atom]][orb][0])])/3] = sym2powerlist[str(basis_data[Z_table[atom]][orb][0])]
            if sum(sym2powerlist[str(basis_data[Z_table[atom]][orb][0])]) == 0:
                CartAng.append(sym2powerlist[str(basis_data[Z_table[atom]][orb][0])])
            else:
                CartAng.extend(sym2powerlist[str(basis_data[Z_table[atom]][orb][0])])
            for ang in range(int(size(sym2powerlist[str(basis_data[Z_table[atom]][orb][0])])/3)):
                orbexp = []
                contrcoef = []
                Center.append(R[N3_count])
                for coef in range(len(basis_data[Z_table[atom]][orb][1])):
                    a,d  = basis_data[Z_table[atom]][orb][1][coef]
                    orbexp.append(a)
                    contrcoef.append(d)
                    #A[N2_count,coef], D[N2_count,coef] = basis_data[Z_table[atom]][orb][1][coef]
                A.append(orbexp)
                D.append(contrcoef)
                N2_count += 1
            N_count += 1
        lenAtomN.append(len(CartAng))
        N3_count += 1


    N = len(CartAng)
    n_electron = int(sum(Z)) #- charge
    print(Z)
    #A = asarray(A,dtype=float)
    #D = asarray(D,dtype=float)
    CartAng = asarray(CartAng,dtype=int)
    Center = asarray(Center,dtype=float)
    n = len(Z)
    R = asarray(l,dtype=float).reshape(n,3)
    lenAtomN = diff(asarray(lenAtomN,dtype=int))
    lenL = []
    locL = [0]
    [lenL.append(len(ap)) for ap in A]
    tmp = 0
    for a in A[0:-1]:
        tmp += len(a)
        locL.append(tmp)

    if units == 'Angstroms':
        R = R*1.889726
        Center = Center*1.889726        
    else:
        pass

    mass = []
    mass_table = mass_tablelist()
    for m in Z:
        mass.append(mass_table[m])

    '''  
    print(locL)
    print(lenL)
    print(D)
    print(Z)
    print(N)
    print(R)
    print(n_electron)
    print(Center)
    print(CartAng)
    '''
    return N,R,Z,A,D,CartAng,Center,n_electron,n,lenL,locL,lenAtomN,atoms,mass

