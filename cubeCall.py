from numpy import arange, asarray,double,intc
import numpy.ctypeslib as npc
import ctypes as ct

def cubeCall(R,C,P,A,D,Center,CartAng,lenL,locL,Z,fnCube,N,n,nos_initial,nos_final,plot_type,libCubePath,logger):    

    import functools

    nptsx = 40
    nptsy = 40
    nptsz = 40
    if any(n < 0 for n in R[:,0]):
        xmin = -7. + min(R[:,0])
    else:
        xmin = -7.
    if any(n < 0 for n in R[:,1]):
        ymin = -7. + min(R[:,1])
    else:
        ymin = -7.
    if any(n < 0 for n in R[:,2]):
        zmin = -7. + min(R[:,2])
    else:
        zmin = -7.
    xmax = 7. + max(R[:,0])
    ymax = 7. + max(R[:,1])
    zmax = 7. + max(R[:,2])

    voxelx = (xmax - xmin)/(nptsx-1)
    voxely = (ymax - ymin)/(nptsy-1)
    voxelz = (zmax - zmin) /(nptsz-1)

    x_plot = arange(xmin,xmax+voxelx,voxelx)
    y_plot = arange(ymin,ymax+voxely,voxely)
    z_plot = arange(zmin,zmax+voxelz,voxelz)


    # CHOOSE: 0 - MOLECULAR ORBITAL; 1 - PROBABILITY DENSITY FUNCTION
    if plot_type == "Molecular Orbital Cube":
        plot_type_index = 0
        CPtranspinmemory = -C.T.copy()
    elif plot_type =="Electron Density Cube":
        plot_type_index = 1
        CPtranspinmemory = P.copy()


    _cube = npc.load_library(libCubePath,'.')   # . load path
    _cube.computecube.argtypes = [npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype=double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=1),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = double,ndim=2),npc.ndpointer(dtype = intc,ndim=2),npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = intc,ndim=1),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_double,ct.c_double,ct.c_double,ct.c_double,ct.c_double,ct.c_double,ct.c_int,ct.c_int,npc.ndpointer(dtype = intc,ndim=1),npc.ndpointer(dtype = double,ndim=2),ct.c_int,ct.c_char_p]
    _cube.computecube.restype = ct.c_void_p
    #c = A.astype(double)
    CartAng2 = CartAng.astype(intc)

    #inputNameCube = fnCube.split('.')[:-1][0]+'MO/'+'_'+basisset+'_'+str(charge)+'.cube'
    Aobj = asarray(functools.reduce(list.__add__,A,[]),dtype=double)
    Dobj = asarray(functools.reduce(list.__add__,D,[]),dtype=double)
    lenLarray = asarray(lenL,dtype=intc)
    locLarray = asarray(locL,dtype=intc)
    Zintc = asarray(Z,dtype=intc)
    inputNameCube = str(fnCube)
    inputNameCube_bytes = str.encode(inputNameCube, 'utf-8')
    inputNameCube_c = ct.create_string_buffer(inputNameCube_bytes)
#print(x_plot,y_plot,z_plot,Aobj,Dobj,CPtranspinmemory,Center,CartAng2,lenLarray,locLarray,N,n,nptsx,nptsy,nptsz,xmin,ymin,zmin,voxelx,voxely,voxelz,nos_initial,nos_final,Zintc,R,plot_type)

    _cube.computecube(x_plot,y_plot,z_plot,Aobj,Dobj,CPtranspinmemory,Center,CartAng2,lenLarray,locLarray,N,n,nptsx,nptsy,nptsz,xmin,ymin,zmin,voxelx,voxely,voxelz,int(nos_initial),int(nos_final),Zintc,R,plot_type_index,inputNameCube_c)


