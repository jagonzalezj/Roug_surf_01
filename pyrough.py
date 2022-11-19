# coding=utf-8
import ran_utile as ru
import numpy as np
import sys, os, subprocess 

import scipy.special as sp
from matplotlib import cm
from Icosahedron import Icosahedron
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from stl import mesh


#print '#######################################################'
#print '#           Running Random Surface script             #'
#print '#     Jonathan Amodeo & Javier Gonzalez....2018       #'  
#print '#######################################################'
 
#  LOADING DEFAULT PARAMETER FILE
par=np.loadtxt('param.in')
 
#comands = ['RS', 'atoms', 'slab', 'wire', 'sphere', 'material', 'dimension']

# 3.1 ON FLY USER DEFINE RS PARAMETERS
if 'RS' in sys.argv:
#    print "==> Generating user defined random surface"
    # put parameters from parameter_file.txt
else:
#    print "==> loading default parameters for random surface"

    if 'b' in sys.argv:	
    	B = float(sys.argv[(sys.argv.index('b'))+1])
    else:
    	B =par[0,0]

    if 'C' in sys.argv:	
    	C1 = float(sys.argv[(sys.argv.index('C'))+1])
    else:
    	C1 =par[0,1]

    if 'N' in sys.argv:	
    	N = int(sys.argv[(sys.argv.index('N'))+1])
    else:
    	N =par[0,2]

    if 'M' in sys.argv:	
    	M = int(sys.argv[(sys.argv.index('M'))+1])
    else:
    	M=par[0,3]



# USER CREATION OF ATOMIC SAMPLE (ATOMSK CALLING)  
# MATERIAL
if 'material' in sys.argv: #Material receive only 3 arguments (FCC, 3.165, [100])
    lat_struct=sys.argv[(sys.argv.index('material'))+1]
    lat_par=sys.argv[(sys.argv.index('material'))+2]
    lat_orient=sys.argv[(sys.argv.index('intercept'))+3]

# DIMENSION  all is considering 100
if 'dimension' in sys.argv: #Material receive only 3 arguments (see param.txt)
    d1=sys.argv[(sys.argv.index('dimension'))+1] # lenght    radii      radii
    d2=sys.argv[(sys.argv.index('dimension'))+2] # wide      height      NULL
    d3=sys.argv[(sys.argv.index('dimension'))+3] # height     NULL       NULL
    if 'slab' in sys.argv:

        dup1= d1/lat_par
        dup2= d2/lat_par
        dup3= d3/lat_par
    elif 'wire' in sys.argv:
        pass
    elif 'sphere' in sys.argv:
        pass

    subprocess.call(['atomsk', '--create', ])





# 3.  INTERCEPTIG OPTION (ASSUMES A .LMP FILE)
if 'intercept' in sys.argv:	
	the_file = sys.argv[(sys.argv.index('intercept'))+1]




# ATOMS or FEM
if 'atoms' in sys.argv[:]:

    # load .lmp file
    data, posicion = ru.lampload(the_file) 

    # critic points in .lmp file
    # grid is non-redundant proj of atom, pnt is the nb of element in grid
    x_grid, y_grid, z_grid, xmax, ymax, zmax, n_pntx, n_pnty, n_pntz = ru.puntos()

    # RANDOM NUMBER GENERATION                                    
    sfrN=np.linspace(-N,N,2*N+1)  # if N=M then step = 2N+1
    sfrM=np.linspace(-M,M,2*M+1)
    m=0+1*np.random.randn(len(sfrM),len(sfrN));   # normal/gaussian for the amplitude
    n=-np.pi/2+(np.pi/2+np.pi/2)*np.random.rand(len(sfrM),len(sfrN)); # uniform for the phase

    if 'slab' in sys.argv:

    	xv,yv = np.meshgrid(x_grid, y_grid, sparse=False, indexing='ij')
    	z=ru.rdnsurf( m, n, B, xv, yv, sfrM, sfrN)
    	ru.stat_analysis(z, N, C1, B)
    	z=C1*z  

        xv=xv*xmax      # rescale to the original sample size
        yv=yv*ymax
        z=z+(zmax-z.max())         # random surface adapted o the maximun height to intercept
        
        # to avoid the arrayinarray problem, javier export the data to a file !
        ru.save_surf(xv,yv,z)      # saving surface data [do not comment!!! (#)]

        datarem = ru.rem_plane()  # intercepting 
        ru.lampsave(the_file, posicion)  # generating new .lmp file  may be put out of loop
        
#        print ' Slab done !'

    elif 'wire' in sys.argv:

    	xv,yv = np.meshgrid(x_grid, z_grid, sparse=False, indexing='ij') # make X-Y axis distributed as .lmp sample #   y_grid change to z_grid in case mistmach
    	z=ru.rdnsurf( m, n, B, xv, yv, sfrM, sfrN)
    	ru.stat_analysis(z, N, C1, B)  
    	z=C1*z
		
        # coordenate trasnformation to cylindric
        # z is the axis of the zire
        r=xmax/4 #-2*xmin 
        z,xv,yv=ru.cilindro(z,z_grid,r)  # change y by z_grid
        
        # if np.sqrt(xv**2+yv**2).max() > xmax/2. :
        #     print 'ERROR: size of radii of cylinder greather than atom sample'
        #     quit()
		
        # translation to center of .lmp file
        xv=xv+xmax/2
        yv=yv+ymax/2
        z=z*zmax  

        datarem=ru.rem_cylinder(data,xv,yv,z) # intercepting
        ru.lampsave(the_file,posicion)  #may be put out of loop 

        print '  wire done !'

    elif 'sphere' in sys.argv:
        radii=xmax/2.5        
        ico = Icosahedron(refine=4, radius=radii)       # Sphere creation + meshin
        nbPoint = len(ico.vertex)

        ### JA: these will not have to be hardcoded
        C1=10000*6
        # JA : name of the paper with details on shape of sph. harmonics
        N_s = 9         # starting degree of the spherical harmonics
        N_e = 15        # ending degree of the spherical harmonics
        l = 3    # degree
        m = 0 

        r = np.zeros(nbPoint)

        for degree in range(N_s,N_e+1,1):           # EQUATION
            print("degree: {}".format(degree))
            _r_amplitude = 0+1*np.random.randn(nbPoint)  
            _r_phase = -np.pi/2.+np.pi*np.random.rand(nbPoint);    
            mod = degree**(-B/2)
            for i, [theta, phi] in enumerate(ico.vertex_tp):
                _phase = sp.sph_harm(0, degree, ico.THETA-theta, ico.PHI-phi).real 
                _phase = 2 * _phase / _phase.ptp()
                r += _r_amplitude[i] * mod * np.cos(_phase + _r_phase)
        
        C2 = 1./nbPoint/N/2.

        ru.stat_sphere(r,C1, C2)

        X = (C1*C2*r+ico.R) * np.sin(ico.PHI) * np.cos(ico.THETA)               
        Y = (C1*C2*r+ico.R) * np.sin(ico.PHI) * np.sin(ico.THETA)
        Z = (C1*C2*r+ico.R) * np.cos(ico.PHI)

        new_vertex = np.array([X, Y, Z]).T
        cube = mesh.Mesh(np.zeros(len(ico.faces), dtype=mesh.Mesh.dtype))

        for i, [a,b,c] in enumerate(ico.faces):
            for j in range(3):
                cube.vectors[i] = [new_vertex[a,:], new_vertex[b,:], new_vertex[c,:]]
       
        cube.save('cube.stl')  # Write the mesh to file "cube.stl"

        subprocess.call(['atomsk','Cu1.lmp',          # intercepting
            '-select','stl', 'center', 'cube.stl',
             '-select', 'invert', '-rmatom', 'select', 
             'sphere.lmp'])
        # atomsk Al_supercell.xsf -select stl center Dog_lowpoly_flowalistik.STL -select invert -rmatom select Dog.cfg
        print '  sphere done !'






# clean .out files  comment if requiered to save the .out
print '****> CLEANING....'
mypath=os.getcwd()
for file in os.listdir(mypath):
    if file.endswith(".out"):
        os.remove(file)

print 'JOB DONE!'

#  GRPHIC THE GEOMETRY
if 'graph' in sys.argv:
    ru.grafica(z,xv,yv)
