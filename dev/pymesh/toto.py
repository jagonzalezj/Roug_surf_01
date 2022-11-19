# Based on : https://pymesh.readthedocs.io/en/latest/api_procedural_mesh_generation.html#cylinder-generation

import pymesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d  import Axes3D 
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

#def filtering():
	
#########                                               
#         RANDOM NUMBER GENERATION (NORMAL / UNIFORM)
#########                                               
def random_normal(M,N):
	m=0+1*np.random.randn(M,N)
	np.savetxt('m.txt', m, fmt='%10.8f')
	return m

def random_uniform(M,N):
	n=-np.pi/2+(np.pi/2+np.pi/2)*np.random.rand(M,N)
	np.savetxt('n.txt', n, fmt='%10.8f')
  	return n
#########---------------------------------------------------------

#########                       
#          PRINCIPAL EQUATION
#########                       
def rdnsurf(m, n, B, xv, yv, sfrM, sfrN):
	print '****> Creatting random surface....'
	z=0.
	for h in range(1,len(sfrM)+1):                        # h=1:length(sfrM)      
	    for i in range(1,len(sfrN)+1):                    # i=1:length(sfrN)
	        if sfrM.item(h-1)==0 or sfrN.item(i-1)==0:    # or(sfrM(h),sfrN(i)) == 0
	            continue
	        else:           
	            mod=(sfrM.item(h-1)**2+sfrN.item(i-1)**2)**(-B/2.0)
	            z=z+(m.item(h-1,i-1)*mod*np.cos(2*np.pi*(sfrM.item(h-1)*xv+sfrN.item(i-1)*yv)+n.item(h-1,i-1)));
	return z
#########---------------------------------------------------------


def plot_3d(x,y,z):        # will be deleted
	fig =plt.figure()
	ax = Axes3D(fig)
	ax.scatter(x,y,z,marker = 'o')
	plt.show()
	return

def grafica(z,xv,yv):      # will be deleted
	fig = plt.figure()
	ax = Axes3D(fig)
	surf=ax.plot_surface(xv, yv,z, rstride=1, linewidth=0, cstride=1, cmap='jet' , antialiased=False) #color ='b')##  #cmap='jet')
	ax.set_aspect('equal')
	ax.set_xlabel("X (nm)")
	ax.set_ylabel("Y (nm)")
	ax.set_zlabel("Height (nm)")
	ax.grid(False)
	#plt.axis('off')
	plt.savefig('figure.jpg', dpi = 300) 
	plt.show()
	plt.close(fig)
	return

def rho(x,y):
	return np.sqrt(x**2+y**2)

def tetha(x,y):
#	return np.arctan(y/x)	
	return np.arctan2(y,x)	

#################################################################	

if __name__ == '__main__':

	# 1 - Make a cylinder
	cyl = pymesh.generate_cylinder([0,0,0], [0,0,5], 1, 1, num_segments=36)

	# 2 - Refine: Make z-slices
	cyl = pymesh.tetrahedralize(cyl, 2, radius_edge_ratio=2, facet_distance=-1.0, feature_angle=120, engine='auto', with_timing=False)

	#  Nodes
	vertices = cyl.vertices

	#plot_3d(vertices[:,0],vertices[:,1],vertices[:,2])

	# 3 - APPLY UNROOL
	#------------------------------------------------------------------------
	nodenumber = range(0,len(vertices))           # generate the nodenumbers column
	vertices = np.insert(vertices,3,nodenumber,1)   # insert nodenumbers column in vertices matrix
	
	# FILTERING (take only nodes in the surface) we assume a perfect centered cylinder so the points 
	# at the surface are those where norm(x,y) is maximum (and moreover has the same value) !!!!!!
	#indexes = np.where( rho(vertices[:,0],vertices[:,1]) == np.max(rho(vertices[:,0],vertices[:,1])))
	#print np.array(indexes)
	stay =[]
	for x in range(0,len(vertices)):
		if rho(vertices[x,0],vertices[x,1]) > 0.9:
			stay.append(x)

	no_need = np.delete(nodenumber,stay)          # delete from nodenumbers the ones in the surface
	nodesurf = np.delete(vertices,no_need,0)      # delete from vertices the ones not at surface
	
	# 3.1 CONVERT TO CYLINDRICS (Rho Theta Z nodenumb) :)
 	cy_nodesurf = np.array([rho(nodesurf[:,0],nodesurf[:,1]), tetha(nodesurf[:,1],nodesurf[:,0]), nodesurf[:,2] ,nodesurf[:,3] ]).T

 	# 3.2 SORT BY Z THEN BY TETHA
 	data  = cy_nodesurf[np.lexsort((cy_nodesurf[:,1],cy_nodesurf[:,2]))][::-1]

 	# 3.3 GENERATE CONECTION MATRIX ------------
    
    # find the uniques elements in z
    # how many time the uniques are repited (regular mesh ==> all uniques are repeated
    # the same number of times ===> so.. runs for only one unique  :)  )
	z_levels = np.unique(cy_nodesurf[:,2])  # how many levels 
	
	fil = len(z_levels)   # number of levels  (Y' = Z )
	col = list(cy_nodesurf[:,2]).count(z_levels[1]) # how many times is repeated (X' = points in perimeter)
	matrix = data[:,3].reshape(fil, col)

	# 4 APPLY EQUATION OF ROUGHNESS
	#------------------------------------------------------------------------

	B=2.8    # parameters
	C1=0.05
	N=30
	M=N
	
	sfrN=np.linspace(-N,N,2*N+1)  # creating vectors for M and N (elements ins the double sum)
	sfrM=np.linspace(-M,M,2*M+1)
	m=0+1*np.random.randn(len(sfrM),len(sfrN));   # normal/gaussian for the amplitude
	n=-np.pi/2+(np.pi/2+np.pi/2)*np.random.rand(len(sfrM),len(sfrN)); # uniform for the phase

	y=np.linspace(0,1,col)        # creating X and Y vector  (Conection matrix size)
	x=np.linspace(0,1,fil)

	xv,yv = np.meshgrid(x, y, sparse=False, indexing='ij')  # making the mesh-grid
	z=rdnsurf( m, n, B, xv, yv, sfrM, sfrN)	   ### MAIN EQUATION 
	z=C1*z       
	print np.shape(z)
#	grafica(z,xv,yv)
	
	# 4.1 displace the dz to avoid negative numbers and mesh nodes conflicts
	min_dz = np.min(z)
	z = z + min_dz 

	# 5 APPLY RE-ROLL
	#------------------------------------------------------------------------

	#clasicall way
	for i in xrange(0,fil):
	 	for j in xrange(0,col):
	 	 	delta_z  = z[i,j]    # take element one by one
	 	 	index = matrix[i,j]  # find the corresponding index in node matrix
	 	 	poss = np.where(vertices[:,3]==index)  # find index location in initial data 
	 	 	phi_p = tetha(vertices[poss,0],vertices[poss,1]) #find the corresponding tetha angle   
	 	 	# subsitute X and Y by X+dX and Y+dY
	 	 	dx = delta_z*np.cos(phi_p)   #conversion to dx !!!!! rad or degree ?????
	 	 	dy = delta_z*np.cos(phi_p)
	 	 	vertices[poss,0] = vertices[poss,0] + dx #subsitute X by x+dx
	 	 	vertices[poss,1] = vertices[poss,1] + dy #subsitute X by x+dx


	plot_3d(vertices[:,0],vertices[:,1],vertices[:,2])
	
	

	#print matrix
 	#np.savetxt('out.txt',data, fmt='%4.f')
 	#copy paste the rest of functions








# apply dz
# apply reroll
# replace (x,y,z) of current node by new x,y,z






