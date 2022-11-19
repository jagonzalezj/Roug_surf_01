# Based on : https://pymesh.readthedocs.io/en/latest/api_procedural_mesh_generation.html#cylinder-generation

import pymesh
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
import numpy as np


#mesh = pymesh.load_mesh("cube.stl")

# 1 - Make a cylinder
# .generate_cylinder(botcentercoord, topcentercoord, botrad, toprad, num_segments=nbsegcrosssec)
cyl = pymesh.generate_cylinder([0,0,0], [0,0,5], 1, 1, num_segments=36)
#print("{}".format(cyl.num_faces))
#cyl = pymesh.generate_tube([0,0,0], [0,0,10], 2, 2, 0, 0, with_quad=True)   

# Refine: Make z-slices
# pymesh.tetrahedralize(cyl, 0.1 ****), not below 0.1 else z-slices become curvy
# ?: find the relation between size, numofseg and this 0.1 variable
cyl = pymesh.tetrahedralize(cyl, 2, radius_edge_ratio=2, facet_distance=-1.0, feature_angle=120, engine='auto', with_timing=False)

# Unroll
    # Get surface nodes
#vertices = cyl.vertices

vertices = cyl.vertices
#print(vertices)
print(type(vertices))

x, y, z = vertices[:,0], vertices[:,1], vertices[:,2]
rho, theta = np.sqrt(x*x+y*y), np.arctan(y/x)

#print("rhos is {}".format(theta))


surfnodes = []
surfnodescyl = []
#verticesindexcyl = [rho, theta, z]
#print(verticesindexcyl)

#quit()

verticesindexcyl = []

for i in range(len(cyl.vertices)):
    verticesindex = np.append(vertices[i], int(i))
#    verticesindexcyl[i] = (rho[i], theta[i], z[i], int(i) ) 
    verticesindexcyl = np.append(verticesindexcyl, [rho[i], theta[i], z[i], int(i)])
    if rho[i]>0.9 :
        #print(vertices[i])
        surfnodes.append(verticesindex)
#        surfnodescyl.append(verticesindexcyl)

#quit()

surfnodes = np.asarray(surfnodes)
#surfnodes = surfnodes[surfnodes[:,2].argsort()]

#surfnodescyl = np.asarray(surfnodescyl)
surfnodescyl = surfnodescyl[surfnodescyl[:,2].argsort()]
#print(surfnodes)
#print(surfnodescyl)

#for i in range(len(surfnodes)):

#print(len(surfnodes))
#print(surfnodes[1][0])

# apply unroll



# apply dz
# apply reroll
# replace (x,y,z) of current node by new x,y,z


#print(vertices)

# Export the cylinder (binary = paraview readable)
#pymesh.save_mesh("cylinder.stl", cyl, ascii=True)
#pymesh.save_mesh("cylindeTri.stl", cylsurfmesh)



#pymesh.save_mesh("cylinder.stl", mesh, ascii=True)

# To be opened with gmsh
#pymesh.save_mesh("cylinder.msh", cyl)

