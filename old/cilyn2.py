import numpy as np
from matplotlib import pyplot as plt

# delclaration of parameters
B=1.8
C1=0.1
N=30
M=N

# Random number generation
m=0+1*np.random.randn(M,N)
n=-np.pi/2+(np.pi/2+np.pi/2)*np.random.rand(M,N)
## Printing histogram
#plt.subplot(121)
#plt.hist(np.hstack(m), 50, normed=1, facecolor='g', alpha=0.75)
#plt.subplot(122)
#plt.hist(np.hstack(n), 50, normed=1, facecolor='g', alpha=0.75)
#plt.show()

# Saving random numbers
###................. incomplete

# Creating vectors  X and Y axis, making mesh
x=np.linspace(0,1,100)
y=np.linspace(0,1,100)
xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')


#f1,col1=np.shape(n)
#f2,col2=np.shape(m)

z=0
for h in range(1,M+1):
	for i in range(1,N+1):	
         mod=((h**2)+(i**2))**(-B/2.0)
         z=z+(m.item(h-1,i-1)*mod*np.cos(2*np.pi*(h*xv+i*yv)+n.item(h-1,i-1)))


z=z*C1
#print z
##### PLANE TO CYLINDER-------------------------------
fil,col=np.shape(z)
v=y
r=1      #radii of cilinder
t=np.linspace(0,2*np.pi,col)
tc,vc=np.meshgrid(t,v, sparse=False, indexing='ij')
xc=(r+z)*np.cos(tc)
yc=(r+z)*np.sin(tc)



from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
surf=ax.plot_surface(xc, yc,v, rstride=1, linewidth=0, cstride=1, cmap='jet')
#surf=ax.plot_surface(xv, yv, C1*z, rstride=1, linewidth=0, cstride=1, cmap='jet')

ax.set_aspect('equal')
plt.show()

