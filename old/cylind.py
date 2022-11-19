import numpy as np
from matplotlib import pyplot as plt

# delclaration of parameters
B=0.7
C1=0.01
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


# Creating vectors  X and Y axis, making mesh
s1=np.linspace(0,1,50)
s2=np.linspace(0,1,50)
#xv, yv = np.meshgrid(s1, s2, sparse=False, indexing='ij')


#f1,col1=np.shape(n)
#f2,col2=np.shape(m)

x1=0
y1=0
for h in range(1,M+1):
	for i in range(1,N+1):	
         mod=((h**2)+(i**2))**(-B/2.0)
         #x1=x1+(m.item(h-1,i-1)*mod*np.cos(2*np.pi*(h*xv+i*yv)+n.item(h-1,i-1)))
         #y2=y1+(m.item(h-1,i-1)*mod*np.cos(2*np.pi*(h*xv+i*yv)+n.item(h-1,i-1)))
         x1=x1+(m.item(h-1,i-1)*mod*np.cos(2*np.pi*(h*s1+i*s2)+n.item(h-1,i-1)))
         y1=y1+(m.item(h-1,i-1)*mod*np.cos(2*np.pi*(h*s1+i*s2)+n.item(h-1,i-1)))

x=np.cos(2*np.pi*s1)*(1+C1*x1)
y=np.sin(2*np.pi*s1)*(1+C1*y1)
z=s2*2*np.pi



#import pdb
#pdb.set_trace()

#print x1, y1

plt.plot(x,y, linewidth=2.0)
#plt.axis('equal')
plt.show()

#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = Axes3D(fig)
#surf=ax.plot_surface(xv, yv, C1*z, rstride=1, linewidth=0, cstride=1, cmap='jet')
#ax.set_aspect('equal')
#plt.show()