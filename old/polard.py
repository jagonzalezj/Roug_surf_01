import numpy as np
from matplotlib import pyplot as plt

# delclaration of parameters
C1=0.01
N=10
a=0.2
b=0.2

# Random number generation
m=0+1*np.random.randn(N+3,1)
n=-np.pi/2+(np.pi/2+np.pi/2)*np.random.rand(N,1)
## Printing histogram
#plt.subplot(121)
#plt.hist(np.hstack(m), 50, normed=1, facecolor='g', alpha=0.75)
#plt.subplot(122)
#plt.hist(np.hstack(n), 50, normed=1, facecolor='g', alpha=0.75)
#plt.show()


# Saving random numbers
###................. incomplete

# Creating vectors  X axis
s=np.linspace(0,1,100)

z=0
for i in range(1,N+1):
	#z=z+((a**n)*np.cos(2*np.pi*(b**n)*x))	     # w1

	#z=z+(np.sin(np.pi*(i**b)*x))/(np.pi*(i**b))  # w2

	z=z+((i**(-b/2.0)) * m.item(i-1) * np.cos((2*np.pi*i*s) + n.item(i-1)))   # w3


# polar coordibate
x=np.cos(2*np.pi*s)*(1+C1*z)
y=np.sin(2*np.pi*s)*(1+C1*z)

#print s

plt.plot(x,y, linewidth=2.0)
plt.axis('equal')
plt.show()