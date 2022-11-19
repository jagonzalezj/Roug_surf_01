# coding=utf-8
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from shapely.geometry import Point, Polygon, MultiPoint, LineString


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
	#print '****> Creatting random surface....'
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



#########
#          COOREDENATES TRASFORMATION (CYLINDER / SPHERE)
#########
def cilindro(z,y,r):
	fil,col=np.shape(z)
	t = np.linspace(0,2*np.pi,col)
	tc, vc = np.meshgrid(t,y, sparse=False, indexing='ij')
	xv = (r+z)*np.cos(tc)
	yv = (r+z)*np.sin(tc)
	ce= y + np.zeros((fil,col))   # adapting one column z data to the grid
	z=ce
	return z, xv, yv

def esfera(z,r):
	fil,col=np.shape(z)
	t=np.linspace(0,2*np.pi,col)
	v=np.linspace(0,np.pi,fil)
	tc,vc=np.meshgrid(t,v,)# sparse=False, indexing='ij')
	xc=(r+z)*np.cos(tc)*np.sin(vc)
	yc=(r+z)*np.sin(tc)*np.sin(vc)
	xv=xc
	yv=yc
	z=np.cos(vc)*(r)
	return z, xv, yv
#########---------------------------------------------------------



#########
#          GRAPHICS  Surf
#########
def grafica(z,xv,yv):
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
#########---------------------------------------------------------

#########
#          Write final lmp file with atoms removed
#########
def lampsave(lafila,posicion):
	data = np.loadtxt("removed.out")
	A=len(data)

	print '****> Creatting file ',lafila[0:len(lafila)-4]+'_out.lmp'
	fw=open(lafila[0:len(lafila)-4]+'_out.lmp','w')

	fp = open(lafila,'r')
	for i, line in enumerate(fp):
		if  i==2:
			fw.write('%12d %s\n' % ( A,' atoms'))
		elif i==posicion+2:   # p=11 for old .lmp format 15 for new .lmp format
			break
		else:
			fw.write(line)
	fp.close()

	for x in xrange(0,A):
		fw.write('%8d %4d %16.8f %16.8f %16.8f\n' % (data[x,0], data[x,1], data[x,2], data[x,3], data[x,4]))
	fw.close()
	return
#########---------------------------------------------------------




#########
#         Find the 2dn minimum value per column [latt parameter]  {remember 0}
#########
def secmin(vec): 
	xx=0
	while vec[xx]==0:
		xx=xx+1
	else:
		minimos=vec[xx]
	return minimos
#########---------------------------------------------------------



#########
#          SAVE SURFACE DATA
#########
def save_surf(xv,yv,z):
	np.savetxt('xv.out', xv, fmt='%10.8f')
	np.savetxt('yv.out', yv, fmt='%10.8f')
	np.savetxt('z.out', z, fmt='%10.8f')
	return
#########---------------------------------------------------------


#########
#          Remove_Plane
#########
def rem_plane():
	data = np.loadtxt('matrix.out')    # remove this part and correct the error of format 
	data=data[np.lexsort((data[:,4][::-1],))] ##  TWO TIMES SAME SCRIP FOR ORIENTED SAMPLES, ONLY ONE DO NOT ORGANICE WELL
	data=data[np.lexsort((data[:,4][::-1],))] ##  THE LAST COLUMN INDUCE FAILURE IN 'WHILE' LOOP,(1 OR 2 LAST LAYERS NOT PROCESSED) 
	xv=np.loadtxt('xv.out')   ######## ???????
	yv=np.loadtxt('yv.out')
	z=np.loadtxt('z.out')
	datarem=data
	toremove=[]
	x=0
	az=data[x,4]
	while az>=z.min():
 		ax=data[x,2]
 		ay=data[x,3]
 		az=data[x,4]
 		indx=np.where((xv[:,0])==(ax))
 		indy=np.where((yv[0,:])==(ay))
		if az>=z[indx[0],indy[0]]:
			toremove.append(x)
		x+=1
	datarem=np.delete(datarem,toremove,0)

####### matrix preparation...  DONE!!
	datarem=datarem[np.lexsort((datarem[:, 0][::1], ))]
	datarem=np.delete(datarem,0,1)
	addcol=range(1,len(datarem)+1)
	datarem=np.insert(datarem,0,addcol,1)
	np.savetxt('removed.out', datarem)
	return datarem
#########---------------------------------------------------------



#########
#          Remove_cylinder
#########
def rem_cylinder(data,xv,yv,z):
	data = np.loadtxt('matrix.out')    # remove this part and correct the error of format 	
	data=data[np.lexsort((data[:, 4][::-1], ))]  ### organized data
	az=data[:,4]                                 ### last column (Z)

	singles = np.sort(list(set(az)))                ### repeated numbers
	singles=singles[np.lexsort((singles[::-1], ))]  ### organized inversed

	final=[]   # matrix to keep the points inside 

	for i in xrange(0,len(singles)):
		h=[]
		g=[]
		level=singles[i]                          # levels (z plane)    
		indx=np.where(az==level)                  # porsitio of z planes, how many times 
	
		x=xv[:,i]            # data perimeter cylinder
		y=yv[:,i]
		for dp in xrange(0,len(y)):
			c=x[dp],y[dp]
			h.append(c)	

		ax=data[np.min(indx):np.max(indx)+1,2]        # list of points
 		ay=data[np.min(indx):np.max(indx)+1,3]
		for dp in xrange(0,len(ay)):
			d=ax[dp],ay[dp]
			g.append(d)

		polygon=Polygon(h)
		matprob=data[np.min(indx):np.max(indx)+1,:]

		for f in xrange(0,len(g)):
			b=Point(g[f]).within(polygon)  ### poner codigo de remover data aqui!!!!		
			if b==True:
				final.append(matprob[f])	

	np.savetxt('final_tempo.out', final,fmt='%10.8f')   

	data = np.loadtxt("final_tempo.out")
	data=data[np.lexsort((data[:,0][::1], ))]    #####
	data=np.delete(data,0,1)
	addcol=range(1,len(data)+1)
	data=np.insert(data,0,addcol,1)

	np.savetxt('removed.out', data,fmt='%10.8f')   # change final by removed

	return data
#########---------------------------------------------------------





# MORE

def lampload(lafila):
	print '****> Reading file :', lafila,' ....'
	flag = False
	arr = []
	p=0
	with open(lafila) as f:
		for line in f:			
			if(flag == False and 'Atoms' in line):
				posicion=p # 'position' of word 'Atoms' for diferenciation between old and new lmp format exported to be used in 'lampsave'
				flag = True
				continue
			if(flag == True and len(line.strip())):
				arr.append(list(map(lambda x: float(x), line.split())))
			p+=1
	np.savetxt('matrix.out', arr,fmt='%16.8f')   
	return arr, posicion


def puntos():
	data = np.loadtxt("matrix.out")
	x_pos=data[:,[2]]
	y_pos=data[:,[3]] 
	z_pos=data[:,[4]]
	n_pntx=len(np.unique(x_pos))   # not used (for now!)
	n_pnty=len(np.unique(y_pos))
	n_pntz=len(np.unique(z_pos))
	xmax=x_pos.max()
	ymax=y_pos.max()
	zmax=z_pos.max()
	x_grid=np.unique(x_pos)/xmax
	y_grid=np.unique(y_pos)/ymax
	z_grid=np.unique(z_pos)/zmax
	return x_grid, y_grid, z_grid, xmax, ymax, zmax, n_pntx, n_pnty, n_pntz


def stat_analysis(z, N, C1, B):
	z_an=np.reshape(z,-1)
	print ''
	print '------------ Random Surface Parameters-----------------'
	print '         N=M = ',N,'  C1 = ',C1, '  b = ',B
	print 'No. points = ', len(z_an)
	print 'Mean_Value = ', np.mean(z_an)
	print ' Stand_dev = ', np.std(z_an)
	print '       RMS = ', np.sqrt(np.sum(np.square(z_an))/len(z_an))
	print '  Skewness = ', np.sum(np.power((z_an-np.mean(z_an)),3)/len(z_an))/np.power(np.std(z_an),3)
	print '  Kurtosis = ', np.sum(np.power((z_an-np.mean(z_an)),4)/len(z_an))/np.power(np.std(z_an),4)
	print '--------------------------------------------------------'
	return



# SPHERE UTILES

#colormap = cm.coolwarm
# fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_trisurf(ico.X, ico.Y, ico.Z, triangles=ico.faces, linewidth=0.4, alpha=0.2, edgecolor='r')
        # ax.scatter(ico.X, ico.Y, ico.Z, c='r', marker='o', s=0.2)
        # plt.show()



# fig = plt.figure()                                                      # out
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_trisurf(X, Y, Z, triangles=ico.faces, linewidth=0.3, alpha=0.1, edgecolor='r')

        # color = colormap((r - r.min())/r.ptp())
        # ax.scatter(X, Y, Z, c=color[...,:3], marker='o', s=1.1)
        # plt.show()

def stat_sphere(r,C1, C2):
    print("ave: {}".format(np.mean(r)))                                    # out
    print("RMS: {}".format(np.std(C1*C2*r)))
    print("RMS: {}".format(np.sqrt(C1*C2*np.sum(r*r)/len(r))))
    return
# def sphere(r,C1, C2, ico.R, ico.PHI, ico.THETA):
#     X = (C1*C2*r+ico.R) * np.sin(ico.PHI) * np.cos(ico.THETA)               # out
#     Y = (C1*C2*r+ico.R) * np.sin(ico.PHI) * np.sin(ico.THETA)
#     Z = (C1*C2*r+ico.R) * np.cos(ico.PHI)
#     return
