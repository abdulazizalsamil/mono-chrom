
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
from math import sqrt
#plt.style.use('seaborn-white')

import sys
  

def plots(T):

	filename='./analysis/files/U'+ str(int(T)) + '.txt'

	file=open(filename,'r')


	F = np.genfromtxt(file) 


	xpt=int(F[0][0])
	ypt=int(F[0][1])
	zpt=int(F[0][2])
	n  =F[0][3]


	dt =F[0][0]
	lx =F[0][1]
	ly =F[0][2]
	lz =F[0][3]

	Ux =F[0][0]
	Vy =F[0][1]
	Wz =F[0][2]
	Ps =F[0][3]

	Ts =F[0][0]
	rho=F[0][1]
	nu =F[0][2]
	Re =F[0][3]





	F = np.delete(F, 0, 0)
	F = np.delete(F, 0, 0)
	F = np.delete(F, 0, 0)
	F = np.delete(F, 0, 0)


	X=np.zeros((xpt,ypt,zpt))
	Y=np.zeros((xpt,ypt,zpt))
	Z=np.zeros((xpt,ypt,zpt))

	U=np.zeros((xpt,ypt,zpt))
	#V=np.zeros((xpt,ypt,zpt))
	#W=np.zeros((xpt,ypt,zpt))
	#P=np.zeros((xpt,ypt,zpt))



	fig,ax = plt.subplots()

	count=0

	for i in range(xpt):
		for j in range(ypt):
			for k in range(zpt):
				X[i][j][k]=i/xpt
				Y[i][j][k]=j/ypt
				Z[i][j][k]=-1*k/zpt

				U[i][j][k]=F[count][0]
				#V[i][j][k]=F[count][1]
				#W[i][j][k]=F[count][2]
				#P[i][j][k]=F[count][2]

				count+= 1
			
 


	x=np.zeros((xpt,zpt))
	y=np.zeros((xpt,zpt))
	z=np.zeros((xpt,zpt))

	u=np.zeros((xpt,zpt))
	#v=np.zeros((xpt,zpt))
	#w=np.zeros((xpt,zpt))

	J=int(ypt/2)
	#print(U[xpt-1][ypt-1][0])

	for i in range(xpt):
		for k in range(zpt):
			x[i][k]=X[i][J][k]
			y[i][k]=Y[i][J][k]
			z[i][k]=Z[i][J][k]
			#u[i][k]=sqrt(U[i][J][k]*U[i][J][k]+V[i][J][k]*V[i][J][k]+W[i][J][k]*W[i][J][k])
			u[i][k]=U[i][J][k]
			#v[i][k]=V[i][J][k]
			#w[i][k]=W[i][J][k]



	ax.clear()
	ax.contourf(x, z, u, cmap=cm.jet)
	plt.axis('off')
	

	#print(F.shape)


	plt.xlabel('x - axis')
	plt.ylabel('z - axis')
	timename='time = '+ str(T*n*dt/100)+' s'
	plt.title(timename)
	#plt.colorbar()
 

	plt.show()
	#return pt



T=int(sys.argv[1])

plots(T)



