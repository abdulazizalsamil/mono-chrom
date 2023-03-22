import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from datetime import datetime
#from plots2d import plots
#import plots2d as pl
from math import sqrt

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
	V=np.zeros((xpt,ypt,zpt))
	W=np.zeros((xpt,ypt,zpt))
	P=np.zeros((xpt,ypt,zpt))



	fig,ax = plt.subplots()

	count=0

	for i in range(xpt):
		for j in range(ypt):
			for k in range(zpt):
				X[i][j][k]=i/xpt
				Y[i][j][k]=j/ypt
				Z[i][j][k]=k/zpt

				U[i][j][k]=F[count][0]
				V[i][j][k]=F[count][1]
				W[i][j][k]=F[count][2]
				P[i][j][k]=F[count][2]

				count+= 1
			
 


	x=np.zeros((xpt,zpt))
	y=np.zeros((xpt,zpt))
	z=np.zeros((xpt,zpt))

	u=np.zeros((xpt,zpt))
	v=np.zeros((xpt,zpt))
	w=np.zeros((xpt,zpt))

	J=int(ypt/2)
	#print(U[xpt-1][ypt-1][0])

	for i in range(xpt):
		for k in range(zpt):
			x[i][k]=X[i][J][k]
			y[i][k]=Y[i][J][k]
			z[i][k]=Z[i][J][k]
			u[i][k]=sqrt(U[i][J][k]*U[i][J][k]+V[i][J][k]*V[i][J][k]+W[i][J][k]*W[i][J][k])
			u[i][k]=U[i][J][k]
			v[i][k]=V[i][J][k]
			w[i][k]=W[i][J][k]


	z1=np.zeros(zpt)
	w1=np.zeros(zpt)
	
	Ix=int(xpt/2)
	for k in range(zpt):
		z1[k]=Z[Ix][J][k]
		w1[k]=W[Ix][J][k]




	plt.plot(w1,z1)


	plt.xlabel('w')
	plt.ylabel('z')
	timename='time = '+ str(T*n*dt/100)+' s'
	plt.title(timename)
	#plt.colorbar()
 

	plt.show()
	#return pt




fig,ax = plt.subplots()

now = datetime.now()
current_time = now.strftime("%Y%m%d%H%M%S")



plnam='./analysis/plots/p'+current_time+'.gif'

interval = 1#in seconds     
ani = animation.FuncAnimation(fig,plots,100,interval=interval*1e+3,blit=False)
ani.save(plnam, writer='PillowWriter', fps=3)
plt.show()


