import numpy as np
from volumes2 import *
from coeficientes import *
from fontes import *
from math import *
#class vel(coefs):
class vel(volume):
	def __init__(self,Nx,Ny,dx,dy,dt,A,uqml,vqml,uold,vold, pqml,bp):

		u_e = np.zeros((Nx*Ny,1), dtype=np.float64)
		mp_p = volume(dx,dy,0.0,0.0).mass
		mp_e = volume(dx,dy,0.0,0.0).mass

		#diag = [ A[i][i] for i in range(len(A)) ]
		#A_p = np.zeros((Nx*Ny-4,Nx*Ny-4), dtype=np.float64)
		#for i in range(len(diag)):
	#		A_p[i][i]=diag[i] 
	#	A = A - A_p



		for i in range(1,Nx-1):
			for j in range(1,Ny-1):
				p = i + (j)*Nx
				pe = p + 1
				pw = p -1
				pn = p + Nx
				ps = p - Nx
				p2 = p+1
				pe2 = p2 + 1
				pw2 = p2 -1
				pn2 = p2 + Nx
				ps2 = p2 - Nx

				ux = np.array([uqml[pw],uqml[pe],uqml[pn],uqml[ps]])
				ux2 = np.array([uqml[pw2],uqml[pe2],uqml[pn2],uqml[ps2]])
				ux.shape=(4,1)
				ux2.shape = (4,1)
				#u_e[p] = (A[p][1: ].dot(np.transpose(ux))+A[p+1][1: ].dot(np.transpose(ux2)))+(mp_p+mp_e)/dt*uold[p]-2.0*(pqml[pe]-pqml[p])*dy
				u_e[p] = (A[p][1: ].dot(ux)+A[p+1][1: ].dot(ux2))+(mp_p+mp_e)/dt*uold[p]-2.0*(pqml[pe]-pqml[p])*dy
				aux = A[p][0]+A[p+1][0]
				u_e[p] = u_e[p]/aux

		self.faceleste = u_e
		#return(u_e)

	#def v_n(self,Nx,Ny,dx,dy,dt,A,vqml, vold, pqml):
		v_n = np.zeros((Nx*Ny,1), dtype=np.float64)
		#mp_p = volume(dx,dy,0.0,0.0).mass
		#mp_e = volume(dx,dy,0.0,0.0).mass

		#diag = [ A[i][i] for i in range(len(A)) ]
		#A_p = np.zeros((Nx*Ny-4,Nx*Ny-4), dtype=np.float64)

		for i in range(1,Nx-1):
			for j in range(1,Ny-1):
				p = i + j*Nx
				pn = p + Nx
				pe = p + 1
				pw = p -1
				ps = p - Nx
				p2 = p+1
				pe2 = p2 + 1
				pw2 = p2 -1
				pn2 = p2 + Nx
				ps2 = p2 - Nx

				ux = np.array([vqml[pw],vqml[pe],vqml[pn],vqml[ps]])
				ux2 = np.array([vqml[pw2],vqml[pe2],vqml[pn2],vqml[ps2]])
				ux.shape=(4,1)
				ux2.shape = (4,1)

				v_n[p] = (A[p][1: ].dot(ux)+A[p+1][1: ].dot(ux2))+(mp_p+mp_e)/dt*vold[p]-2.0*(pqml[pn]-pqml[p])*dy
				v_n[p] = v_n[p] + (bp[p]+bp[pn])*dx*dy 
				aux = A[p][0]+A[pn][0]
				v_n[p] = v_n[p]/aux

		self.facenorte = v_n
		#return(v_n)



