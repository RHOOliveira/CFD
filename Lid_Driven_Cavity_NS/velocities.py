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
		Re = np.float64(1.0)

		#uold aqui é a velocidade na face em t-dt
		for i in range(1,Nx-2):
			for j in range(1,Ny-1):  #(1,Ny-1)
				p = i + (j)*Nx
				pe = p + 1
				pw = p -1
				pn = p + Nx
				ps = p - Nx
				pee = pe + 1
				pse = ps + 1
				pne = pn + 1
				#p2 = p+1
				#pe2 = p2 + 1
				#pw2 = p2 -1
				#pn2 = p2 + Nx
				#ps2 = p2 - Nx

				ux = np.array([uqml[pw],uqml[pe],uqml[pn],uqml[ps]])
				#ux2 = np.array([uqml[pw2],uqml[pe2],uqml[pn2],uqml[ps2]])
				ux2 = np.array([uqml[p],uqml[pee],uqml[pne],uqml[pse]])
				ux.shape=(4,1)
				ux2.shape = (4,1)
				#u_e[p] = (A[p][1: ].dot(np.transpose(ux))+A[p+1][1: ].dot(np.transpose(ux2)))+(mp_p+mp_e)/dt*uold[p]-2.0*(pqml[pe]-pqml[p])*dy
				u_e[p] = (A[p][1: ].dot(ux)+A[p+1][1: ].dot(ux2))+(mp_p+mp_e)/dt*uold[p]-2.0*(pqml[pe]-pqml[p])*dy
				aux = A[p][0]+A[p+1][0]
				u_e[p] = u_e[p]/aux

		self.faceleste = u_e
		#return(u_e)


		v_n = np.zeros((Nx*Ny,1), dtype=np.float64)

		Bs = np.zeros((Nx*Ny,1), dtype = np.float64)
		#fon shi 
		for i in range(1,Nx-1):
			for j in range(1,Ny-1):
				p = i + j*Nx
				#pn =  p + Nx
				#ps =  p - Nx
				x = np.float64((i - 0.5)*dx) 
				y = np.float64((j - 0.5)*dy)

				f = np.float64(x**4.0-2.0*x**3.0+x**2.0)
				F0 = np.float64(0.2*x**5.0-0.5*x**4.0+x**3/3.0)
				d1f = np.float64(4.0*x**3.0-6.0*x**2.0+2.0*x)
				d2f = np.float64(12.0*x**2.0-12.0*x+2.0)
				d3f = np.float64(24.0*x-12.0)
				F1 = np.float64(f*d2f - d1f**2.0)
				F2 = np.float64(0.5*f**2.0)
				g = np.float64(y**4.0 - y**2.0)
				d1g = np.float64(4.0*y**3.0-2.0*y)
				d2g = np.float64(12.0*y**2.0-2.0)
				d3g = np.float64(24.0*y)
				G1 = np.float64(g*d3g-d1g*d2g) 

				Bs[p] = -(8.0/Re*(24.0*F0+2.0*d1f*d2g+d3f*g)+64.0*(F2*G1-g*d1g*F1))




		for i in range(1,Nx-1):
			for j in range(1,Ny-2):
				p = i + j*Nx
				pn = p + Nx
				pe = p + 1
				pw = p -1
				ps = p - Nx
				pnw = pn - 1
				pne = pn + 1
				pnn = pn + Nx
				#p2 = p+1
				#pe2 = p2 + 1
				#pw2 = p2 -1
				#pn2 = p2 + Nx
				#ps2 = p2 - Nx

				#ux = np.array([vqml[pw],vqml[pe],vqml[pn],vqml[ps]])
				#ux2 = np.array([vqml[pw2],vqml[pe2],vqml[pn2],vqml[ps2]])
				ux = np.array([vqml[pw],vqml[pe],vqml[pn],vqml[ps]])
				ux2 = np.array([vqml[pnw],vqml[pne],vqml[pnn],vqml[p]])
				ux.shape=(4,1)
				ux2.shape = (4,1)

				v_n[p] = (A[p][1: ].dot(ux)+A[pn][1: ].dot(ux2))+(mp_p+mp_e)/dt*vold[p]-2.0*(pqml[pn]-pqml[p])*dx
				v_n[p] = (v_n[p] -(Bs[p]+Bs[pn])*dx*dy)/(A[p][0]+A[pn][0]) 

				#analisar sinal de Bs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		self.facenorte = v_n
		#return(v_n)



