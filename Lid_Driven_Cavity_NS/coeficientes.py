#coding=latin-1
import numpy as np
from volumes2 import *
from velocities import *
class coefs(volume):
	def __init__(self,Nx,Ny,dt,dx,dy,u_e, v_n):
  
		rho = np.float64(1.0)
		mi = np.float64(1.0)
		
		a = np.zeros((Nx*Ny,5), dtype=np.float64) # 0 p 1 w 2 e 3 n 4 s
		alphae = np.float64(None)
		alphaw = np.float64(None)
		alphan = np.float64(None)
		alphas = np.float64(None)


		for i in range(1,(Nx-1)) :
			for j in range(1,(Ny-1)):
				p = int(i + (j)*Nx)
				pn = int(p + Nx)
				pe = int(p + 1)
				pw = int(p - 1)
				ps = int(p - Nx)
				alphae = 0.5*np.sign(u_e[p])
				alphaw = 0.5*np.sign(u_e[pw])
				alphan = 0.5*np.sign(v_n[p])
				alphas = 0.5*np.sign(v_n[ps])

				massa = volume(dx,dy,0.0,0.0).mass

				vol =volume(dx,dy,u_e[pw],v_n[p])
				a[p][1]= vol.mdot_e*(0.5+alphaw) + mi*dy/dx

				vol =volume(dx,dy,u_e[p],v_n[p])
				a[p][3]= -vol.mdot_n*(0.5-alphan) + mi*dx/dy

				vol =volume(dx,dy,u_e[p],v_n[p])
				a[p][2]= -vol.mdot_e*(0.5-alphae) + mi*dy/dx

				vol =volume(dx,dy,u_e[p],v_n[ps])
				a[p][4]= vol.mdot_n*(0.5+alphas) + mi*dx/dy

				a[p][0] =(a[p][1]+a[p][2]+a[p][3]+a[p][4]+massa/dt)  #<===========================

		self.reais_qml = a

		a = np.zeros((Nx*Ny,5), dtype=np.float64)
		#contorno sul
		j = int(0)
		for i in range(1,(Nx-1)) :
				p = int(i + (j)*Nx)
				pn = int(p + Nx)
				pe = int(p + 1)
				pw = int(p - 1)
				ps = int(p - Nx)

				a[p][3]= -1.0
				a[p][0] =1.0
		#contorno norte
		j = int(Ny-1)
		for i in range(1,(Nx-1)) :
				p = int(i + (j)*Nx)
				pn = int(p + Nx)
				pe = int(p + 1)
				pw = int(p - 1)
				ps = int(p - Nx)

				a[p][4]= -1.0
				a[p][0] =  1.0 
		#contorno leste	
		i = int(Nx-1)
		for j in range(1,(Ny-1)) :
				p = int(i + (j)*Nx)
				pn = int(p + Nx)
				pe = int(p + 1)
				pw = int(p - 1)
				ps = int(p - Nx)

				a[p][1]= -1.0
				a[p][0] =  1.0 
		#oontorno oeste
		i = int(0)
		for j in range(1,(Ny-1)) :
				p = int(i + (j)*Nx)
				pn = int(p + Nx)
				pe = int(p + 1)
				pw = int(p - 1)
				ps = int(p - Nx)

				a[p][2]= -1.0
				a[p][0] =  1.0

		self.ficticios_qml = a


class simplec():
	def __init__(self, Nx,Ny,rho,dt, dx,dy):

		dup = np.zeros((Nx*Ny,1), dtype=np.float64)
		dvp = np.zeros((Nx*Ny,1), dtype=np.float64)
		de = np.zeros((Nx*Ny,1), dtype=np.float64)
		dn = np.zeros((Nx*Ny,1), dtype=np.float64)

		for i in range(1,(Nx-1)) :
			for j in range(1,(Ny-1)):
				p = i + j*Nx
				dup[p] = dt/(rho*dx)
				dvp[p] = dt/(rho*dy)

		for i in range(1,(Nx-1)) :
			for j in range(1,(Ny-1)):
				p = i + j*Nx
				pn = p + Nx
				de[p] = 0.5*(dup[p]+dup[p+1])
				dn[p] = 0.5*(dvp[p]+dvp[pn])	

		self.du = dup
		self.dv = dvp
		self.de = de
		self.dn = dn

class eq_pressao(simplec):
	def __init__ (self, Nx,Ny,rho,dt, dx,dy, u_e,v_n, plin):
		a = np.zeros((Nx*Ny,5), dtype = np.float64)
		bpp = np.zeros((Nx*Ny,1), dtype = np.float64)
		for i in range(1,Nx-1):
			for j in range(1, Ny-1):
				p = i + j*Nx
				pw = p - 1
				pn = p + Nx
				ps = p - Nx
				pe = p + 1
				a[p][1] = simplec( Nx,Ny,rho,dt, dx,dy).de[pw]*dy
				a[p][2] = simplec( Nx,Ny,rho,dt, dx,dy).de[p]*dy
				a[p][3] = simplec( Nx,Ny,rho,dt, dx,dy).dn[p]*dx
				a[p][4] = simplec( Nx,Ny,rho,dt, dx,dy).dn[ps]*dx
				a[p][0] = (a[p][1]+a[p][2]+a[p][3]+a[p][4])  

				# ghosts cell
		#print(App)
		#contorno leste
		i = Ny-1
		for j in range(1,Nx-1):
			p = i + j*Nx
			pw = p -1
			pww = pw - 1
			a[p][1] = 1.0
			a[p][0] = 1.0                #<==================
			bpp[p] = plin[pw]-plin[pww]

		#contorno oeste
		i = 0
		for j in range(1,Nx-1):
			p = i + j*Nx
			pe = p +1
			pee = pe + 1
			a[p][2] = 1.0
			a[p][0] =  1.0 
			bpp[p] = plin[pe]-plin[pee]

		#contorno norte
		j = Ny-1
		for i in range(1,Nx-1):
			p = i + j*Nx
			ps = p - Nx
			pss = ps - Nx
			a[p][4] = 1.0
			bpp[p] = plin[ps]-plin[pss]
			a[p][0] = 1.0
		#contorno sul
		j = 0
		for i in range(1,Nx-1):
			p = i + j*Nx
			pn = p + Nx
			pnn = pn + Nx
			a[p][3] = 1.0
			bpp[p] = plin[pn]-plin[pnn]
			a[p][0] = 1.0

		#App = -App



		for i in range(1,Nx-1):
			for j in range(1, Ny-1):
				p = i + j*Nx
				p2 = p - 2
				p2w = p - 1
				p2s = p - Nx
				bpp[p] = (u_e[pw]-u_e[p])*dy+(v_n[ps]-v_n[p])*dx  # <---------- verificar indices

		self.fontes = bpp
		self.coeficientes = a 





				 



















			





