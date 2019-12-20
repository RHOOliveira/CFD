#coding=latin-1
import numpy as np
import volumes2
#class coefs():
def coef_reais_qml_u_v(self,Nx,Ny,dx,dy,u_e, v_n):
		#import numpy as np
		#import volumes2
		#import volume  
		rho = np.float32(1.0)
		mi = np.float32(1.0)
		A = np.zeros((Nx*Ny,Ny*Nx), dtype=np.float32)
		#Areais = ( a for a in A if a <= 1.0) 
		for i in range(1,Nx) :
			for j in range(1,Ny):
				vol =volume(dx,dy,u_e[i-1][j],v_n[i][j-1])
				A[i][j-1]= vol.mdot_e + mi*dy/dx
				vol =volume(dx,dy,u_e[i-1][j],v_n[i][j-1])
				A[i][j+Nx]= vol.mdot_n + mi*dy/dx
				#A[i][j]=volume.m_p(dx,dy,u_e[i-1][j],v_n[i][j-1])
		B = self.A
		return(B)

def fluxos(self, dx,dy, u_e, v_n):
		import numpy as np
		Nx = 13
		Ny = 13
		rho = np.float32(1.0)
		mi = np.float32(1.0)
		dz = np.float32(1.0)

		m_p = rho*dx*dy*dz
		mdot_e = rho*u_e*dy
		mdot_n = rho*v_n*dx
		#self.mdot_e = mdot_e
		#self.mdot_n = mdot_n

		return(mdot_e, mdot_n)
