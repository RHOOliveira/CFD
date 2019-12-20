#coding=latin-1
import numpy as np
from volumes2 import *
class coefs(volume):
	def coef_reais_qml_u_v(self,Nx,Ny,dx,dy,u_e, v_n):
		#import numpy as np
		#import volumes2
		#import volume  
		rho = np.float32(1.0)
		mi = np.float32(1.0)
		A = np.zeros((Nx*Ny,Ny*Nx), dtype=np.float32)
		#Areais = ( a for a in A if a <= 1.0) 
		for i in range(1,(Nx)) :
			for j in range(1,(Ny)):
				p = int(i + (j-1)*Nx)
				pn = int(p + Nx)
				pe = int(p + 1)
				pw = int(p - 1)
				ps = int(p - Nx)
				if (ps >= 0 & pn >=0):
					vol =volume(dx,dy,u_e[(i-1)%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pw]= vol.mdot_e + mi*dy/dx

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pn]= vol.mdot_n + mi*dx/dy

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pe]= vol.mdot_e + mi*dy/dx

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][(j-Nx)%Nx])
					A[p][ps]= vol.mdot_n + mi*dy/dx

					A[p][p] =A[p][pw]+A[p][pn]+A[p][pe]+A[p][ps] 

				elif (ps < 0):

					vol =volume(dx,dy,u_e[(i-1)%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pw]= vol.mdot_e + mi*dy/dx

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pn]= vol.mdot_n + mi*dx/dy

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pe]= vol.mdot_e + mi*dy/dx

						#vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][(j-Nx)%Nx])
					A[p][ps]= 0.0 #vol.mdot_n + mi*dy/dx

					A[p][p] =A[p][pw]+A[p][pn]+A[p][pe]+A[p][ps] 

				else:
					vol =volume(dx,dy,u_e[(i-1)%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pw]= vol.mdot_e + mi*dy/dx

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][j%Nx])
					A[p][pe]= vol.mdot_e + mi*dy/dx

					vol =volume(dx,dy,u_e[i%Ny][j%Nx],v_n[i%Ny][(j-Nx)%Nx])
					A[p][ps]= vol.mdot_n + mi*dy/dx

					A[p][p] =A[p][pw]+A[p][pn]+A[p][pe]+A[p][ps] 

		#B = self.matriz_coef
		B=A
		return(B)
