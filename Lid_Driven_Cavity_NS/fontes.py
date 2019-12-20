import numpy as np
from volumes2 import *
from velocities import *

class fonte(volume):

	def __init__(self, Nx,Ny,dx,dy,dt,p_qml,uold, vold):
		bp_u = np.zeros((Nx*Ny,1), dtype=np.float64)

		for i in range(1,Nx-1):
			for j in range(1,Ny-1):
				p = i + j*Nx
				pe = p + 1
				pw = p -1
				mp = volume(dx,dy,0.0,0.0)
				bp_u[p] = -0.5*(p_qml[pe]-p_qml[pw])*dy+mp.mass/dt*uold[p]
				

		# ghost cells
		j = Ny
		for i in range(1,Nx-1):
			xp = (i - 0.5)*dx
			p = i + (j-1)*Nx
			aux = 16.0*(xp**4.0-2.0*xp**3.0+xp**2.0)
			bp_u[p] = aux #aux

		self.qml_u = bp_u
		#return(bp_u)

	#def b_qml_v(self, Nx,Ny,dx,dy,dt,p_qml,vold):
		bp_v = np.zeros((Nx*Ny,1), dtype=np.float64)
		dz = np.float64(1.0)
		Re = np.float64(1.0)
		mp = volume(dx,dy,0.0,0.0).mass
		for i in range(1,Nx-1):
			for j in range(1,Ny-1):
				p = i + j*Nx
				pn =  p + Nx
				ps =  p - Nx
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

				aux = 8.0/Re*(24.0*F0+2.0*d1f*d2g+d3f*g)+64.0*(F2*G1-g*d1g*F1)

				bp_v[p] = -aux*dx*dy-0.5*(p_qml[pn]-p_qml[ps])*dx+mp/dt*vold[p]
		self.qml_v=bp_v
		#return(bp_v)





