#coding=latin-1
class volume():
	def __init__ (self,dx,dy, u_e, v_n):
		import numpy as np

		rho = np.float64(1.0)
		mi = np.float64(1.0)
		dz = np.float64(1.0)

		m_p = rho*dx*dy*dz
		mdot_e = rho*u_e*dy
		mdot_n = rho*v_n*dx
		self.mdot_e = mdot_e
		self.mdot_n = mdot_n
		self.mass = m_p 

