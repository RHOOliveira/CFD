#coding=latin-1
class volume():
	def __init__ (self, dx,dy, u_e, v_n):
		import numpy as np
		rho = np.float32(1.0)
		mi = np.float32(1.0)
		dz = np.float32(1.0) 
		m_p = rho*dx*dy*dz
		mdot_e = rho*u_e*dy
		mdot_n = rho*v_n*dx
		self.mdot_e = mdot_e
		self.mdot_n = mdot_n

#def coef_reais_qml_u_v(mi,dx,dy,dt):
#	import numpy 
	#chamar subrotina das mdot
def veloc_faces():
	import numpy as np
	np.zeros((2,), dtype=float32)
