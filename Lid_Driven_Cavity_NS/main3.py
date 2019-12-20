#coding=latin-1
from math import *
import numpy as np 
import string
#from volumes import *
from coeficientes import *
from volumes2 import *
from fontes import *
from solvers import *


np.set_printoptions(precision=15)
x = np.float32(None)
x = np.matrix(None, dtype=np.float64)

# alterar em volumes tbm
Nx = 13
Ny = 13

Nxy = Nx*Ny
dx = np.float64(1.0/(Nx-2.0))
dy = np.float64(1.0/(Ny-2.0))
dt = np.float64(0.1)
rho = np.float64(1.0) 
u_e = np.ones((Nx*Ny,1), dtype=np.float64)
v_n = np.ones((Nx*Ny,1), dtype=np.float64)
u_qml = np.zeros((Nx*Ny,1), dtype=np.float64)
v_qml = np.zeros((Nx*Ny,1), dtype=np.float64)
uold = np.ones((Nx*Ny,1), dtype=np.float64)
vold = np.ones((Nx*Ny,1), dtype=np.float64)
p_qml = np.zeros((Nx*Ny,1), dtype = np.float64)
pl = np.zeros((Nx*Ny,1), dtype = np.float64)
itt = 10
for t in range(itt):
		coeficientes = coefs(Nx,Ny,dt,dx,dy,u_e,v_n)
		ar= coeficientes.reais_qml
		af= coeficientes.ficticios_qml
		A = ar + af
		#print(ar+af)
		fontes = fonte(Nx,Ny,dx,dy,dt,p_qml,uold, vold)
		bp_u = fontes.qml_u
		bp_v = fontes.qml_v

		
		tol = 1.0e-6
		itmax = int(5)
		u_qml = gseidel(Nx,Ny,A,bp_u,u_qml,tol,itmax)
		v_qml = gseidel(Nx,Ny,A,bp_v,v_qml,tol,itmax)


		velocidade = vel(Nx,Ny,dx,dy,dt,A,u_qml,v_qml,uold,vold,p_qml, bp_v) 
		v_n = velocidade.facenorte
		u_e = velocidade.faceleste

		correcao_press = eq_pressao( Nx,Ny,rho,dt, dx,dy, u_e,v_n, pl)
		coefp = correcao_press.coeficientes
		fontp = correcao_press.fontes

		itmax = int(10)
		pl = gseidel(Nx,Ny,coefp,fontp,pl,tol,itmax)
		print(pl)

		p_qml = p_qml + pl  

		pqml_ficticios = extrapolacao(Nx,Ny,p_qml)

		up, vp = acoplamento( Nx,Ny,rho,dt, dx,dy, u_qml, v_qml, pl, bp_u, bp_v)

		u_e, v_n = acoplamento_faces( Nx,Ny,rho,dt, dx,dy, u_e, v_n, pl)

		uold  = up
		vold = vp
		u_qml = up
		v_qml = vp

#w = open("pressao.txt", "w")
#w.write(p_qml)
#w.close()

#w = open("u_qml.txt", "w")
#w.write(u_qml)
#w.close()

print(u_qml)
#print(p_qml)
