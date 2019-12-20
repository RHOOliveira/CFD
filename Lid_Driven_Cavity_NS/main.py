#coding=latin-1
from math import *
import numpy as np 
import string
#from volumes import *
from coeficientes import *
from volumes2 import *
from fontes import *
from solvers import *
#from modulos import *     # FUNCIONOU
x = np.float32(None)
x = np.matrix(None, dtype=np.float32)

#x[1][1]=np.float32(2.3)

Nx = 4
Ny = 4
Nxy = Nx*Ny
dx = np.float32(1.0/(Nx-2.0))
dy = np.float32(1.0/(Ny-2.0))
dt = np.float32(1.0)

pl = np.zeros((Nx,Ny), dtype=np.float32)
p_qml = np.zeros((Nx,Ny), dtype=np.float32)

#u_qml = np.zeros((Nx,Ny), dtype=np.float32)
#v_qml = np.zeros((Nx,Ny), dtype=np.float32)

u_e = np.ones((Nx,Ny), dtype=np.float32)
v_n = np.ones((Nx,Ny), dtype=np.float32)

u_e = np.ones((Nx*Ny,1), dtype=np.float32)
v_n = np.ones((Nx*Ny,1), dtype=np.float32)

u_qml = np.zeros((Nx*Ny,1), dtype=np.float32)
v_qml = np.zeros((Nx*Ny,1), dtype=np.float32)
uold = 0.5*np.ones((Nx*Ny,1), dtype=np.float32)
vold = np.ones((Nx*Ny,1), dtype=np.float32)




#A = np.array([], dtype=np.float32)
A = np.zeros((Nx,Ny), dtype=np.float32)
B = np.zeros((Nx,Ny), dtype=np.float32)

#for i 
#x = 1.0
#x = volume(dx,dy,u_e[2][2],v_n[2][2])#.(dx,dy,1.0,2.0)

#print((x.mdot_e))
print((6-Nx)%Nx)
A = coefs(dx,dy,u_e,v_n)
B = coefs(dx,dy,u_e,v_n)

#print(A)
#print(A.coef_reais_qml_u_v(Nx,Ny,dx,dy,u_e,v_n))
#print(B.coef_fic_qml_u_v(Nx,Ny))
A=(A.coef_reais_qml_u_v(Nx,Ny,dt,dx,dy,u_e,v_n))
B=(B.coef_fic_qml_u_v(Nx,Ny))
A = A + B
#A = np.matrix(A) + np.matrix(B) 
print(A)

#print(A)
#print(B)
#cr = A.coef_reais_qml_u_v(Nx,Ny,dx,dy,u_e,v_n)
#print(cr)

#escreve = open("coeficientes.txt", "w")
#escreve.write(str(np.matrix(A)))
#escreve.close()

bp = np.zeros((Nx*Ny,1), dtype=np.float32)
#print(bp)
mp = volume(dx,dy,5.0,3.0).mass
#print(mp)
bp = fonte(dx,dy,u_e,v_n).b_qml_u(Nx,Ny,dx,dy,dt,p_qml,uold)

#pqml = 0.0
#print(bp.b_qml_u(Nx,Ny,dx,dy,dt,p_qml,uold))

ufaces = vel(dx,dy,u_e,v_n)
#print(ufaces.u_e(Nx,Ny,dx,dy,dt,A,u_qml, uold, p_qml))

#pl[0][1]=1.0
#print(pl)
aum = np.ones((Nx*Nx,1), dtype=np.float32) 	
#print(np.multiply(A, aum))
#print(A.dot(aum))
#print(A.dot(aum)[6])
tol = float(1e-5)
it = 1

#xx = gauss_seidel0(A,bp,uold,tol,it)

#print(A)











#print(A.diagonal())
diag = [ A[i][i] for i in range(len(A)) ]
diag2 = np.zeros((Nx*Ny,Nx*Ny), dtype=np.float32)
for i in range(len(diag)):
	diag2[i][i]=diag[i] 
#print((diag2))

Up = np.zeros((Nx*Ny,Ny*Nx), dtype = np.float32)
Lp = Up
Dp = Up
for i in range(0,Nx*Ny):
	for j in range(0,Ny*Nx):
		if (j>i):
			Up[i][j] = A[i][j]
		elif (j<i):
			Lp[i][j] = A[i][j]
		else:
			Dp[i][j] = A[i][j]
#print(Dp)
#print(diag2)
#print(Dp==diag2)
#DECOMPOSIÇÃO LU
D = np.diag(np.diag(A))  
L = np.tril(A)-D  
U=np.triu(A)-D
#print(D+L+U==A)
print(L)
print(D)
#print(bp)