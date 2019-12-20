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

np.set_printoptions(precision=15)
x = np.float32(None)
x = np.matrix(None, dtype=np.float64)


Nx = 4
Ny = 4
Nxy = Nx*Ny
dx = np.float64(1.0/(Nx-2.0))
dy = np.float64(1.0/(Ny-2.0))
dt = np.float64(1.0)
rho = np.float64(1.0) 
pl = np.zeros((Nx,Ny), dtype=np.float64)
p_qml = np.zeros((Nx*Ny,1), dtype=np.float64)

#u_qml = np.zeros((Nx,Ny), dtype=np.float32)
#v_qml = np.zeros((Nx,Ny), dtype=np.float32)

u_e = np.ones((Nx,Ny), dtype=np.float64)
v_n = np.ones((Nx,Ny), dtype=np.float64)

u_e = np.ones((Nx*Ny,1), dtype=np.float64)
v_n = np.ones((Nx*Ny,1), dtype=np.float64)

u_qml = np.zeros((Nx*Ny,1), dtype=np.float64)
v_qml = np.zeros((Nx*Ny,1), dtype=np.float64)
uold = np.ones((Nx*Ny,1), dtype=np.float64)
vold = np.ones((Nx*Ny,1), dtype=np.float64)




#A = np.array([], dtype=np.float32)
A = np.zeros((Nx,Ny), dtype=np.float64)
B = np.zeros((Nx,Ny), dtype=np.float64)

#for i 
#x = 1.0
#x = volume(dx,dy,u_e[2][2],v_n[2][2])#.(dx,dy,1.0,2.0)

#print((x.mdot_e))
print((6-Nx)%Nx)
ap = coefs(dx,dy,u_e,v_n)
B = coefs(dx,dy,u_e,v_n)

#print(A)
#print(A.coef_reais_qml_u_v(Nx,Ny,dx,dy,u_e,v_n))
#print(B.coef_fic_qml_u_v(Nx,Ny))
ap=(ap.coef_reais_qml_u_v(Nx,Ny,dt,dx,dy,u_e,v_n))
B=(B.coef_fic_qml_u_v(Nx,Ny))
#A = A + B

A = np.zeros((Nx*Ny-4,Nx*Ny-4), dtype = np.float64)
k = Nx-2
for i in range(1,Nx-2):
    for j in range(1,Ny-2):
                    p = i + j*Nx
                    pw = p - 1
                    pe = p + 1
                    pn = p + Nx
                    ps = p - Nx
                    A[k][pw] = ap[p][1]
                    A[k][pe] = ap[p][2]
                    A[k][pn] = ap[p][3]
                    A[k][ps] = ap[p][4]
                    A[k][p] = ap[p][0]
                    k += 1
k = 0          
j = 0            
for i in range(1,Nx-1):
                    p = i + j*Nx
                    pw = p - 1
                    pe = p + 1
                    pn = p + Nx
                    ps = p - Nx
                    A[k][pw] = ap[p][1]
                    A[k][pe] = ap[p][2]
                    A[k][pn] = ap[p][3]
                    A[k][ps] = ap[p][4]
                    A[k][p] = ap[p][0]
                    k += 1

j = Ny - 1
k = j*Nx
for i in range(1,Nx-3):
                    p = i + j*Nx
                    pw = p - 1
                    pe = p + 1
                    pn = p + Nx
                    ps = p - Nx
                    A[k][pw] = ap[p][1]
                    A[k][pe] = ap[p][2]
                    A[k][pn] = ap[p][3]
                    A[k][ps] = ap[p][4]
                    A[k][p] = ap[p][0]
                    k += 1
print(A)

#A = np.delete(A, [0,Nx-1,-1,-1-Nx], axis=0)
#A = np.delete(A, [0, Nx-1, Nx*Ny], axis=0)

#extraindo os cantos 

#A = np.delete(A, 0, axis=0)
#A = np.delete(A, (Nx*(Nx)-2), axis=0)
#A = np.delete(A, (Nx*(Nx)-5), axis=0)
#A = np.delete(A, Nx-2, axis=0)
#A = np.delete(A, 0, axis=1)
#A = np.delete(A, (Nx*(Nx)-2), axis=1)
#A = np.delete(A, (Nx*(Nx)-5), axis=1)
#A = np.delete(A, Nx-2, axis=1)

#A = np.delete(A, [i for i in range(0,Nx*Nx) if A[i][i]==0 ], axis= 0)

#A = np.delete(A, [i for i in range(len(A)) if i%(Nx-1)==0 ], axis= 1)

#A = np.delete(A, [i*(Nx-1) for i in range(len(A)) ], axis= 1)
#A = np.delete(A, [i*(Nx-1) for i in range(len(A)) ], axis= 1)
print(len(A))

#A = np.delete(A, [ ], axis= 0)
#A = np.matrix(A) + np.matrix(B) 
print(A)

#print(A)
#print(B)
#cr = A.coef_reais_qml_u_v(Nx,Ny,dx,dy,u_e,v_n)
#print(cr)

#escreve = open("coeficientes.txt", "w")
#escreve.write(str(np.matrix(A)))
#escreve.close()

bp = np.zeros((Nx*Ny,1), dtype=np.float64)

mp = volume(dx,dy,5.0,3.0).mass
#print(mp)
bp = fonte(dx,dy,u_e,v_n).b_qml_u(Nx,Ny,dx,dy,dt,p_qml,uold)
bp_v = fonte(dx,dy,u_e,v_n).b_qml_v(Nx,Ny,dx,dy,dt,p_qml,vold)
print(bp_v)


#bp = np.delete(bp, 0, axis=0)
#bp = np.delete(bp, (Nx*(Nx)-2), axis=0)
#bp = np.delete(bp, (Nx*(Nx)-5), axis=0)
#bp = np.delete(bp, Nx-2, axis=0)

#pqml = 0.0
#print(bp.b_qml_u(Nx,Ny,dx,dy,dt,p_qml,uold))

ufaces = vel(dx,dy,u_e,v_n)
#print(ufaces.u_e(Nx,Ny,dx,dy,dt,A,u_qml, uold, p_qml))

#pl[0][1]=1.0
#print(pl)
aum = np.ones((Nx*Nx,1), dtype=np.float64) 	
#print(np.multiply(A, aum))
#print(A.dot(aum))
#print(A.dot(aum)[6])

tol = np.float64(1e-10)
it = 2

uold = 0.5*np.ones((Nx*Ny,1), dtype = np.float64)
vold = 0.5*np.ones((Nx*Ny,1), dtype = np.float64)
p_qml = np.zeros((Nx*Ny,1), dtype = np.float64)


A = coefs(dx,dy,u_e,v_n)
B = coefs(dx,dy,u_e,v_n)
A=(A.coef_reais_qml_u_v(Nx,Ny,dt,dx,dy,u_e,v_n))
B=(B.coef_fic_qml_u_v(Nx,Ny))
A = A + B
N = 10
#ver = gseidel(A,bp)
#uqml = np.float64(gauss_seidel0(A,bp,uold,tol,it))
#vqml = np.float64(gauss_seidel0(A,bp_v,vold,tol,it))
print(uqml)
print(vqml)

u_e = vel(dx,dy,u_e,v_n).u_e(Nx,Ny,dx,dy,dt,A,uqml, uold, p_qml) 
v_n = vel(dx,dy,u_e,v_n).u_e(Nx,Ny,dx,dy,dt,A,vqml, vold, p_qml) 
print(np.around(u_e, decimals=10))
#print(v_n)


d_u = simplec( Nx,Ny,rho,dt, dx,dy).du
print(d_u)
plin = np.ones(Nx*Ny, dtype = np.float64)
bpp = eq_pressao( Nx,Ny,rho,dt, dx,dy, u_e,v_n, plin).fontes
App = eq_pressao( Nx,Ny,rho,dt, dx,dy, u_e,v_n, plin).coef

plin = np.ones((Nx*Ny-4,1), dtype = np.float64)
#plin = np.float32(gauss_seidel0(App,bpp,plin,tol,it))

print(App[:][:])

u_e = np.zeros((Nx*Ny-4,1), dtype = np.float64)
v_e = np.zeros((Nx*Ny-4,1), dtype = np.float64)








