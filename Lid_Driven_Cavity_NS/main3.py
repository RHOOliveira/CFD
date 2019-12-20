#coding=latin-1
from math import *
import numpy as np 
import string
#from volumes import *
from coeficientes import *
from volumes2 import *
from fontes import *
from solvers import *
from posprocess import *
import pandas as pd
from graficos import *  
from pandas import DataFrame 

np.set_printoptions(precision=11)

Nx = 13
Ny = 13

Nxy = Nx*Ny
dx = np.float64(1.0/(Nx-2.0))
dy = np.float64(1.0/(Ny-2.0))

dt = np.float64(0.1) #0.01
rho = np.float64(1.0)
mi = np.float64(1.0)

u_e = np.zeros((Nx*Ny,1), dtype=np.float64)
v_n = np.zeros((Nx*Ny,1), dtype=np.float64)
u_qml = np.zeros((Nx*Ny,1), dtype=np.float64)
v_qml = np.zeros((Nx*Ny,1), dtype=np.float64)
uold = np.zeros((Nx*Ny,1), dtype=np.float64)
vold = np.zeros((Nx*Ny,1), dtype=np.float64)
p_qml = np.zeros((Nx*Ny,1), dtype = np.float64)
pl = np.zeros((Nx*Ny,1), dtype = np.float64)

#initial guess
exata = solexata(Nx,Ny,dx,dy)
u_qml = exata.u
v_qml = exata.v
p_qml = exata.pressao
correct0 = correction(Nx,Ny,dx,dy,p_qml,pl)
p_qml = correct0.press 
#print(u_qml)
uiter = np.array([], dtype=np.float64)
viter = np.array([], dtype=np.float64)
pliter = np.array([], dtype=np.float64)


itt = 50
for t in range(itt):
		coeficientes = coefs(Nx,Ny,dt,dx,dy,u_e,v_n)
		ar= coeficientes.reais_qml
		af= coeficientes.ficticios_qml
		A = ar + af
		#print(ar+af)

		fontes = fonte(Nx,Ny,dx,dy,dt,p_qml,uold, vold)
		bp_u = fontes.qml_u
		bp_v = fontes.qml_v
		
		
		tol = 1.0e-8
		Itv = int(5)
		u_qml = gseidel(Nx,Ny,A,bp_u,u_qml,tol,Itv)
		v_qml = gseidel(Nx,Ny,A,bp_v,v_qml,tol,Itv)

		#print(u_qml)

		velocidade = vel(Nx,Ny,dx,dy,dt,A,u_qml,v_qml,u_e,v_n,p_qml, bp_v) 
		v_n = velocidade.facenorte
		u_e = velocidade.faceleste

		for ciclomassa in range(1):
			comp_press = eq_pressao( Nx,Ny,rho,dt, dx,dy, u_e,v_n,pl)
			coefp = comp_press.coeficientes
			fontp = comp_press.fontes
			#print(fontp)
		

			tol = 1e-8
			itmax = int(10)
			pl = gseidel(Nx,Ny,coefp,fontp,pl,tol,itmax)
			
			correct = correction(Nx,Ny,dx,dy,p_qml,pl)
			p_qml = correct.press 


			up, vp = acoplamento( Nx,Ny,rho,dt, dx,dy, u_qml, v_qml, pl, bp_u, bp_v)

			ue, vn = acoplamento_faces( Nx,Ny,rho,dt, dx,dy, u_e, v_n, pl)

			u_e = ue
			v_n = vn

		uold  = up
		vold = vp
		u_qml = up
		v_qml = vp

		uiter = dados_itera(Nx,Ny,up,uiter)
		viter = dados_itera(Nx,Ny,vp,viter)
		pliter = dados_itera(Nx,Ny,pl,pliter)
		#para vizualizar a evolucao
		#if True : #(t%10==0):
		#	up2 = contornos(Nx,Ny,up)
		#	vp2 = contornos(Nx,Ny,vp)
		#	matrizu = dados2d(Nx,Ny,up2)
		#	matrizv = dados2d(Nx,Ny,vp2)
		#	plot_field(t,Nx,Ny,matrizu,matrizv)

#para pegar o campo melhor convergido		
#t = itt
#up2 = contornos(Nx,Ny,up)
#vp2 = contornos(Nx,Ny,vp)
#matrizu = dados2d(Nx,Ny,up2)
#matrizv = dados2d(Nx,Ny,vp2)
#plot_field(t,Nx,Ny,matrizu,matrizv)


force, forceex, errof = force_tampa(Nx,Ny,dx,dy,mi,rho,u_qml)
flux_m, flux_mex, erroflux = mdot0(Nx, dx, rho)

varsec = { 'fluxo': [flux_m],
'fluxo analitico': [flux_mex],
'erro' : [erroflux]}

varsec2 = { 'forca': [force],
'forca analitica': [forceex],
'erro' : [errof]}

print(varsec)
# print(varsec2)

df = pd.DataFrame(varsec, columns= ['fluxo', 'fluxo analitico', 'erro'])
export_csv = df.to_csv (r'fluxo.csv', index = None, header=True)
df = pd.DataFrame(varsec2, columns= ['forca', 'forca analitica', 'erro'])
export_csv = df.to_csv (r'forca.csv', index = None, header=True)



#========================================================

#GERAÇÃO DE ARQUIVOS DE DADOS PARA O POS PROCESSAMENTO 

#print(u_qml)
#print(v_qml)
#print(p_qml)

#ARRUMAR CONTORNOS PARA PLOTAR
uplot = contornos(Nx,Ny,u_qml)
vplot = contornos(Nx,Ny,v_qml)

xpos, ypos = positions(Nx, Ny, dx,dy)
xp = np.zeros((Nx), dtype=np.float64)



#matx = dados2d(Nx,Ny,xpos)
#maty = dados2d(Nx,Ny,ypos)

for i in range(0,Nx):
	p = i + (ceil(Nx/2.0)-1)*Nx
	xp[i]= xpos[p]
#xp = [xpos[p] for p in aux]
#xp.shape = (1,Nx)

#VELOCIDADES NO PERFIL
eixo = 1 #vertical
uc = central_node(eixo,Nx,uplot)
eixo = 0 # horizontal
vc = central_node(eixo,Nx,vplot)
dados0 ={'position': xp,
			'vely': vc,
			'velx': uc
		}
#	'velx': uc
#print(dados)
#item2
#todas as velocidades
uplot.shape = (Nx*Ny)
vplot.shape = (Ny*Nx)
dados1 = {'veli': uplot,
			'velj': vplot
		} 

#print(uplot)
df2 = pd.DataFrame(dados0, columns= ['position', 'vely'])
df3 = pd.DataFrame(dados0, columns= ['position', 'velx'])

export_csv = df2.to_csv (r'x_vely.csv', index = None, header=True)
export_csv = df3.to_csv (r'y_velx.csv', index = None, header=True)

df4 = pd.DataFrame(dados1, columns= ['veli', 'velj'])
export_csv = df4.to_csv (r'velx_vely.csv', index = None, header=True)

#solucoes exatas
uex = exata.u #np.array(exata.u, dtype = np.float64)
vex = exata.v #np.array(exata.v, dtype = np.float64)

uex.shape = (Nx*Ny)
vex.shape = (Ny*Nx)
dados_exata= {'uexata': uex,
				'vexata': vex
				}
df4 = pd.DataFrame(dados_exata, columns=['uexata', 'vexata'])
export_csv = df4.to_csv (r'velx_vely_exatas.csv', index = None, header=True)

#err_u, err_v = 
num_errors(Nx,Ny,uex,vex, uplot, vplot)


#SOLUCOES EXATAS NO PERFIL
eixo = 1 
uc_ex = central_node(eixo,Nx,exata.u)
eixo = 0
vc_ex = central_node(eixo,Nx,exata.v)
dados0 ={'vely': vc_ex,
			'velx': uc_ex
		}

df4 = pd.DataFrame(dados0, columns=['vely', 'velx'])
export_csv = df4.to_csv (r'vuperfil_exatas.csv', index = None, header=True)

dados = {'uit': uiter,
			'vit': viter,
			'plit': pliter
		 } 
df = pd.DataFrame(dados, columns=['uit','vit','plit'])
export_to_csv = df.to_csv(r'data_iteracao.csv', index = None, header=True)
# df = pd.DataFrame(dados1, columns=['erro'])
# export_to_csv = df.to_csv(r'errov.csv', index = None, header=True)
#print (df2)
arquivo = 'item2.png'
eixo = 0
plotar(arquivo,eixo,x,vely,vely_ex)

arquivo = 'item4.png'
eixo = 1
plotar(arquivo,eixo,x,velx,velx_ex)

iteracao = itt
plot_field(iteracao,Nx,Ny,matvx,matvy)

#itt = 50
plotar_iter(itt)
