#coding:utf-8
import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy as np 
from posprocess import *
np.set_printoptions(precision=10)
dados = pd.read_csv('x_vely.csv')
dados2 = pd.read_csv('y_velx.csv')
dados3 = pd.read_csv('velx_vely.csv')
dados_exatas = pd.read_csv('velx_vely_exatas.csv')
dados_perfil_exata = pd.read_csv('vuperfil_exatas.csv')
#print(dados.head())
x = np.float64(dados['position'])
vely = np.float64(dados['vely'])
velx = np.float64(dados2['velx'])
vely_ex = np.float64(dados_perfil_exata['vely'])
velx_ex = np.float64(dados_perfil_exata['velx'])

veli = dados3['veli']
velj = dados3['velj']

Nx = 13# len(veli)
Ny = 13# len(velj)
dx = np.float64(1.0/(Nx-2.0))
dy = np.float64(1.0/(Ny-2.0))
xpos, ypos = positions(Nx, Ny, dx,dy)
matx = dados2d(Nx,Ny,xpos)
maty = dados2d(Nx,Ny,ypos)
matvx = dados2d(Nx,Ny,veli)
matvy = dados2d(Nx,Ny,velj)

def plotar(arquivo,eixo,x,vely,vely_ex):
	fig, ax = plt.subplots()
	#plt.subplot(2, 1, 1)
	if eixo==0:
		xlabel = 'x [m]'
		ylabel = 'Vertical Velocity V [m/s]'
		title = 'Perfil de velocidade em y = 1/2 [m]'

		ax.plot(x, vely, marker='D', color='purple', linewidth=1)
		ax.plot(x, vely_ex, marker='.', color='blue', linewidth=1)



	elif eixo ==1:
		xlabel = 'Horizontal Velocity U [m/s]'
		ylabel = 'y[m]'
		title = 'Perfil de velocidade em x = 1/2 [m]'

		ax.plot(velx ,x, marker='D', color='purple', linewidth=1)
		ax.plot(velx_ex,x, marker='.', color='blue', linewidth=1)

	else: 
		#if (eixo!=3):
		print('Type 0 to the x axis or 1 to the y axis') # or 3 to the iteration')
		#else:
			#xlabel = 'Iteração'
			#ylabel = ' '
			#title = 'Convergencia de u, v e da correção da pressão'

			#ax.plot(vely ,x, marker='D', color='purple', linewidth=1) #iteracao
			#ax.plot(vely_ex,x, marker='.', color='blue', linewidth=1)


	#ax.plot(x, vely, marker='D', color='purple', linewidth=1)
	#ax.plot(x, vely_ex, marker='.', color='blue', linewidth=1)
	ax.legend(['Numérico', 'Analítico'])
	ax.set(xlabel= xlabel, ylabel= ylabel,
	       title=title)
	ax.grid()
	fig.savefig(arquivo)
	#else:
	#			print('Type 0 to the x axis and 1 to the y axis')
# plt.show()
# arquivo = 'item2.png'
# eixo = 0
# plotar(arquivo,eixo,x,vely,vely_ex)

# arquivo = 'item4.png'
# eixo = 1
# plotar(arquivo,eixo,x,velx,velx_ex)

def plot_field(iteracao,Nx,Ny,matvx,matvy):

		norma_vel = np.zeros((Nx,Ny), dtype = np.float64)
		for i in range(Nx):#
			for j in range(Ny):
				norma_vel[i][j] = np.sqrt(matvx[i][j]**2.0+matvy[i][j]**2.0)
		X = np.linspace(0,1, Nx)
		Y = np.linspace(0,1, Ny)
		#ax2 = fig.add_subplot(gs[1, 0])
		fig2, ax2 = plt.subplots()

		#  Varying line width along a streamline
		#lw = 4*norma_vel/norma_vel.max() #5*speed / speed.max()
		#lw = 2
		#ax2.streamplot(X, Y, matvx, matvy, density=0.6, color='k', linewidth=lw)
		#ax2.set_title('Varying Line Width')

		#QUIVER mais bonito para campo vetorial 
		#q = ax2.quiver(Y, X, matvy/norma_vel.max() , matvx/norma_vel.max())
		q = ax2.quiver(Y, X, matvy, matvx)
		#ax2.quiverkey(q, X=0.9, Y=0.9, U=lw,
        #     label='Quiver key, length = 10', labelpos='E')
		#ax2.quiverkey(q, X=0.9, Y=0.9, U=lw,
        #     label='Quiver key, length = 10', labelpos='E')
		#ax.set_aspect('equal')
		plt.show()
		arquivo = 'campo'+str(iteracao)+'.png'
		fig2.savefig(arquivo)

# iteracao = 101
# plot_field(iteracao,Nx,Ny,matvx,matvy)


def plotar_iter(itt):
	#i = int(np.ceil(Nx/2.0))
	#j = int(np.ceil(Ny/2.0))
	#pc = i + j*Nx
	dados = pd.read_csv('data_iteracao.csv')
	uit = np.abs(dados['uit'])
	vit = np.abs(dados['vit'])
	plit = np.abs(dados['plit'])
	#itt = len(plit)
	it = np.linspace(1,itt, num=itt)
	fig, ax = plt.subplots()
	#plt.subplot(2, 1, 1)
	plt.yscale('log')
	xlabel = 'Iteração'
	ylabel = 'Módulo de u, v e p\' em escala Logarítmica'
	title = 'Convergencia de u, v e p\' '
	ax.plot(it, uit, marker='p', color='red', linewidth=1)
	ax.plot(it, vit, marker='.', color='blue', linewidth=1)
	ax.plot(it, plit, marker='x', color='green', linewidth=1)
	ax.legend(['u(1/2,1/2)','v(1/2,1/2)', 'p\'(1/2,1/2)'])
	ax.set(xlabel= xlabel, ylabel= ylabel,
	title=title)
	ax.grid()
	plt.show()
	fig.savefig('iteracao.png')
	return()

# itt = 50
# plotar_iter(itt)



