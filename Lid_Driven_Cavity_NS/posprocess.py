#coding=latin-1
import numpy as np 
import pandas as pd

#pos processamento
def contornos(Nx,Ny,u):
		u1 = np.zeros((Nx*Ny), dtype = np.float64)
		#u1 = np.array((u), dtype = np.float64)
		u1 = u 
		#contorno norte 
		j = Ny-1
		for i in range(1,Nx-1):
			p = i + j*Nx
			ps = p - Nx			           
			u1[p] = 0.5*(u[ps]+u[p])
		#contorno sul
		j = 0
		for i in range(1,Nx-1):
			p = i + j*Nx
			pn = p + Nx			            
			u1[p] = 0.5*(u[pn]+u[p])
		#contorno leste
		i= Nx-1
		for j in range(1,Nx-1):
			p = i + j*Nx
			pw = p -1			            
			u1[p] = 0.5*(u[pw]+u[p])

		#contorno oeste
		i= 0
		for j in range(1,Nx-1):
			p = i + j*Nx
			pe = p +1
			u1[p] = 0.5*(u[pe]+u[p])


		return(u1)

def positions(Nx, Ny, dx,dy):
	#dx = np.float64(1.0/(Nx-2))
	x = np.zeros((Nx*Ny,1), dtype = np.float64)
	y = np.zeros((Nx*Ny,1), dtype = np.float64)
	for i in range( 1, Nx-1):
		for j in range(1, Ny-1):
			p = i + j*Nx
			x[p] = np.float64(dx*(i-0.5))
			y[p] = np.float64(dx*(j-0.5)) 

	j = Ny-1
	for i in range(1,Nx-1):
			p = i + j*Nx
					           
			y[p] = np.float64(1.0)
			x[p] = np.float64(dx*(i-0.5)) 
	#contorno sul
	j = 0
	for i in range(1,Nx-1):
			p = i + j*Nx
			pn = p + Nx			            
			y[p] = np.float64(0.0)
			x[p] = np.float64(dx*(i-0.5)) 
	#contorno leste
	i= Nx-1
	for j in range(1,Nx-1):
			p = i + j*Nx
			pw = p -1			            
			y[p] = np.float64(dx*(j-0.5))
			x[p] = np.float64(1.0) 

	#contorno oeste
	i= 0
	for j in range(1,Nx-1):
			p = i + j*Nx
			pe = p +1
			y[p] = np.float64(dx*(j-0.5))
			x[p] = np.float64(0.0) 


	return( x, y)

def central_node(eixo,Nx, u):
	ucentral = np.zeros((Nx), dtype=np.float64)
	#ucentral.shape = (1,Nx)
	jc = int(np.ceil(np.float64(Nx/2.0))-1)
	#jc = np.ceil(np.float64(Nx/2.0))
	if eixo == 0:
		for i in range(0,Nx): 
			p = int(i + jc*Nx)
			ucentral[i] = u[p]
	elif eixo == 1:
		for i in range(0,Nx): 
			p = int(jc + i*Nx)
			ucentral[i] = u[p]


	return( ucentral)

def positions2d(Nx,Ny,x):
	matx = np.zeros((Nx,Ny), dtype = np.float64)
	for i in range(Nx):
		for j in range(Ny):
			p = i + j*Nx
			matx[i][j] = x[p]

	return(matx)


def dados2d(Nx,Ny,u):
	matriz = np.zeros((Nx,Ny), dtype = np.float64)
	for i in range(Nx):
		for j in range(Ny):
			p = int(i + j*Nx)
			matriz[i][j] = u[p]



	return(matriz)

def num_errors(Nx,Ny,uexata,vexata, unum, vnum):
	errosu = np.zeros(Nx*Ny, dtype = np.float64)
	errosv = np.zeros(Nx*Ny, dtype = np.float64)
	eruc = np.zeros(Nx, dtype = np.float64)
	ervc = np.zeros(Nx, dtype = np.float64)
	errosu = uexata - unum
	errosv = vexata - vnum
#	dados0 = {'erro': errosu}  #
#	dados1 = {'erro': errosv}

	#ERROS NO CORTE DO PERFIL
	j = np.ceil(Nx/2.0)-1
	for i in range(Nx):
		#for j in range(Ny):
		p = int(i + j*Nx)
		p0 = int(j + i*Nx)
		eruc[i] = errosu[p0]
		ervc[i] = errosv[p]

	dados0 = {'erro': eruc}  
	dados1 = {'erro': ervc}



	df = pd.DataFrame(dados0, columns=['erro'])
	export_to_csv = df.to_csv(r'errou.csv', index = None, header=True)
	df = pd.DataFrame(dados1, columns=['erro'])
	export_to_csv = df.to_csv(r'errov.csv', index = None, header=True)

	#return( errosu, errosv)


def dados_itera(Nx,Ny,u,uiter):
	i = int(np.ceil(Nx/2.0)-1)
	j = int(np.ceil(Ny/2.0)-1)
	pc = i + j*Nx
	#uiter = np.array([], dtype=np.float64)
	uiter = np.append(uiter, [u[pc]])
	return(uiter)


def mdot0(Nx, dx, rho):
	np.set_printoptions(precision=11)
	dados = pd.read_csv('velx_vely.csv')
	u = dados['veli']
	v = dados['velj']
	jc = int(np.ceil(Nx/2.0)-1) 
	vaux = [v[i+jc*Nx] for i in range(0,int(np.ceil(Nx/2.0)-1))]
	vaux = np.sum(vaux, dtype = np.float64)+ 0.5*v[int(np.ceil(Nx/2.0)-1)]

	mdot = rho*dx*vaux
	

# exata 
	mdotexata = np.float64(3.0/32.0)
	erro_mdot = mdotexata - mdot 


	return(mdot, mdotexata, erro_mdot)



def force_tampa(Nx,Ny,dx,dy,mi,rho,uqml):
	np.set_printoptions(precision=11)
	#dados = pd.read_csv('velx_vely.csv')
	#u = dados['veli']
	j = Ny-1
	F = np.float64(0.0)
	du = np.float64(0.0)

	for i in range(1,Nx-1):
		p = i +j*Nx
		ps = p - Nx 
		du += uqml[p]-uqml[ps]
	F = np.float64(rho*dx*du/dy)

	Fexata = np.float64(8.0/3.0*mi)
	Erro = Fexata - F 
	
	return(F, Fexata, Erro)

	



















