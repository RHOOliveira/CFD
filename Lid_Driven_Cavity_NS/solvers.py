from math import *
import numpy as np
from coeficientes import *

def gseidel(Nx,Ny,a,b,u0,tol,itmax):

    k = 0
    while (k < itmax):
        u = np.zeros((Nx*Ny,1), dtype = np.float64)
        for i in range(1,Nx-1):
            for j in range(1,Ny-1):
                p = i + j*Nx
                pe = p +1
                pn = p + Nx
                ps = p - Nx
                pw = p-1
                u[p] = (a[p][1]*u0[pw]+a[p][2]*u0[pe]+a[p][3]*u0[pn]+a[p][4]*u0[ps]+b[p])/a[p][0]

        #contorno sul
        j = int(0)
        for i in range(1,(Nx-1)) :
                    p = int(i + (j)*Nx)
                    pn = int(p + Nx)


                    u[p]= -u0[pn] + b[p]

        #contorno norte
        j = int(Ny-1)
        for i in range(1,(Nx-1)) :
                    p = int(i + (j)*Nx)
                    ps = int(p - Nx)

                    u[p]= -u0[ps] + b[p]

        #contorno leste 
        i = int(Nx-1)
        for j in range(1,(Ny-1)) :
                    p = int(i + (j)*Nx)
                    pw = int(p - 1)
                    #ps = int(p - Nx)

                    u[p]= -u0[pw] + b[p] 
        #oontorno oeste
        i = int(0)
        for j in range(1,(Ny-1)) :
                    p = int(i + (j)*Nx)
                    pe = int(p + 1)
 
                    u[p]= -u0[pe] + b[p]

        if (np.linalg.norm(u-u0,np.inf) < tol):
                print('tolerancia atingida com', k, 'iteracoes')
                k = itmax+1
                
                #return u
                

        else:
            k += 1
            u0 = np.copy(u)

    return(u) 

def extrapolacao(Nx,Ny,pqml):
    y = np.zeros((Nx*Ny,1), dtype = np.float64)
    #contorno norte
    j = Ny-1
    for i in range(1,Nx-1):
        p = i + j*Nx
        ps = p - Nx
        pss= ps - Nx
        #y[p] = 2.0*pqml[ps]-pqml[pss]
        y[p] = pqml[pss - Nx] -3.0*pqml[pss]+3.0*pqml[ps]
    #contorno sul
    j = 0
    for i in range(1,Nx-1):
        p = i + j*Nx
        pn = p + Nx
        pnn= pn + Nx
        #y[p] = 2.0*pqml[pn]-pqml[pnn]
        y[p] = pqml[pnn + Nx] -3.0*pqml[pnn]+3.0*pqml[pn]

    #contorno leste
    i= Nx-1
    for j in range(1,Nx-1):
        p = i + j*Nx
        pw = p -1
        pww= pw -1
        #y[p] = 2.0*pqml[pw]-pqml[pww]
        y[p] = pqml[pww -1] -3.0*pqml[pww]+3.0*pqml[pw]

    #contorno oeste
    i= 0
    for j in range(1,Nx-1):
        p = i + j*Nx
        pe = p +1
        pee= pe +1
        #y[p] = 2.0*pqml[pe]-pqml[pee]
        y[p] = pqml[pee + 1] -3.0*pqml[pee]+3.0*pqml[pe]

    return(y)

def acoplamento( Nx,Ny,rho,dt, dx,dy, uqml, vqml, pl, bp_u, bp_v):
    up = np.zeros((Nx*Ny,1), dtype = np.float64)
    vp = np.zeros((Nx*Ny,1), dtype = np.float64)
    csimplec = simplec(Nx,Ny,rho,dt, dx,dy)
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            p = i + (j)*Nx
            pe = p +  1 
            pw = p -1
            pn = p + Nx
            ps = p - Nx
            up[p] = uqml[p] - csimplec.du[p]*(pl[pe]-pl[pw])/2.0 
            vp[p] = vqml[p] - csimplec.dv[p]*(pl[pn]-pl[ps])/2.0

    #contorno sul
    j = int(0)
    for i in range(1,(Nx-1)) :
            p = int(i + (j)*Nx)
            pn = int(p + Nx)


            up[p]= -up[pn] + bp_u[p]
            vp[p]= -vp[pn] + bp_v[p]

    #contorno norte
    j = int(Ny-1)
    for i in range(1,(Nx-1)) :
            p = int(i + (j)*Nx)
            ps = int(p - Nx)

            up[p]= -up[ps] + bp_u[p]
            vp[p]= -vp[ps] + bp_v[p]

    #contorno leste 
    i = int(Nx-1)
    for j in range(1,(Ny-1)) :
            p = int(i + (j)*Nx)
            pw = int(p - 1)
            

            up[p]= -up[pw] + bp_u[p] 
            vp[p]= -vp[pw] + bp_v[p] 
    #oontorno oeste
    i = int(0)
    for j in range(1,(Ny-1)) :
            p = int(i + (j)*Nx)
            pe = int(p + 1)

            up[p]= -up[pe] + bp_u[p]
            vp[p]= -vp[pe] + bp_v[p] 







    return(up, vp)

def acoplamento_faces( Nx,Ny,rho,dt, dx,dy, u_e, v_n, pl):
    ue = np.zeros((Nx*Ny,1), dtype = np.float64)
    vn = np.zeros((Nx*Ny,1), dtype = np.float64)
    csimplec = simplec(Nx,Ny,rho,dt, dx,dy)

    for i in range(1,Nx-2):
        for j in range(1,Ny-1):   #(1,Ny-1)
            p = i + (j)*Nx
            pe = p +  1 
            pw = p -1
            pn = p + Nx
            ps = p - Nx
            ue[p] = u_e[p] - csimplec.de[p]*(pl[pe]-pl[pw])/2.0 
            #vn[p] = v_n[p] - csimplec.dn[p]*(pl[pn]-pl[ps])/2.0 

    for i in range(1,Nx-1):
        for j in range(1,Ny-2):   #(1,Ny-1)
            p = i + (j)*Nx
            pe = p +  1 
            pw = p -1
            pn = p + Nx
            ps = p - Nx
            #ue[p] = u_e[p] - csimplec.de[p]*(pl[pe]-pl[pw])/2.0 
            vn[p] = v_n[p] - csimplec.dn[p]*(pl[pn]-pl[ps])/2.0 



    return(ue, vn)

class solexata():
    def __init__(self,Nx,Ny,dx,dy):
        uexata = np.zeros((Nx*Ny,1), dtype = np.float64)
        vexata = np.zeros((Nx*Ny,1), dtype = np.float64)
        pexata = np.zeros((Nx*Ny,1), dtype = np.float64)
        Bs = np.zeros((Nx*Ny,1), dtype = np.float64)
        
        Re = np.float64(1.0)
        #mp = volume(dx,dy,0.0,0.0).mass
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

                Bs[p] = -(8.0/Re*(24.0*F0+2.0*d1f*d2g+d3f*g)+64.0*(F2*G1-g*d1g*F1))
                uexata[p] = 8.0*f*d1g
                vexata[p] = -8.0*d1f*g
                pexata[p] = (8.0/Re*(F0*d3g+d1f*d1g)+64.0*F2*(g*d2g-d1g**2.0))

        j = Ny
        for i in range(1,Nx-1):
            x = np.float64((i - 0.5)*dx)
            p = i + (j-1)*Nx
            uexata[p] = (16.0*(x**4.0-2.0*x**3.0+x**2.0))

        self.u = uexata
        self.v = vexata
        self.pressao = pexata
        self.Bs = Bs

class correction():
    def __init__ (self, Nx,Ny,dx,dy,pqml,pl):

        pressao = np.zeros((Nx*Ny,1), dtype = np.float64)
        pressao = pqml+pl


        j = Ny-1
        for i in range(1,Nx-1):
            p = i + j*Nx
            ps = p - Nx
            pss= ps - Nx
            
            #pressao[p] = pqml[pss - Nx] -3.0*pqml[pss]+3.0*pqml[ps]
            pressao[p] = 2.0*pqml[ps]-pqml[pss]
        #contorno sul
        j = 0
        for i in range(1,Nx-1):
            p = i + j*Nx
            pn = p + Nx
            pnn= pn + Nx
            
            #pressao[p] = pqml[pnn + Nx] -3.0*pqml[pnn]+3.0*pqml[pn]
            pressao[p] = 2.0*pqml[pn]-pqml[pnn]
        #contorno leste
        i= Nx-1
        for j in range(1,Nx-1):
            p = i + j*Nx
            pw = p -1
            pww= pw -1
            
            #pressao[p] = pqml[pww -1] -3.0*pqml[pww]+3.0*pqml[pw]
            pressao[p] = 2.0*pqml[pw]-pqml[pww]

        #contorno oeste
        i= 0
        for j in range(1,Nx-1):
            p = i + j*Nx
            pe = p +1
            pee= pe +1
            
            #pressao[p] = pqml[pee + 1] -3.0*pqml[pee]+3.0*pqml[pe]
            pressao[p] = 2.0*pqml[pe]-pqml[pee]
        self.press = pressao 































