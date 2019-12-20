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


                    u[p]= -u0[pn] + 2.0*b[p]

        #contorno norte
        j = int(Ny-1)
        for i in range(1,(Nx-1)) :
                    p = int(i + (j)*Nx)
                    ps = int(p - Nx)

                    u[p]= -u0[ps] + 2.0*b[p]

        #contorno leste 
        i = int(Nx-1)
        for j in range(1,(Ny-1)) :
                    p = int(i + (j)*Nx)
                    pw = int(p - 1)
                    #ps = int(p - Nx)

                    u[p]= -u0[pw] + 2.0*b[p] 
        #oontorno oeste
        i = int(0)
        for j in range(1,(Ny-1)) :
                    p = int(i + (j)*Nx)
                    pe = int(p + 1)
 
                    u[p]= -u0[pe] + 2.0*b[p]

        if (np.linalg.norm(u-u0,np.inf) < tol):
                print('tolerancia atingida com', k, 'iteracoes')
                k = itmax+1
                
                #return u
                

        else:
            k += 1
            u0 = np.copy(u)

    return(u) 

def extrapolacao(Nx,Ny,pqml):
    #contorno norte
    j = Ny-1
    for i in range(1,Nx-1):
        p = i + j*Nx
        ps = p - Nx
        pss= ps - Nx
        pqml[p] = 2.0*pqml[ps]-pqml[pss]

    #contorno sul
    j = 0
    for i in range(1,Nx-1):
        p = i + j*Nx
        pn = p + Nx
        pnn= pn + Nx
        pqml[p] = 2.0*pqml[pn]-pqml[pnn]

    #contorno leste
    i= Nx-1
    for j in range(1,Nx-1):
        p = i + j*Nx
        pw = p -1
        pww= pw -1
        pqml[p] = 2.0*pqml[pw]-pqml[pww]

    #contorno oeste
    i= 0
    for j in range(1,Nx-1):
        p = i + j*Nx
        pe = p +1
        pee= pe +1
        pqml[p] = 2.0*pqml[pe]-pqml[pee]

    return(pqml)

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


            up[p]= -up[pn] + 2.0*bp_u[p]
            vp[p]= -vp[pn] + 2.0*bp_v[p]

    #contorno norte
    j = int(Ny-1)
    for i in range(1,(Nx-1)) :
            p = int(i + (j)*Nx)
            ps = int(p - Nx)

            up[p]= -up[ps] + 2.0*bp_u[p]
            vp[p]= -vp[ps] + 2.0*bp_v[p]

    #contorno leste 
    i = int(Nx-1)
    for j in range(1,(Ny-1)) :
            p = int(i + (j)*Nx)
            pw = int(p - 1)
            

            up[p]= -up[pw] + 2.0*bp_u[p] 
            vp[p]= -vp[pw] + 2.0*bp_v[p] 
    #oontorno oeste
    i = int(0)
    for j in range(1,(Ny-1)) :
            p = int(i + (j)*Nx)
            pe = int(p + 1)

            up[p]= -up[pe] + 2.0*bp_u[p]
            vp[p]= -vp[pe] + 2.0*bp_v[p] 







    return(up, vp)

def acoplamento_faces( Nx,Ny,rho,dt, dx,dy, u_e, v_n, pl):
    ue = np.zeros((Nx*Ny,1), dtype = np.float64)
    vn = np.zeros((Nx*Ny,1), dtype = np.float64)
    csimplec = simplec(Nx,Ny,rho,dt, dx,dy)
    for i in range(1,Nx-2):
        for j in range(1,Ny-2):
            p = i + (j)*Nx
            pe = p +  1 
            pw = p -1
            pn = p + Nx
            ps = p - Nx
            ue[p] = u_e[p] - csimplec.de[p]*(pl[pe]-pl[pw])/2.0 
            vn[p] = v_n[p] - csimplec.dn[p]*(pl[pn]-pl[ps])/2.0 

    return(ue, vn)



















