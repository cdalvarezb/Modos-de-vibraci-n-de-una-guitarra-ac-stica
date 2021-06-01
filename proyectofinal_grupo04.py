# -*- coding: utf-8 -*-
"""
Created on Fri May 28 04:06:36 2021

@author: cdalv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
    


POS = pd.read_excel('malla_guitarra_python_569_nodes.xlsx',dtype={'ejex':float,
'ejey':float,'ejez':float},sheet_name='POS') #  Lee el excel con los puntos de posición
LINES = pd.read_excel('malla_guitarra_python_569_nodes.xlsx',dtype={'a':float,
'b':float,'c':float},sheet_name='LINES')
TRIANGLES = pd.read_excel('malla_guitarra_python_569_nodes.xlsx',dtype={'d':float,
'e':float,'f':float,'g':float},sheet_name='TRIANGLES') #  Lee el excel con los elementos de la discretización


a = pd.DataFrame(POS)
b = pd.DataFrame(LINES)
c = pd.DataFrame(TRIANGLES)


mshPOS = a.to_numpy()  #  Posición de los puntos
mshLINES = b.to_numpy()  #  Puntos de la superficie
mshQUADS = c.to_numpy()  #  Matriz de conectividad
mshPNT = np.array(len(mshPOS))  #  Número de nodos




'''
                            PARÁMETROS DEL MATERIAL
'''

rho = 570  #  Densidad básica
#Elasticity = 94*(0.01)**2  #  Módulo de elasticidad

Elasticity = 8.9*10**8  #  Módulo de elasticidad



'''
                      LECTOR DE LA TAPA DE LA GUITARRA
'''

#def reading_membrane(net,n,jnb,xe,nee,nbc,nee1,n1):
net = len(mshQUADS)
n = len(mshPOS)
n1 = len(mshLINES)
ni1 = 2
xe = np.zeros((3*n,1))
ni = 3
nee = np.zeros((net,ni))
nee1 = np.zeros((n1,ni1))

ii = 0

for i in range(n):
    ii= ii + 1
    xe[3*ii - 3] = mshPOS[i, 0]
    xe[3*ii - 2] = mshPOS[i, 1]
    xe[3*ii - 1] = mshPOS[i, 1]

for i in range(net):
    nee[i, 0:ni] = mshQUADS[i, 0:ni]

for i in range(n1):
    nee1[i, 0:ni1] = mshLINES[i, 0:ni1]


##########  DIBUJA TODOS LOS NODOS Y LOS TRIANGULOS DE LA MALLA  #################

xee = np.zeros((ni,2))

p1 = []
k1 = []

for i in range(net):
    k2 = []
    p2 = []
    for j in range(ni):

        p1.append(xe[int(3*nee[i, j] - 3)])
        k1.append(xe[int(3*nee[i, j] - 2)])
        p2.append(xe[int(3*nee[i, j] - 3)])
        k2.append(xe[int(3*nee[i, j] - 2)])
        xee[j, 0] = xe[int(3*nee[i, j]- 3)]
        xee[j, 1] = xe[int(3*nee[i, j] - 2)]
    p2.append(xe[int(3*nee[i, 0] - 3)])
    k2.append(xe[int(3*nee[i, 0] - 2)])

    plt.plot(p2,k2,color="black")
    plt.title('visualización de la malla en python')
    plt.xlabel('Eje x [mm]')
    plt.ylabel('Eje y [mm]')

plt.scatter(p1,k1)


#######   DIBUJA DE COLOR VERDE LOS NODOS DE LOS BORDES    ###########
p2 = []
k2 = []

for i in range(n1):
    for j in range(ni1):

        p2.append(xe[int(3*nee1[i,j]-3)])
        k2.append(xe[int(3*nee1[i,j]-2)])

plt.scatter(p2,k2, color = '#88c999')
plt.savefig('malla_en_python.pdf')
plt.show()

jnb = np.size(mshPNT)
nbc = np.zeros((jnb, 1))



'''
        FINITE ELEMENTS METHOD (FEM) PARA HALLLAR LOS MODOS DE FRECUENCIA

'''

nn = mshPNT

u = np.zeros((nn, 1))

nbc = np.zeros((nn, 1))
nbp = np.zeros((nn, 1))
q = np.zeros((nn, 1))

jnb = 0

for i in range(n1):
    jnb = jnb + 1
    nbc[jnb, 0] = nee1[i, 0]
    for j in range(ni1):
        u[ int(nee1[i,j]) ] = int(0)


kee = np.zeros((nn, nn))
qee = np.zeros((nn, 1))
Mee = np.zeros((nn, nn))
Bigke = np.zeros((nn, nn))
BigMe = np.zeros((nn, nn))

nf = 3

for i in range(net):

    keu = np.zeros((nf, nf))
    Meu = np.zeros((nf, nf))
    qeu = np.zeros((nf, 1))
    dfdx = np.zeros((nf, 2))
    f = np.zeros((nf, 1))

    ni = 3
    xi = np.zeros((ni, 2))
    wi = np.ones((ni, 2))
    wi[:,:] = 1/3

    xi[0, 0] = 0.5
    xi[0, 1] = 0
    xi[1, 0] = 0
    xi[1, 1] = 0.5
    xi[2, 0] = 0.5
    xi[2, 1] = 0.5

    Jacb = np.zeros((2,2))

    for k1 in range(nf):

        int_kuu = 0

        for ii in range(ni):

            for jj in range(ni):

                xii = 0.5*( xi[ii, 0] + 1)
                yii = 0.5*( xi[jj, 1] + 1)

                f[1] = yii
                dfdx[1, 0] = 0.0
                dfdx[1, 1] = 1.0

                f[0] = 1 - yii - xii
                dfdx[0, 0] = -1
                dfdx[0, 1] = -1

                f[2] = xii
                dfdx[2, 0] = 1.0
                dfdx[2, 1] = 0.0

                for j1 in range(1,3):

                    jnx = 0
                    jny = 0

                    for i1 in range(nf):

                        jnx = jnx + xe[int(3*nee[i, i1]) - 3]*dfdx[i1, j1-1]
                        jny = jny + xe[int(3*nee[i, i1]) - 2]*dfdx[i1, j1-1]

                    Jacb[j1-1, 0] = jnx
                    Jacb[j1-1, 1] = jny

                jd = np.linalg.det(Jacb)
                int_kuu =  int_kuu + jd*wi[ii-1, 0]*wi[jj-1, 1]*(q[int(nee[i, k1]-1), 0]*f[k1-1])

        qeu[k1, 0] = int_kuu
        ii = nee[i, k1]
        qee[int(ii-1)][0] = qee[int(ii-1), 0] + qeu[int(k1-1), 0]

    #############################################################

    ni = 3
    xi = np.zeros((ni, 2))
    wi = np.ones((ni, 2))
    wi[:,:] = 1/3

    xi[0, 0] = 0.5
    xi[0, 1] = 0
    xi[1, 0] = 0
    xi[1, 1] = 0.5
    xi[2, 0] = 0.5
    xi[2, 1] = 0.5

    Jacb = np.zeros((2, 2))

    for k1 in range(nf):
        for k2 in range(nf):

            int_kuu = 0
            int_kee = 0

            for ii in range(ni):
                for jj in range(ni):

                    xii = 0.5*( xi[ii, 0] + 1)
                    yii = 0.5*( xi[jj, 1] + 1)

                    f[1] = yii

                    f[1] = yii
                    dfdx[1, 0] = 0.0
                    dfdx[1, 1] = 1.0

                    f[0] = 1 - yii - xii
                    dfdx[0, 0] = -1
                    dfdx[0, 1] = -1

                    f[2] = xii
                    dfdx[2, 0] = 1.0
                    dfdx[2, 1] = 0.0

                    for j1 in range(1,3):
                        jnx = 0
                        jny = 0

                        for i1 in range(nf):

                            jnx = jnx + xe[int(3*nee[i, i1]) - 3]*dfdx[i1, j1-1]
                            jny = jny + xe[int(3*nee[i, i1]) - 2]*dfdx[i1, j1-1]

                        Jacb[j1-1, 0] = jnx
                        Jacb[j1-1, 1] = jny

                    jd = np.linalg.det(Jacb)
                    iJacob = np.linalg.inv(Jacb)

                    dxx = (iJacob[0, 0]*dfdx[k1, 0] + iJacob[0, 1]*dfdx[k1, 1])*(iJacob[0, 0]*dfdx[k2, 0] + iJacob[0, 1]*dfdx[k2, 1])
                    dyy = (iJacob[1, 0]*dfdx[k1, 0] + iJacob[1, 1]*dfdx[k1, 1])*(iJacob[1, 0]*dfdx[k2, 0] + iJacob[1, 1]*dfdx[k2, 1] )

                    int_kuu = int_kuu + jd*wi[ii, 0]*wi[jj, 1]*(dxx + dyy)
                    int_kee =  int_kee + jd*wi[ii, 0]*wi[jj, 1]*f[k1]*f[k2]

            keu[k1, k2] = int_kuu
            Meu[k1, k2] = int_kee
            ii = nee[i, k1]
            jj = nee[i, k2]
            kee[int(ii-1)][int(jj-1)] =  float(kee[int(ii-1), int(jj-1)]) + Elasticity*float(keu[k1, k2])
            Mee[int(ii-1)][int(jj-1)] = float(Mee[int(ii-1), int(jj-1)]) + rho*float(Meu[k1, k2])




for i in range(jnb):

    ii = nbc[i, 0]
    kee[int(ii)][:] = 0
    kee[int(ii)][int(ii)] = 1
    qee[int(ii), 0] = u[int(ii), 0]


nnl = nn - jnb

Bigkee = np.zeros((nnl, nnl))
BigMee = np.zeros((nnl, nnl))
Ckee = np.zeros((nnl, nnl))

Bigkee[:, :] = kee[ jnb + 0:nn , jnb + 0:nn ]
BigMee[:, :] = Mee[ jnb + 0:nn , jnb + 0:nn ]
Ckee[:, :] = np.matmul(np.linalg.inv(BigMee),Bigkee)
b_bonito = np.linalg.inv(BigMee)

# lambda son los valores propios, V son los vectores propios

lambda1,V = np.linalg.eig(Ckee)
u1 = np.zeros((nnl, 1))


for jj in range(nnl):
    
    ine = len(mshLINES)

    lambda1[jj]
    u1[:, 0] = np.real(V[:, jj])
    u[jnb + 0:nn ] =  u1

    uii = np.zeros((ni, 1))

    triag = tri.Triangulation(mshPOS[:,0], mshPOS[:,1], mshQUADS[:,:3]-1)
    #v = np.zeros(71)
    v = np.zeros(mshPNT)
    #v[30:] = V[jj]
    v[ine:] = V[jj]
    plt.tricontourf(triag, v)

    p1 = []
    k1 = []

    for i in range(net):
        k2 = []
        p2 = []
        for j in range(ni):
            p2.append(xe[int(3*nee[i, j] - 3)])
            k2.append(xe[int(3*nee[i, j] - 2)])
            p2.append(xe[int(3*nee[i, 0] - 3)])
            k2.append(xe[int(3*nee[i, 0] - 2)])
            
            aaa = round(float(float(lambda1[jj])), 4)
            plt.plot(p2,k2,color="black",linewidth=.3)
            plt.title(f'Madera caoba modo {jj+1}')
            plt.xlabel('Eje x [mm]')
            plt.ylabel('Eje y [mm]')

            #plt.savefig('x.pdf')

    print('Modo de frecuencia:', round(np.sqrt(lambda1[jj])/(2*np.pi),7), 'Hz')
    cbar = plt.colorbar()
    cbar.set_label('Deformación eje z [mm]')

    plt.show()

