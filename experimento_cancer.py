#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:36:03 2019

@author: laura
"""



#bibliotecas

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#funcao f

global s1,s2

global g1,g2

global alfa,beta

global eta

global phi, phi1, phi2

global mi1, mi2

global S0, R0

global teta1, teta2
#####parametros######
s1,s2,gama1,gama2,alfa,eta,phi,mi1,mi2,phi1=0.192,0.084,0.15,0.02,1./3.,0.02,0.1,0.1,0.15,0.38

beta=2./3. - alfa
S0,R0=10.*10**(-6),2.*10**(-6)
teta1,teta2=2.,0.02


def F1(X):
    S=X[0,0]
    R=X[1,0]
    C=X[2,0]
    f1=-s1*S*np.log((R+S)/C)-gama1*S+gama2*R
    f2=-s2*R*np.log((R+S)/C)+gama1*S-gama2*R
    f3=(R+S)**((2./3.)-alfa)*C**(1.-beta)*(1-((R+S)**alfa)*(C**beta))-eta*C
    saida=np.matrix([[f1],[f2],[f3]])
    return saida



def F2(X,U):
    S=X[0,0]
    R=X[1,0]
    C=X[2,0]
    Q=X[3,0]#D1
    A=X[4,0]#D2
    u=U[3,0]
    v=U[4,0]
#    v=1.
#    u=1.
    f1=-s1*S*np.log((R+S)/C)-(gama1+phi*Q)*S+gama2*R
    f2=-s2*R*np.log((R+S)/C)+gama1*S-gama2*R
    f3=(R+S)**((2./3.)-alfa)*C**(1.-beta)*(1.-((R+S)**alfa)*(C**beta))-(eta+mi1*Q+mi2*v)*C
    f4=u-phi1*Q
    f5=v
    saida=np.matrix([[f1],[f2],[f3],[f4],[f5]])
    return saida

def G(X,U,L):
    S=X[0,0]
    R=X[1,0]
    C=X[2,0]
    Q=X[3,0]#D1
    A=X[4,0]#D2
    l1=L[0,0]
    l2=L[1,0]
    l3=L[2,0]
    l4=L[3,0]
    l5=L[4,0]
    u=U[3,0]
    v=U[4,0]
    g1=l1*(s1*(np.log((R+S)/C)+(S/(R+S)))+gama1-phi*Q)+l2*(s2*(R/(R+S))-gama1)+l3*((-C**(1.-beta))*((2./3.)-alfa)*(R+S)**(-alfa-1./3.)+(2./3.)*C*(R+S)**(-1./3.))-1./S0
    g2=l1*(s1*S/(R+S)-gama2)+l2*(s2*(np.log((R+S)/C)+(R/(S+R)))+gama2)+l3*((-C**(1.-beta))*(2./3.-alfa)*(R+S)**(-alfa-1./3.)+(2./3.)*C*(R+S)**(1./2.))-1./R0
    g3=-l1*s1*S/C-l2*s2*R/C+l3*(-((R+S)**(2./3.-alfa))*(1.-beta)*C+(R+S)**(2./3.)+eta+mi1*Q+mi2*v)
    g4=l1*phi*S+l3*mi1*C+l4*phi1
    g5=0.
    saida=np.matrix([[g1],[g2],[g3],[g4],[g5]])
    return saida



#condicao inicial

S,R,C,Q,A,u,v=S0,R0,15.*10**(-6),0.,100.,0.,0.#A=0.002
k=1.#mudar essa cte dps
l1,l2,l3,l4,l5=1.,1.,0.,0.,k
X0=np.zeros((5,1))
X0[0,0],X0[1,0],X0[2,0],X0[3,0],X0[4,0]=S,R,C,Q,A
L0=np.zeros((5,1))
L0[0,0],L0[1,0],L0[2,0],L0[3,0],L0[4,0]=l1,l2,l3,l4,l5
U0=np.zeros((5,1))
U0[0,0],U0[1,0],U0[2,0],U0[3,0],U0[4,0]=0.,0.,0.,0.,0.
#variaveis auxiliares

num=100 #numero de iteracoes

#####inicio sem controle#####

Y0=np.zeros((3,1))
Y0[0,0],Y0[1,0],Y0[2,0]=S,R,C


#variaveis auxiliares

num=100 #numero de iteracoes

Y1=np.zeros((3,1)) #vetor do passo seguinte

dt=1. #taxa de variacao do tempo
t=np.arange(0,num,dt)
Ry=np.zeros((3,num)) #matriz de resultados

Ry[:,0]=Y0[:,0]

#aplicandoo o metodo

for k in range(1,num):
    K1=F1(Y0)
    K2=F1(Y0+dt*K1*0.5)
    K3=F1(Y0+dt*K2*0.5)
    K4=F1(Y0+dt*K3)
    Y1=Y0+(dt/6)*(K1+2*K2+2*K3+K4)
    
    Ry[0,k]=Y1[0,0]
    Ry[1,k]=Y1[1,0]
    Ry[2,k]=Y1[2,0]
    Y0=Y1

#####fim sem controle######


#####inicio com controle#####


X1=np.zeros((5,1)) #vetor do passo seguinte
L1=np.zeros((5,1))
U1=np.zeros((5,1))

dt=1#taxa de variacao do tempo
t=np.arange(0,num,dt)
R1=np.zeros((5,num)) #matriz de resultados
R2=np.zeros((5,num))
R1[:,0]=X0[:,0]
L2=np.zeros((5,num))
umax=1./2.
vmax=75.
Qmax=2./10.
Amax=300.
delta,test=0.001,-1
#aplicando o metodo
i=0
while ((test < 0)and(X1[3,0]>=0)and(X1[3,0]<=Qmax)and(X1[4,0]<=Amax)and(X1[4,0]>=0)):
    
    oldu = U0
    oldx = X0
    oldlambda = L0    
    for k in range(1,num):
        K1=F2(X0,U0)#######################################
        K2=F2(X0+dt*K1*0.5,(U0+U1)*0.5)
        K3=F2(X0+dt*K2*0.5,(U0+U1)*0.5)
        K4=F2(X0+dt*K3,U1)
        X1=X0+(dt/6)*(K1+2*K2+2*K3+K4)
        
        R1[0,k]=X1[0,0]
        R1[1,k]=X1[1,0]
        R1[2,k]=X1[2,0]
        R1[3,k]=X1[3,0]
        R1[4,k]=X1[4,0]
        X0=X1
        
    
    for k in range(1,num):
        j=num-k
        X1[0,0]=R1[0,j]
        X1[1,0]=R1[1,j]
        X1[2,0]=R1[2,j]
        X1[3,0]=R1[3,j]
        X1[4,0]=R1[4,j]
        X0[0,0]=R1[0,j-1]
        X0[1,0]=R1[1,j-1]
        X0[2,0]=R1[2,j-1]
        X0[3,0]=R1[3,j-1]
        X0[4,0]=R1[4,j-1]
        
       
        K1=G(X1,U0,L0)
        K2=G(0.5*(X0+X1),U0,L0-dt*K1*0.5)
        
        K3=G(0.5*(X0+X1),U0,L0-dt*K2*0.5)
        K4=G(X0,U0,L0-dt*K3)
        L1=L0-(dt/6)*(K1+2*K2+2*K3+K4)
        L0=L1
        R2[1,k]=L2[1,0]
        R2[2,k]=L2[2,0]
        R2[3,k]=L2[3,0]
        R2[4,k]=L2[4,0]
    i=i+1
    c=R1[2,i]
    u1=(-1/teta1)*L0[3,0]
    v1=(-1/teta2)*(L0[4,0]-mi2*c*L0[2,0])
    if (u1>=umax):
        U1[3,0]=umax
    elif(u1<=0):
        U1[3,0]=0.
    else:
        U1[3,0] = (-1/teta1)*L0[3,0]    
    if(v1>=vmax):
        U1[4,0]=vmax
    elif(v1<=0):
        U1[4,0]=0.
    else:
        U1[4,0] = (-1/teta2)*(L0[4,0]-mi2*c*L0[2,0])
    U0 = 0.5*(U1 + oldu)
    
    R2[0,i]=U0[0,0]
    R2[1,i]=U0[1,0]
    R2[2,i]=U0[2,0]
    R2[3,i]=U0[3,0]
    R2[4,i]=U0[4,0]
    
    #%Convergence Test
    temp1 = delta*np.sum(np.linalg.norm(U0)) - np.sum(np.linalg.norm(oldu-U0))
    temp2 = delta*np.sum(np.linalg.norm(X0)) - np.sum(np.linalg.norm(oldx-X0))
    temp3 = delta*np.sum(np.linalg.norm(L0)) - np.sum(np.linalg.norm(oldlambda-L0))
    print("temps = ",temp1, temp2, temp3,"i=",i)
    
    test = np.min([temp1, temp2,temp3])

#####fim com controle#####

######grafico 3d

#####com controle#####

xr=R1[0,:]
yr=R1[1,:]
zr=R1[2,:]
kr=R1[3,:]
wr=R1[4,:]
ur=R2[3,:]
vr=R2[4,:]

lu=L2[3,:]
lv=L2[3,:]


#####sem controle####
xyr=Ry[0,:]
yyr=Ry[1,:]
zyr=Ry[2,:]

#fig0y=plt.figure()
#ax0y=fig0y.add_subplot(111,projection='3d')
#ax0y.scatter(xyr,yyr,zyr)
#plt.show()

fig0=plt.figure()
#fig1=plt.figure()
ax0=fig0.add_subplot(111,projection='3d')
#ax1=fig1.add_subplot(111,projection='3d')
#ax0.plot(zr,zr,zr)
ax0.plot(xr,yr,zr,color='red')
ax0.scatter(xyr,yyr,zyr)
plt.show()
####grafico 2d
#fig1=plt.figure()
#ax1=fig1.add_subplot(3,1,1)
#ax2=fig1.add_subplot(3,1,2)
#ax3=fig1.add_subplot(3,1,3)
linha11,=plt.plot( t, xr, label="Células sensíveis com controle", color='green')
linha12,=plt.plot(t,xyr,label="Células sensíveis sem controle",color='green',ls='-.')
plt.xlabel("t")
plt.ylabel("S")
plt.legend(handles=[linha11,linha12],loc='best')
plt.show()

linha21,=plt.plot( t, yr, label="Células resistentes com controle", color='blue')
linha22,=plt.plot(t,yyr, label="Células resistentes sem controle",color='blue',ls='-.')
plt.xlabel("t")
plt.ylabel("R")
plt.legend(handles=[linha21,linha22],loc='best')
plt.show()

linha31,=plt.plot( t, zr, label="Capacidade de carga com controle", color='red')
linha32,=plt.plot(t,zyr, label="Capacidade de carga sem controle",color='red',ls='-.')
plt.xlabel("t")
plt.ylabel("C")
plt.legend(handles=[linha31,linha32],loc='best')
plt.show()



linha4,=plt.plot( xr+yr, zr, label="Capacidade de carga com controle", color='pink')
#linha32,=plt.plot(t,zyr, label="Capacidade de carga sem controle",color='red',ls='-.')
plt.xlabel("R+S")
plt.ylabel("C")
plt.legend(handles=[linha4],loc='best')
plt.show()


#ax1.plot(t,xr,c='b')
#ax2.plot(t,yr,c='m')
#ax3.plot(t,zr,c='r')


linha4,=plt.plot( t, ur, label="controle u", color='purple')

plt.xlabel("t")
plt.ylabel("U")
plt.legend(handles=[linha4])
plt.show()

linha5,=plt.plot( t, vr, label="controle v", color='brown')

plt.xlabel("t")
plt.ylabel("V")
plt.legend(handles=[linha5])
plt.show()

#Xr=xr+yr
#
#linha6,=plt.plot( lu, ur, label="ul", color='brown')
#
#plt.xlabel("Lu")
#plt.ylabel("u")
#plt.legend(handles=[linha6])
#plt.show()
#
#
#
#linha7,=plt.plot( lv, vr, label="vl", color='gray')
#
#plt.xlabel("Lv")
#plt.ylabel("v")
#plt.legend(handles=[linha7])
#plt.show()
#
#
#linha8,=plt.plot( zr, ur, label="u", color='gray')
#
#plt.xlabel("C")
#plt.ylabel("u")
#plt.legend(handles=[linha8])
#plt.show()
#


linha9,=plt.plot( t, kr, label="Quimioterapia", color='gray')

plt.xlabel("t")
plt.ylabel("Q")
plt.legend(handles=[linha9])
plt.show()


linha10,=plt.plot( t, wr, label="Antiangiogênese", color='orange')

plt.xlabel("t")
plt.ylabel("A")
plt.legend(handles=[linha10])
plt.show()


iss=0
maiors=0
for i in range(num):
    if R1[0,i]>=maiors:
        maiors=R1[0,i]
        iss=i
        

maiorr=0
ir=0
for i in range (num):
    if R1[1,i]>=maiorr:
        maiorr=R1[1,i]
        ir=i

ic=0
maiorc=0
for i in range (num):
    if R1[2,i]>=maiorc:
        maiorc=R1[2,i]
        ic=i

iu=0
maioru=0
for i in range (num):
    if R2[3,i]>=maioru:
        maioru=R2[3,i]
        iu=i

iv=0
maiorv=0
for i in range (num):
    if R2[4,i]>=maiorv:
        maiorv=R2[4,i]
        iv=i




print ("pico s",maiors,"indice s",iss)

print ("pico r",maiorr,"indice r",ir)

print ("pico c",maiorc,"indice c",ic)

print ("pico u",maioru,"indice u",iu)

print ("pico v",maiorv,"indice v",iv)

print("resultado",R1[:,99])

print("resultado",R2[:,99])
#plt.scatter(xr,yr)
