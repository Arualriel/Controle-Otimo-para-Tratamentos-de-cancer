#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 21:51:06 2020

@author: laura
"""


#bibliotecas

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#funcao f

global p0,q0,r0,s,m,gama,alfa,beta,delta,xi,b,d,mi,phi1,phi2,phi3,u0,u1,umax,A,B,C
#####parametros######
p0,q0,r0,s,m,gama=.20,.30,0.0001,10.,1.,0.01 ##s=10. ou 0.01
alfa,beta,delta,teta=0.0529,0.00291,0.3743,1.
xi,b,d,mi,phi1=0.0347,5,0.0667,0.0,0.005
phi2,phi3,u0,u1,umax=0.06,0.02,0.0015,0.0015,100.
######artigo######
A,B,C=0.9068,0.4216,0.0029
##################
def F(X,u):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    f1=-xi*p*(np.log(p/q))-teta*p*r-phi1*u
    f2=b*p-(mi+d*p**(2./3.))*q-phi2*u
    f3=alfa*(p-beta*p**2)*r+gama-delta*r+phi3*u
    saida=np.matrix([[f1],[f2],[f3]])
    return saida

def G(X,u,L):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    l1=L[0,0]
    l2=L[1,0]
    l3=L[2,0]
    e1=l1*(xi*(1.+np.log(p/q))+teta*r+phi1*u)-l2*(b-(2./3.)*d*q*p**(-1./3.))
    e2=-l3*(alfa*(1.-2.*beta*p)*r)
    g1=e1+e2
    g2=-l1*xi*(p/q)+l2*(mi+d*p**(2./3.)+phi2*u)
    g3=l1*teta*p-l3*(alfa*(p-beta*p**2.)-delta+phi3*u)
    saida=np.matrix([[g1],[g2],[g3]])
    return saida

def colchete(X):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    c1=(phi2-phi1)*xi*p+phi3*teta*r
    c2=(phi1-phi2)*b*p-(2./3.)*phi1*d*q*p**(2./3.)
    c3=phi1*alfa*(p-2*beta*p**2)*r+phi3*gama
    saida=np.matrix([[c1],[c2],[c3]])
    return saida

def ge(X):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    ge1=-phi1*p
    ge2=-phi2*q
    ge3=phi3*r
    saida=np.matrix([[ge1],[ge2],[ge3]])
    return saida

def colchete1(X):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    a1=(xi**2.)*(phi2-phi1)*p+xi*phi3*r*teta*p
    a2=phi3*teta*p*alfa*(p-beta*p**2.)*r-phi3*teta*p*delta*r
    a3=phi3*teta*p*gama-xi*(p/q)*(phi1-phi2)*b*p
    a4=xi*(p/q)*(2./3.)*phi1*d*q*p**(2./3.)
    a5=teta*p*phi1*alfa*(p-2.*beta*p**2.)*r+teta*p*phi3*gama
    b1=-(phi1-phi2)*b*xi*p*np.log(p/q)+(4./9.)*phi1*d*q*xi*(np.log(p/q))*p**(2./3.)
    b2=(4./9.)*teta*p*r*phi1*d*q*p**(-1./3.)-(2./3.)*phi1*b*d*p**(5./3.)
    b3=(2./3.)*phi1*d*q*(mi+d*p**(2./3.))*p**(2./3.)
    b4=-phi3*r*teta*p*b+phi3*r*teta*(2./3.)*d*q*p**(2./3.)
    b5=(mi+d*p**(2./3.))*(phi1-phi2)*b*p-(mi+d*p**(2./3.))*(2./3.)*phi1*d*q*p**(2./3.)
    b6=-teta*p*r*(phi1-phi2)*b-(phi2-phi1)*p*xi*b+(phi2-phi1)*xi*(2./3.)*q*d*p**(2./3.)
    d1=-phi1*alfa*(1.-4.*beta*p)*r*xi*p*np.log(p/q)-teta*p*r*phi1*alfa*(1.-4.*beta*p)*r
    d2=phi1*alfa*(p-2.*beta*p**2.)*(p-beta*p**2.)*r+phi1*alfa*(p-2.*beta*p**2.)*gama
    d3=-phi1*alfa*(p-2.*beta*p**2.)*delta*r-alfa*(1.-2.*beta*p)*r*(phi2-phi1)*p*xi
    d4=-alfa*(1.-2.*beta*p)*r*phi3*r*teta*p
    d5=-alfa*(p-beta*p**2.)*phi3*gama+delta*phi1*alfa*(p-2*beta*p**2.)*r+delta*phi3*gama
    d6=-alfa*(p-beta*p**2.)*phi1*alfa*(p-2*beta*p**2)*r
    c11=a1+a2+a3+a4+a5
    c12=b1+b2+b3+b4+b5+b6
    c13=d1+d2+d3+d4+d5+d6
    saida=np.matrix([[c11],[c12],[c13]])
    return saida

def colchete2(X):
    p=X[0,0]
    q=X[1,0]
    r=X[2,0]
    c21=(phi3**2.)*r*teta*p
    c22=(phi1**2.-phi2**2.)*b*p+(4./9.)*(phi1**2.)*d*q*p**(2./3.)
    c23=-(phi1**2.)*alfa*(p-4*beta*p**2)*r-gama*phi3**2.
    saida=np.matrix([[c21],[c22],[c23]])
    return saida

    



#condicao inicial

X0=np.zeros((3,1))
X0[0,0],X0[1,0],X0[2,0]=p0,q0,r0

#condicao final

L0=np.zeros((3,1))
L0[0,0],L0[1,0],L0[2,0]=A,B,-C

#variaveis auxiliares

num=49 #numero de iteracoes

X1=np.zeros((3,1)) #vetor do passo seguinte
L1=np.zeros((3,1))

dt=1. #taxa de variacao do tempo
t=np.arange(0,num,dt)
Rx=np.zeros((3,num)) #matriz de resultados
Rl=np.zeros((3,num))
Ru=np.zeros((1,num))
g=np.zeros((3,0))

Rx[:,0]=X0[:,0]
Rl[:,0]=L0[:,0]
Ru[0,0]=u0

delta,test=0.001,-1

#aplicandoo o metodo


M=np.zeros((3,3))





i=1
while (test < 0) and (i<num):
    
    oldu = u0
    oldx = X0
    oldlambda = L0    
    for k in range(1,num):
        K1=F(X0,u0)#######################################
        K2=F(X0+dt*K1*0.5,(u0+u1)*0.5)
        K3=F(X0+dt*K2*0.5,(u0+u1)*0.5)
        K4=F(X0+dt*K3,u1)
        X1=X0+(dt/6.0)*(K1+2.0*K2+2.0*K3+K4)
        
        Rx[0,k]=np.float(X1[0,0])
        Rx[1,k]=np.float(X1[1,0])
        Rx[2,k]=np.float(X1[2,0])
        X0=X1
        
    
    for k in range(1,num):
        j=num-k
        X1[0,0]=Rx[0,j]
        X1[1,0]=Rx[1,j]
        X1[2,0]=Rx[2,j]
        X0[0,0]=Rx[0,j-1]
        X0[1,0]=Rx[1,j-1]
        X0[2,0]=Rx[2,j-1]
        #c=R1[2,j-1]
        K1=G(X1,u0,L0)
        K2=G(0.5*(X0+X1),u0,L0-dt*K1*0.5)
        
        K3=G(0.5*(X0+X1),u0,L0-dt*K2*0.5)
        K4=G(X0,u0,L0-dt*K3)
        L1=L0-(dt/6.0)*(K1+2.0*K2+2.0*K3+K4)
        L0=L1
        Rl[0,k]=L1[0,0]
        Rl[1,k]=L1[1,0]
        Rl[2,k]=L1[2,0]
        
    
    i=i+1
    g=ge(X0)
    Col=colchete(X0)
    M[0,0]=X0[0,0]
    M[1,0]=X0[1,0]
    M[2,0]=X0[2,0]
    M[0,1]=g[0,0]
    M[1,1]=g[1,0]
    M[2,1]=g[2,0]
    M[0,2]=Col[0,0]
    M[1,2]=Col[1,0]
    M[2,2]=Col[2,0]
    N1=colchete1(X0)
    N2=colchete2(X0)
    sigma=np.linalg.solve(M,N1)
    ro=np.linalg.solve(M,N2)
    u=np.float(-(sigma[0]+sigma[1])/(ro[0]+ro[1]))
    
    
    
    if (u>=umax):
        u1=umax
    elif(u<=0):
        u1=0.
    else:
        u1 = u    
    
    u0 = 0.5*(u1 + oldu)
    
    Ru[0,i-1]=u0

    
    
    
    
    #%Convergence Test
    #temp1 = delta*np.sum(np.linalg.norm(U0)) - np.sum(np.linalg.norm(oldu-U0))
    temp1 = delta*np.sum(np.abs(u0)) - np.sum(np.abs(oldu-u0))
    
    temp2 = delta*np.sum(np.linalg.norm(X0)) - np.sum(np.linalg.norm(oldx-X0))
#    temp2 = delta*np.sum(np.linalg.norm(X0)) - np.sum(np.linalg.norm(oldx-X0))
    temp3 = delta*np.sum(np.linalg.norm(L0)) - np.sum(np.linalg.norm(oldlambda-L0))
    
    test = np.min([temp1, temp2,temp3])
    #print(test)






######grafico 3d
P=Rx[0,:]
Q=Rx[1,:]
R=Rx[2,:]
U=Ru[0,:]
print(U)
L1=Rl[0,:]
L2=Rl[1,:]
L3=Rl[2,:]
fig0=plt.figure()
#fig1=plt.figure()
ax0=fig0.add_subplot(111,projection='3d')
#ax1=fig1.add_subplot(111,projection='3d')
#ax0.plot(zr,zr,zr)
ax0.plot(P,Q,R)
plt.show()
#fig1=plt.figure()
#ax1=fig1.add_subplot(111,projection='3d')
#ax1.plot(P,Q,Ru)
#plt.show()
n=7
x = range(n)
y = range(n)
X,Y = np.meshgrid(x,y)
Zinter = np.reshape(U,(n,n))
Z = np.zeros((n,n))

  
for i in range(n):
    for j in range(n):
        Z[i,j]=Zinter[i,j]

fig = plt.figure(figsize=plt.figaspect(0.5))

ax1 = fig.add_subplot(111, projection='3d')
for i in range(n):
	for j in range(n):
		x = X[i,j]
		y = Y[i,j]
		z = Z[i,j]
		
		ax0.scatter(x,y,z)

ax1.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,	linewidth=0.2, antialiased=False)

plt.show()



####grafico 2d
#fig1=plt.figure()
#ax1=fig1.add_subplot(3,1,1)
#ax2=fig1.add_subplot(3,1,2)
#ax3=fig1.add_subplot(3,1,3)
linha1,=plt.plot( t, P, label="Volume tumoral", color='red')

plt.xlabel("t")
plt.ylabel("P")
plt.legend(handles=[linha1])
plt.show()

linha2,=plt.plot( t, Q, label="Capacidade de carga tumoral", color='navy')

plt.xlabel("t")
plt.ylabel("Q")
plt.legend(handles=[linha2])
plt.show()

linha3,=plt.plot( t, R, label="Densidade de células imunocompetentes", color='green')

plt.xlabel("t")
plt.ylabel("R")
plt.legend(handles=[linha3])
plt.show()



#ax1.plot(t,xr,c='b')
#ax2.plot(t,yr,c='m')
#ax3.plot(t,zr,c='r')


linha4,=plt.plot( t, U, label="Controle u", color='purple')

plt.xlabel("t")
plt.ylabel("U")
plt.legend(handles=[linha4])
plt.show()

linha5,=plt.plot( P, Q, label="-----", color='orange')

plt.xlabel("P")
plt.ylabel("Q")
plt.legend(handles=[linha5])
plt.show()

linha6,=plt.plot( P, R, label="------", color='brown')

plt.xlabel("P")
plt.ylabel("R")
plt.legend(handles=[linha6])
plt.show()

linha7,=plt.plot( Q, R, label="-------", color='pink')

plt.xlabel("Q")
plt.ylabel("R")
plt.legend(handles=[linha7])
plt.show()

ip=0
maiorp=0
for i in range(num):
    if Rx[0,i]>=maiorp:
        maiorp=Rx[0,i]
        ip=i
        

maiorq=0
iq=0
for i in range (num):
    if Rx[1,i]>=maiorq:
        maiorq=Rx[1,i]
        iq=i

ir=0
maiorr=0
for i in range (num):
    if Rx[2,i]>=maiorr:
        maiorr=Rx[2,i]
        ir=i

iu=0
maioru=0
for i in range (num):
    if Ru[0,i]>=maioru:
        maioru=Ru[0,i]
        iu=i



print ("pico p",maiorp,"indice p",ip)

print ("pico q",maiorq,"indice q",iq)

print ("pico r",maiorr,"indice r",ir)

print ("pico u",maioru,"indice u",iu)

print("resultado",Rx[:,num-1])

print("resultado",Ru[:,num-1])
#plt.scatter(xr,yr)





