#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:44:00 2019

@author: laura
"""
import numpy as np
import matplotlib.pyplot as plt

global xii,pinf,teta,alfa,beta,gamma,delta,p0,r0,v
xii,pinf,teta,alfa,beta,gamma,delta,v=0.5618,780.*10.**(6),1.*10.**(7),0.00484,0.00264,0.1181,0.37451,1.
p0,r0=600.*10.**(6),0.10


def F(X):
    p=X[0,0]
    r=X[1,0]
    f1=-xii*p*(1-p/pinf)**v-teta*p*r
    f2=alfa*r*(p-beta*p**2)+gamma-delta*r
    saida=np.matrix([[f1],[f2]])
    return saida

X0=np.zeros((2,1))
X1=np.zeros((2,1))
X0[0,0],X0[1,0]=p0,r0

num=1000
dt=1

t=np.arange(0,num,dt)
R=np.zeros((2,num))
R[:,0]=X0[:,0]


for k in range(1,num):
    K1=F(X0)
    K2=F(X0+dt*K1*0.5)
    K3=F(X0+dt*K2*0.5)
    K4=F(X0+dt*K3)
    X1=X0+(dt/6.)*(K1+2.*K2+2.*K3+K4)
    
    R[0,k]=X1[0,0]
    R[1,k]=X1[1,0]
    X0=X1
    
    
pr=R[0,:]
rr=R[1,:]

linha1,=plt.plot( t, pr, label="Volume de celulas tumorais", color='red')
plt.xlabel("t")
plt.ylabel("p")
plt.legend(handles=[linha1],loc='best')
plt.show()

linha2,=plt.plot( t, rr, label="Células imunológicas", color='navy')
plt.xlabel("t")
plt.ylabel("r")
plt.legend(handles=[linha2],loc='best')
plt.show()

linha3,=plt.plot( pr, rr, label="Retrato de fase", color='purple')
plt.xlabel("p")
plt.ylabel("r")
plt.legend(handles=[linha3],loc='best')
plt.show()

