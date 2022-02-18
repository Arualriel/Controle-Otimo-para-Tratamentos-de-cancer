import numpy as np

test = -1
delta = 0.01
N = 4000
t = np.linspace(1,5,N)
h = 1.0/N
h2 = h/2.0

u = np.zeros(N)
x = np.zeros(N)
x[0] = 2

lamb = np.zeros(N)
while (test < 0):
    oldu = u
    oldx = x
    oldlambda = lamb

    #%RK4 Forward For x
    for i in range(N-1):
        k1 = x[i] + u[i]
        k2 = x[i] + h2*k1 + 0.5*(u[i] + u[i+1])
        k3 = x[i] + h2*k2 + 0.5*(u[i] + u[i+1])
        k4 = x[i] + h*k3 + u[i+1]
        x[i+1] = x[i] + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
    #print("x=",x)
    #%RK4 Backwards For lambda   
    for i in range(1,N):
        j = N  - i
        #print("j=",j)
        k1 = -u[j] + 2*x[j] - lamb[j]
        k2 = -0.5*(u[j] + u[j-1]) - (lamb[j] - h2*k1) + (x[j] + x[j-1])
        k3 = -0.5*(u[j] + u[j-1]) - (lamb[j] - h2*k2) + (x[j] + x[j-1])
        k4 = -u[j-1] - (lamb[j] - h*k3) + 2*x[j-1]
        lamb[j-1] = lamb[j] - (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
    
    #print("lambda=",lamb)

    #%Update u
    u1 = (1/2)*(x + lamb)
    u = 0.5*(u1 + oldu)

    #%Convergence Test
    temp1 = delta*np.sum(np.abs(u)) - np.sum(np.abs(oldu - u))
    temp2 = delta*np.sum(np.abs(x)) - np.sum(np.abs(oldx - x))
    temp3 = delta*np.sum(np.abs(lamb)) - np.sum(np.abs(oldlambda - lamb))
    #print("temps = ",temp1, temp2, temp3)
    
    test = np.min([temp1, temp2,temp3])

import matplotlib.pyplot as plt

plt.plot(x,u)
plt.show()
