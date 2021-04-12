# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 16:27:31 2021

@author: Joyce
"""

"""
Created on Mon Mar 15 18:11:20 2021

@author: Lucas and Joyce
"""
from __future__ import print_function
import numpy as np
import sympy 
import math
from scipy.linalg import eigh
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.animation import FuncAnimation
from time import sleep

### Modelagem do sistema 

m_roda = 7 # kg
m_carro = 230
g = 9.81
L0 = 2.385 #[m]
L1 = (3/5)*1.7
L2 = 1.7-L1
L3 = 1.4

Ix_carro = 47
Iz_carro = 20

k2a = 150*9.81/0.01
k3a = 150*9.81/0.01
k4a = 150*9.81/0.01
k5a = 150*9.81/0.01
k2b = k2a/8
k3b = k3a/8
k4b = 1.5*k2b
k5b = 1.5*k3b

#### Declarando variaveis

y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1, y2, y3, y4, y5, gama, teta = sympy.symbols('y1_dot, y2_dot, y3_dot, y4_dot, y5_dot, gama_dot, teta_dot, y1, y2, y3, y4, y5, gama, teta')      
var = [y2, y3, y4, y5, y1, gama, teta]
var_dot = [y2_dot, y3_dot, y4_dot, y5_dot, y1_dot, gama_dot, teta_dot]

#### energia cinetica e potencial 

Ec = m_carro*(y1_dot**2)/2 + Ix_carro*(teta_dot**2)/2 + Iz_carro*(gama_dot**2)/2 + m_roda*((y2_dot**2)/2 + (y3_dot**2)/2 + (y4_dot**2)/2 + (y5_dot**2)/2) 
Ep = ((y1 + (gama*L1) + teta*(L3/2) - y3)**2)*k3b/2 + ((y1 - (gama*L2) + teta*(L3/2)- y2)**2)*k2b/2  + ((y1 + (gama*L2) - teta*(L3/2) - y4)**2)*k4b/2 + ((y1 - (gama*L2) - teta*(L3/2) - y5)**2)*k5b/2  + (y2**2)*k2a/2 + (y3**2)*k3a/2 + (y4**2)*k4a/2 + (y5**2)*k5a/2 + m_roda*g*(y1 + y2 + y3 + y4 + y5)


##### Lagrangeano 
L = Ec - Ep
L = sympy.simplify(L)

#### Derivadas em relacao aos deslocamentos
#### Derivadas em relacao as velocidades
#### Equecoes do movimento para cada grau de liberdade 
dL_dy = []
dL_dy_dot = []
graus = []

for i in range(len(var)):
    dL_dy.append(0)
    dL_dy_dot.append(0)
    graus.append(0)
    
    dL_dy[i] = sympy.diff(L, var[i])
    dL_dy_dot[i] = sympy.diff(L, var_dot[i])
    graus[i] = dL_dy_dot[i] + dL_dy[i]
    print(graus[i],'\n')
    
### Montando Matriz de Rigidez
K = []

for i in range(len(var)):
    k = []
    for j in range(len(var)):
        k.append(0)
    K.append(k)
    
print(K)

for i in range(len(var)):
    for j in range(len(var)):
       
       K[i][j] = graus[i].coeff(var[j],1)

for i in range(len(var)):
    for j in range(len(var)):
        if K[i][j] == K[j][i]:
           pass
        else:
           print('Erro: Matriz de rigidez assimétrica')
           break
print('A matriz de rigidez é simétrica') 

### Montando Matriz de Massa

M = []

for i in range(len(var)):
    m = []
    for j in range(len(var)):
        m.append(0)
    M.append(m)
    
print(M)

for i in range(len(var)):
    for j in range(len(var)):
       
       M[i][j] = graus[i].coeff(var_dot[j],1)
print(M)
for i in range(len(var)):
    for j in range(len(var)):
        if M[i][j] == M[j][i]:
           pass
        else:
           print('Erro: Matriz de massa assimétrica')
           break
print('A matriz de massa é simétrica')

##### Montando a Matriz de amortecimento 

C = []

for i in range(len(var)):
    c = []
    for j in range(len(var)):
        c.append(0)
    C.append(c)
    
print(C)

for i in range(len(var)):
    for j in range(len(var)):
       
       C[i][j] = graus[i].coeff(var_dot[j],1)
print(C)
for i in range(len(var)):
    for j in range(len(var)):
        if C[i][j] == C[j][i]:
           pass
        else:
           print('Erro: Matriz de amortecimento é assimétrica')
           break
print('A matriz de amortecimento é simétrica')


##### Analise Modal 

##### autovalores e autovetores

for i in range(len(var)):
    for j in range(len(var)):
        K[i][j] = float(K[i][j])
        M[i][j] = float(M[i][j])
        
M = np.array(M)
K = np.array(K)        
[wn, B] = eigh(K,M)
wn = np.sqrt(wn**2)
print(wn)
print(B) 

freq_n = np.sqrt(((wn)**2)**(1/2))
print(freq_n)

meff = []
for i in range(len(var)):
    meff.append(0)
    meff[i] = np.dot(np.dot(B[:,i].T,M),B[:,i])
    
meff = np.diag(meff)
for i in range(len(var)):
    for j in range(len(var)):
        if meff[i][j] == meff[j][i]:
           pass
        else:
           print('Erro: Matriz de rigidez assimétrica')
           break
print('A matriz é simétrica')

 #### Condições iniciais
# [roda2 roda3 roda4 roda5 chassi1 rolagemX6 rolagemY7]
#x0 = [0,0.5,0,0,0,0,0]
#x0 = [0,0.3,0.3,0,0,0,0]
#x0 = [0,0,0.2,0,0,0.1,0]
#0 = [0,0,0.3,0,0,0,0]
x0 = [0.3,0.3,0.3,0.3,0,0,0]
v0 = [0,0,0,0,0,0,0]
y0 = np.dot(np.dot(B.T,M),x0)
y0dot = np.dot(np.dot(B.T,M),v0)

t_inc0 = 0.005 #segundos
t_end0 = 30#segundos    
y = []
t = np.arange(0, t_end0, t_inc0)
for r in range(len(var)):
    y.append(0)
    y[r] = y0[r]*np.cos(freq_n[r]*t) + y0dot[r]/(freq_n[r]*np.sin(freq_n[r]*t))
    
x = np.dot(B,y)
x1 = x[0][:]
x2 = x[1][:]
x3 = x[2][:]
x4 = x[3][:]
x5 = x[4][:]
x6 = x[5][:]
x7 = x[6][:]

np.savetxt('7GDL_Resposta4', np.transpose([x1,x2,x3,x4,x5,x6,x7]), fmt='%1.5f')

#print(x)

### PLOTAGENS BABILÔNICAS

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


### fixando eixos do gráfico
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')
ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')
ax.set_zlim3d([-1.0, 1.0])
ax.set_zlabel('Z')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_ = 0.5 * np.outer(np.cos(u), np.sin(v)) 
y = 0.5 * np.outer(np.sin(u), np.sin(v))
z = 0.5 * np.outer(np.ones(np.size(u)), np.cos(v)) + x[0][0]
surf = ax.plot_surface(x_, y, z)

def animate(frame):
    # Make data
    z = 0.5 * np.outer(np.ones(np.size(u)), np.cos(v)) + x[0][frame]
    # Plot the surface
    surf(x_,y,z)
    
    
    
anim = FuncAnimation(fig, animate, frames=400, interval=20)
