#!/usr/bin/env python
# encoding: utf-8
"""
Creado por: Gustavo Quintero
email: gustavoquintero9105@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt

N = input("Escreve N: ")
item = input("Qual teste deseja fazer (1: Item a errado, 2: Item a certo, 3: Item b ou 4: Item c) ")

dx = 1.0 / N

lam = 0.25; dt = lam*(dx**2); T = 1; M = 1 / dt; M = int(np.ceil(M))

x = np.linspace(0.0, 1.0, N+1); t = np.linspace(0.0, 1.0, M+1)

uold = np.zeros(N+1); A = np.zeros((11,N+1)); aux = np.zeros(N+1)


"============================================================="
"DEFINICAO DAS FUNCOES A SEREM USADAS"
"============================================================="

def r(t):
	return 10000.0 * (1.0 - (2.0 * (t**2)))

def g(x,aux,h,p):

	for i in range(len(x)):
		if ((p - (0.5 * h)) <= x[i]) & (x[i] <= (p + (0.5 * h))):
			aux[i] = 1.0 / h		
	else:
		aux[i] = 0.0
	return aux

def f(t,x,item):
	if item == 1:
		return 10*(x**2)*(x-1) - 60*x*t + 20*t
	elif item == 2:
		return (10.0 * np.cos(10.0 * t) * (x**2) * ((1.0 - x)**2)) - ((1.0 + np.sin(10.0 * t)) * ((12.0 * (x**2)) - 12.0*x + 2.0))
	elif item == 3:
		return np.exp(t - x) * ((25.0 * (t**2) * np.cos(5.0 * t * x)) + (10.0 * t * np.sin(5.0 * t * x)) - (5.0 * x * np.sin(5.0 * t * x)))
	else:
		return r(t) * g(x, aux, dx, 0.25)

def u0(x,item):
	if item == 1:
		return 0.0
	elif item == 2:
		return ((x**2) * ((1.0 - x)**2))
	elif item == 3:
		return np.exp(-x)
	else:
		return 0.0

def g1(t,item):
	if item == 1:
		return 0.0
	elif item == 2:
		return 0.0
	elif item == 3:
		return np.exp(t)
	else:
		return 0.0

def g2(t,item):
	if item == 1:
		return 0.0
	elif item == 2:
		return 0.0
	elif item == 3:
		return np.exp(t - 1.0) * np.cos(5.0 * t)
	else:
		return 0.0

def sol(t,x,item):
	if item == 1:
		return 10.0 * t * (x**2) * (x - 1.0)
	elif item == 2:
		return (1.0 + np.sin(10.0 * t)) * (x**2) * ((1.0 - x)**2)
	else:
		return np.exp(t - x) * np.cos(5.0 * t * x)

def heat(uold, f, dt, dx):
	u = uold.copy()
	u[1:-1] = u[1:-1] + dt*(u[:-2] - 2*u[1:-1] + u[2:])/(dx**2) + dt*f[:]
	return u

"============================================================="
"SOLUCAO"
"============================================================="

uold[:] = u0(x,item); A[0,:] = uold.copy(); ind = 1

for k in range(M):
	faux = f(t[k], x, item)
	uold[0] = g1(t[k], item); uold[N] = g2(t[k], item)
	unew = heat(uold, faux[1:-1], dt, dx)
	uold = unew.copy(); unew[:] = 0.0 

	if k == ind*int((M)/10):
		A[ind,:] = uold.copy()
		ind = ind + 1

	if k == M-1:
		A[10,:] = uold.copy()

for i in range(11):
	tk = 0.1*i
	plt.plot(x,A[i,:],label='t = ' + str(tk))

if item == 4:
	exata = uold
else:
	exata = sol(t[-1], x, item)

erro = np.max(np.abs(uold - exata))

#plt.plot(x,exata,color='black',label='Exata')
plt.xlabel('Barra')
plt.ylabel('Temperatura')
plt.legend()
plt.show()

print(erro)




