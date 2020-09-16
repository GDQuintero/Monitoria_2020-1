
import numpy as np
import matplotlib.pyplot as plt

N = int(input("Escreve N: "))
metodo = int(input("Qual metodo deseja usar (1: Euler implicito, 2: Crank-Nikolson)? "))
item = int(input("Qual teste deseja fazer (1: Item a errado, 2: Item a certo, 3: Item b ou 4: Item c) "))
p = float(input("Digite p: "))
T = 1.0; dx = 1.0 / N; dt = dx; lam = N; M = N

barra = np.linspace(0.0, 1.0, N+1); t = np.linspace(0.0, 1.0, M+1)

uold = np.zeros(N+1); A = np.zeros((11,N+1)); aux = np.zeros(N+1)

Adiag = np.zeros(N-1); Asub = np.zeros(N-2); b = np.zeros(N-1)

A = np.zeros((11,N+1)); z = np.zeros(N-1); y = np.zeros(N-1); x = np.zeros(N-1)

graf = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

"============================================================="
"DEFINICAO DAS FUNCOES A SEREM USADAS"
"============================================================="

def r(t):
	return 10000.0 * (1.0 - 2.0*(t**2))

def g(x,aux,h,p):

	for i in range(len(x)):
		if ((p - (0.5 * h)) <= x[i]) & (x[i] <= (p + (0.5 * h))):
			aux[i] = 1.0 / h		
	else:
		aux[i] = 0.0
	return aux

def f(t,x,item):
	if item == 1:
		return 10.0*(x**2)*(x-1) - 60.0*x*t + 20.0*t
	elif item == 2:
		return (10.0 * np.cos(10.0 * t) * (x**2) * ((1.0 - x)**2)) - ((1.0 + np.sin(10.0 * t)) * ((12.0 * (x**2)) - 12.0*x + 2.0))
	elif item == 3:
		return np.exp(t - x) * ((25.0 * (t**2) * np.cos(5.0 * t * x)) - (10.0 * t * np.sin(5.0 * t * x)) - (5.0 * x * np.sin(5.0 * t * x)))
	else:
		return r(t) * g(x, aux, dx, p)

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

def fatoracao_ldl(Adiag,Asub,n):
	D = np.zeros(n); L = np.zeros(n-1)
	D[0] = Adiag[0]; L[0] = Asub[0] / D[0]

	for i in range(1,n-1):
		D[i] = Adiag[i] - D[i-1] * (L[i-1]**2)
		L[i] = Asub[i] / D[i]

	D[n-1] = Adiag[n-1] - D[n-2] * (L[n-2]**2)

	return L, D

def solver(L,D,b,x,y,z,n):
	
	z[0] = b[0]

	for i in range(1,n):
		z[i] = b[i] - L[i-1] * z[i-1]

	for i in range(n):
		y[i] = z[i] / D[i]

	x[n-1] = y[n-1]

	for i in range(n-2,-1,-1):
		x[i] = y[i] - L[i] * x[i+1]
		
	return x

"============================================================="
"SOLUCAO"
"============================================================="

if metodo == 1:
	Adiag[:] = 1.0 + (2.0 * lam) 
	Asub[:] = -lam
else:
	Adiag[:] = 1.0 + lam 
	Asub[:] = -0.5 * lam

L, D = fatoracao_ldl(Adiag,Asub,N-1)

#print(D)
#print()
#print(L)

uold[:] = u0(barra,item); A[0,:] = uold.copy(); ind = 1

for k in range(M):
	uold[0] = g1(t[k],item); uold[N] = g2(t[k],item)

	if metodo == 1:
		fk = f(t[k+1], barra, item)
		b[0] = uold[1] + dt * fk[1] + (lam * g1(t[k+1],item))
		b[1:-1] = uold[2:-2] + dt*fk[2:-2]
		b[-1] = uold[-2] + dt * fk[-2] + (lam * g2(t[k+1],item))
		
	else: 
		fk = f(t[k], barra, item); fkk = f(t[k+1], barra, item)
		b[0] = (0.5 * lam * (g1(t[k],item) + g1(t[k+1],item) + uold[2] - 2.0*uold[1])) + uold[1] + ((0.5 * dt) * (fk[1] + fkk[1]))
		b[1:-1] = (0.5 * lam * (uold[1:-3] -2.0*uold[2:-2] + uold[3:-1])) + uold[2:-2] + (0.5 * dt * (fk[2:-2] + fkk[2:-2]))
		b[N-2] = (0.5 * lam * (g2(t[k],item) + g2(t[k+1],item) + uold[-3] - 2.0*uold[-2])) + uold[-2] + ((0.5 * dt) * (fk[-2] + fkk[-2]))

	unew = solver(L,D,b,x,y,z,N-1)
	uold[1:-1] = unew.copy()
	x[:] = 0.0; y[:] = 0.0; z[:] = 0.0

	if k == ind*int((M)/10):
		A[ind,:] = uold.copy()
		ind = ind + 1

	if k == M-1:
		A[10,:] = uold.copy()

for i in range(11):
	tk = 0.1*i
	plt.plot(barra,A[i,:],label='t = ' + str(graf[i]))

if item == 4:
	exata = uold
else:
	exata = sol(t[-1], barra, item)

erro = np.max(np.abs(uold - exata))
plt.plot(barra,A[10,:],label='t = ' + str(graf[10]))
plt.xlabel('Barra')
plt.ylabel('Temperatura')
plt.legend()
plt.show()

#print(D)
#print(erro)





