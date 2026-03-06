import numpy as np
import math as  mt
import random as rm
import pandas as pd
from numpy import linalg as LA
import matplotlib.pyplot as plt 
import math as  mt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np
import matplotlib.pyplot as plt
import math as mt

def repartitionDensite(X,Nmc):
	a = -0.75
	b = 0.75
	Nx = 100
	delta = (b-a)/Nx 
	x = [0]*Nx
	proba = [0]*Nx
	densite = [0]*Nx
	for i in range(Nx):
		x[i]=a+delta*i
		compteur = 0
		for j in range (Nmc):
			if (X[j]>=x[i] and X[j] <= x[i]+delta):
				compteur += 1
		proba[i] = compteur/Nmc
		densite[i]=proba[i]/delta
	return [x,densite]

def prixHeston(N, k, eta, S0,theta,rho, r,v0, reduction = False,N1 = [], N2 = []):
	T = 0.5
	deltaT = T / (N+1)
	V = [v0 for i in range(N + 1)]
	S = [S0 for i in range(N + 1)]
	Vs = [v0 for i in range(N + 1)]
	Ss = [S0 for i in range(N + 1)]
	for i in range(N):
		V[i+1] = V[i] + k*(theta-V[i])*deltaT +	eta*mt.sqrt(abs(V[i]))*mt.sqrt(deltaT)*N1[i] +1/4*(eta**2)*deltaT*((N1[i]**2)-1)
		S[i+1] = S[i]*mt.exp((r-V[i]/2)*deltaT +mt.sqrt(abs(V[i]))*(rho*mt.sqrt(deltaT)*N1[i] + mt.sqrt(1-rho**2)*mt.sqrt(deltaT)*N2[i]))
		Vs[i+1] = Vs[i] + k * (theta - Vs[i]) * deltaT + eta*mt.sqrt(abs(Vs[i])) * mt.sqrt(deltaT) * (-N1[i]) + 1/4*(eta**2)*deltaT *((-N1[i]**2) -1)
		Ss[i+1] = Ss[i] * mt.exp( (r- Vs[i]/2) * deltaT +	mt.sqrt(abs(Vs[i])) * ( rho * mt.sqrt(deltaT) *(- N1[i]) + mt.sqrt(1-rho **2) * mt.sqrt(deltaT) * (-N2[i]) ))
	if reduction:
		return V, S, Vs, Ss
	else :
		return V, S

def Echantillons():
	simuS = [[] for i in range(5)]
	simuV = [[] for i in range(5)]
	N1 = [[np.random.randn() for j in range(100)] for i in range(5)]
	N2 = [[np.random.randn() for j in range(100)] for i in range(5)]
	for i in range(5):
		simuV, simuS = prixHeston(100, 2, 0.3, 1,0.04,0.9,0.01,0.04,False,N1[i], N2[i])
		plt.plot(simuV)
	plt.title("Courbes d’évolution de volatilité")
	plt.xlabel("N")
	plt.ylabel("volatilité")
	plt.show()

	for i in range(5):
		simuV, simuS =prixHeston(100, 2, 0.3, 1,0.04,0.9,0.01,0.04,False,N1[i], N2[i])
		plt.plot(simuS)
	plt.title("Courbes d’évolution de l'actif")
	plt.xlabel("N")
	plt.ylabel("Prix de l'actif S")
	plt.show()


def simulation(Nmc, rho):
	simuS = [[] for i in range(Nmc)]
	simuV = [[] for i in range(Nmc)]
	R = [[] for i in range(Nmc)]
	N = 500
	N1 = [[np.random.randn() for j in range(N)] for i in range(Nmc)]
	N2 = [[np.random.randn() for j in range(N)] for i in range(Nmc)]
	for i in range(Nmc):
		simuV, simuS = prixHeston(N, 2, 0.3, 1,0.04,rho,0.01,0.04,False,N1[i], N2[i])
		R[i] = np.log(simuS[N]/1)
	x=repartitionDensite(R,Nmc)
	plt.rcParams["figure.figsize"]=[16,9]
	plt.plot(x[0], x[1], label = f'Densité empirique où rho = {rho}')
	return simuS, simuV, R

#Nmc = 10000
#Echantillons()
#sim1, sim2 , R = simulation(Nmc, 0.9)
#sim1, sim2 , R = simulation(Nmc, 0)
#sim1, sim2 , R = simulation(Nmc, -0.9)
#plt.legend()
#plt.title("Densités empiriques selon la valeur de rho")
#plt.show()

def prixEstimateur(Nmc, N, T, v0,k, eta, S0, K,theta,rho, r,reduction,N1 = [], N2 = []):
	somme = 0
	for i in range(10):
		if reduction == True: 
			V, S, Vs, Ss = prixHeston(N,k, eta, S0, theta,rho, r,v0,reduction,N1[i], N2[i])		
			somme += (max(S[N]-K,0) + max(Ss[N]-K,0)) / 2
		else :
			V, S = prixHeston(N,k, eta, S0, theta,rho, r,v0,reduction,N1[i], N2[i])		
			somme += (max(S[N]-K,0))
	
	return mt.exp(-r*T)*somme/10

def estimateur():
	Ns = 500
	Nmc = 10
	N=100
	estimateur_v, estimateur_v_reduit = [i for i in range(Ns)], [i for i in range(Ns)]
	for i in range(Ns):
		N1 = [[np.random.randn() for j in range(N)] for i in range(Nmc)]
		N2 = [[np.random.randn() for j in range(N)] for i in range(Nmc)]

		estimateur_v[i] = prixEstimateur(Nmc, N,0.5,0.04,2, 0.3,i/Ns*20,10,0.04,0.5,0.01,False, N1, N2) 
		estimateur_v_reduit[i] = prixEstimateur(Nmc, N,0.5,0.04,2, 0.3,i/Ns*20,10,0.04,0.5,0.01,True, N1, N2) 
	t = [i/Ns*20 for i in range(Ns)]
	plt.plot(t,estimateur_v, label="estimateur 1")
	plt.plot(t, estimateur_v_reduit,label="estimateur 2 (réduction de la variance)")
	plt.legend()
	plt.title("Graphe du Call pour différents estimateurs V")
	plt.xlabel("S0")
	plt.ylabel("V(S0,0)")
	plt.show()
#estimateur()	

"""def grecTheta(N, k, eta, S0,theta,rho, r,v0, reduction = False,N1 = None, N2 = None):
	N = 100
	Nmc = 10000
	T = 0.5
	deltaT = T / (N+1)

	if N1 is None:
		N1 = np.random.normal(0, 1, (N, Nmc))
     
	if N2 is None:
		N2 = np.random.normal(0, 1, (N, Nmc))
 
	V = [v0 for i in range(N + 1)]
	S = [S0 for i in range(N + 1)]
	Vs = [v0 for i in range(N + 1)]
	Ss = [S0 for i in range(N + 1)]
	for i in range(N):
		V[i+1] = V[i] + k*(theta-V[i])*deltaT +	eta*mt.sqrt(abs(V[i]))*mt.sqrt(deltaT)*N1[i] +1/4*(eta**2)*deltaT*((N1[i]**2)-1)
		S[i+1] = S[i]*mt.exp((r-V[i]/2)*deltaT +mt.sqrt(abs(V[i]))*(rho*mt.sqrt(deltaT)*N1[i] + mt.sqrt(1-rho**2)*mt.sqrt(deltaT)*N2[i]))
		Vs[i+1] = Vs[i] + k * (theta - Vs[i]) * deltaT + eta*mt.sqrt(abs(Vs[i])) * mt.sqrt(deltaT) * (-N1[i]) + 1/4*(eta**2)*deltaT *((-N1[i]**2) -1)
		Ss[i+1] = Ss[i] * mt.exp( (r- Vs[i]/2) * deltaT +	mt.sqrt(abs(Vs[i])) * ( rho * mt.sqrt(deltaT) *(- N1[i]) + mt.sqrt(1-rho **2) * mt.sqrt(deltaT) * (-N2[i]) ))
	if reduction:
		return V, S, Vs, Ss, N1, N2
	else :
		return V, S, N1, N2"""

Nmc=1000
Nc = 50
N = 100
h = 0.1
K_list = [i*20/Nc for i in range(Nc)]
N1 = np.random.normal(0, 1, (N, Nmc))
N2 = np.random.normal(0, 1, (N, Nmc))


def GrecTheta(Nmc, h, N, T, v0,k ,eta, S0, K,theta,rho, r, reduction,N1, N2):    
	A = prixEstimateur(Nmc, N, T, v0,k, eta, S0, K,theta+h,rho, r,reduction,N1, N2)
	B = prixEstimateur(Nmc, N, T, v0,k, eta, S0, K,theta-h,rho, r,reduction,N1, N2)
	return (A-B)/(2*h)

def GrecEta(Nmc, h, N, T, v0,k ,eta, S0, K,theta,rho, r, reduction,N1, N2):
	A = prixEstimateur(Nmc, N, T, v0,k, eta+h, S0, K,theta,rho, r,reduction,N1, N2)
	B = prixEstimateur(Nmc, N, T, v0,k, eta-h, S0, K,theta,rho, r,reduction,N1, N2)
	return (A-B)/(2*h)

theta = [GrecTheta(Nmc, 0.1, 100, 0.5, 0.04, 3, 0.5, 10, i, 0.2, 0.5, 0.1, True, N1, N2) for i in K_list]
eta = [GrecEta(Nmc, 0.1, 100, 0.5, 0.04, 3, 0.5, 10, i, 0.2, 0.5, 0.1, True, N1, N2) for i in K_list]

# Création des subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 8))  

# Graphique de Theta
axs[0].plot(K_list, theta, label="Theta", color="blue")
axs[0].set_title("Theta en fonction de K")
axs[0].set_xlabel("K")
axs[0].set_ylabel("Theta")
axs[0].grid(False)

# Graphique de Eta
axs[1].plot(K_list, eta, label="Eta", color="orange")
axs[1].set_title("Eta en fonction de K")
axs[1].set_xlabel("K")
axs[1].set_ylabel("Eta")
axs[1].grid(False)

#plt.show()

#Calibration de Heston avec volatilité stochastique
N = 100
Nmc = 1000
def LevenbergMarquart():
	N1 = [[np.random.randn() for j in range(N)] for i in range(Nmc)]
	N2 = [[np.random.randn() for j in range(N)] for i in range(Nmc)]
	Res, Heston = [0 for i in range(21)], [0 for i in range(21)]
	Kp_market = [i*(4/10) + 8 for i in range(21)]
	Vp_market=[2.0944,1.7488,1.4266,1.1456,0.8919,0.7068,0.5641,0.4187,0.3166,0.2425,0.1860,0.1370,0.0967,0.0715,0.0547,0.0381,0.0306,0.0239,0.0163,0.0139,0.086]
	lambdas = 0.01
	tetha0, eta0 = 0.2, 0.5
	eps = 0.0001
	J = np.array([np.zeros(2) for i in range(21)])
	compteur = 0
	d = [tetha0, eta0]
	beta = d
	while (LA.norm(d)>eps):
		for p in range(21):
			Heston[p] = prixEstimateur(Nmc, N,0.5,0.04,3, beta[1],10,Kp_market[p],beta[0],0.5,0.01,True, N1, N2) 
			Res[p]= Vp_market[p]- Heston[p]
			J[p][0] = - GrecTheta(Nmc, 0.1, N, 0.5, 0.04,3,beta[1],10,Kp_market[p],beta[0],0.5, 0.01, True,N1, N2)
			J[p][1] = - GrecEta(Nmc, 0.1, N, 0.5, 0.04,3,beta[1],10,Kp_market[p],beta[0],0.5, 0.01, True,N1, N2)
		d = np.dot(-LA.inv(np.dot(J.transpose(), J) + lambdas*np.identity(2)), np.dot(J.transpose(), Res))
		beta = beta + d.transpose()
		if beta[0] > 1 or beta[0] < 0:
			beta[0] = tetha0
		if beta[1] > 1 or beta[1] < 0:
			beta[1] = eta0
		compteur += 1
		print(compteur, LA.norm(d))
	print('Theta =',beta[0])
	print('Eta =',beta[1])    

	plt.plot(Kp_market,Heston,label = 'V Heston calibré')
	plt.plot(Kp_market,Vp_market, '*', label = 'V marché')
	plt.title("Prix de marché - prix de Heston calibré")
	plt.xlabel("K marché")
	plt.ylabel("Prix de l'option V")
	plt.legend()
	plt.show()
LevenbergMarquart()


