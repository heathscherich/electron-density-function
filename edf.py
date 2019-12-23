from math import factorial, sqrt, radians, sin, cos, pi
from cmath import exp
import pdb

import numpy as np
import matplotlib.pyplot as plt

bohr_radius = 1	# in terms of bohr radii
j = complex(0, 1)

# Generalized Laguerre polynomials using the recursion formula
# https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
def Laguerre_of_rho(alpha, k):
	L = []
	L.append(1)
	L.append(1 + alpha - rho)
	
	for i in range(k):
		# k = i+1
		next = (2*(i+1) + 1 + alpha + rho)*L[i+1] - (i+1 + alpha)*L[i]
		next = next/(i+1 + 1)
		L.append(next)
		
	return L[-1]
		
# Spherical harmonics from the following table
# https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
def spherical_harmonic(l, m, theta, phi):
	if l == 0 and m == 0:
		return .5*sqrt(1/pi)
		
	if l == 1:
		if m == -1:
			return .5*sqrt(3/(2*pi))*exp(-j*phi)*sin(theta)
		if m == 0:
			return .5*sqrt(3/pi)*cos(theta)
		if m == 1:
			return -.5*sqrt(3/(2*pi))*exp(j*phi)*sin(theta)
			
	if l == 2:
		if m == -2:
			return .25*sqrt(15/(2*pi))*exp(-2*j*phi)*sin(theta)**2
		if m == -1:
			return .5*sqrt(15/(2*pi))*exp(-j*phi)*sin(theta)*cos(theta)
		if m == 0:
			return .25*sqrt(5/pi)*(3*cos(theta)**2 - 1)
		if m == 1:
			return -.5*sqrt(15/(2*pi))*exp(j*phi)*sin(theta)*cos(theta)
		if m == 2:
			return .25*sqrt(15/(2*pi))*exp(2*j*phi)*sin(theta)**2

num_r = 25
num_t = 100
num_p = 3
n_range = range(1, 4)
r_space = list(range(0, num_r))
r_space[:] = [ r*10*bohr_radius/num_r for r in r_space ]	
theta_space = list(range(0, num_t))
theta_space[:] = [ t*360/num_t for t in theta_space ]
phi_space = list(range(0, num_p))
phi_space[:] = [ p*180/num_p for p in phi_space ]

wavefunction = []
for n in n_range:
	wavefunction.append([])
	for l in range(n):
		wavefunction[n-1].append([])
		for m in range(-l, l+1):
			wavefunction[n-1][l].append([])
			for r in range(len(r_space)):
				r = r_space[r]
				rho = 2*r/(n*bohr_radius)
				
				for theta in range(len(theta_space)):
					theta = radians(theta_space[theta])
					for phi in range(len(phi_space)):
						# Equation from the solution of the Schrodinger Equation
						# https://en.wikipedia.org/wiki/Hydrogen_atom#Wavefunction
						phi = radians(phi_space[phi])
						psi = sqrt(((2/(n*bohr_radius))**3)*(factorial(n-l-1)/(2*n*factorial(n+l))))
						psi = psi*exp(-rho/2)*rho**l
						psi = psi*Laguerre_of_rho(2*l+1, n-l-1)
						psi = psi*spherical_harmonic(l, m, theta, phi)
			
						wavefunction[n-1][l][m+l].append([-psi, r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)])
						
# Graphs each electron density function in it's own window
xs = []
ys = []
max_psi = []
i=1
for n in n_range:
	for l in range(n):
		for m in range(-l, l+1):
			xs.append([])
			ys.append([])
			max_psi.append(0)
			
			for b in range(len(wavefunction[n-1][l][m+l])):
				wf_mag = sqrt(wavefunction[n-1][l][m+l][b][0].real**2 + wavefunction[n-1][l][m+l][b][0].imag**2)
				if wf_mag > max_psi[i-1]:
					max_psi[i-1] = wf_mag
				xs[i-1].append(wavefunction[n-1][l][m+l][b][1])
				ys[i-1].append(wavefunction[n-1][l][m+l][b][2])
				
			rgba = np.zeros((25, len(xs[i-1]), 4))
			rgba[i, :, 0] = 1.0
			print(len(xs[i-1]))
			for a in range(len(xs[i-1])):
				op = sqrt(wavefunction[n-1][l][m+l][a][0].real**2 + wavefunction[n-1][l][m+l][a][0].imag**2)/max_psi[i-1]
				rgba[i, a, 3] = op
		
			plt.figure()
			plt.suptitle("n: " + str(n) + " l: " + str(l) + " m: " + str(m))
			plt.scatter(xs[i-1], ys[i-1], color=rgba[i])
			i+=1
plt.show()