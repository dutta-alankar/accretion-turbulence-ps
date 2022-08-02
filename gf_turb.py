# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 18:22:29 2022

@author: Alankar
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# constants 
mu = 0.61
mp = 1.676e-24
Msun = 2e33
yr = 365*24*60**2
pc = 3.0857e18
kB = 1.38e-16
G  = 6.7e-8 
km = 1e5
s  = 1

kpc = 1e3*pc
mu = 1.0
gamma = 0.088*(1.5/1.5)**-3 #5/3.

UNIT_DENSITY  = mp
UNIT_LENGTH   = pc
UNIT_VELOCITY = km/s
UNIT_TIME     = UNIT_LENGTH/UNIT_VELOCITY

#parameters
R = 40*kpc
H = 330*pc
ndens1 = 0.1
ndens2 = 5.5e-2*ndens1
rho1 = ndens1*mu*mp
rho2 = ndens2*mu*mp
alpha1 = ndens1/(ndens1+ndens2)
alpha2 = 1-alpha1
vrel = 40*km/s
vrot = 158*km/s
B0 = 1e-6 #Gauss
kmag = 2*np.pi/(200*pc)
Sigma = 1*Msun*(pc**-2)
g_rot = vrot**2/R * (H/R)
g_dis = 2*np.pi*Sigma*G
g_tot = g_rot+g_dis
Mdot  = 2*np.pi * 2*H * R * rho2 * vrel
print("Mdot = %.2f Msun/yr"%(Mdot/(Msun/yr)))

def ns(k):
    kh   = alpha1*alpha2*(k*vrel)**2
    mag  = (B0**2/(2*np.pi*(rho1+rho2))) * (k**3/kmag)
    grav = g_tot*k*(alpha1-alpha2)
    val  = kh - mag - grav
    #val  = np.piecewise(val, [val>0,], [lambda x:x, 0.]) 
    #if (val<0): print('Problem!')
    return np.sqrt(val)

def ns2(k):
    kh   = alpha1*alpha2*(k*vrel)**2
    mag  = (B0**2/(2*np.pi*(rho1+rho2))) * (k**3/kmag)
    grav = g_tot*k*(alpha1-alpha2)
    val  = kh - mag - grav
    return val

def ns_derv(k):
    kh   = 2*alpha1*alpha2*k*vrel**2
    mag  = 3*(B0**2/(2*np.pi*(rho1+rho2))) * (k**2/kmag)
    grav = g_tot*k*(alpha1-alpha2)
    val  = kh - mag - grav
    val  = val/(2*ns(k))
    return val

def nc_derv(k, nc): 
    num   = (2-np.sqrt(gamma*nc/ns(k)))*ns_derv(k) + (2*gamma/k)*nc
    denom = gamma * (3-np.sqrt(ns(k)/(gamma*nc))) 
    return (num/denom)

def roots():
    A = 3*B0**2/(4*np.pi*(rho1+rho2)*kmag)
    B = -alpha1*alpha2*vrel**2
    C = 0.5*g_tot*(alpha1-alpha2)
    discrim = np.sqrt(B**2-4*A*C)
    return np.array([-B+discrim, -B-discrim])/(2*A)

k_start = 2*np.pi/(50*kpc)
k_stop  = 2*np.pi/(200*pc)
k_vals = np.logspace(np.log10(k_start), np.log10(k_stop), 100)
plt.plot((2*np.pi/k_vals)/kpc, ns_derv(k_vals)-2*ns(k_vals)/k_vals)
plt.xlabel(r'$L=2\pi/k \ [kpc]$')
plt.ylabel(r'$\frac{d}{dk} \left(\frac{n_s}{k^2}\right)$')
plt.grid()
plt.savefig('plots/L0-driving.png', transparent=True)
plt.show()
plt.close()

k_start = 2*np.pi/(500*pc)
k_stop  = 2*np.pi/(10*pc)
k_vals = np.logspace(np.log10(k_start), np.log10(k_stop), 100)
plt.plot((2*np.pi/k_vals)/pc, ns2(k_vals), color='tab:blue')
plt.plot((2*np.pi/k_vals)/pc, -ns2(k_vals), color='tab:blue')
plt.xlabel(r'$L=2\pi/k \ [pc]$')
plt.ylabel(r'$n_s(k)$')
plt.yscale('log')
plt.grid()
plt.savefig('plots/Lmin-dissp.png', transparent=True)
plt.show()
plt.close()

#k0 and kf from previous plots
k0 = 2*np.pi/(9.0*kpc)
kf = 2*np.pi/(240*pc)

k_vals = np.logspace(np.log10(k0), np.log10(kf), 1000)
plt.plot(k_vals/k0, ns(k_vals)/ns(k0))   
plt.xlabel(r'$k/k_0$')
plt.ylabel(r'$n_s(k)/n_s(k_0)$')
plt.yscale('log')
plt.xscale('log')
plt.xticks([1,2,5,10,20,50],[1,2,5,10,20,50])
plt.yticks([0.5,1.0,5,10,20],[0.5,1.0,5,10,20])
plt.grid()
plt.savefig('plots/turb-growth-rate.png', transparent=True)
plt.show()
plt.close()

eps = 1e-8*k0
init = ns(k0)/gamma
k_evals = k_vals = np.logspace(np.log10(k0+eps), np.log10(kf-eps), 1000)
res  = solve_ivp(nc_derv, np.array([k0, kf]), 
                 np.array([init,]), t_eval=k_evals)
nc = res.y.flatten()
pspec = -(np.sqrt(gamma*nc) + \
          np.sqrt(ns(k_evals)))**2*\
          np.gradient(nc/k_evals**2, k_evals)
          
plt.plot(k_evals/k0, pspec*((ns(k0)/k0)**-2)*k0*gamma )   
plt.xlabel(r'$k/k_0$')
plt.ylabel(r'$F(k)\times (n_s(k_0)/k_0)^{-2}k_0 \gamma$')
plt.yscale('log')
plt.xscale('log')
plt.xticks([1,2,5,10,20,50],[1,2,5,10,20,50])
plt.yticks([0.05,0.1,0.5,1.0,5,10],[0.05,0.1,0.5,1.0,5,10])
plt.grid()
plt.savefig('plots/power-spec-norm.png', transparent=True)
plt.show()    
plt.close()

plt.plot((2*np.pi/k_evals)/kpc, pspec*((ns(k0)/k0)**-2)*k0*gamma )   
plt.xlabel(r'$L\ [kpc]$')
plt.ylabel(r'$F(k)\times (n_s(k_0)/k_0)^{-2}k_0 \gamma$')
plt.yscale('log')
plt.xscale('log')
plt.xticks([0.5,1,2,5,10],[0.5,1,2,5,10])
plt.yticks([0.05,0.1,0.5,1.0,5,10],[0.05,0.1,0.5,1.0,5,10])
plt.grid()
plt.show()    
