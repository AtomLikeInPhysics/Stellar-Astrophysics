import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

G = 6.67e-11
rho0 = 1e5
Ohm = 3.5e-8
C = 1e-11
Msun = 1e30
Rsun = 7e8




class Star(object):
    def __init__(self,M,R):
        self.M = M
        self.R = R

    def __call__(self,init,t):
        def rho(self,r):
            return(rho0 * (1 - (r/self.R)**2))
        def m(self,r):
            return(4/3 * np.pi * r**3 * rho(r))
        def dPdr(self,r):
            return(-G * m(r)*rho(r) / r**2)
        radius = init[0]
        velocity = init[1]
        d_radius = velocity
        d_velocity = (-G*self.M/radius**2) + (-G * 4/3 * np.pi * radius**3 * (rho0 * (1 - (radius/self.R)**2))**2 / radius**2) + C*(np.cos(Ohm*t))
        return np.array([d_radius, d_velocity], float)

    
star = Star(1.18*Msun,300*Rsun)

tstart = 0
tstop = 10000
N = 10000
h = (tstop-tstart)/N

xpnt = []
ypnt = []
t = np.linspace(tstart,tstop,N)
init = np.array([300*Rsun*1.1,0],float)

for i in t:
    xpnt.append(init[0])
    ypnt.append(init[1])
    k1 = h*star(init,i)
    k2 = h*star(init + 0.5*k1, i + 0.5*h)
    k3 = h*star(init + 0.5*k2, i + 0.5*h)
    k4 = h*star(init + k3, i + h)
    init += 1/6 * (k1 + 2*k2 + 2*k3 + k4)

# plt.plot(t,np.array(xpnt))
plt.plot(t,-2.5*np.log10(np.array(xpnt)**2 * 5.67e-8 * 5800**4) + 71.2)
plt.xlim(0,1000)
plt.ylim(-1,30)
# plt.plot(t,ypnt)
plt.ylabel('dRadius')
plt.show()