import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import sys

n = int(sys.argv[1])
s = 1
rh = 1

def f(r,n): 
	return 1 - (rh/r)**(n+1)

def deriv(p,r,om,l,n):
	ReR, ImR, ReS, ImS = p
	h = f(r,n)
	ang = l*(l+1)/(h*r**2)
	der = 2/r
	img = 2*om/h
	return [ReS, ImS, ang*ReR - der*ReS - img*ImS, ang*ImR - der*ImS + img*ReS]

###########################################################
# initialising the solver
###########################################################
steps = 10**6
A = 1
eps = 10**(-5)
rpoints = np.linspace(1+eps, 10**4, steps)

def grey(om,l,n):
	###########################################################
	# solving the ODE
	###########################################################
	omeg = l*(l+1)/(2*om)
	p0 = [A,0,0,A*omeg]
	sol = odeint(deriv, p0, rpoints, args=(om,l,n))

	###########################################################
	# extracting the reflection coefficient
	###########################################################
	ReF = sol[-1,0]
	ImF = sol[-1,1]
	B = (ReF**2 + ImF**2)

	return (2*l+1)*A**2/B # returns the absorption coefficient, not greybody factor

#n = 0
#om = 0.36
#l = 1
#omeg = l*(l+1)/(2*om)
#p0 = [A,0,0,A*omeg]
#sol = odeint(deriv, p0, rpoints, args=(om,l,n))
#plt.plot(rpoints, sol[:,0])
#plt.plot(rpoints, sol[:,1])
#plt.show()
#quit()

lmax = 30
lmin = 1
xpoints, xmin, xmax = 200, 0.01, 5
gs = []
gbs = []
xs = np.power(10, np.linspace(np.log10(xmin),np.log10(xmax),xpoints))
count = 0
outfile = "gb" + str(n) + str(s) + ".txt"
with open(outfile,'w') as myfile:
	for om in xs:
		count += 1
		g = 0
		for l in range(lmin,lmax):
			print(str(count)+"/"+str(xpoints)+' '+str(l) , end='\r', flush=True)
			g += grey(om,l,n)
		gs.append(g)
		gbs.append(g/(rh**2*om**2))
		myfile.write(str(g)+'\n')

#gbs = []
#with open(outfile, 'r') as myfile:
#	for line, om in zip(myfile, xs):
#		gb = float(line[0:-1])
#		gbs.append(gbb/(rh**2*om**2))

plt.plot(xs, gbs)
plt.show()
