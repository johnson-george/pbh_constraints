import numpy as np
import math
import matplotlib.pyplot as plt
import sys

# Once the greybody tables are produced, this script checks that they give me what I expect. I just need to change the value of s and n immediately below, so it reads in the correct files. Checking this for all values of n and s will confirm these can be safely fed into primary.c.

n = int(sys.argv[1])
s = float(sys.argv[2])
if s != 0.5: s = int(s)
rh = 1

def temp(rh, n):
	return (n+1)/(4*math.pi*rh)
T = temp(rh, n)

filename = str(n) + "spin_" + str(s) + ".txt"
with open(filename, "r") as f:
	lines = f.readlines()
	d = lines[1]
	da = d.split()
	dat = [float(x) for x in da] 
data = dat[1:]

filename = str(n) + "spin_" + str(s) + "_fits.txt"
with open(filename, "r") as g:
	lines = g.readlines()
	e = lines[1]
	ea = e.split()
	a = [float(x) for x in ea] 

def distinv(E, T, s):
	zeta = int(2*s)
	return np.exp(E/T) - (-1)**zeta

def cross(g, om, n):
	if s < 1.5:
		return math.pi*g/om**2
	else:
		return 2**n*math.pi**((n+1)/2)*math.gamma((n+1)/2)*(n+1)*g/om**(n+2)

xs = np.power(10, np.linspace(np.log10(0.01),np.log10(5),200))
gb = [cross(data[i], xs[i], n)*distinv(xs[i],T,s) for i in range(len(xs))]
plt.plot(xs, gb)

#oga = [0, 9.92034e-001, -5.20053e-001, -5.13437e+000, 7.45615e-001, -1.40633e-002, -3.06190e-002,  5.42759e+000] # the scalar fit parameters from the original code
#mya = [0, 1.00000e+000, -4.97150e-001, -5.45751e+000, 8.29304e-001, 0.00000e+000, 0.00000e+000, 2.00000e+000] # my improved fit parameters

def lowQ(x, a):
	loglowQ = a[1]*np.log10(x) + a[2]
	loglowQ -= np.log10(distinv(x,T,s))
	return np.power(10, loglowQ)

def highQ(x, a):
	loghighQ = a[3]*x + a[4] + a[5]*np.cos(a[7]*x) + a[6]*np.sin(a[7]*x)
	loghighQ += a[7]*np.log10(x) 
	return np.power(10, loghighQ)

xslow = np.power(10, np.linspace(np.log10(0.001),np.log10(0.01),200))
xshigh = np.power(10, np.linspace(np.log10(5),np.log10(10),200))
gblow = [cross(lowQ(x, a),x,n)*distinv(x,T,s) for x in xslow]
gbhigh = [cross(highQ(x, a),x,n)*distinv(x,T,s) for x in xshigh]
plt.plot(xslow, gblow)
plt.plot(xshigh, gbhigh)
plt.ylim([0, 40])

plt.show()

#xs2 = [1.00000e-002,  1.03172e-002,   1.06445e-002,   1.09822e-002,   1.13305e-002,   1.16900e-002,   1.20608e-002,   1.24434e-002,   1.28381e-002,   1.32454e-002,   1.36655e-002,   1.40990e-002,   1.45463e-002,   1.50077e-002,   1.54838e-002,   1.59750e-002,   1.64817e-002,   1.70046e-002,   1.75440e-002,   1.81005e-002,   1.86747e-002,   1.92671e-002,   1.98783e-002,   2.05088e-002,   2.11594e-002,   2.18306e-002,   2.25232e-002,   2.32376e-002,   2.39748e-002,   2.47353e-002,   2.55200e-002,   2.63295e-002,   2.71647e-002,   2.80264e-002,   2.89155e-002,   2.98327e-002,   3.07791e-002,   3.17555e-002,   3.27628e-002,   3.38021e-002,   3.48744e-002,   3.59807e-002,   3.71220e-002,   3.82996e-002,   3.95146e-002,   4.07680e-002,   4.20613e-002,   4.33955e-002,   4.47721e-002,   4.61924e-002,   4.76577e-002,   4.91695e-002,   5.07293e-002,   5.23385e-002,   5.39988e-002,   5.57117e-002,   5.74790e-002,   5.93023e-002,   6.11835e-002,   6.31244e-002,   6.51268e-002,   6.71928e-002,   6.93242e-002,   7.15233e-002,   7.37922e-002,   7.61330e-002,   7.85481e-002,   8.10398e-002,   8.36106e-002,   8.62628e-002,   8.89993e-002,   9.18225e-002,   9.47353e-002,   9.77405e-002,   1.00841e-001,   1.04040e-001,   1.07340e-001,   1.10745e-001,   1.14258e-001,   1.17883e-001,   1.21622e-001,   1.25480e-001,   1.29461e-001,   1.33568e-001,   1.37805e-001,   1.42176e-001,   1.46686e-001,   1.51339e-001,   1.56140e-001,   1.61093e-001,   1.66203e-001,   1.71476e-001,   1.76915e-001,   1.82527e-001,   1.88317e-001,   1.94291e-001,   2.00454e-001,   2.06813e-001,   2.13374e-001,   2.20142e-001,   2.27126e-001,   2.34331e-001,   2.41764e-001,   2.49433e-001,   2.57346e-001,   2.65509e-001,   2.73932e-001,   2.82621e-001,   2.91587e-001,   3.00836e-001,   3.10379e-001,   3.20225e-001,   3.30383e-001,   3.40864e-001,   3.51677e-001,   3.62833e-001,   3.74342e-001,   3.86217e-001,   3.98469e-001,   4.11109e-001,   4.24150e-001,   4.37605e-001,   4.51487e-001,   4.65809e-001,   4.80585e-001,   4.95830e-001,   5.11559e-001,   5.27786e-001,   5.44529e-001,   5.61802e-001,   5.79624e-001,   5.98011e-001,   6.16981e-001,   6.36552e-001,   6.56745e-001,   6.77578e-001,   6.99072e-001,   7.21248e-001,   7.44128e-001,   7.67733e-001,   7.92087e-001,   8.17213e-001,   8.43137e-001,   8.69883e-001,   8.97477e-001,   9.25947e-001,   9.55320e-001,   9.85625e-001,   1.01689e+000,   1.04915e+000,   1.08243e+000,   1.11677e+000,   1.15219e+000,   1.18874e+000,   1.22645e+000,   1.26536e+000,   1.30550e+000,   1.34691e+000,   1.38964e+000,   1.43372e+000,   1.47920e+000,   1.52612e+000,   1.57453e+000,   1.62448e+000,   1.67601e+000,   1.72918e+000,   1.78403e+000,   1.84062e+000,   1.89901e+000,   1.95925e+000,   2.02140e+000,   2.08553e+000,   2.15168e+000,   2.21994e+000,   2.29036e+000,   2.36301e+000,   2.43797e+000,   2.51531e+000,   2.59510e+000,   2.67742e+000,   2.76235e+000,   2.84998e+000,   2.94039e+000,   3.03366e+000,   3.12990e+000,   3.22918e+000,   3.33162e+000,   3.43730e+000,   3.54634e+000,   3.65884e+000,   3.77491e+000,   3.89465e+000,   4.01820e+000,   4.14566e+000,   4.27717e+000,   4.41285e+000,   4.55284e+000,   4.69726e+000,   4.84627e+000,   5.00000e+000]
