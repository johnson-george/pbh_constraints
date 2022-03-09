import math
import numpy as np
import sys

n = int(sys.argv[1])
s = float(sys.argv[2])
if s != 0.5: s = int(s)
rh = 1

def temp(rh, n):
	return (n+1)/(4*math.pi*rh)
T = temp(rh, n)

def dist(E, T, s):
	zeta = int(2*s)
	return 1/(np.exp(E/T) - (-1)**zeta)

xs = np.power(10, np.linspace(np.log10(0.01),np.log10(5),200))
a = [0.00000e+000, 1.00000e-003, 1.93070e-003, 3.72759e-003, 7.19686e-003, 1.38950e-002, 2.68270e-002, 5.17947e-002, 1.00000e-001, 1.36364e-001, 1.72727e-001, 2.09091e-001, 2.45455e-001, 2.81818e-001, 3.18182e-001, 3.54545e-001, 3.90909e-001, 4.27273e-001, 4.63636e-001, 5.00000e-001, 5.36364e-001, 5.72727e-001, 6.09091e-001, 6.45455e-001, 6.81818e-001, 7.18182e-001, 7.54545e-001, 7.90909e-001, 8.27273e-001, 8.63636e-001, 9.00000e-001, 9.30481e-001, 9.51671e-001, 9.66402e-001, 9.76643e-001, 9.83762e-001, 9.88712e-001, 9.92152e-001, 9.94544e-001, 9.96207e-001, 9.97363e-001, 9.98167e-001, 9.98726e-001, 9.99114e-001, 9.99384e-001, 9.99572e-001, 9.99702e-001, 9.99793e-001, 9.99856e-001, 9.99900e-001]

coeffs = [] # this will be a list of 200 absorption coefficients
infile = "gb" + str(n) + str(s) + ".txt" # taking the files produced by gb.py
with open(infile) as g:
	for line in g:
		out = line[0:-1]
		coeffs.append(float(out))
gammas = [coeffs[i]*dist(xs[i], T, s) for i in range(len(xs))] # turning coeffs into gammas

def eformat(f, prec, exp_digits): # this function just gets the scientific notation correct
    s = "%.*e"%(prec, f)
    mantissa, exp = s.split('e')
    return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))

###########################################################
# producing the files for greybody data between x = 0.01 and x = 5
###########################################################

pad = 15
filename = str(n) + "spin_" + str(s) + ".txt"
with open(filename, "w") as f:
	f.write("a/x".rjust(pad))
	for x in xs:
		y = eformat(x, 5, 3)
		f.write(y.rjust(pad))
	f.write("\n")
	for m in range(50):
		b = eformat(a[m], 5, 3)
		f.write(b.rjust(pad))
		for g in gammas:
			h = eformat(g, 5, 3)
			f.write(h.rjust(pad))
		f.write("\n")

###########################################################
# producing the fit files
###########################################################

# limiting cross sections for brane emission
rc = np.power((n+3)/2, 1/(n+1))*np.sqrt((n+3)/(n+1))
sig = math.pi*rc**2 # the high energy asymptotic value of the greybody factor (all spins)
sigs = 4*math.pi # the low energy scalar greybody factor
sigf = math.pi*np.power(2, (3*n-1)/(n+1)) # the low energy fermion greybody factor
# below are the constants in sig = sigg*x^p for spin 1 particles
sigg = 16*math.pi*(math.gamma(1/(n+1))*math.gamma(2/(n+1))/math.gamma(3/(n+1)))**2/(3*(n+1)**2)
p = 2

# higher dimensional cross sections --- sgn is cross section, gn is absorption coefficient (sans om) 
sgn = 2*math.pi**(n/2+1)*((n+3)/2)**((n+2)/(n+1))*((n+3)/(n+1))**(n/2+1)/((n+2)*math.gamma(n/2+1))
gn = sgn/(2**n*math.pi**((n+1)/2)*math.gamma((n+1)/2)*(n+1)) # = sig/pi for n = 0

# the seven parameter fit lists
highfits = [-np.log10(np.exp(1))/T, np.log10(sig/math.pi), 0, 0, 2]
lowfits0 = [2, np.log10(sigs/math.pi)] 
lowfits05 = [2, np.log10(sigf/math.pi)] 
lowfits1 = [2+p, np.log10(sigg/math.pi)] 
lowfits = [lowfits0, lowfits05, lowfits1]
highfits2 = [-np.log10(np.exp(1))/T, np.log10(gn), 0, 0, n+2]
lowfits2 = [n+6,0]

if s < 1.5:
	fits = lowfits[int(2*s)] + highfits
else:
	fits = lowfits2 + highfits2

heads = ["a1", "b1", "a2", "b2", "c2", "d2", "e2"]
filename2 = str(n) + "spin_" + str(s) + "_fits.txt"
with open(filename2, "w") as f:
	f.write("a/fits".rjust(pad))
	for head in heads:
		f.write(head.rjust(pad))
	f.write("\n")
	for m in range(50):
		b = eformat(a[m], 5, 3)
		f.write(b.rjust(pad))
		for fit in fits:
			h = eformat(fit, 5, 3)
			f.write(h.rjust(pad))
		f.write("\n")
