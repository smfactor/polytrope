import numpy as np
import matplotlib.pyplot as plt
import polytrope as poly
import scipy.interpolate as intp

data = poly.poly(2.9,0.001,0.00001,'poly2p9')

z=data[0]
th=data[1]
mz2dthdz=data[2]
m3ozdthdz=data[3]

z2p9 = intp.interp1d(z[-2:-1],th[-2:-1])(0.)
print z2p9
