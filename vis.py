import numpy as np
import matplotlib.pyplot as plt

# xmin, xmax= -2*np.pi, 2*np.pi
# ymin, ymax = -2*np.pi, 2*np.pi
# Ngrid = 10
#
# X = np.linspace(xmin, xmax, Ngrid)
# Y = np.linspace(ymin, ymax, Ngrid)
#
# Y,X= np.mgrid[xmin:xmax:Ngrid*1j, ymin:ymax:Ngrid*1j]
# VX = Y
# VY = -X
# print(X,Y)
# ##########field line calculation#######
# r=[np.pi,0]
# dl=0.0001
#
#
#
# plt.streamplot(X,Y,VX,VY, density=3)
# plt.xlim(xmin,xmax)
# plt.ylim(ymin,ymax)
# plt.show()

B4LINE = np.loadtxt("B4LINE.dat")
print(B4LINE)

z = B4LINE[:,0]
r = B4LINE[:,1]
vr = B4LINE[:,2]
vz = B4LINE[:,3]

plt.quiver(r,z,vr,vz)
plt.savefig("inputfile.png")
plt.show()
