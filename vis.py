import numpy as np
import matplotlib.pyplot as plt

def weight(x,y):
    return 1/(x*x+y*y)

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

def localAction(point, radius=0.0001):

    r_distances=r-point[0]
    z_distances=z-point[1]
    indices_in_radius = (r_distances**2+z_distances**2<radius**2)
    number_points_inside_radius=np.count_nonzero(indices_in_radius)
    # print(indices_in_radius)
    print(number_points_inside_radius)

    # fig = plt.gcf()
    # plt.quiver(r,z,vr,vz, label="Input data")
    # plt.scatter(r,z, label="Input data", alpha=0.05)
    # plt.scatter(r[indices_in_radius],z[indices_in_radius],
    #     label="Input data in circle", alpha=0.5)
    # plt.plot(point[0], point[1], "bo-", label="plotted point" )
    # plt.xlabel("r")
    # plt.ylabel("z")
    # fig.gca().add_artist(plt.Circle(point, radius, color='g', alpha=0.5))
    # plt.grid()
    # plt.legend()
    # plt.axes().set_aspect('equal', 'datalim')
    # plt.xlim(point[0]-2*radius, point[0]+2*radius)
    # plt.ylim(point[1]-2*radius, point[1]+2*radius)
    # plt.show()

    iterations=0

    if(iterations>10):
        return radius
    if(number_points_inside_radius>8):
        localAction(point, radius*0.9)
        iterations+=1
    elif(number_points_inside_radius<8):
        localAction(point, radius*1.055)
        iterations+=1
    else:
        return radius


#read data
B4LINE = np.loadtxt("B4LINE.dat")
#print(B4LINE)
#read input data and
z = B4LINE[:,0]
r = B4LINE[:,1]
vr = B4LINE[:,2]
vz = B4LINE[:,3]
plt.quiver(r,z,vr,vz, label="Input data")
plt.scatter(r,z, label="Input data", alpha=0.05)

point = np.array([0.00068850, 0.003625])

# StartingZ=np.linspace(min(z),max(z),100)
# #StartingZ=np.hstack((StartingZ, np.zeros_like(StartingZ)))
# plt.plot(np.zeros_like(StartingZ), StartingZ, label="Starting data")



plt.xlabel("r")
plt.ylabel("z")
plt.grid()
plt.legend()
# plt.show()
plt.clf()

localAction(point)
