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

def localAction(point, radius=0.0001, iterations=0):

    r_distances=r-point[0]
    z_distances=z-point[1]
    indices_in_radius = (r_distances**2+z_distances**2<radius**2)
    number_points_inside_radius=np.count_nonzero(indices_in_radius)
    # print(indices_in_radius)
    # print(number_points_inside_radius)

    # fig = plt.gcf()
    # plt.quiver(r[indices_in_radius],z[indices_in_radius],vr[indices_in_radius],vz[indices_in_radius], label="Input data")
    # plt.scatter(r[indices_in_radius],z[indices_in_radius], label="Input data", alpha=0.05)
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

    if(iterations>30):
        # print("Next!")
        r_distances_inside=r_distances[indices_in_radius]
        z_distances_inside=z_distances[indices_in_radius]
        r_velocities_inside=vr[indices_in_radius]
        z_velocities_inside=vz[indices_in_radius]
        weights=weight(r_distances_inside, z_distances_inside)
        weight_sum=np.sum(weights)
        vr_interpolated=np.sum(weights*r_velocities_inside)/weight_sum
        vz_interpolated=np.sum(weights*z_velocities_inside)/weight_sum
        return np.array([vr_interpolated, vz_interpolated]), radius
    if(number_points_inside_radius>8):
        # density = number_points_inside_radius/np.pi/radius**2
        # radius=np.sqrt(density/8/np.pi)
        iterations+=1
        # return localAction(point, radius)
        return localAction(point, radius*0.9, iterations)
    elif(number_points_inside_radius<8):
        iterations+=1
        return localAction(point, radius*1.055, iterations)
    else:
        r_distances_inside=r_distances[indices_in_radius]
        z_distances_inside=z_distances[indices_in_radius]
        r_velocities_inside=vr[indices_in_radius]
        z_velocities_inside=vz[indices_in_radius]
        weights=weight(r_distances_inside, z_distances_inside)
        weight_sum=np.sum(weights)
        vr_interpolated=np.sum(weights*r_velocities_inside)/weight_sum
        vz_interpolated=np.sum(weights*z_velocities_inside)/weight_sum
        # print(vr_interpolated, vz_interpolated)
        # print(r_velocities_inside)
        # print(z_velocities_inside)
        # fig = plt.gcf()
        # plt.quiver(r[indices_in_radius],z[indices_in_radius],vr[indices_in_radius],
        #     vz[indices_in_radius], label="Input data")
        # plt.quiver(point[0],point[1],vr_interpolated,vz_interpolated,
        #     label="Found point")
        # #plt.scatter(r,z, label="Input data", alpha=0.05)
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
        # print (np.array([vr_interpolated, vz_interpolated]), radius)
        # print("REACHED")
        # return (vr_interpolated, vz_interpolated, radius)
        return np.array([vr_interpolated, vz_interpolated]), radius

def step_to_next_point(points):
    last_point=points[-1]
    interpolation = localAction(last_point)
    # print("last point", last_point)
    # print("this",interpolation)
    local_velocities, step=localAction(last_point)
    shift=step*local_velocities/np.sqrt(np.sum(local_velocities**2))
    new_point = last_point + shift
    # print("diff", new_point-last_point)
    points=np.vstack((points,new_point))
    # print(new_point)
    return points

#read data
B4LINE = np.loadtxt("B4LINE.dat")
#print(B4LINE)
#read input data and
z = B4LINE[:,0]
r = B4LINE[:,1]
vr = B4LINE[:,2]
vz = B4LINE[:,3]

fig = plt.gcf()
plt.quiver(r,z,vr,vz, label="Input data")

# StartingZ=np.linspace(min(z)*0.99,max(z)*0.99, 20)#max(z)*0.99,100)
# StartingR=np.linspace(min(r),0.000895139,20)
StartingZ=np.linspace(min(z)*0.99,min(z)*0.985, 2)#max(z)*0.99,100)
StartingR=np.linspace(min(r),0.0035,100)
# StartingZ=np.vstack(((np.ones_like(StartingZ)*(max(r)-min(r))/10+min(r)).T, StartingZ.T)).T
# StartingZ=np.vstack((StartingR.T, StartingZ.T)).T
def fieldline(point):
    plt.plot(point[0], point[1], "go", alpha=0.05)
    points=np.array([[point[0], point[1]]])
    print(points[-1])
    i=0
    while(points[-1,0] < max(r) and points[-1,0]>min(r) and points[-1,1]<max(z) and points[-1,1]>min(z)):
        i+=1
        # print(str(i) +"\t", end="")
        points=step_to_next_point(points)
    plt.plot(points[-1,0], points[-1,1], "rx")
    plt.plot(points[:,0], points[:,1], "b-")
#plt.scatter(r,z, label="Input data", alpha=0.05)
#plt.scatter(r,z,
#    label="Input data in circle", alpha=0.5)
# fieldline(np.array([0.000812395, 0.00130208]))
# for linenumber in range(len(StartingZ)):
#     point=np.array([StartingZ[linenumber,0], StartingZ[linenumber,1]])
#     fieldline(point)
for i in StartingZ:
    for j in StartingR:
        point=np.array([j,i])
        fieldline(point)


# fieldline(np.array([0.00052395, 0.0000208]))
# # fieldline(np.array([[0.005, -0.003]]))
# fieldline(np.array([0.001, 0.002]))
# # fieldline(np.array([[0.001, 0.0002]]))
# # fieldline(np.array([[0.001, 0.0001]]))
# # fieldline(np.array([[0.001, 0.001]]))
# # fieldline(np.array([[0.001, -0.0001]]))
# fieldline(np.array([0.001, 0.004]))
# fieldline(np.array([0.001, 0.003]))
# fieldline(np.array([0.001, -0.004]))
# fieldline(np.array([0.001, -0.005]))
# fieldline(np.array([0.001, -0.003]))
# fieldline(np.array([0.001, -0.002]))
# fieldline(np.array([[0.001, -0.001]]))
plt.xlabel("r")
plt.ylabel("z")
plt.grid()
plt.legend()
plt.xlim(min(r),max(r))
plt.ylim(min(z),max(z))
plt.savefig("line4grid.png")
plt.show()
