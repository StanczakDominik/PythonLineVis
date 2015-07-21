import numpy as np
import matplotlib.pyplot as plt

def weight(x,y):
    return 1/(x*x+y*y)

def localAction(point, radius=0.0001, iterations=0):

    r_distances=r-point[1]
    z_distances=z-point[0]
    indices_in_radius = (r_distances**2+z_distances**2<radius**2)
    number_points_inside_radius=np.count_nonzero(indices_in_radius)
##    print(indices_in_radius)

    #DEBUG\GRAPHICS
##    print(number_points_inside_radius)
##    fig = plt.gcf()
##    plt.quiver(z[indices_in_radius],r[indices_in_radius],vz[indices_in_radius],vr[indices_in_radius], label="Input data", color="g")
##    plt.scatter(z,r, label="Input data", alpha=0.05)
##    plt.scatter(z[indices_in_radius],r[indices_in_radius],
##        label="Input data in circle", alpha=0.5)
##    plt.plot(point[0], point[1], "bo-", label="plotted point" )
##    plt.xlabel("z")
##    plt.ylabel("r")
##    fig.gca().add_artist(plt.Circle(point, radius, color='g', alpha=0.5))
##    plt.grid()
##    plt.legend()
##    plt.axes().set_aspect('equal', 'datalim')
##    plt.xlim(point[0]-2*radius, point[0]+2*radius)
##    plt.ylim(point[1]-2*radius, point[1]+2*radius)
##    plt.show()

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
        return np.array([vz_interpolated, vr_interpolated]), radius
    if(number_points_inside_radius>8):
        #TODO: faster scaling 
        # density = number_points_inside_radius/np.pi/radius**2
        # radius=np.sqrt(density/8/np.pi)
        iterations+=1
        return localAction(point, radius*0.8, iterations)
    elif(number_points_inside_radius<8):
        iterations+=1
        return localAction(point, radius*1.03, iterations)
    else:
        r_distances_inside=r_distances[indices_in_radius]
        z_distances_inside=z_distances[indices_in_radius]
        r_velocities_inside=vr[indices_in_radius]
        z_velocities_inside=vz[indices_in_radius]
        weights=weight(r_distances_inside, z_distances_inside)
        weight_sum=np.sum(weights)
        vr_interpolated=np.sum(weights*r_velocities_inside)/weight_sum
        vz_interpolated=np.sum(weights*z_velocities_inside)/weight_sum
        return np.array([vz_interpolated, vr_interpolated]), radius

def step_to_next_point(points):
    last_point=points[-1]
    # print("last point", last_point)
    # print("this",interpolation)
    local_velocities, step=localAction(last_point)
    shift=step*local_velocities/np.sqrt(np.sum(local_velocities**2))/5
    new_point = last_point + shift
    # print("diff", new_point-last_point)
    points=np.vstack((points,new_point))
    print(new_point)
    return points

#read data
B4LINE = np.loadtxt("B4LINE.dat")
#print(B4LINE)
#read input data and
z = B4LINE[:,0]
r = B4LINE[:,1]
vz = B4LINE[:,3]
vr = B4LINE[:,2]



##======ZOOOOOOMING=========
##fig = plt.gcf()
##plt.quiver(z,r,vz,vr, label="Input data")
##plt.xlabel("z")
##plt.ylabel("r")
##plt.grid()
##plt.legend()
##plt.xlim(-0.004,-0.002)
##plt.ylim(0.002,0.0035)
####plt.xlim(min(z),max(z))
####plt.ylim(min(r),max(r))
##plt.show()
##=================

fig = plt.gcf()
plt.quiver(z,r,vz,vr, label="Input data")

# StartingZ=np.linspace(min(z)*0.99,max(z)*0.99, 20)#max(z)*0.99,100)
# StartingR=np.linspace(min(r),0.000895139,20)
##StartingZ=np.linspace(min(z)*0.95,min(z)*0.98, 1)#max(z)*0.99,100)
StartingZ=np.linspace(min(z)*0.95,max(z)*0.95, 2)#max(z)*0.99,100)
StartingR=np.linspace(max(r)*0.01,max(r)*0.9,5)
# StartingZ=np.vstack(((np.ones_like(StartingZ)*(max(r)-min(r))/10+min(r)).T, StartingZ.T)).T
# StartingZ=np.vstack((StartingR.T, StartingZ.T)).T
def fieldline(point):
    plt.plot(point[0], point[1], "go", alpha=0.5)
    points=np.array([[point[0], point[1]]])
    print(points[-1])
    i=0
    while(points[-1,0] < max(z) and points[-1,0]>min(z) and points[-1,1]<max(r) and points[-1,1]>min(r)):
        i+=1
        #print(str(i) +"\t", end="")
        points=step_to_next_point(points)
    plt.plot(points[-1,0], points[-1,1], "rx")
    plt.plot(points[:,0], points[:,1], "-")
#plt.scatter(r,z, label="Input data", alpha=0.05)
#plt.scatter(r,z,
#    label="Input data in circle", alpha=0.5)
# fieldline(np.array([0.000812395, 0.00130208]))
# for linenumber in range(len(StartingZ)):
#     point=np.array([StartingZ[linenumber,0], StartingZ[linenumber,1]])
#     fieldline(point)
for i in StartingZ:
    for j in StartingR:
        point=np.array([i,j])
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
plt.xlabel("z")
plt.ylabel("r")
plt.grid()
plt.legend()
plt.xlim(min(z),max(z))
plt.ylim(min(r),max(r))
plt.savefig("local_line4grid.png")
plt.show()
