import numpy as np
import matplotlib.pyplot as plt

B4LINE = np.loadtxt("B4LINE.dat")
X = B4LINE[:,0]
Y = B4LINE[:,1]
xmin=min(X)
xmax=max(X)
ymin=min(Y)
ymax=max(Y)
VX = B4LINE[:,3]
VY = B4LINE[:,2]
fig=plt.figure(figsize=(11,7),dpi=100)
plt.xlabel("z")
plt.ylabel("r")
plt.quiver(X[::5],Y[::5],VX[::5],VY[::5], alpha=1)

def weight(x,y):
    r=x*x+y*y
    r3=r*r*r*r*r
    return 1/r3

def field(r, step_length):
    radius=3*step_length
    #TODO: time this part
    x_distances=X-r[0]
    y_distances=Y-r[1]
    indices_in_radius = x_distances**2+y_distances**2<radius**2
    number_points_inside_radius=np.count_nonzero(indices_in_radius)
    #print(r, number_points_inside_radius, radius)
    if(number_points_inside_radius<5):
        return field(r,step_length*2)
    x_distances_inside=x_distances[indices_in_radius]
    y_distances_inside=x_distances[indices_in_radius]
    x_velocities_inside=VX[indices_in_radius]
    y_velocities_inside=VY[indices_in_radius]
    #TODO: use linear interpolation
    weights=weight(x_distances_inside, y_distances_inside)
    weight_sum=np.sum(weights)
    vx_interpolated=np.sum(weights*x_velocities_inside)/weight_sum
    vy_interpolated=np.sum(weights*y_velocities_inside)/weight_sum
    v_vector=np.array([vx_interpolated, vy_interpolated])

    v_vector/= np.linalg.norm(v_vector)

    #print(r, number_points_inside_radius, radius)


    return v_vector, step_length

acceptable_error = 0.00001 #relative to step length
alpha_coefficient = 0.25 #

def step(r,step_length):
    v, step_length = field(r, step_length)

    rk1 = r + v*step_length
    vk1,dummy = field(rk1, step_length)
    rk = r+(v+vk1)/2 * step_length #the first approximation

    rk1half = r + v*step_length/2
    vk1half,dummy = field(rk1half, step_length)
    rkhalf = r+(v+vk1half)/2 * step_length/2

    vkhalf,dummy = field(rkhalf, step_length)
    rkstar1 = rkhalf + vkhalf * step_length/2
    vkstar1,dummy = field(rkstar1, step_length)
    rkstar = rkhalf + (vkhalf+vkstar1)/2 * step_length/2 #the second approximation
    vkstar,dummy = field(rkstar, step_length)

    relative_difference=np.linalg.norm(rk-rkstar)/step_length
    if(relative_difference>acceptable_error):
        step_length/=2
    elif(relative_difference<alpha_coefficient*acceptable_error):
        step_length*=2

    return rkstar, vkstar, step_length

def line(starting_r):
    checking_if_closed = False
    closed = False
    out_of_bounds = False
    r=starting_r
    step_length = 0.00001

    r_array=np.copy(starting_r)
    #TODO: FIELD INTERPOLATION
    v_array,step_length=field(starting_r, step_length)
    iterations=0
    while not closed and not out_of_bounds:

        distance=np.linalg.norm(r-starting_r)
        if(not checking_if_closed and distance>step_length*3):
            checking_if_closed = True
        if (checking_if_closed and distance<0.5*step_length):
            closed=True
        if (r[0]>xmax or r[0] < xmin or r[1]>ymax or r[1]<ymin):
            out_of_bounds = True
        iterations+=1
        r,v, step_length=step(r, step_length)
        r_array=np.vstack((r_array,r))
        v_array=np.vstack((v_array,v))

        if(iterations>10000):
            print("Does not converge from r = " + str(starting_r))
            print("Step size was " + str(step_length))
            plt.plot(r_array[:,0], r_array[:,1], "r-")

            empty_result=np.zeros_like(r_array)
            return empty_result, empty_result
    print("Line starting at " + str(starting_r) + " closed in " + str(iterations) + " iterations. Final step length was " + str(step_length))
    return r_array, v_array

rki = [np.array([-0.0049, 0.0002]),
        np.array([-0.0049, 0.0008]),
        np.array([-0.0049, 0.0005]),
        np.array([-0.0049, 0.0010]),
        np.array([-0.0049, 0.0015]),
        np.array([-0.0049, 0.0020]),
        np.array([-0.0049, 0.0021]),
        np.array([-0.0049, 0.0022]),
        np.array([-0.0049, 0.0023]),
        np.array([-0.0049, 0.0024]),
        np.array([-0.0049, 0.0025]),
        np.array([-0.0049, 0.0030]),
        np.array([-0.0049, 0.0035]),
        np.array([-0.0049, 0.00367806])]
for r in rki:
    r_array, v_array = line(r)
    plt.plot(r_array[:,0], r_array[:,1], "-")
    plt.quiver(r_array[:,0], r_array[:,1], v_array[:,0], v_array[:,1], alpha=1, color="blue")

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.title("i dont even")
plt.savefig("output.png")
plt.show()
