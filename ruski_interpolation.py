import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1,1,500)
y = np.linspace(-1,1,500)
X,Y = np.meshgrid(x,y)
VX = -Y
VY = X
fig=plt.figure(figsize=(11,7),dpi=200)
plt.quiver(X[::20, ::20],Y[::20, ::20],VX[::20, ::20],VY[::20, ::20], alpha=1)

def weight(x,y):
    return 1/(x*x+y*y)

def field(r, step_length):
    radius=step_length
    #TODO: time this part
    x_distances=X-r[0]
    y_distances=Y-r[1]
    indices_in_radius = x_distances**2+y_distances**2<radius**2
    number_points_inside_radius=np.count_nonzero(indices_in_radius)
    if(number_points_inside_radius<10):
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


    return np.array([vx_interpolated, vy_interpolated])

acceptable_error = 0.001 #relative to step length
alpha_coefficient = 0.25 #

def step(r,step_length):
    v = field(r, step_length)

    rk1 = r + v*step_length
    vk1 = field(rk1, step_length)
    rk = r+(v+vk1)/2 * step_length #the first approximation

    rk1half = r + v*step_length/2
    vk1half = field(rk1half, step_length)
    rkhalf = r+(v+vk1half)/2 * step_length/2

    vkhalf = field(rkhalf, step_length)
    rkstar1 = rkhalf + vkhalf * step_length/2
    vkstar1 = field(rkstar1, step_length)
    rkstar = rkhalf + (vkhalf+vkstar1)/2 * step_length/2 #the second approximation
    vkstar = field(rkstar, step_length)

    relative_difference=np.linalg.norm(rk-rkstar)/step_length
    if(relative_difference>acceptable_error):
        step_length/=2
    elif(relative_difference<alpha_coefficient*acceptable_error):
        step_length*=2

    return rkstar, vkstar, step_length

def line(starting_r):
    checking_if_closed = False
    closed = False
    r=starting_r
    step_length = 0.1

    r_array=np.copy(starting_r)
    #TODO: FIELD INTERPOLATION
    v_array=field(starting_r, step_length)
    iterations=0
    while not closed:

        distance=np.linalg.norm(r-starting_r)
        if(not checking_if_closed and distance>step_length*3):
            checking_if_closed = True
        if (checking_if_closed and distance<0.5*step_length):
            closed=True

        iterations+=1
        r,v, step_length=step(r, step_length)
        r_array=np.vstack((r_array,r))
        v_array=np.vstack((v_array,v))

        if(iterations>1000):
            print("Does not converge from r = " + str(starting_r))
            print("Step size was " + str(step_length))
            plt.plot(r_array[:,0], r_array[:,1], "r-")

            empty_result=np.zeros_like(r_array)
            return empty_result, empty_result
    print("Closed in " + str(iterations) + " iterations. Final step length was " + str(step_length))
    return r_array, v_array

rki = np.linspace(0.2,0.9,10)
for r in rki:
    r = np.array([r,0])
    r_array, v_array = line(r)
    plt.plot(r_array[:,0], r_array[:,1], "b-")
    plt.quiver(r_array[:,0], r_array[:,1], v_array[:,0], v_array[:,1], alpha=0.1, color="blue")

plt.xlim(-1,1)
plt.ylim(-1,1)
plt.show()
