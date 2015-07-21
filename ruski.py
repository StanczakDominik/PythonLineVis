import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1,1,100)
y = np.linspace(-1,1,100)
X,Y = np.meshgrid(x,y)
VX = -Y
VY = X
fig=plt.figure(figsize=(11,7),dpi=200)
plt.quiver(X,Y,VX,VY, alpha=1)

def field(r):
    return np.array([-r[1], r[0]])


acceptable_error = 0.001 #relative to step length
alpha_coefficient = 0.25 #

def step(r,step_length):
    v = field(r)

    rk1 = r + v*step_length
    vk1 = field(rk1)
    rk = r+(v+vk1)/2 * step_length #the first approximation

    rk1half = r + v*step_length/2
    vk1half = field(rk1half)
    rkhalf = r+(v+vk1half)/2 * step_length/2

    vkhalf = field(rkhalf)
    rkstar1 = rkhalf + vkhalf * step_length/2
    vkstar1 = field(rkstar1)
    rkstar = rkhalf + (vkhalf+vkstar1)/2 * step_length/2 #the second approximation
    vkstar = field(rkstar)

    relative_difference=np.linalg.norm(rk-rkstar)/step_length
    if(relative_difference>acceptable_error):
        step_length/=2
        # print(str(relative_difference) + " Decreasing")
    elif(relative_difference<alpha_coefficient*acceptable_error):
        step_length*=2
        # print(str(relative_difference) + " Increasing")

    return rkstar, vkstar, step_length

def line(starting_r):
    checking_if_closed = False
    closed = False
    r=starting_r
    step_length = 0.1

    r_array=np.copy(starting_r)
    v_array=field(starting_r)
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
    print("Closed in " + str(iterations) + " iterations")
    print("Final step length was " + str(step_length))
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
