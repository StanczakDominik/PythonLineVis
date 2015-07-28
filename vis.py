#Wizualizacja linii pola
#Dominik Stańczak
#wersja 24.07.2015

#Wymaga bibliotek numpy, matplotlib, scipy
#Najprościej dobrać się do nich poprzez dystrybucję WinPython
#http://winpython.github.io/
#aktualna wersja do znalezienia na
#https://github.com/StanczakDominik/PythonLineVis

#W razie problemów proszę o wysłanie pliku debug.log i ewentualnie
# plików z danymi na mojego maila.

CZY_ZAPISYWAC_OBRAZKI=0

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import logging
logging.basicConfig(filename='debug.log',filemode='w',level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.captureWarnings(True)

import time
timestr = time.strftime("%Y_%m_%d-%H_%M_%S") #czas uruchomienia pliku

fig=plt.figure(figsize=(11,7),dpi=100)
plt.grid()
plt.xlabel("z")
plt.ylabel("r")

print("""Podaj nazwę pliku z danymi brzegu (sformatowanymi w ten sposób: z, r, v_z, v_r)
    Może to również być cała ścieżka do pliku, dopuszczalne spacje itd. Przykład:
    A11LINE.dat
    Wciśnij enter bez wpisywania nazwy aby pominąć wizualizację brzegu""")
data_file_name = input()            #czyta nazwę pliku jako string
if data_file_name:
    print(" Wizualizuję brzeg.")
    data_file = np.loadtxt(data_file_name)  #czyta dane liczbowe
    print(data_file)
    print(" koniec pliku")

    Xbrzeg = data_file[:,0]             #tworzy jednowymiarowe tablice...
    Ybrzeg = data_file[:,1]             #danych zczytanych z pliku
    VXbrzeg = data_file[:,2]
    VYbrzeg = data_file[:,3]
    print(Xbrzeg.shape)
    #normalizacja prędkości na brzegu do 1
##    Vbrzeg=np.sqrt(VXbrzeg**2+VYbrzeg**2)
##    VXbrzeg/=Vbrzeg
##    VYbrzeg/=Vbrzeg

    #sortowanie punktów brzegowych według położeń na osi x (z)
    indices = np.argsort(Xbrzeg)   #zwraca indeksy - kolejność posortowanych danych)
    Xbrzeg=Xbrzeg[indices]
    Ybrzeg=Ybrzeg[indices]
    VXbrzeg=VXbrzeg[indices]
    VYbrzeg=VYbrzeg[indices]

    #interpolacja spline'm oraz obliczanie wektorów normalnych
    print(Xbrzeg.shape)
    print(Ybrzeg.shape)
    tck=interpolate.splrep(Xbrzeg,Ybrzeg,k=5,s=0.0000001)
    print(tck)
    xspline=np.linspace(min(Xbrzeg),max(Xbrzeg),5000)
    
    spline=interpolate.splev(xspline,tck,der=0)
    spline_gradient=interpolate.splev(Xbrzeg,tck,der=1)
    spline_normal_delta_x=1/np.sqrt(1+spline_gradient**2)
    spline_normal_delta_y=spline_normal_delta_x*spline_gradient
    print(" Zaznaczam geometryczne (z dopasowania spline'a) wektory normalne na niebiesko, zaś z danych - na zielono.")
    plt.plot(Xbrzeg, Ybrzeg, "g-")
    plt.plot(xspline, spline, "b-")
    plt.quiver(Xbrzeg, Ybrzeg, VXbrzeg, VYbrzeg, alpha=1, angles='xy', scale_units='xy', color="green") #brzeg jako strzałki
    plt.quiver(Xbrzeg,Ybrzeg, -spline_normal_delta_y, spline_normal_delta_x, alpha=1, angles='xy', scale_units='xy', color="blue")

print("""Czy wyświetlać sąsiadów? Jeśli nie, pozostaw puste. Jeśli tak,
    podaj nazwę pliku z ich danymi. Przykład: sasiad.dat. Format taki, jaki generuje załączony program fortranowy sasiad.for""")
neighbor_file_name=input()
if neighbor_file_name:
    data = np.loadtxt(neighbor_file_name)
    number_of_batches=int(len(data)/9)
    for i in range(number_of_batches):
        batch=data[i*9:9*(i+1),2:4]

        #NEIGHBORS
        RA = batch[:,0]    #R position of particles
        ZA = batch[:,1]    #Z position of particles
        center_z = ZA[4]
        center_r = RA[4]
        center_r_array=np.ones(9)*center_r
        center_z_array=np.ones(9)*center_z

        r_array=np.vstack((center_r_array,RA))
        z_array=np.vstack((center_z_array,ZA))
        plt.plot(z_array,r_array, "k-")
        plt.plot(center_z, center_r, "o")

if CZY_ZAPISYWAC_OBRAZKI and (neighbor_file_name or data_file_name):
    brzeg_file_name="brzeg_sasiedzi"+timestr+".png"
    plt.savefig(brzeg_file_name)
    print(" Brzeg zapisany do pliku " + brzeg_file_name)




print("Czy wyświetlić w tej chwili dane brzegu? Pozostaw puste, jeśli nie.")
show_brzeg_data=input()
if show_brzeg_data:
    plt.show()
    fig=plt.figure(figsize=(11,7),dpi=100)
    plt.grid()
    plt.xlabel("z")
    plt.ylabel("r")



##SYMULACJA




print("""Czy interpolować linie? Jeśli tak, podaj nazwę pliku z danymi całej symulacji.
    Jeśli nie, pozostaw puste.
    Formatowanie z, r, vz, vr. Przykład: B4LINE.dat""")
simulation_data_file_name = input()            #czyta nazwę pliku jako string
if simulation_data_file_name:
    simulation_data_file = np.loadtxt(simulation_data_file_name)
    X = simulation_data_file[:,0]
    Y = simulation_data_file[:,1]
    xmin=min(X)
    xmax=max(X)
    ymin=min(Y)
    ymax=max(Y)
    VX = simulation_data_file[:,3]
    VY = simulation_data_file[:,2]
    plt.xlabel("z")
    plt.ylabel("r")

    try_to_cutoff_lines=False
    print("""Czy próbować obcinać interpolowane linie w próżni?
    Domyślnie - nie, wciśnij enter. Aby je obcinać, wpisz dowolny znak.""")
    try_to_cutoff_lines=input()


    plt.quiver(X[::1],Y[::1],VX[::1],VY[::1], alpha=1, angles='xy', scale_units='xy')
    #rozwiązanie skali strzałek dzięki uprzejmości miłych ludzi z
    #http://stackoverflow.com/questions/12079842/quiver-plot-arrow-aspect-ratio


    def weight(x,y):
        #zastosowana waga punktów: 1/r^3
        r=x*x+y*y
        r3=r*r*r
        return 1/r3


    def field(r, step_length):
        cutoff=False
        radius=3*step_length
        x_distances=X-r[0]
        y_distances=Y-r[1]
        distances_squared = x_distances**2+y_distances**2

        if try_to_cutoff_lines and np.min(distances_squared)>0.00005**2:
            cutoff=True

        indices_in_radius = distances_squared<radius**2
        number_points_inside_radius=np.count_nonzero(indices_in_radius)
        if(number_points_inside_radius<15):
            return field(r,step_length*2)
        x_distances_inside=x_distances[indices_in_radius]
        y_distances_inside=y_distances[indices_in_radius]
        x_velocities_inside=VX[indices_in_radius]
        y_velocities_inside=VY[indices_in_radius]
        weights=weight(x_distances_inside, y_distances_inside)
        weight_sum=np.sum(weights)
        vx_interpolated=np.sum(weights*x_velocities_inside)/weight_sum
        vy_interpolated=np.sum(weights*y_velocities_inside)/weight_sum
        v_vector=np.array([vx_interpolated, vy_interpolated])

        v_vector/= np.linalg.norm(v_vector)

        return v_vector, step_length, cutoff

    acceptable_error = 0.0000001 #relative to step length
    alpha_coefficient = 0.5 #

    def step(r,step_length):
        v, step_length, cutoff = field(r, step_length)

        rk1 = r + v*step_length
        vk1,dummy, dummy2 = field(rk1, step_length)
        rk = r+(v+vk1)/2 * step_length #the first approximation

        rk1half = r + v*step_length/2
        vk1half,dummy, dummy2 = field(rk1half, step_length)
        rkhalf = r+(v+vk1half)/2 * step_length/2

        vkhalf,dummy, dummy2 = field(rkhalf, step_length)
        rkstar1 = rkhalf + vkhalf * step_length/2
        vkstar1,dummy, dummy2 = field(rkstar1, step_length)
        rkstar = rkhalf + (vkhalf+vkstar1)/2 * step_length/2 #the second approximation
        vkstar,dummy, dummy2 = field(rkstar, step_length)

        relative_difference=np.linalg.norm(rk-rkstar)/step_length
        if(relative_difference>acceptable_error):
            step_length/=2
        elif(relative_difference<alpha_coefficient*acceptable_error):
            step_length*=2

        return rkstar, vkstar, step_length, cutoff

    def line(starting_r):
        checking_if_closed = False
        closed = False
        out_of_bounds = False
        r=starting_r
        step_length = 0.000001
        cutoff=False
        r_array=np.copy(starting_r)
        v_array,step_length, cutoff=field(starting_r, step_length)
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
            r,v, step_length, cutoff=step(r, step_length)
            r_array=np.vstack((r_array,r))
            v_array=np.vstack((v_array,v))

            if cutoff:
                print("Obcięto linię zaczynającą się w " + str(starting_r) + " zapętla się.")
                break

            if(iterations>10000):
                print("Linia zaczynająca się w " + str(starting_r) + " zapętla się.")
                plt.plot(r_array[:,0], r_array[:,1], "r-")
                empty_result=np.zeros_like(r_array)
                return empty_result, empty_result
        print("Linia zaczynająca się w " + str(starting_r) + " zamyka się lub wychodzi poza obszar w " + str(iterations) + " iteracji.")
        return r_array, v_array

    print("""Podaj nazwę pliku z danymi punktów początkowych linii.
    Domyślnie vis_punkty_poczatkowe.dat (aby tego nie zmieniać, pozostaw puste)
    Dane z których program zaczyna interpolację powinny być sformatowane w ten sposób: z, r""")
    starting_point_file_name=input()
    if not starting_point_file_name:
        starting_point_file_name = "vis_punkty_poczatkowe.dat"

    print("Zaczynam interpolację.")
    rki = np.loadtxt(starting_point_file_name)
    #
    # rki = [np.array([-0.0049, 0.0002]),
    #         np.array([-0.0049, 0.0008]),
    #         np.array([-0.0049, 0.0005]),
    #         np.array([-0.0049, 0.0010]),
    #         np.array([-0.0049, 0.0011]),
    #         np.array([-0.0049, 0.0012]),
    #         np.array([-0.0049, 0.0013]),
    #         np.array([-0.0049, 0.0014]),
    #         np.array([-0.0049, 0.0015]),
    #         np.array([-0.0049, 0.0015]),
    #         np.array([-0.0049, 0.0016]),
    #         np.array([-0.0049, 0.0017]),
    #         np.array([-0.0049, 0.0018]),
    #         np.array([-0.0049, 0.0019]),
    #         np.array([-0.0049, 0.0020]),
    #         np.array([-0.0049, 0.0021]),
    #         np.array([-0.0049, 0.0022]),
    #         np.array([-0.0049, 0.0023]),
    #         np.array([-0.0049, 0.0024]),
    #         np.array([-0.0049, 0.0025]),
    #         np.array([-0.0049, 0.0026]),
    #         np.array([-0.0049, 0.0027]),
    #         np.array([-0.0049, 0.0028]),
    #         np.array([-0.0049, 0.0029]),
    #         np.array([-0.0049, 0.0030]),
    #         np.array([-0.0049, 0.0031]),
    #         np.array([-0.0049, 0.0032]),
    #         np.array([-0.0049, 0.0033]),
    #         np.array([-0.0049, 0.0034]),
    #         np.array([-0.0049, 0.0035]),
    # ##        np.array([-0.0049, 0.00367806]),
    # ##        np.array([-0.0027527, 0.00226012]),
    # ##        np.array([-0.00250663, 0.00197791]),
    # ##        np.array([-0.00152236, 0.00144696]), #bugs out
    # ##        np.array([-0.00121477, 0.00118866]), #bugs out
    # ##        np.array([-0.000661117, 0.000935142]), #bugs
    # ##        np.array([1.55715e-5, 0.000916008]) #bugs out as well
    #         ]
    for r in rki:
        r_array, v_array = line(r)
        plt.plot(r_array[:,0], r_array[:,1], "-")
        plt.quiver(r_array[:,0], r_array[:,1], v_array[:,0], v_array[:,1], alpha=0, color="blue",
         angles='xy', scale_units='xy')

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.grid()
    if CZY_ZAPISYWAC_OBRAZKI:
        visualization_file_name=simulation_data_file_name[:-4] + "wizualizacja"+timestr+".png"
        plt.savefig(visualization_file_name)
        print("Wizualizacja wyników zapisana do pliku " + visualization_file_name)

    plt.show()
print("Koniec działania programu.")
