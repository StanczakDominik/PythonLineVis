#!/usr/bin/env python3
import time
import tkinter as Tkinter
import os
from tkinter.filedialog import askopenfilename
import tkinter.messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import logging
import subprocess
logging.basicConfig(filename='debug.log',filemode='w',level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.captureWarnings(True)

L=0
K=0
WhileRunning=False
def plot_format():
    plt.grid()
    plt.xlabel("z")
    plt.ylabel("r")

class MyDialog:
    def __init__(self, parent):
        top = self.top = Tkinter.Toplevel(parent)
        Tkinter.Label(top, text="Indeksy L oraz K\nPrzed wpisaniem zamknij wykres!").pack()
        self.L = Tkinter.Entry(top)
        self.L.pack(padx=5)

        self.K = Tkinter.Entry(top)
        self.K.pack(padx=5)

        b = Tkinter.Button(top, text="OK", command=self.ok)
        b.pack()
        quitbutton = Tkinter.Button(top, text="Koniec", command=self.quitc)
        quitbutton.pack()

    def ok(self):
        global L, K, WhileRunning
        L=int(self.L.get())
        K=int(self.K.get())
        print ("value is", L, K)
        WhileRunning=True
        self.top.destroy()
    def quitc(self):
        global WhileRunning
        WhileRunning=False
        self.top.destroy()

class window(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent=parent
        self.initialize()
    def initialize(self):
        self.grid()
        self.SeparatePlots=True
        self.WhileRunning=False
        #Brzeg
        self.Brzeg = Tkinter.Button(self, text=u"Pokaż wektory normalne na brzegu", command = self.OnBrzegClick)
        self.Brzeg.grid(column=0,row=0,sticky='NSEW')

        self.BrzegLabel = Tkinter.Label(self, text=u"Dane w formacie Z, R, Vz, Vr (brzeg)")
        self.BrzegLabel.grid(column=1,row=0,sticky='NSEW')


        #Sasiedzi na brzegu
        self.Sasiedzi = Tkinter.Button(self, text=u"Pokaż sąsiadów na brzegu", command = self.OnSasiedziClick)
        self.Sasiedzi.grid(column=0,row=1,sticky='NSEW')
        self.SasiedziLabel = Tkinter.Label(self, text=u"Dane w formacie LI, KI, RA, ZA")
        self.SasiedziLabel.grid(column=1,row=1,sticky='NSEW')

        #Sąsiedzi wewnątrz
        self.SasiedziInside = Tkinter.Button(self,text=u"Pokaż sąsiadów wewnątrz", command = self.OnSasiedziInsideClick)
        self.SasiedziInside.grid(column=0,row=2,sticky='NSEW')

        self.SasiedziInsideLabel = Tkinter.Label(self, text=u"Dane z pliku preparowanego przez sasiad.exe")
        self.SasiedziInsideLabel.grid(column=1,row=2,sticky='NSEW')

        #Interpolacja linii
        self.Linie = Tkinter.Button(self, text=u"Interpoluj linie", command = self.OnLinieClick)
        self.Linie.grid(column=0,row=3,sticky='NSEW')

        self.LinieLabel = Tkinter.Label(self, text=u"Dane w formacie Z, R, Vz, Vr (brzeg")
        self.LinieLabel.grid(column=1,row=3,sticky='NSEW')

        #Spline
        self.Spline = Tkinter.Button(self, text=u"Fituj spline'a", command = self.OnSplineClick)
        self.Spline.grid(column=0,row=4,sticky='NSEW')

        self.SplineLabel = Tkinter.Label(self, text=u"Dane w formacie Z, R, Vz, Vr (brzeg)")
        self.SplineLabel.grid(column=1,row=4,sticky='NSEW')

        #Łączenie wyników
        # self.Show = Tkinter.Button(self,bg="lightblue", text=u"Łącz wyniki", command = self.OnShowClick)
        # self.Show.grid(column=0,row=5,sticky='NSEW')

        # self.ShowLabel = Tkinter.Label(self, text=u"Wchodzi w tryb łączenia danych.")
        # self.ShowLabel.grid(column=1,row=5,sticky='NSEW')


        #Komunikaty
        # self.Messenger = Tkinter.Label(self, text=u"Komunikaty", anchor="c")
        # self.Messenger.grid(column=0,row=6,sticky='NSEW', columnspan=2)

        # self.bind('<Return>', self.OnEnterPress)

        self.grid_columnconfigure(0,weight=1)
        self.resizable(0,0)
        self.update()
        self.geometry(self.geometry())

        self.wm_attributes('-topmost', 1)
    def messenger_print(self, message):
        # self.Messenger.configure(text=message)
        pass

    def OnSasiedziClick(self):
        neighbor_file_name = askopenfilename()
        print(neighbor_file_name)
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
        plot_format()
        if(self.SeparatePlots):
            plt.axes().set_aspect('equal', 'datalim')
            plt.show()
    def OnSasiedziInsideClick(self):
        global WhileRunning
        subprocess.call("sasiad.exe")
        file_name = "dane2sasiedzi"
        print("Czytam dane z plików...")
        plik = np.loadtxt(file_name)
        fort15 =plik[:,:2].astype(int)
        fort16=plik[:,2:]
##        fort15 = np.loadtxt(indices_file_name, dtype=int)
##        fort16 = np.loadtxt(positions_file_name)
        WhileRunning=True
        print("hm")
        d = MyDialog(self)
        self.wait_window(d.top)

        while(WhileRunning):
            try:
                number_of_batches=int(len(fort15)/9)
                print(L,K)
                first_index_correct=fort15[:,0]==L
                second_index_correct=fort15[:,1]==K
                print(first_index_correct, second_index_correct)
                resulting_array=first_index_correct*second_index_correct
                resulting_indices=np.nonzero(resulting_array)
                print(resulting_indices)

                found_values=fort16[resulting_indices]
                Z = found_values[:,1]
                R = found_values[:,0]

                center_z = Z[4]
                center_r = R[4]
                center_r_array=np.ones(9)*center_r
                center_z_array=np.ones(9)*center_z

                r_array=np.vstack((center_r_array,R))
                z_array=np.vstack((center_z_array,Z))
                plt.plot(fort16[:,1], fort16[:,0], "k.")
                plt.plot(z_array,r_array, "r-")
                plt.plot(center_z, center_r, "ro")
                plot_format()
                #plt.plot(Z,R, 'bo')
                #plt.plot(Z[4], R[4], 'ro')
                print("wtf")
                if(self.SeparatePlots):
                    plt.axes().set_aspect('equal', 'datalim')
                    plt.show(block=False)
                print("boo")
            except IndexError:
                tkinter.messagebox.showwarning("Nietrafiony indeks", "Nie ma takiego indeksu. Spróbuj z innym.")
            d = MyDialog(self)
            self.wait_window(d.top)
    def OnBrzegClick(self):
        data_file_name = askopenfilename()            #czyta nazwę pliku jako string
        self.messenger_print(data_file_name)
        self.messenger_print(" Wizualizuję brzeg.")
        data_file = np.loadtxt(data_file_name)  #czyta dane liczbowe

        Xbrzeg = data_file[:,0]             #tworzy jednowymiarowe tablice...
        Ybrzeg = data_file[:,1]             #danych zczytanych z pliku
        VXbrzeg = data_file[:,2]
        VYbrzeg = data_file[:,3]
        print(Xbrzeg.shape)
        indices = np.argsort(Xbrzeg)   #zwraca indeksy - kolejność posortowanych danych)
        Xbrzeg=Xbrzeg[indices]
        Ybrzeg=Ybrzeg[indices]
        VXbrzeg=VXbrzeg[indices]
        VYbrzeg=VYbrzeg[indices]
        plt.plot(Xbrzeg, Ybrzeg, "g-")
        plt.quiver(Xbrzeg, Ybrzeg, VXbrzeg, VYbrzeg, alpha=1, angles='xy', scale_units='xy', color="green") #brzeg jako strzałki
        plot_format()
        if(self.SeparatePlots):
            plt.axes().set_aspect('equal', 'datalim')
            plt.show()
    def OnLinieClick(self):
        simulation_data_file_name = askopenfilename()            #czyta nazwę pliku jako string
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
        # print("""Czy próbować obcinać interpolowane linie w próżni?
        # Domyślnie - nie, wciśnij enter. Aby je obcinać, wpisz dowolny znak.""")
        # try_to_cutoff_lines=input()


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
                    self.messenger_print("Obcięto linię  " + str(starting_r))
                    break

                if(iterations>10000):
                    self.messenger_print("Linia " + str(starting_r) + " zapętla się.")
                    plt.plot(r_array[:,0], r_array[:,1], "r-")
                    empty_result=np.zeros_like(r_array)
                    return empty_result, empty_result
            self.messenger_print("Linia " + str(starting_r) + " zamyka się lub wychodzi z obszaru w " + str(iterations) + " kroku.")
            return r_array, v_array

        starting_point_file_name = "vis_punkty_poczatkowe.dat"

        print("Zaczynam interpolację.")
        rki = np.loadtxt(starting_point_file_name)

        for r in rki:
            r_array, v_array = line(r)
            plt.plot(r_array[:,0], r_array[:,1], "-")
            plt.quiver(r_array[:,0], r_array[:,1], v_array[:,0], v_array[:,1], alpha=0, color="blue",
             angles='xy', scale_units='xy')

        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        plt.grid()
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()
    def OnSplineClick(self):
        data_file_name = askopenfilename()            #czyta nazwę pliku jako string
        self.messenger_print(data_file_name)
        self.messenger_print(" Wizualizuję brzeg.")
        data_file = np.loadtxt(data_file_name)  #czyta dane liczbowe
        print(data_file)
        print(" koniec pliku")

        Xbrzeg = data_file[:,0]             #tworzy jednowymiarowe tablice...
        Ybrzeg = data_file[:,1]             #danych zczytanych z pliku
        VXbrzeg = data_file[:,2]
        VYbrzeg = data_file[:,3]
        print(Xbrzeg.shape)
        indices = np.argsort(Xbrzeg)   #zwraca indeksy - kolejność posortowanych danych)
        Xbrzeg=Xbrzeg[indices]
        Ybrzeg=Ybrzeg[indices]
        VXbrzeg=VXbrzeg[indices]
        VYbrzeg=VYbrzeg[indices]
        print(Xbrzeg.shape)
        print(Ybrzeg.shape)
        tck=interpolate.splrep(Xbrzeg,Ybrzeg,k=5,s=0.0000001)
        print(tck)
        xspline=np.linspace(min(Xbrzeg),max(Xbrzeg),5000)

        spline=interpolate.splev(xspline,tck,der=0)
        spline_gradient=interpolate.splev(Xbrzeg,tck,der=1)
        spline_normal_delta_x=1/np.sqrt(1+spline_gradient**2)
        spline_normal_delta_y=spline_normal_delta_x*spline_gradient
        self.messenger_print("Zaznaczam geometryczne (z dopasowania spline'a) wektory normalne na niebiesko, zaś z danych - na zielono.")
        plt.plot(Xbrzeg, Ybrzeg, "g-")
        plt.plot(xspline, spline, "b-")
        plt.quiver(Xbrzeg, Ybrzeg, VXbrzeg, VYbrzeg, alpha=1, angles='xy', scale_units='xy', color="green") #brzeg jako strzałki
        plt.quiver(Xbrzeg,Ybrzeg, -spline_normal_delta_y, spline_normal_delta_x, alpha=1, angles='xy', scale_units='xy', color="blue")
        plot_format()
        if(self.SeparatePlots):
            plt.axes().set_aspect('equal', 'datalim')
            plt.show()
    def OnShowClick(self):
        pass

if __name__ == "__main__":
    app = window(None)
    app.title("Wizualizacja danych")
    app.mainloop()
