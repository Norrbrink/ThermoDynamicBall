"""
A module for creating simulations
"""

import matplotlib.pyplot as plt
import scipy as sp
import random

from thermodynamic_ball import Ball, Container

class Simulation():
    """
    A composition of the Ball and Container classes
    This class is used to create Simulations of Ideal Gases inside a container
    and investigate certain Physics parameters of the gases given certain initial conditions
    """
    def __init__(self, mass=10.0, radius=1.0, numberofballs=10, \
                 containerradius=10.0, containermass=1e20, sigma=20): 
        #initialing parameters of the simulation
        self._container = Container(containermass, containerradius) #creating the container
        self.balls = [self._container] #list of the balls
        self._N = numberofballs 
        
        #creating a guassion distribution for the initial velocties of the balls
        dist = sp.random.normal(0, sigma, numberofballs*2) 
        
        #creating an inscribed square grid of the container for the intial postions of the balls
        positionx = sp.arange(-(containerradius-radius-0.1)/sp.sqrt(2), \
            (containerradius-radius-0.1)/sp.sqrt(2), 2*radius + 0.5)
        positiony = sp.arange(-(containerradius-radius-0.1)/sp.sqrt(2), \
            (containerradius-radius-0.1)/sp.sqrt(2), 2*radius + 0.5)
        position = []
        for i in range(0, len(positionx)):
            for j in range(0, len(positiony)):
                position.append([positionx[i], positiony[j]])
        random.shuffle(position) #randomly shuffling this list to start the balls at random postions
        
        #creating the balls
        for i in range(0, self._N):
            self.balls.append(Ball(mass, radius, position[i], [dist[i], \
                              dist[numberofballs + i]]))
        
        #initialising varibles of the gas
        self._t = 0.0 #time elapsed of simulation
        self._t_list = [] #list of times at which collisions occur
        self._kinetic_energy_total = [] #list of total kinetic energy of the gas throughout time
        self._temperature = [] #list of temperature of the gas throughout time
        self._pressure_evolution = [[], []] #list of pressure caused by a collision, list of time of the collision  
        self._total_momentum = [] #list of net momentum of the gas throughout time
        
    def next_collision(self): 
        #function that checks when the next collsion will occur and causes it
        #finding which balls will collide next
        deltat = 1e10 
        firstball = 0
        secondball = 0
        for i in range(0, len(self.balls)): #checking how long until each ball collides
            for j in range(0, len(self.balls)):
                if i != j: 
                    currentdt = (self.balls[i].time_to_collision(self.balls[j]))
                    if currentdt < deltat and currentdt > 1e-10: #finding the minimum time until collsion
                        deltat = currentdt
                        firstball = i
                        secondball = j
        
        #performing the collision
        self._t += deltat #setting the time of the collsion
        self._t_list.append(self._t)
        self.move_all(deltat) #moving all the balls to the point of collision
        self.balls[firstball].collide(self.balls[secondball]) #causing the collsion
        
        #evaluating variables at time of collision
        if self.balls[firstball].isContainer: #calculating pressure
            self._pressure_evolution[0].append(self._container._pressuretime/ \
                                              (self._t - \
                                              self.balls[secondball].bounce_time))
            self._pressure_evolution[1].append(self._t)
            self.balls[secondball].bounce_time = self._t
        elif self.balls[secondball].isContainer:
            self._pressure_evolution[0].append(self._container._pressuretime/ \
                                              (self._t - \
                                              self.balls[firstball].bounce_time))
            self._pressure_evolution[1].append(self._t)
            self.balls[firstball].bounce_time = self._t
        self.kinetic_energy_total() 
        self.temperature()
        self.total_momentum()


    def run(self, num_frames, animate=False): 
        #a function for progressing the simulation along a for a given number of collsions (frames)
        #and animating the frames
        if animate:
            f = plt.figure()
            ax = plt.axes(xlim=(-self._container._R, self._container._R), \
                          ylim=(-self._container._R, self._container._R))
            ax.add_artist(self._container.patch)
            for i in range(1, len(self.balls)):
                ax.add_patch(self.balls[i].patch) 
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                ax.set_title("Time=%s" %(self._t))
                plt.pause(0.1)
        if animate:
            plt.show()
        
    def kinetic_energy_total(self): 
        #a function for calculating the total kinetic energy of the gas
        KE_total = 0.0
        for i in range(0, len(self.balls)):
            KE_total += self.balls[i].kinetic_energy()
        self._kinetic_energy_total.append(KE_total)

    def temperature(self): 
        #a function for calculating the temperature of the gas
        self._temperature.append(self._kinetic_energy_total[-1]/ \
                                (3*self._N*1.38e-23))      
    
    def total_momentum(self): 
        #a function for calculating the total momentum of the gas
        tot_time_momentum = 0
        for i in range(len(self.balls)):
            tot_time_momentum += self.balls[i].linear_momentum()
        self._total_momentum.append(tot_time_momentum)
    
    def move_all(self, dt): 
        #a function for moving all the balls in the simulation
        for i in range(0, len(self.balls)):
            self.balls[i].move(dt)
    
    def pressure(self): 
        #a function for calculating the pressure of the gas 
        self._pressure_mean_data = []
        for i in range(2, len(self._pressure_evolution[0])): 
            self._pressure_mean_data.append(self._pressure_evolution[0][i])
        self._pressure = sp.mean(self._pressure_mean_data)
        return self._pressure
    
    def distance_from_container(self): 
        #a function for creating a histogram of the distance between each ball and the container
        distancefromcontainer = []
        for i in range(1, self._N + 1):
            distancefromcontainer.append(\
                sp.sqrt(abs(sp.dot(self.balls[0].pos() - \
                self.balls[i].pos(), (self.balls[0].pos() - \
                self.balls[i].pos())))))
        binslistcon = sp.arange(0.0, self.balls[0]._R + 0.5, 1.0)
        plt.hist(distancefromcontainer, bins=binslistcon)
        plt.title("Separation between balls and container")
        plt.xlabel("Distance between each balls and the container (m)")
        plt.ylabel("Number of balls that are this distance away")
        #plt.savefig("Fig.1.png") #a command to save the figure
        plt.show()
   
    def distance_between_balls(self):
        #a function for creating a histogram of the distance between each ball
        distancebetweenballs = []
        for i in range(1, self._N + 1):
            for j in range(1, self._N + 1):
                if i != j:
                    distancebetweenballs.append(\
                        sp.sqrt(abs(sp.dot((self.balls[i].pos() - \
                        self.balls[j].pos()), (self.balls[i].pos() - \
                        self.balls[j].pos())))))
        binslistballs = sp.arange(0.0, 2*self.balls[0]._R + 0.5, 1.0)
        plt.hist(distancebetweenballs, bins=binslistballs)
        plt.title("Separation between balls")
        plt.xlabel("Distance between each ball (m)")
        plt.ylabel("Number of balls that are this distance away")
        #plt.savefig("Fig.2.png") #a command to save the figure
        plt.show()
        
    def velocity_distribution(self, plotmaxboltz=False):
        #a function for creating a histogram of the speed of the molecules in the gas
        #and if plotmaxboltz==True comparing this to the Maxwell-Boltzman distribution
        velocities = []
        mass = []
        pdf = []
        a = 3.5
        b = -3
        for i in range(1, self._N + 1):
            velocities.append(sp.sqrt(sp.dot(self.balls[i].vel(), \
                self.balls[i].vel())))
            mass.append(self.balls[i]._m)
        maxwellboltzman = sp.linspace(0.0, max(velocities), len(velocities))
        for i in range(0, len(velocities)):
            pdf.append(a*(maxwellboltzman[i])*sp.exp(((\
                       -0.5*mass[i]*(maxwellboltzman[i] + b)**2)/ \
                       (1.38e-23*self._temperature[-1]))))
        binslistvel = sp.arange(0.0, max(velocities), 1.5)
        plt.hist(velocities, bins=binslistvel)
        if plotmaxboltz:
            plt.plot(maxwellboltzman, pdf, c='orange')
        plt.title("Distribution of Velocities of the balls")
        plt.xlabel("Velocities ($m\cdot s^{-1}$)")
        plt.ylabel("Number of balls that have this velocity away")
        #plt.savefig("Fig.7.png") #a command to save the figure
        plt.show()
    
    def kinetic_energy_conservation(self):
        # a function for ploting kinetic energy through time of the simulation
        zeropoint = [0.0]
        plt.plot(self._t_list, self._kinetic_energy_total)
        plt.scatter(zeropoint, zeropoint, color='white')
        plt.title("Total Kinetic Energy of the Gas throughout time")
        plt.xlabel("Time (s)")
        plt.ylabel("Total Kinetic Energy of System (J)" )
        #plt.savefig("Fig.3.png") #a command to save the figure
        plt.show()
        
    def momentum_conservation(self):
        # a function for ploting the components of momentum through time of the simulation
        zeropoint = [0.0]
        xcomponent = []
        ycomponent = []
        for i in range(0, len(self._total_momentum)):
            xcomponent.append(self._total_momentum[i][0])
            ycomponent.append(self._total_momentum[i][1])
        plt.plot(self._t_list, xcomponent, c='red', \
                 label='X-Component of Momentum')
        plt.plot(self._t_list, ycomponent, c='blue', \
                 label='Y-Component of Momentum')
        plt.scatter(zeropoint, zeropoint, color='white')
        plt.title("Total Momentum of the Gas throughout time")
        plt.xlabel("Time (s)")
        plt.ylabel("Total Momentum of System ($kg m s^{-1}$)" )
        plt.legend()
        #plt.savefig("Fig.4.png") #a command to save the figure
        plt.show()
