"""
A module for the class Ball and its derived class Container
"""

import numpy
import scipy as sp
import matplotlib.pyplot as plt

class Ball():
    """
    A Class for creating Ball objects that will work as a molecule in a gas
    """
    def __init__(self, mass=10.0, radius=1.0, position=[3.0, 4.0], \
                 velocity=[1.0, 0.0]):
        #initialing parameters of the ball
        self._m = mass 
        self._R = radius
        self._r = sp.array(position)
        self._v = sp.array(velocity)
        self.patch = plt.Circle(self._r, self._R, fc='r', ec='black') #visual patch of ball
        self.isContainer = False
        self.bounce_time = 0.0 #used for calculating pressure
        
    def pos(self): 
        #a function for obtaining the hidden variable of the position of the particle
        return self._r
    
    def vel(self): 
        #a function for obtaining the hidden variable of the veloctiy of the particle
        return self._v
    
    def move(self, dt): 
        #a function for moving a particle due to its velocity for a time dt
        self._r += self._v * dt
    
    def time_to_collision(self, other): 
        #a function for determining the time until the given balls collide
        r_relative = self._r - other._r 
        v_relative = self._v - other._v
        if self.isContainer or other.isContainer: #checks for if either of the balls are the container
            R_relative = self._R - other._R
        else:
            R_relative = self._R + other._R
        dotrr = sp.dot(r_relative, r_relative)
        dotrv = sp.dot(r_relative, v_relative)
        dotvv = sp.dot(v_relative, v_relative)
        #δt=(-v.r±sqrt((r.v)^2-v.v*(r.r-R^2)))/v.v
        deltatplus = (-dotrv + sp.sqrt(dotrv**2 - \
            dotvv*(dotrr - R_relative**2)))/(dotvv)
        #δt=(-v.r±sqrt((r.v)^2-v.v*(r.r-R^2)))/v.v
        deltatminus = (-dotrv - sp.sqrt(dotrv**2 - \
            dotvv*(dotrr - R_relative**2)))/(dotvv)  
        if deltatminus > 1e-14 and type(deltatminus) != numpy.complex128:
            return deltatminus  #return the smallest change in time given that it is positive and real
        elif deltatplus > 0.0 and type(deltatplus) != numpy.complex128:
            return deltatplus #return the change in time given that it is positive and real
        else: 
            return 1e10 #returning a large number as in next_collision the smallest 
                        #change in time is being looked for
        
    def collide(self, other): 
        #a function to cause the changes that the collision would change
        r_relative = (self._r - other._r)
        r_relative_norm = r_relative/sp.sqrt(sp.dot(r_relative, r_relative)) 
        r_relative_perp_norm = sp.array([-r_relative_norm[1], \
            r_relative_norm[0]])
        self.v_parr = sp.dot(self._v, r_relative_norm)
        self.v_perp = sp.dot(self._v, r_relative_perp_norm)
        other.v_parr = sp.dot(other._v, r_relative_norm)
        other.v_perp = sp.dot(other._v, r_relative_perp_norm)
        self.k_initial = self.kinetic_energy() #used for checking Energy conservation
        other.k_initial = other.kinetic_energy()
        self.momentum_inital = self.linear_momentum() #used for calculating pressure
        other.momentum_inital = other.linear_momentum() 
        self._v = (((self._m - other._m)/(self._m + other._m))*self.v_parr + \
            ((2*other._m)/(self._m + other._m))*other.v_parr)* \
            r_relative_norm + self.v_perp*r_relative_perp_norm 
        other._v = (((2*self._m)/(self._m + other._m))*self.v_parr + \
            ((other._m - self._m)/(self._m + other._m))*other.v_parr)* \
            r_relative_norm + other.v_perp*r_relative_perp_norm
        self.kinetic_energy()
        if (self.k_initial + other.k_initial) - (self.kinetic_energy() + \
            other.kinetic_energy()) > (self.k_initial + other.k_initial)*1e-3: #checking for energy conservation
           raise Exception("Oh no energy is not conserved!")
        if self.isContainer: #checking if either ball is the container to calculate pressure
            self.collided(other)
        if other.isContainer:
            other.collided(self)
                
    def kinetic_energy(self): 
        #a function that returns the kinetic energy of a ball
        return 0.5*self._m*sp.dot(self._v, self._v)
       
    def linear_momentum(self): 
        #a function that returns the linear momentum of a ball
        return self._m*self._v
        
class Container(Ball):
    """
    A Class derived from Ball for creating the container object that will work as the container of the gas
    """
    def __init__(self, m=100000.0, r=10.0): 
        #initialing parameters of the container
        Ball.__init__(self, mass=m, radius=r , position=[0.0, 0.0], \
                      velocity=[0.0, 0.0])
        self.patch = plt.Circle(self._r, self._R, ec='b', fc='None') #visual patch of container
        self.isContainer = True
        self._pressuretime = 0.0
        
    def collided(self, other): 
        #used for calculating pressure
        self.delta_linear_momentum = other.linear_momentum() - \
                other.momentum_inital
        self._pressuretime = sp.sqrt(sp.dot(self.delta_linear_momentum, \
                                self.delta_linear_momentum))/(2*sp.pi*self._R)
          
        
            