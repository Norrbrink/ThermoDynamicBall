"""
A module for testing parameters of the gas and producing plots
"""
#%% importing modules
import matplotlib.pyplot as plt
import scipy as sp
from thermodynamic_simulation import Simulation
#%% #Parameters for plotting graphs
params = {
            'axes.labelsize': 24,
            'font.size': 24,
            'legend.fontsize': 22,
            'xtick.labelsize': 22,
            'ytick.labelsize': 22,
            'figure.figsize': [11, 11/1.618],
            'font.family': 'Times New Roman'
            } 
plt.rcParams.update(params)
 
#%% #An example simulation, used for most of the important figures
#Run this block to be able to make figures 1, 2, 3, 4 and 7 
#(can experiment with values)
Sim = Simulation(1.0, 1.0, 500, 40.0)
Sim.run(1000, False)

#%% Task 11, Investigating Pressure as temperature is changed
#Used for figure 5 and 8, but takes several hours to run
stds = [1, 5, 10, 15, 20]
Molecules_stds = []
for i in range(0, len(stds)):
    Molecules_stds.append(Simulation(1.0, 1.0, 100, 20.0, sigma=stds[i]))
for i in range(0, len(stds)):
    Molecules_stds[i].run(1000, False)

stds_pressure = []
stds_temperature = []
for i in range(0, len(Molecules_stds)):
    stds_pressure.append(Molecules_stds[i].pressure())
    stds_temperature.append(Molecules_stds[i]._temperature[-1])




#%% Task 11 Testing how the volume of the container affects the simulation
#Does not need to be run
radius = [35, 40, 45, 50]
Molecules_radius = []
for i in range(0, len(radius)):
    Molecules_radius.append(Simulation(1.0, 1.0, 300, radius[i], 40.0))
for i in range(0, len(radius)):
    Molecules_radius[i].run(1000, False)

radius_pressure = []
radius_temperature = []
for i in range(0, len(Molecules_radius)):
    radius_pressure.append(Molecules_radius[i].pressure())
    radius_temperature.append(Molecules_radius[i]._temperature[-1])

#%% Task 
#plt.plot(radius, radius_pressure)
#plt.show()
#%% 
#plt.plot(radius, radius_temperature)
#plt.show()


#%% Task 11 Testing how the number of molecules affects the simulation
#does not need to be run
molecules = [100, 200, 300, 400, 500]
Molecules = []
for i in range(0, len(molecules)):
    Molecules.append(Simulation(1.0, 1.0, molecules[i], 40.0))
for i in range(0, len(Molecules)):
    Molecules[i].run(1000, False)

molecule_pressure = []
molecule_temperature = []
for i in range(0, len(Molecules)):
    molecule_pressure.append(Molecules[i].pressure())
    molecule_temperature.append(Molecules[i]._temperature[-1])

#%%
#plt.plot(molecules, molecule_pressure)
#plt.show()
#%%
#plt.plot(molecules, molecule_temperature)
#plt.show()

#%% Task 9 Fig. 1
Sim.distance_from_container()

#%% Task 9 Fig. 2
Sim.distance_between_balls()

#%% Task 11 Fig 3
Sim.kinetic_energy_conservation()

#%% Task 11 Fig 4
Sim.momentum_conservation()

#%% Task 12 Fig. 5
TempPressLine = sp.polyfit(stds_temperature, stds_pressure, 1)
x = sp.linspace(0, max(stds_temperature), 1000)
y = TempPressLine[0]*x + TempPressLine[1]
plt.scatter(stds_temperature, stds_pressure, c='r', label= 'Measured points')
plt.plot(x, y, label='Line of best fit, Gradient: %s' %TempPressLine[0]) 
plt.title("Pressure dependence on varying Temperature")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (Pa)" )
plt.legend()
#plt.savefig("Fig.5.png") #a command to save the figure
plt.show()

#%% Task 12 Testing how the relative ball size affects the equation of state  
#This block will take several hours to run as it is running 9 simulations 
#of 300 balls for 1000 collisions
#used for Figure 6
stds = [1, 15, 20]
radius_stds = [1, 1e-2, 1e-4]
Molecules_radius_stds = [[], [], []]
for i in range(0, len(radius_stds)):
    for j in range(0, len(stds)):
        Molecules_radius_stds[i].append(Simulation(1.0, radius_stds[i], 300, \
            40.0, sigma=stds[j]))
for i in range(0, len(radius_stds)):
    for j in range(0, len(stds)):
        Molecules_radius_stds[i][j].run(1000, False)

radius_stds_pressure = [[], [], []]
radius_stds_temperature = [[], [], []]
for i in range(0, len(radius_stds)):
    for j in range(0, len(stds)):
        radius_stds_pressure[i].append(Molecules_radius_stds[i][j].pressure())
        radius_stds_temperature[i].append( \
            Molecules_radius_stds[i][j]._temperature[-1])

#%% Task 12 Fig 6
radius_stds_line0 = sp.polyfit(radius_stds_temperature[0], \
    radius_stds_pressure[0], 1)
radius_stds_line1 = sp.polyfit(radius_stds_temperature[1], \
    radius_stds_pressure[1], 1)
radius_stds_line2 = sp.polyfit(radius_stds_temperature[2], \
    radius_stds_pressure[2], 1)
x0 = sp.linspace(0, max(radius_stds_temperature[0]), 1000)
x1 = sp.linspace(0, max(radius_stds_temperature[1]), 1000)
x2 = sp.linspace(0, max(radius_stds_temperature[2]), 1000)
y0 = radius_stds_line0[0]*x0 + radius_stds_line0[1]
y1 = radius_stds_line1[0]*x1 + radius_stds_line1[1]
y2 = radius_stds_line2[0]*x2 + radius_stds_line2[1]
plt.scatter(radius_stds_temperature[0], radius_stds_pressure[0], c='r', \
    label= 'Ball radius: %s m' %radius_stds[0])
plt.scatter(radius_stds_temperature[1], radius_stds_pressure[1], c='b', \
    label= 'Ball radius: %s m' %radius_stds[1])
plt.scatter(radius_stds_temperature[2], radius_stds_pressure[2], c='g', \
    label= 'Ball radius: %s m' %radius_stds[2])
plt.plot(x0, y0, label='$k_B$: %s $J K^{-1}$' %radius_stds_line0[0], c='r')
plt.plot(x1, y1, label='$k_B$: %s $J K^{-1}$' %radius_stds_line1[0], c='b')
plt.plot(x2, y2, label='$k_B$: %s $J K^{-1}$' %radius_stds_line2[0], c='g')
plt.title("Equation of State Dependence on Ball radius")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (Pa)" )
plt.legend()
#plt.savefig("Fig.6.png") #a command to save the figure
plt.show()




#%% Task 13 Fig 7
Sim.velocity_distribution()

#%% Task 14 Fig 8
# a = 0
# P(V-Nb)=NKbT
# P=(NKb/V)T+N/Vb
TempPressLine = sp.polyfit(stds_temperature, stds_pressure, 1)
x = sp.linspace(0, max(stds_temperature), 1000)
y = TempPressLine[0]*x + TempPressLine[1]
b = TempPressLine[1]*((sp.pi*Molecules_stds[0]._container._R**2)/ \
    (Molecules_stds[0]._N))
plt.scatter(stds_temperature, stds_pressure, c='r', label= 'Measured points')
plt.plot(x, y, label='Line of best fit, b is equal to: %s' %b)
plt.title("Wan der Vaals equation fit")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (Pa)" )
plt.legend()
#plt.savefig("Fig.8.png") #a command to save the figure
plt.show()
