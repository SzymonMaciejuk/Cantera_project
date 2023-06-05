import cantera as ct
import numpy as np

import matplotlib.pyplot as plt

#########################################################################
# Input Parameters
#########################################################################

# reaction mechanism, kinetics type and compositions
reaction_mechanism = 'nDodecane_Reitz.yaml'
phase_name = 'nDodecane_IG'
comp_air = 'o2:1, n2:3.76'
comp_fuel = 'c6h12:1' #here we enter fuel type we want to compare

f = 600. / 60.  # engine speed [1/s] (3000 rpm)
V_H = .5e-3  # displaced volume [m**3]
epsilon = 20.  # compression ratio [-]
d_piston = 0.083  # piston diameter [m]

# turbocharger temperature, pressure, and composition
T_inlet = 300.  # K
p_inlet = 1.3e5  # Pa
comp_inlet = comp_air

# outlet pressure
p_outlet = 1.2e5  # Pa

# fuel properties (gaseous!)
T_injector = 300.  # K
p_injector = 1600e5  # Pa
comp_injector = comp_fuel

# ambient properties
T_ambient = 300.  # K
p_ambient = 1e5  # Pa
comp_ambient = comp_air

# Inlet valve friction coefficient, open and close timings
inlet_valve_coeff = 1.e-6
inlet_open = -18. / 180. * np.pi
inlet_close = 198. / 180. * np.pi

# Outlet valve friction coefficient, open and close timings
outlet_valve_coeff = 1.e-6
outlet_open = 522. / 180 * np.pi
outlet_close = 18. / 180. * np.pi

# Fuel mass, injector open and close timings
injector_open = 170. / 180. * np.pi
injector_close = 185. / 180. * np.pi
injector_mass = 3.2e-5  # kg

# Simulation time and parameters
sim_n_revolutions = 2
delta_T_max = 20.
rtol = 1.e-12
atol = 1.e-16

#####################################################################
# Set up IC engine Parameters and Functions
#####################################################################

V_oT = V_H / (epsilon - 1.)
A_piston = .25 * np.pi * d_piston ** 2
stroke = V_H / A_piston


def crank_angle(t):
    """Convert time to crank angle"""
    return np.remainder(2 * np.pi * f * t, 8 * np.pi)


def piston_speed(t):
    """Approximate piston speed with sinusoidal velocity profile"""
    return - stroke / 2 * 2 * np.pi * f * np.sin(crank_angle(t))


#####################################################################
# Set up Reactor Network
#####################################################################

# load reaction mechanism
gas = ct.Solution(reaction_mechanism, phase_name)

# define initial state and set up reactor
gas.TPX = T_inlet, p_inlet, comp_inlet
cyl = ct.IdealGasReactor(gas)
cyl.volume = V_oT

# define inlet state
gas.TPX = T_inlet, p_inlet, comp_inlet
inlet = ct.Reservoir(gas)

# inlet valve
inlet_valve = ct.Valve(inlet, cyl)
inlet_delta = np.mod(inlet_close - inlet_open, 4 * np.pi)
inlet_valve.valve_coeff = inlet_valve_coeff
inlet_valve.set_time_function(
    lambda t: np.mod(crank_angle(t) - inlet_open, 4 * np.pi) < inlet_delta)

# define injector state (gaseous!)
gas.TPX = T_injector, p_injector, comp_injector
injector = ct.Reservoir(gas)

# injector is modeled as a mass flow controller
injector_mfc = ct.MassFlowController(injector, cyl)
injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
injector_t_open = (injector_close - injector_open) / 2. / np.pi / f
injector_mfc.mass_flow_coeff = injector_mass / injector_t_open
injector_mfc.set_time_function(
    lambda t: np.mod(crank_angle(t) - injector_open, 4 * np.pi) < injector_delta)

# define outlet pressure (temperature and composition don't matter)
gas.TPX = T_ambient, p_outlet, comp_ambient
outlet = ct.Reservoir(gas)

# outlet valve
outlet_valve = ct.Valve(cyl, outlet)
outlet_delta = np.mod(outlet_close - outlet_open, 4 * np.pi)
outlet_valve.valve_coeff = outlet_valve_coeff
outlet_valve.set_time_function(
    lambda t: np.mod(crank_angle(t) - outlet_open, 4 * np.pi) < outlet_delta)

# define ambient pressure (temperature and composition don't matter)
gas.TPX = T_ambient, p_ambient, comp_ambient
ambient_air = ct.Reservoir(gas)

# piston is modeled as a moving wall
piston = ct.Wall(ambient_air, cyl)
piston.area = A_piston
piston.set_velocity(piston_speed)

# create a reactor network containing the cylinder and limit advance step
sim = ct.ReactorNet([cyl])
sim.rtol, sim.atol = rtol, atol
cyl.set_advance_limit('temperature', delta_T_max)

#####################################################################
# Run Simulation
#####################################################################

# set up output data arrays
states = ct.SolutionArray(
    cyl.thermo,
    extra=('t', 'ca', 'V', 'dWv_dt'),
)

# simulate with a maximum resolution of 1 deg crank angle
dt = 1. / (360 * f)
t_stop = sim_n_revolutions / f
while sim.time < t_stop:

    # perform time integration
    sim.advance(sim.time + dt)

    # calculate results to be stored
    dWv_dt = - (cyl.thermo.P - ambient_air.thermo.P) * A_piston * \
        piston_speed(sim.time)

    # append output data
    states.append(cyl.thermo.state,
                  t=sim.time, ca=crank_angle(sim.time),
                  V=cyl.volume, dWv_dt=dWv_dt)


#######################################################################
# Plot Results in matplotlib
#######################################################################

def ca_ticks(t):
    """Helper function converts time to rounded crank angle."""
    return np.round(crank_angle(t) * 180 / np.pi, decimals=1)


t = states.t

# pressure and temperature diagram
xlim1 = 330
xlim2 = xlim1 + 60
xticks = np.arange(300, 420, 5)
xminor = np.arange(300, 420, 0.5)
fig, ax = plt.subplots(nrows=2)
ax[0].plot(ca_ticks(t), states.P / 1.e5)
ax[0].set_title('$C_6H_1$$_2$')
ax[0].margins(x=-0.4)
ax[0].set_ylabel('$p$ [bar]')
ax[0].set_xlabel(r'$\phi$ $[\degree]$')
ax[0].set_xticks(xticks)
ax[0].set_xticklabels(xticks)
ax[0].set_xticks(xminor,minor=True)
ax[0].set_yticks(np.arange(0,300,50))
ax[0].set_yticks(np.arange(0,300,10),minor=True)
ax[0].grid(which='major',linewidth=1)
ax[0].grid(which='minor',linewidth=0.3)
ax[0].set_xlim([xlim1,xlim2])
ax[0].set_ylim([0,250])
ax[1].plot(ca_ticks(t), states.T)
ax[1].margins(x=-0.4)
ax[1].set_ylabel('$T$ [K]')
ax[1].set_xlabel(r'$\phi$ $[\degree]$')
ax[1].set_xticks(xticks)
ax[1].set_xticklabels(xticks)
ax[1].set_xticks(xminor,minor=True)
ax[1].set_yticks(np.arange(0,3500,500))
ax[1].set_yticks(np.arange(0,3500,100),minor=True)
ax[1].grid(which='major',linewidth=1)
ax[1].grid(which='minor',linewidth=0.3)
ax[1].set_xlim([xlim1,xlim2])
ax[1].set_ylim([500,3000])
plt.show()

## p-V diagram
#fig, ax = plt.subplots()
#ax.plot(states.V[t > 0.04] * 1000, states.P[t > 0.04] / 1.e5)
#ax.set_xlabel('$V$ [l]')
#ax.set_ylabel('$p$ [bar]')
#plt.show()

## phi_V diagram
#xticks = np.arange(0.05, 0.155, 0.005)
#fig, ax = plt.subplots()
#ax.plot(t, states.V * 1000)
#ax.margins(x=-0.2)
#ax.set_ylabel('$V$ [l]')
#ax.set_xlabel(r'$\phi$ [deg]')
#ax.set_xticks(xticks)
#ax.set_xticklabels(ca_ticks(xticks))
#plt.show()
