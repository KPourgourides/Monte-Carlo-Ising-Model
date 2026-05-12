"""
TU Delft Computational Physics - Project 2: Ising Model
Authors: Kyproula Mitsidi, Konstantinos Pourgourides
April - May 2026

~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

This module hosts all the functions that are necessary for the Monte Carlo simulation.
The functions are used to initialize the spin grid, and to make
the time evolution of the system according to the Metropolis condition.
"""

import observables as obs
import numpy as np
import matplotlib.pyplot as plt

def declare_simulation_constants(side_length, magnetic_field, temperature, ferromagnetism=True):

    """
    Initializes and assigns the global simulation constants.

    Parameters
    ----------
    side_length : int
        The length of the spin grid.
    magnetic_field : float
        External magnetic field applied to the system.
    temperature : float
        Temperature of the system in reduced units.
    ferromagnetism : bool
        Determines whether the interaction strength J will be positive (True) or negative (False - Antiferromagnetism)

    Returns
    -------
    None
        Sets the global simulation parameters and equilibration step count.
    """

    global L, J, B, BETA, K_B, TEMPER, TEMPER_C, N0
    K_B = 1
    if ferromagnetism:
        J = 1
    else:
        J = -1
    L = side_length
    B = magnetic_field
    BETA = 1/(K_B*temperature)
    TEMPER = temperature
    TEMPER_C = 2.269

    if L <= 16:
            N0 = 3000
    elif L <= 32 and L > 16:
            N0 = 6000
    elif L > 32 and L <= 64:
            N0 = 15000
    else:
            N0 = 20000


def nearest_neighbors(spin_grid, i, j):

    """
    Computes the sum of the nearest neighboring spins for a lattice site. Applies periodic boundary conditions.

    Parameters
    ----------
    spin_grid : np.ndarray
        Two-dimensional array representing the spin configuration.
    i : int
        Row index of the selected spin grid site.
    j : int
        Column index of the selected spin grid site.

    Returns
    -------
    nn : int
        Sum of the four nearest neighboring spins using periodic
        boundary conditions.
    """
    
    nn = (spin_grid[(i + 1) % L, j] 
        + spin_grid[i, (j + 1) % L] 
        + spin_grid[i - 1, j] 
        + spin_grid[i, j - 1])
    return nn
    

def simulate(cycles, initial_configuration, spin_grid_input=0):

    """
    Evolution of the spin grid using a Monte Carlo simulation which utilizes the Metropolis algorithm (Standard algorithm).

    Parameters
    ----------
    cycles : int
        Number of Monte Carlo cycles to simulate.
    initial_configuration : str
        Specifies the starting spin grid configuration. Options are:
        - 'random' : initializes a random spin grid with values -1 and 1
        - 'aligned' : initializes all spins to +1
        - 'custom' : uses the provided spin_grid_input
    spin_grid_input : np.ndarray, optional
        initial spin grid provided by user when initial_configuration is 'custom'.

    Returns
    -------
    magnetization : np.ndarray
        Array of magnetization values at each cycle.
    energy : np.ndarray
        Array of energy values at each cycle.
    time_grid : np.ndarray
        Array of cycle indices corresponding to the output observables.
    """

    #initialization
    if initial_configuration=='custom':
        spin_grid = spin_grid_input
    elif initial_configuration=='random':
        spin_grid = np.random.choice([-1, 1], size=(L, L))
    elif initial_configuration=='aligned':
        spin_grid = np.ones((L, L))
    else:
        raise TypeError("initial_configuration must be 'random' or 'aligned' or 'custom' ")
        
    N = spin_grid.size
    
    magnetization = np.zeros((cycles + 1))
    energy = np.zeros((cycles + 1))
    
    magnetization[0] = obs.get_magnetization(spin_grid)
    energy[0] = obs.get_energy(spin_grid)

    for c in range(cycles):
        
        for step in range(N):

            #randomly choosing spin
            i,j = np.random.randint(0, L, size=2)
            delta_H = 2*spin_grid[i, j]*(J * nearest_neighbors(spin_grid, i, j) + B)
    
            if delta_H < 0 or np.random.random() < np.exp(-BETA * delta_H):
                spin_grid[i, j] *= -1
            
        magnetization[c + 1] = obs.get_magnetization(spin_grid)
        energy[c + 1] = obs.get_energy(spin_grid)


    time_grid = np.arange(0, cycles + 1, 1)

    return magnetization, energy, time_grid


def nearest_neighbors_numpyfied(spin_grid):

    """
    Computes the sum of nearest neighboring spins for every site in the spin grid using a vectorized 
    implementation (utilizes the NumPy library). Applies periodic boundary conditions.

    Parameters
    ----------
    spin_grid : np.ndarray
        Two-dimensional array representing the spin configuration.

    Returns
    -------
    nn_tot : np.ndarray
        Array of the same shape as spin_grid containing the sum of the
        four nearest neighbors for each site.
    """
    
    nn_horizontal = np.roll(spin_grid, -1, axis=0) + np.roll(spin_grid, 1, axis=0)
    nn_vertical = np.roll(spin_grid, -1, axis=1) + np.roll(spin_grid, 1, axis=1)
    nn_tot = nn_horizontal + nn_vertical
    
    return nn_tot

    
def simulate_numpyfied(cycles, initial_configuration, spin_grid_input=0):
    """
    Evolution of the spin grid using a Monte Carlo simulation which utilizes the Metropolis algorithm. 
    The spin grid is updated in a fully NumPy-vectorized manner (odd/even checkerboard decomposition).

    Parameters
    ----------
    cycles : int
        Number of Monte Carlo cycles to simulate.
    initial_configuration : str
        Specifies the starting spin grid configuration. Options are:
        - 'random' : initializes a random spin grid with values -1 and 1
        - 'aligned' : initializes all spins to +1
        - 'custom' : uses the provided spin_grid_input
    spin_grid_input : np.ndarray, optional
        initial spin grid provided by user when initial_configuration is 'custom'.

    Returns
    -------
    magnetization : np.ndarray
        Array of magnetization values at each cycle.
    energy : np.ndarray
        Array of energy values at each cycle.
    time_grid : np.ndarray
        Array of cycle indices corresponding to the output observables.
    """

    if cycles < N0 + 2000:
        raise ValueError(rf"cycles must be at least {N0 + 2000}")
        
    #initialization
    if initial_configuration == 'custom':
        spin_grid = spin_grid_input.copy()
    elif initial_configuration == 'random':
        spin_grid = np.random.choice([-1, 1], size=(L, L))
    elif initial_configuration == 'aligned':
        spin_grid = np.ones((L, L))
    else:
        raise TypeError("initial_configuration must be 'random' or 'aligned' or 'custom' ")
        
    N = spin_grid.size

    magnetization = np.zeros((cycles + 1))
    energy = np.zeros((cycles + 1))
    
    magnetization[0] = obs.get_magnetization(spin_grid)
    energy[0] = obs.get_energy(spin_grid)

    for c in range(cycles):
        
        for oddness in [0, 1]:

            nearest_neighbors = nearest_neighbors_numpyfied(spin_grid)
            delta_H = 2*spin_grid*(J * nearest_neighbors + B)
            
            boltzmann_matrix = np.exp(- BETA * delta_H)
            random_matrix = np.random.rand(L,L)
            
            flip_matrix = (delta_H < 0) | (random_matrix < boltzmann_matrix)
            checkboard_parity = np.sum(np.indices((L, L)),axis = 0)%2 == oddness
            spin_grid[checkboard_parity & flip_matrix] *= -1
        
        magnetization[c + 1] = obs.get_magnetization(spin_grid)
        energy[c + 1] = obs.get_energy(spin_grid)

    time_grid = np.arange(0, cycles + 1, 1)

    return magnetization, energy, time_grid

    
def plot_equilibrium(magnetization, time_grid):

    """
    Plots the time evolution of the absolute magnetization and indicates
    the estimated equilibration point of the system.

    Parameters
    ----------
    magnetization : np.ndarray
        Array of magnetization values over simulation time.
    time_grid : np.ndarray
        Array of Monte Carlo cycle indices corresponding to magnetization values.

    Returns
    -------
    None
        Displays a matplotlib plot.
    """
    plt.figure(figsize=(10, 3),dpi = 100)
    plt.plot(time_grid, abs(magnetization), color='blue', linewidth=0.5, alpha=0.5)
    plt.plot(N0, abs(magnetization[N0]), 'ko', label = '$n_0$')
    plt.ylim(-0.1, 1.1)
    plt.xlabel('MC cycle')
    plt.ylabel(r'$\langle |m| \rangle$')
    plt.legend()
    plt.tight_layout()
    plt.show(close = True)
