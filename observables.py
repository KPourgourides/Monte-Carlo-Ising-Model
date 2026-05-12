"""
TU Delft Computational Physics - Project 2: Ising Model
Authors: Kyproula Mitsidi, Konstantinos Pourgourides
April - May 2026

~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

This module hosts all the functions that are necessary for the calculation and plotting of
observables, as well as their errors (when applicable).
"""


import ising_model as ismo
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import random as rand


def block_bootstrap(A, n, tau=False, resample=True):
    """
    Implements the block bootsrtap method.
    
    Parameters
    ----------
    A : np.array
        Values of observable over time
    n : int
        Repetitions of resampling method
    tau : float
        Exponential fit parameter of the autocorrelation function, if known (default False)
    resample : bool
         Determines whether the method is performed for 1 or n repetitions
         
    Returns
    --------
    if resample is False:
    
    avgA_per_n[i] : float
        Average of observable for 1 resampling
    avgAsqr_per_n[i] : float
        Average of squared observable for 1 resampling
        
    ---------
    
    if resample is True:
    
    avgA : float
        Average of observable over n resamplings
    sigma : float
        Estimated error of the observable
    """

    if tau == False:
        tau = int(autocorrelation_function(magnetization[ismo.N0:], fit_until=100)[2])
    A = A[::tau]
    avgA_per_n = np.zeros((n))
    avgAsqr_per_n = np.zeros((n))
    
    for i in range(n):
        sample = np.array(rand.choices(A, k=len(A)))
        avgA_per_n[i]  = np.mean(sample)
        avgAsqr_per_n[i] = np.mean(sample ** 2)
        if not resample: 

            return avgA_per_n[i], avgAsqr_per_n[i]
                
    avgA = np.mean(avgA_per_n)
    avgAsqr = np.mean(avgAsqr_per_n)
    sigma = np.sqrt(avgAsqr - avgA ** 2)

    return avgA, sigma


def autocorrelation_function(A, fit_until, show_plot=False, save_plot=False):
    """
    Applies the autocorrelation function method to calculate errors of time-averages.

    Parameters
    ----------
    A : np.ndarray
        Any time-series we want to apply the autocorrelation function to
    fit_until : int
        The time step until calculations are carried out and fit is performed (well-behaved part)
    show_plot : bool
        Show autocorrelation function fit to find chi or not (default False)
    save_plot : str,bool
        Name of plot in case of save (default False)

    Returns
    -------
    A_avg : float
        Time average of observable A
    A_error : float
        Error of the time average of observable A
    tau : float
        Exponential fit parameter of the autocorrelation function
    """
    time_grid = np.arange(len(A))
    
    def model(t, tau):
        """
        Returns the fit equation  exp(-t/tau) for autocorrelation function.

        Parameters
        ----------
        t : np.ndarray
            Time grid
        tau : float
            Exponential fit parameter

        Returns
        -------
        equation : np.ndarray
            Exponential fit function
        """
        equation = np.exp(-t / tau)
        return equation
        
    def fit(chi):
        """
        Performs the exponential fit on the autocorrelation function.

        Parameters
        ----------
        chi : np.ndarray
            The autocorrelation function with shape(num_tsteps)

        Returns
        -------
        tau : float
            The exponential fit parameter
        """

        params, covariance = curve_fit(model, time_grid[:fit_until], chi, p0=[20], bounds=(0.9, np.inf))
        if show_plot:
            fig,ax = plt.subplots(figsize=(10,5), dpi=300)
            #plt.title('Autocorrelation Function')
            plt.xlabel('timestep')
            plt.ylabel(r'$\chi_A(x)$')
            plt.plot(time_grid[:fit_until], chi, label=r'$\chi_A(x)$')
            plt.plot(time_grid[:fit_until], model(time_grid[:fit_until],*params),linestyle='--',label=rf'$exp(-t/ \tau) $, $\tau$ ={params[-1]:.2f}')
            plt.legend()
            plt.tight_layout()
            if save_plot!=False:
                fig.savefig(f'journal_plots/{save_plot}.png', dpi=300, bbox_inches='tight')
            plt.show(close=True)
        tau = params[-1]

        return tau

    chi = np.zeros((fit_until))
    N = len(A)
    for i,t in enumerate(time_grid[:fit_until]):
        a = (N-i)*sum(A[0:N-i] * A[i:N])
        b = np.sum(A[0:N-i]) * np.sum(A[i:N])
        c = (N-i) * np.sum(A[0:N-i] ** 2) - (np.sum(A[0:N-i])) ** 2
        d = (N-i) * np.sum(A[i:N] ** 2) - (np.sum(A[i:N])) ** 2
        if np.sqrt(c) * np.sqrt(d) != 0:
            chi[i] = (a - b)/(np.sqrt(c) * np.sqrt(d))
        else:
            chi[i] = 0
    tau = int(np.round(fit(chi), decimals=0))
    A_error = np.sqrt(2 * tau / N) * np.sqrt(np.mean(A ** 2) - (np.mean(A)) ** 2)
    A_avg = np.mean(A)

    return A_avg, A_error, tau


def get_magnetization(spin_grid):
    """
    Computes the average magnetization of a spin lattice.

    Parameters
    ----------
    spin_grid : np.ndarray
        Two-dimensional array representing the spin configuration.

    Returns
    -------
    m : float
        Average magnetization per lattice site.
    """
    m = (1 / spin_grid.size)*np.sum(spin_grid)
    return m

    
def get_avg_abs_magnetization(magnetization):
    """
    Computes the time-averaged absolute magnetization and its statistical error
    using the autocorrelation function method.

    Parameters
    ----------
    magnetization : np.ndarray
        Time series of magnetization values.

    Returns
    -------
    avg_abs_magn : float
        Mean absolute magnetization after equilibration.
    err_avg_abs_magn : float
        Estimated statistical error of the mean absolute magnetization.
    """
    avg_abs_magn, err_avg_abs_magn, tau = autocorrelation_function(abs(magnetization[ismo.N0:]), fit_until=100)
    return avg_abs_magn, err_avg_abs_magn

    
def get_avg_magnetization(magnetization):
    """
    Computes the time-averaged magnetization and its statistical error
    using the autocorrelation function method.

    Parameters
    ----------
    magnetization : np.ndarray
        Time series of magnetization values.

    Returns
    -------
    avg_magn : float
        Mean magnetization after equilibration.
    err_avg_magn : float
        Estimated statistical error of the mean magnetization.
    """
    avg_magn, err_avg_magn, tau = autocorrelation_function(magnetization[ismo.N0:], fit_until=100)
    return avg_magn, err_avg_magn

    
def get_energy(spin_grid):
    """
    Computes the total energy of a 2D Ising spin configuration with periodic
    boundary conditions.

    Parameters
    ----------
    spin_grid : np.ndarray
        Two-dimensional array representing the spin configuration.

    Returns
    -------
    h_tot : float
        Total energy of the spin configuration.
    """
    h_mag = -ismo.B * np.sum(spin_grid)
    h_int = -ismo.J * (np.sum(spin_grid * np.roll(spin_grid, -1, axis=0)) + np.sum(spin_grid * np.roll(spin_grid, -1, axis=1)))
    h_tot = h_int + h_mag
    return h_tot

    
def get_avg_energy(energy):
    """
    Computes the time-averaged energy and its statistical error
    using the autocorrelation function method.

    Parameters
    ----------
    energy : np.ndarray
        Time series of energy measurements from the simulation.

    Returns
    -------
    avg_energy : float
        Mean energy after equilibration.
    err_avg_energy : float
        Estimated statistical error of the mean energy.
    """
    avg_energy, err_avg_energy, tau = autocorrelation_function(energy[ismo.N0:], fit_until=100)
    return avg_energy, err_avg_energy


def get_avg_susceptibility(magnetization):
    """
    Computes the magnetic susceptibility and its statistical error
    using block bootstrap resampling.

    Parameters
    ----------
    magnetization : np.ndarray
        Time series of magnetization values from the simulation.

    Returns
    -------
    avg_chi : float
        Mean magnetic susceptibility estimate.
    error_chi : float
        Standard deviation of susceptibility estimates from bootstrap resampling.
    """
    
    n = 1500
    chi_samples = np.zeros((n))
    tau = int(autocorrelation_function(magnetization[ismo.N0:], fit_until=100)[2])
    
    #Bootstrap
    for i in range(n):
        avgA_per_n,  avgAsqr_per_n = block_bootstrap(magnetization[ismo.N0:], n, tau, resample=False)
        var =  avgAsqr_per_n - avgA_per_n ** 2  
        chi_samples[i] = (ismo.BETA * ismo.L ** 2) * var
    
    avg_chi = np.mean(chi_samples)
    error_chi = np.std(chi_samples)
    
    return avg_chi, error_chi


def get_avg_specific_heat(energy):
    """
    Computes the specific heat and its statistical error using
    block bootstrap resampling.

    Parameters
    ----------
    energy : np.ndarray
        Time series of energy measurements from the simulation.

    Returns
    -------
    avg_c : float
        Mean specific heat estimate.
    error_c : float
        Standard deviation of specific heat estimates from bootstrap resampling.
    """
    n = 1500
    c_samples = np.zeros((n))
    tau = int(autocorrelation_function(energy[ismo.N0:], fit_until =100)[2])

    #Bootstrap
    for i in range(n):
        avgA_per_n,  avgAsqr_per_n = block_bootstrap(energy[ismo.N0:], n, tau, resample=False)
        var =  avgAsqr_per_n - avgA_per_n ** 2  
        c_samples[i] = (ismo.BETA ** 2 / (ismo.L ** 2))*var
    
    avg_c = np.mean(c_samples)
    error_c = np.std(c_samples)
    
    return avg_c, error_c

    
def get_obs_vs_temperature(cycles=17000, L=32, B=0, T_init=2, T_final=3.1, step=0.05):
    """
    Computes thermodynamic observables of a 2D Ising model as a function
    of temperature using Monte Carlo simulations.

    Parameters
    ----------
    cycles : int, optional
        Number of Monte Carlo cycles per temperature point.
    L : int, optional
        Linear size of the spin lattice.
    B : float, optional
        External magnetic field strength.
    T_init : float, optional
        Starting temperature of the sweep.
    T_final : float, optional
        Final temperature of the sweep (exclusive).
    step : float, optional
        Temperature increment.

    Returns
    -------
    temperature_grid : np.ndarray
        Array of temperature values used in the simulation.
    obs_m : list of np.ndarray
        [average absolute magnetization, error on absolute magnetization]
    obs_chi : list of np.ndarray
        [magnetic susceptibility, error on susceptibility]
    obs_c : list of np.ndarray
        [specific heat, error on specific heat]
    obs_e : list of np.ndarray
        [average energy, error on energy]
    """

    temperature_grid = np.arange(T_init, T_final, step)

    absmag, err_absmag = np.zeros((temperature_grid.size)), np.zeros((temperature_grid.size))
    susceptibility, err_susceptibility = np.zeros((temperature_grid.size)), np.zeros((temperature_grid.size))
    specific_heat, err_specific_heat = np.zeros((temperature_grid.size)), np.zeros((temperature_grid.size))
    avg_energy, err_avg_energy = np.zeros((temperature_grid.size)), np.zeros((temperature_grid.size))

    for t,value in enumerate(temperature_grid):
    
        ismo.declare_simulation_constants(side_length=L, magnetic_field=B, temperature=value)
        magnetization, energy, time_grid = ismo.simulate_numpyfied(cycles=cycles, initial_configuration='random')

        abs_magnetization, err_abs_magnetization = get_avg_abs_magnetization(magnetization)
        chi, err_chi = get_avg_susceptibility(magnetization)
        c, err_c = get_avg_specific_heat(energy)
        avg_e, err_avg_e = get_avg_energy(energy)
    
        absmag[t], err_absmag[t] = abs_magnetization, err_abs_magnetization
        susceptibility[t], err_susceptibility[t] = chi, err_chi
        specific_heat[t], err_specific_heat[t] = c, err_c
        avg_energy[t], err_avg_energy[t] = avg_e, err_avg_e 

        obs_m = [absmag, err_absmag]
        obs_chi = [susceptibility, err_susceptibility]
        obs_c = [specific_heat, err_specific_heat]
        obs_e = [avg_energy, err_avg_energy]
        
        print(f"\rStatus: {(100*(t+1)/temperature_grid.size):.1f}%", end="", flush=True)
        
    return temperature_grid, obs_m, obs_chi, obs_c, obs_e


def plot_obs_vs_temperature(temperature_grid, obs_m, obs_chi, obs_c, obs_e, save_plot=False):
    """
    Plots thermodynamic observables of a 2D Ising model as a function of temperature.
    Error bars are included for all observables. A vertical line marks the
    critical temperature of the system.

    Parameters
    ----------
    temperature_grid : np.ndarray
        Array of temperature values used in the simulation.
    obs_m : list of np.ndarray
        [average absolute magnetization, error on absolute magnetization]
    obs_chi : list of np.ndarray
        [magnetic susceptibility, error on susceptibility]
    obs_c : list of np.ndarray
        [specific heat, error on specific heat]
    obs_e : list of np.ndarray
        [average energy, error on energy]
    save_plot : bool
        Saves the plot if True.

    Returns
    -------
    None
        Displays a matplotlib figure and does not return a value.
    """

    absmag, err_absmag = obs_m [0], obs_m[1]
    susceptibility, err_susceptibility = obs_chi[0], obs_chi[1]
    specific_heat, err_specific_heat = obs_c[0], obs_c[1]
    avg_energy, err_avg_energy = obs_e[0], obs_e[1]
    
    fig,ax = plt.subplots(nrows=2, ncols=2, figsize=(10,6), dpi=300)
    
    ax[0,0].errorbar(temperature_grid, absmag, yerr=err_absmag, fmt='o', markersize=3.5, capsize=1, color='red', label=r'$\langle |m| \rangle$')
    ax[0,0].set_ylabel(r'$\langle |m| \rangle$')
    ax[0,0].set_yticks(np.arange(0, 1.2, 0.2))
   
    ax[1,0].errorbar(temperature_grid, susceptibility, yerr=err_susceptibility, fmt='o', markersize=3.5, capsize=1, color='blue', 
                     label=r'$\langle \chi \rangle$')
    ax[1,0].set_xlabel(r'Temperature $(\mathcal{J}/k_B)$')
    ax[1,0].set_ylabel(r'$\langle \chi \rangle$ $(1/\mathcal{J})$')

    ax[1,1].errorbar(temperature_grid, specific_heat, yerr=err_specific_heat, fmt='o', markersize=3.5, capsize=1, color='orange', 
                     label=r'$\langle c \rangle$')
    ax[1,1].set_xlabel(r'Temperature $(\mathcal{J}/k_B)$')
    ax[1,1].set_ylabel(r'$\langle c \rangle$ $(k_B)$')
    
    ax[0,1].errorbar(temperature_grid, avg_energy / ismo.L ** 2, yerr=err_avg_energy / ismo.L ** 2, fmt='o', markersize=3.5, capsize=1, color='magenta',
                     label=r'$\langle E \rangle / N$')
    ax[0,1].set_ylabel(r'$\langle E \rangle$ $(\mathcal{J})$')
    
    for axis in ax.flat:
        axis.axvline(x=ismo.TEMPER_C, linestyle='--', color='gray', label = r'$T_c$', alpha=0.3)
        axis.legend()
        
    fig.suptitle(f'Simulation for grid size {ismo.L}x{ismo.L}')
    plt.tight_layout()
    if save_plot:
        fig.savefig(f'journal_plots/physical_observables_vs_temperature.png', dpi=400, bbox_inches='tight')
    plt.show(close=True)

    
def get_magnetization_vs_field(cycles=17000, L=32, T=1, B_init=-1, B_final=1.01, step=0.01):
    """
    Computes magnetization of a 2D Ising model as a function of external
    magnetic field using Monte Carlo simulations.

    Parameters
    ----------
    cycles : int, optional
        Number of Monte Carlo cycles per field value.
    L : int, optional
        Linear size of the spin lattice.
    T : float, optional
        Temperature of the system.
    B_init : float, optional
        Initial magnetic field value.
    B_final : float, optional
        Final magnetic field value (exclusive).
    step : float, optional
        Increment in magnetic field.

    Returns
    -------
    field_grid : np.ndarray
        Array of magnetic field values.
    absmag : list
        Average absolute magnetization for each field.
    absmagerr : list
        Statistical error of the absolute magnetization.
    mag : list
        Average magnetization for each field.
    magerr : list
        Statistical error of the magnetization.
    """
    
    field_grid = np.arange(B_init, B_final, step)

    mag=[]
    magerr=[]
    absmag=[]
    absmagerr=[]
    
    for i,b in enumerate(field_grid):

        ismo.declare_simulation_constants(side_length=L, magnetic_field=b, temperature = T)
        magnetization, energy, time_grid = ismo.simulate_numpyfied(cycles=cycles, initial_configuration='random')
        
        absmag_temp, absmag_err_temp = get_avg_abs_magnetization(magnetization)
        mag_temp, mag_err_temp = get_avg_magnetization(magnetization)

        absmag.append(absmag_temp)
        absmagerr.append(absmag_err_temp)

        mag.append(mag_temp)
        magerr.append(mag_err_temp)

        print(f"\rStatus: {(100 * (i + 1) / field_grid.size):.1f}%", end="", flush=True)

    return field_grid, absmag, absmagerr, mag, magerr


def plot_magnetization_vs_field(field_grid, absmag, absmagerr, mag, magerr, save_plot=False):
    """
    Plots magnetization and absolute magnetization of a 2D Ising model
    as a function of external magnetic field.

    Parameters
    ----------
    field_grid : np.ndarray
        Array of magnetic field values.
    absmag : array-like
        Average absolute magnetization for each field value.
    absmagerr : array-like
        Statistical error of the absolute magnetization.
    mag : array-like
        Average magnetization for each field value.
    magerr : array-like
        Statistical error of the magnetization.
    save_plot : bool
        Saves plot if True.

    Returns
    -------
    None
        Displays a matplotlib figure.
    """

    fig,ax = plt.subplots(2, 1, figsize=(10, 6), dpi= 300) 
        
    ax[0].set_ylabel(r"$\langle m \rangle$")
    ax[0].set_yticks(np.arange(-1, 1.5, 0.5))
    ax[0].set_xticklabels([])
    ax[0].errorbar(field_grid, mag, yerr = magerr, fmt='ko', capsize=3)
    ax[1].errorbar(field_grid, absmag, yerr = absmagerr, fmt='ko', capsize=3)
    ax[1].set_ylabel(r"$\langle |m| \rangle$")
    ax[1].set_xlabel("H")
    ax[1].set_yticks(np.arange(0, 1.25, 0.25))

    fig.suptitle(f'Simulation for grid size {ismo.L}x{ismo.L} and T = {ismo.TEMPER} J/kB')
    
    plt.tight_layout()
    if save_plot:
        fig.savefig(f'journal_plots/magnetization_vs_magnetic_field.png', dpi=400, bbox_inches='tight')
    plt.show()