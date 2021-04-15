# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 12:26:07 2021

@author: Cameron Cummins
"""

import numpy as np
import matplotlib.pyplot as plt
"""
    The following code is sourced from https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method
"""

# In-line function for second equation
equation_two = lambda theta, mu, yt, dt, sigma, noise : theta * (mu - yt) * dt + sigma * noise
def secondEquation(y_init:float, theta:float, mu:float, sigma:float,
                 noise_center:float = 0.0, num_sims:int = 5, num_points:int = 1000, t_init:float = 0, t_final:float = 5,):
    """
    Parameters
    ----------
    y_init : float
        initial function value
    theta : float, optional
        Value for theta.
    mu : float, optional
        Value for mu.
    sigma : float, optional
        Value for sigma.
    num_sims : int
        number of simulations to run and calculate data for. The default is 5.
    num_points : int
        number of data points to calculate. The default is 1000.
    t_init : float
        initial time vlaue. The default is 0.
    t_final : float
        final time value. The default is 5.
    noise_center : float, optional
        Mean for Gaussian noise generator (numpy.random.normal). The default is 0.0.

    Returns
    -------
    sims : list
        list of simulation results, each containing a dictionary with respective 't' and 'y' values
    """

    # In-line function to generate noise
    noise = lambda : np.random.normal(loc = noise_center, scale = np.sqrt(dt)) #same noise as first equation

    equation_two = lambda theta, mu, yt, dt, sigma, noise : theta * (mu - yt) * dt + sigma * noise

    # Create list for simulations
    sims = []

    # Calculate time delta
    dt = (t_final - t_init) / num_points

    for j in range(num_sims): #going through each sim
        #array for time and y values for this run
        t_vals = []
        y_vals = []

        # Initialize
        t = t_init
        y = y_init
        t_vals.append(t)
        y_vals.append(y)

        #go through each step of the sim
        for i in range(num_points):
            t += dt
            yt = y
            dy = equation_two(theta, mu, yt, dt, sigma, noise()) #current value of y goes into the calculation of the dy
            y = y + dy
            #adding values to the list
            t_vals.append(t)
            y_vals.append(y)
        sims.append({"t_set":t_vals, "y_set":y_vals}) #for each sim, add the t and y values for each to a list of dictionaries
    return sims





def firstEquation(y_init:float, mu:float, sigma:float,
                   num_sims:int = 5, num_points:int = 1000, t_init:float = 0, t_final:float = 5, noise_center:float = 0.0):
    """
    Parameters
    ----------
     y_init : float
        initial function value
    mu : float, optional
        Value for mu.
    sigma : float, optional
        Value for sigma.
    num_sims : int
        number of simulations to run and calculate data for. The default is 5.
    num_points : int
        number of data points to calculate. The default is 1000.
    t_init : float
        initial time vlaue. The default is 0.
    t_final : float
        final time value. The default is 5.
    noise_center : float, optional
        Mean for Gaussian noise generator (numpy.random.normal). The default is 0.0.

    Returns
    -------
    sims : list
        list of simulation results, each containing a dictionary with respective 't' and 'y' values
    """
    # In-line function for first equation
    equation_one = lambda mu, dt, sigma, noise : mu * dt + sigma * noise

    # In-line function to generate noise
    noise = lambda : np.random.normal(loc = noise_center, scale = np.sqrt(dt))

    # Create list for simulations
    sims = []

    # Calculate time delta
    dt = (t_final - t_init) / num_points


    for sim in range(num_sims):
        # Create lists for storing points
        t_values = []
        y_values = []

        # Initialize
        t = t_init
        y = y_init
        t_values.append(t)
        y_values.append(y)

        # Iterate through 'N' to calculate values of 't' and 'y'
        for index in range(num_points):
            t += dt
            y = y + equation_one(mu, dt, sigma, noise())
            t_values.append(t)
            y_values.append(y)
        sims.append({"t_set":t_values, "y_set":y_values})
    return sims
"""
    End of sourced code.
"""

def plotAndOutputResults(simulation_results:list, path:str):
    """
    Parameters
    ----------
    simulation_results : list
        list of results for any number of simulations run
    path : str
        path to output image of figure to

    Returns
    -------
    None.
    """
    figure = plt.figure()
    figure_axis = figure.add_subplot()
    # print("number of sims", len(simulation_results))

    for result in simulation_results:
        figure_axis.plot(result["t_set"], result["y_set"], color = "gray")


    #computing the average of all the simulations
    avg_t = []
    avg_y = []
    #Going through each of the simulaitons y_0 then each y_1 then each y_2... y_N
    #each time computing a mean of all the y_0 and y_1 and y_2...y_N
    for i in range(len(simulation_results[0]["t_set"])):

        avg_t.append(simulation_results[0]["t_set"][i]) #all t_values are the same, so taking the first one abritrarily
        tot_yi = 0
        number_y = len(simulation_results) #total number of sims we will get
        k = 0
        for sim in simulation_results:
            tot_yi += sim["y_set"][i]
        avg_y.append(tot_yi / number_y)

    #plotting the avg graph
    figure_axis.plot(avg_t, avg_y, label = "Average", c = "red")






    figure_axis.set_xlabel("time (s)")
    h = figure_axis.set_ylabel("y")
    h.set_rotation(0)
    figure_axis.legend()


    figure.savefig(path)

    #want to return the avg graph so we can compare
    return avg_t, avg_y

# Example useage using paramters from wikipedia
# secondEquation(y_init, theta, mu, sigma)
# firstEquation(y_init, mu, sigma)
avg_t_y00, avg_y_y00 = plotAndOutputResults(secondEquation(0, 0.7, 1.5, 0.06), "SecondEquationResults.png")
avg_t_y01, avg_y_y01 = plotAndOutputResults(secondEquation(1, 0.7, 1.5, 0.06), "SecondEquationResults1.png")
avg_t_y02, avg_y_y02 = plotAndOutputResults(secondEquation(2, 0.7, 1.5, 0.06), "SecondEquationResults2.png")
plotAndOutputResults(firstEquation(0, 1.5, 0.06), "FirstEquationResults.png")
plotAndOutputResults(firstEquation(1, 1.5, 0.06), "FirstEquationResults1.png")
plotAndOutputResults(firstEquation(2, 1.5, 0.06), "FirstEquationResults2.png")

plt.clf()
figure = plt.figure()
figure_axis = figure.add_subplot()
figure_axis.plot(avg_t_y00, avg_y_y00, label = "y_init = 0")
figure_axis.plot(avg_t_y01, avg_y_y01, label = "y_init = 1")
figure_axis.plot(avg_t_y02, avg_y_y02, label = "y_init = 2")
figure_axis.set_xlabel("time (s)")
h = figure_axis.set_ylabel("y")
h.set_rotation(0)
figure_axis.legend()
figure.savefig("Avg_y_init.png")
