#!/usr/bin/env python
"""
Run a simulation for Shrimp Project
Rebecca Li 2019
"""
import numpy as np
from scipy.integrate import odeint
from timeit import default_timer as timer
import matplotlib.pyplot as plt

from mathUtil import rpm2RadiansPerSecond
from plotShrimp import plotOdeStates, plotOdeOutputs
from shrimpOde import flyerOde, _odeOutputs
from shrimpClasses import w3ShrimpParams
from shrimpVisualizer import drawShrimp


def initialOdeState():
    """Set up nice state for us """
    r_w2b_w = np.array([0, 0, 0])
    vel_w2b_w = np.array([0, 0, 0])
    euler_w2f = np.array([0.1, 0.01, 0]) * np.pi / 180.
    angvel_w2f_f = np.array([0, 0, 0])
    yaw_b2p = 0
    yaw_f2b = 0
    yawDot_f2b = rpm2RadiansPerSecond(10)
    yawDot_b2p = rpm2RadiansPerSecond(500)
    odeState = np.zeros(16)
    odeState[0:3] = r_w2b_w
    odeState[3:6] = vel_w2b_w
    odeState[6:9] = euler_w2f
    odeState[9:12] = angvel_w2f_f
    odeState[12] = yawDot_f2b
    odeState[13] = yawDot_b2p
    odeState[14] = yaw_f2b
    odeState[15] = yaw_b2p
    return odeState


def runSimulation(tf=0.3, plot=False, viz=False, test=False):
    """ Run a sim!"""
    y0 = initialOdeState()
    parameters = w3ShrimpParams()
    dt = 0.01
    t = np.arange(0, tf, dt)
    print('Simulating shrimp...')
    startTime = timer()
    states = odeint(flyerOde, y0, t, args=(parameters,))
    endTime = timer()
    print('Time elapsed for simulation: %f' % (endTime - startTime))
    if plot:
        plotOdeStates(states, t, test)
        plotOdeOutputs(_odeOutputs)
        plt.show()
    if viz:
        drawShrimp(parameters, t, states, autoplay=False, gif=test)
    return parameters, states, t


def testRunSim():
    """Test stuff"""
    runSimulation(tf=0.3, plot=True, test=False)


if __name__ == "__main__":
    # TODO: add args
    runSimulation(tf=0.3, plot=True, viz=False)

