#!/usr/bin/env python
"""
Run a simulation for Shrimp Project
Rebecca Li 2019
"""
import numpy as np
from scipy.integrate import odeint

from mathUtil import rpm2RadiansPerSecond
from plotShrimp import plotOdeStates
from shrimpOde import flyerOde
from shrimpClasses import defaultShrimpParams


def initialOdeState():
    """Set up nice state for us """
    r_w2b_w = np.array([0, 0, 1])
    vel_w2b_w = np.array([0, 0, 0])
    euler_w2f = np.array([0.01, 0.01, 0]) * np.pi / 180.
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


def runSimulation(tf=1.0):
    """ Run a sim!"""
    y0 = initialOdeState()
    parameters = defaultShrimpParams()
    dt = 0.01
    t = np.arange(0, tf, dt)
    states = odeint(flyerOde, y0, t, args=(parameters,))
    plotOdeStates(states, t)
    return (states, t)


def testRunSim():
    """Test stuff"""
    runSimulation(tf=0.5)


if __name__ == "__main__":
    runSimulation()
