#!/usr/bin/env python
"""
Plots for Shrimp Project
Rebecca Li 2019
"""
import matplotlib.pyplot as plt
import numpy as np

def plot4(name, a, b, c, d, ts):
    """Plot 3 values a,b,c, and their associated times"""

    def plotState(state, ylab):
        """ Helper fun """
        plt.plot(ts, state, "-")
        plt.ylabel(ylab)

    plt.figure()
    ax = plt.subplot(411)
    plotState(a[0], a[1])
    plt.title(name)
    plt.subplot(412, sharex=ax)
    plotState(b[0], b[1])
    plt.subplot(413, sharex=ax)
    plotState(c[0], c[1])
    plt.subplot(414, sharex=ax)
    plotState(d[0], d[1])
    plt.xlabel("Time [s]")


def plot3(name, a, b, c, ts):
    """Plot 3 values a,b,c, and their associated times"""

    def plotState(state, ylab):
        """ Helper fun """
        plt.plot(ts, state, "-")
        plt.ylabel(ylab)

    plt.figure()
    ax = plt.subplot(311)
    plotState(a[0], a[1])
    plt.title(name)
    plt.subplot(312, sharex=ax)
    plotState(b[0], b[1])
    plt.subplot(313, sharex=ax)
    plotState(c[0], c[1])
    plt.xlabel("Time [s]")


def plotPositions(odeStates, times):
    """Plots Positions"""
    x_b2w_w = odeStates[:, 0]
    y_b2w_w = odeStates[:, 1]
    z_b2w_w = odeStates[:, 2]

    a = (x_b2w_w, 'x [m]')
    b = (y_b2w_w, 'y [m]')
    c = (z_b2w_w, 'z [m]')
    name = "Position"
    plot3(name, a, b, c, times)


def plotVelocities(odeStates, times):
    """Plots Velocities"""
    vx_b2w_w = odeStates[:, 3]
    vy_b2w_w = odeStates[:, 4]
    vz_b2w_w = odeStates[:, 5]

    a = (vx_b2w_w, 'vel x [m/s]')
    b = (vy_b2w_w, 'vel y [m/s]')
    c = (vz_b2w_w, 'vel z [m/s]')
    name = "Velocity"
    plot3(name, a, b, c, times)


def plotEuler(odeStates, times):
    """Plots euler angles """
    yaw = odeStates[:, 6]
    roll = odeStates[:, 7]
    pitch = odeStates[:, 8]

    a = (roll, 'roll [rad]')
    b = (pitch, 'pitch [rad]')
    c = (yaw, 'yaw [rad]')
    name = "ZYX Euler Angles"
    plot3(name, a, b, c, times)


def plotAngVel(odeStates, times):
    """Plots angular velocities """
    p = odeStates[:, 9]
    q = odeStates[:, 10]
    r = odeStates[:, 11]

    a = (p, 'p [rad/s]')
    b = (q, 'q [rad/s]')
    c = (r, 'r [rad/s]')
    name = "Body Angular Velocities"
    plot3(name, a, b, c, times)


def plotYaws(odeStates, times):
    """Plots angular velocities """
    yawDot_f2b = odeStates[:, 12]
    yawDot_b2p = odeStates[:, 13]
    yaw_f2b = odeStates[:, 14]
    yaw_b2p = odeStates[:, 15]

    a = (yaw_f2b, 'yaw_f2b  [rad/s]')
    b = (yawDot_f2b, 'yawDot_f2b  [rad/s]')
    c = (yaw_b2p, 'yaw_b2p  [rad/s]')
    d = (yawDot_b2p, 'yawDot_b2p  [rad/s]')
    name = "Yaw angular position and velocities"
    plot4(name, a, b, c, d, times)

def plotOdeOutputs(outputs):
    """ Plots the outputs dictionary
    """
    # this is hacky, TODO fix
    keys = list(outputs.keys())
    keys.pop(keys.index('times'))
    times = outputs['times']
    for k in keys:
        val = np.array(outputs[k])
        if len(np.shape(val)) == 1:
            plt.figure()
            plt.plot(times, val)
            plt.xlabel('Time [s]')
            plt.ylabel(k)
        else:
            numItemsToPlot = np.shape(val)[1]
            fig, axs = plt.subplots(numItemsToPlot, 1)
            for i, v in enumerate(val.T):
                axs[i].plot(times, v)
                axs[i].set_xlabel('Time [s]')
                axs[i].set_ylabel('index %d' % i)
            axs[0].set_title(k)





def plotOdeStates(odeStates, times, test=False):
    """ Plot all the things.
        Arguments:
            odeStates (np.ndarray): nx16 array of odeState vectors
            times (np.ndarray): nx1 array of times
    """
    plotPositions(odeStates, times)
    plotVelocities(odeStates, times)
    plotEuler(odeStates, times)
    plotAngVel(odeStates, times)
    plotYaws(odeStates, times)
    if not test:
        plt.show()
