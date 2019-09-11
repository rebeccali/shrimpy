#!/usr/bin/env python
"""
Plots for Shrimp Project
Rebecca Li 2019
"""
from mathUtil import yawIndex, rollIndex, pitchIndex
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

    plt.figure()
    plt.plot(x_b2w_w, y_b2w_w)
    plt.title('Position XY')
    plt.xlabel('X position [m]')
    plt.xlabel('Y position [m]')


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
    eulerAngles = odeStates[:, 6:9]

    yawDeg = np.array(eulerAngles[:, yawIndex]) * 180. / np.pi
    rollDeg = np.array(eulerAngles[:, rollIndex]) * 180. / np.pi
    pitchDeg = np.array(eulerAngles[:, pitchIndex]) * 180. / np.pi

    a = (rollDeg, 'roll [deg]')
    b = (pitchDeg, 'pitch [deg]')
    c = (yawDeg, 'yaw [deg]')
    name = "ZYX Euler Angles"
    plot3(name, a, b, c, times)


def plotAngVel(odeStates, times):
    """Plots angular velocities """
    pDeg = np.array(odeStates[:, 9]) * 180. / np.pi
    qDeg = np.array(odeStates[:, 10]) * 180. / np.pi
    rDeg = np.array(odeStates[:, 11]) * 180. / np.pi

    a = (pDeg, 'p [deg/s]')
    b = (qDeg, 'q [deg/s]')
    c = (rDeg, 'r [deg/s]')
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


def plotXyzType(outputs, typeName):
    """Plot all outputs that contain the xyz typeName, like forces"""
    times = outputs['times']
    keys = list(outputs.keys())
    typeKeys = [k for k in keys if typeName in k]
    axesNames = 'xyz'
    fig, axs = plt.subplots(3, 1)
    for (i, axisName) in enumerate(axesNames):
        for k in typeKeys:
            print('plotting %s' % k)
            axs[i].plot(times, np.array(outputs[k])[:, i])
        axs[i].set_xlabel('Time [s]')
        axs[i].set_ylabel(axisName)
        axs[i].legend(typeKeys)


def plotForces(outputs):
    """Plot all forces on a single plot"""
    plotXyzType(outputs, 'forces')


def plotMoments(outputs):
    """Plot all moments on a single plot"""
    plotXyzType(outputs, 'moments')


def plotOdeOutputs(outputs):
    """ Plots the outputs dictionary
    """
    # this is hacky, TODO fix
    keys = list(outputs.keys())
    keys.pop(keys.index('times'))
    times = outputs['times']
    for k in keys:
        print('plotting %s' % k)
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
    plt.tight_layout()


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
    plt.tight_layout()
# if not test:
# plt.show()
