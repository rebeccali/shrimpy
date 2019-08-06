#!/usr/bin/env python
"""
Blade Element Theory for Shrimp Project
Rebecca Li 2019
"""
import argparse
# import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from shrimpClasses import PropellerParameters

def getParamsFromTaguchiArray(index, pitchRootDeg):
    aspectRatios = [4,4,4,4.5,4.5,4.5,5,5,5]
    cambers = [6,4,6,6,4,6,6,4,6]
    symmetric = [True, False, False, True, False, False,  True, False, False]
    twistDeg = [20,10,0,10,0,20,20,0,10]
    bumpType = [0,1,2,2,0,1,1,2,0]
    underbump = [True, False, False, False, False, False, False, False, False]
    radius = 16  # mm
    maxChord = 8.8
    minChord = 5
    numBlades = 2

    return PropellerParameters()


def main(args):
    """ Process Prop Data """
    filename = '../lib/propdata.csv'
    if args.f is not None:
        filename = args.f
    df = pd.read_csv(filename)

    columnNames = list(df.columns)
    assert columnNames[0] == 'Wing Version', 'Expected "Wing Version" as first column'
    assert columnNames[1] == 'Voltage - Power Supply (V)', 'Expected "Voltage - Power Supply (V)" as first column'




    n = 4  # TODO derive from actual length of each propellor data
    propellerDataframes = [df[i * n:(i + 1) * n] for i in range((len(df) + n - 1) // n)]

    for prop in propellerDataframes:
        propName = prop[columnNames[0]][0]  # First row name
        if '_' in propName:
            # Get angle of attack from propeller name
            i = propName.index('_')
            pitchDeg = float(propName[i+1:i+3])
            pitch = pitchDeg * np.pi/180.
        voltage = prop[columnNames[1]]
        thrust = prop[columnNames[2]]

        # Write an interpolation function that takes in voltage and gives thrust
        getThrustFromVoltage = interp1d(voltage, thrust, fill_value='extrapolate')

        # Plot the interpolation
        vs = np.arange(min(voltage)-0.3, max(voltage)+0.3, 0.1)
        ts = getThrustFromVoltage(vs)
        plt.figure()
        plt.scatter(voltage, thrust)
        plt.plot(vs, ts)
        plt.title(propName)
        break
    plt.show()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Prop Data CSV file to be processed",
                        type=str)
    args = parser.parse_args()
    main(args)
