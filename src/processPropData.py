#!/usr/bin/env python
"""
Process Propeller data for shrimp project
TO USE: run `python processPropData -f myprop.csv`. You can see examples of csv format
in the `lib/wingV6Data.csv` folder.
Also change the propeller parameters in the first few lines (chord at the tip, etc) to
whatever you want
Rebecca Li 2019
"""
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

from bladeDynamics import getPropForceMoment
from mathUtil import rotPerSec2RadiansPerSecond
from shrimpClasses import PropellerType, zeroShrimpState, PropellerParameters, defaultShrimpParams

# Edit these parameters for the particular wing you are processing!!
chordTip = 6e-3  # m
chordRoot = 9.7e-3  # m
radiusTip = 16e-3  # m
radiusRoot = 2e-3  # m
rho = 1.225  # kg/m^3
# chordTip = radiusTip * radiusTip / (chordRoot * aspectRatio)

# Pre-initializing for speeeed
_shrimpParams = defaultShrimpParams()  # Only rho matters
_shrimpParams.rho = rho
_shrimpState = zeroShrimpState()


def getThrustDragProp(propParams, rotorFreqHz):
    # Set up parameters
    _shrimpState.yawDot_b2p = rotPerSec2RadiansPerSecond(rotorFreqHz)

    (forces_b, moments_b) = getPropForceMoment(propParams, _shrimpParams, _shrimpState)
    thrust = forces_b[2]
    drag = moments_b[2]
    return thrust, drag


def getParamsFromTaguchiArray(taguchiIndex, pitchRootDeg):
    """ Returns PropellerParameters from Taguchi Array
        See Complete Guide to Shrimp Wings google doc for an explanation
    """
    # Unused Parameters
    # cambers = [6,4,6,6,4,6,6,4,6]
    # symmetric = [True, False, False, True, False, False,  True, False, False]
    # bumpType = [0,1,2,2,0,1,1,2,0]
    # underbump = [True, False, False, False, False, False, False, False, False]

    aspectRatios = [4, 4, 4, 4.5, 4.5, 4.5, 5, 5, 5]
    twistDegs = [20, 10, 0, 10, 0, 20, 20, 0, 10]

    # Particulars
    aspectRatio = aspectRatios[taguchiIndex - 1]
    twistDeg = twistDegs[taguchiIndex - 1]

    propType = PropellerType.SHAFT
    height_b2p = 0.00001
    numBlades = 2
    radiusRootTip = (radiusRoot, radiusTip)

    pitchRootTip = np.array([pitchRootDeg, twistDeg + pitchRootDeg]) * np.pi / 180.
    chordRootTip = np.array([chordRoot, chordTip])

    return PropellerParameters.fromRootTipParams(numBlades, pitchRootTip, radiusRootTip, chordRootTip,
                                                 height_b2p, propType)


def cleanUpPropDataframe(df):
    """ Remove unused columns and rows """
    # remove empty start column
    df.drop(columns='Unnamed: 0')
    # remove repetitive headers
    dfClean = df[df['Wing Version'] != 'Wing Version']
    return dfClean


def main(args):
    """ Process Prop Data """
    filename = args.f
    dfUnclean = pd.read_csv(filename)
    df = cleanUpPropDataframe(dfUnclean)

    wingColIndex = df.columns[['Wing Version' in c for c in df.columns]]  # This is a pandas.index
    thrustGramsColIndex = df.columns[['Thrust' in c for c in df.columns]]  # This is a pandas.index
    voltageColIndex = df.columns[['Voltage - Arduino' in c for c in df.columns]]  # This is a pandas.index
    doubleFreqColIndex = df.columns[['Frequency - Osc' in c for c in df.columns]]  # This is a pandas.index

    assert not wingColIndex.empty, 'Cannot find "Wing Version" Column'
    assert not thrustGramsColIndex.empty, 'Cannot find "Thrust" Column'
    assert not voltageColIndex.empty, 'Cannot find "Voltage - Arduino" Column'
    assert not doubleFreqColIndex.empty, 'Cannot find "Frequency - Osc" Column'

    n = 4
    propellerDataframes = [df[i * n:(i + 1) * n] for i in range((len(df) + n - 1) // n)]
    figs = []
    clFudgeFactors = []
    pitchRootDegs = []
    pitchTipDegs = []
    for prop in propellerDataframes:

        propName = list(prop[wingColIndex].values.flatten())[0]  # First row name
        print('Processing %s .' % propName)
        # Get angle of attack from propeller name
        assert '_' in propName, 'Expected underscore in propeller name, instead got %s' % propName
        i = propName.index('_')
        pitchRootDeg = float(propName[i + 1:i + 3])
        taguchiIndex = int(propName[i - 1:i])
        params = getParamsFromTaguchiArray(taguchiIndex, pitchRootDeg)
        pitchTipDeg = params.getPitchFromRadius(params.radiusRootTip[1]) * 180. / np.pi
        # Parse the table data
        voltage = [float(i) for i in prop[voltageColIndex].values]
        gramsToNewtons = 0.001 * 9.81
        thrust = [float(i) * gramsToNewtons for i in prop[thrustGramsColIndex].values]
        freqHz = [float(i) * 0.5 for i in prop[doubleFreqColIndex].values]
        pitchRootDegs.append(pitchRootDeg)
        pitchTipDegs.append(pitchTipDeg)
        # Search for the CL fudge factor
        clFudge = 0.
        def fitFreqVsThrust(clFudgeFactor):
            fudgedParams = params
            fudgedParams.clFudgeFactor = clFudgeFactor
            thrustDragFudgedTuple = [getThrustDragProp(fudgedParams, f) for f in freqHz]
            thrustFudged = np.array([t for (t, d) in thrustDragFudgedTuple])
            # return error compared to true thrust
            thrustError = np.linalg.norm(np.array(thrust) - thrustFudged)
            return thrustError

        # Bound the CL search to be reasonable
        result = minimize_scalar(fitFreqVsThrust, bounds=(-3., 3.), method='bounded')
        if result.success:
            clFudge = result.x
            print('Fudge factor: %f' % result.x)
        else:
            raise Exception("Failed to optimize CL: %s" % result.message)

        params.clFudgeFactor = clFudge
        clFudgeFactors.append(clFudge)

        # Compare to the analytical data including the fudge factor
        freqHzSpan = np.linspace(min(freqHz) - 5, max(freqHz) + 5, 10)
        thrustDragTuple = [getThrustDragProp(params, f) for f in freqHzSpan]
        thrustAnalytical = [t for (t, d) in thrustDragTuple]
        dragAnalytical = [d for (t, d) in thrustDragTuple]

        # Write an interpolation function that takes in voltage and gives thrust
        getThrustFromVoltage = interp1d(voltage, thrust, fill_value='extrapolate')

        # Plot the interpolation
        vs = np.arange(min(voltage) - 0.3, max(voltage) + 0.3, 0.1)
        ts = getThrustFromVoltage(vs)
        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(12, 8))
        figs.append(fig)
        ax1.scatter(voltage, thrust, label='Data')
        ax1.plot(vs, ts, label='Interpolation')
        ax1.set_xlabel('Voltage [V]')
        ax1.set_ylabel('Thrust [N]')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        titleText = '%s : (root, tip) Pitch %2.1f, %2.1f [deg]' % (propName, pitchRootDeg, pitchTipDeg)
        ax1.set_title(titleText)

        ax2.scatter(freqHz, thrust, label='Data')
        ax2.plot(freqHzSpan, thrustAnalytical, 'r-', label=('Flat Plate Model, Cl Fudge factor %f' % clFudge))
        ax2.set_xlabel('Frequency [Hz]')
        ax2.set_ylabel('Thrust [N]')
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        ax3.plot(freqHzSpan, dragAnalytical, 'r-', label='Flat Plate Model')
        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_ylabel('Drag Moment [Nm]')
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig = plt.figure()
    plt.scatter(pitchTipDegs, clFudgeFactors, label='Pitch Tip Deg', marker='o')
    plt.scatter(pitchRootDegs, clFudgeFactors, label='Pitch Root Deg', marker='*')
    plt.xlabel('Blade Pitch [Deg]')
    plt.ylabel('Cl Fudge Factor')
    plt.title('Cl fudge facter vs blade pitch')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figs.append(fig)

    # plt.show()
    # save file
    plotFilename = os.path.splitext(filename)[0] + '.pdf'
    pp = PdfPages(plotFilename)
    for f in figs:
        pp.savefig(f, bbox_inches='tight')
    pp.close()
    print('Saved plots to %s' % plotFilename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Prop Data CSV file to be processed", type=str, default='../lib/wingV3Data.csv')
    parser.add_argument("--header", help="Header for every line", action='store_true', default=False)
    _args = parser.parse_args()
    main(_args)
