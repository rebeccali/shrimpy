import numpy as np

from shrimpClasses import PropellerType, PropellerParameters, defaultShrimpParams, ShrimpState


def w3ShaftPropParams():
    """ Parameters for W3 Propeller
    """
    chordTip = 6e-3  # m
    chordRoot = 9.7e-3  # m
    radiusTip = 16e-3  # m
    radiusRoot = 2e-3  # m
    rho = 1.225  # kg/m^3
    numBlades = 2
    pitchRootTip = (20. * np.pi / 180., 20. * np.pi / 180.)
    height_b2p = 10.55 / 1000  # Piccolissimo_V11, positive is UP
    chordRootTip = (chordRoot, chordTip)
    radiusRootTip = (radiusRoot, radiusTip)
    propType = PropellerType.SHAFT
    clFudge = -0.34  # from data
    clFudge = 0.04
    p = PropellerParameters.fromRootTipParams(numBlades, pitchRootTip, radiusRootTip,
                                              chordRootTip, height_b2p, propType)
    p.clFudgeFactor = clFudge
    return p


def w3ShrimpParams():
    """ Shrimp Params with W3"""
    p = defaultShrimpParams()
    p.shaftPropParams = w3ShaftPropParams()
    return p


def zeroOdeState():
    return np.zeros(16)


def zeroShrimpState():
    """ Generates Shrimp state where everything is zero/nominal """
    odeState = np.zeros(16)
    state = ShrimpState.fromOdeState(odeState)
    return state
