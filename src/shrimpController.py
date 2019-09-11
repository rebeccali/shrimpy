#!/usr/bin/env python
"""
Control Scheme for Shrimp Project
Rebecca Li 2019
"""
import numpy as np
from shrimpClasses import defaultShrimpParams, dummyShrimpState


def shrimpController(p, s):
    """ A super simple PWM controller for the shrimp vehicle """
    vel_z = s.vel_w2b_w[2]
    command_vel_z = 0
    error_z = command_vel_z - vel_z
    pwmBias = 0.5
    Kp = 0.1
    pwmCommand = np.clip(pwmBias - error_z * Kp, 0, 1)
    pwmCommand = 1.
    return pwmCommand


def testShrimpController():
    p = defaultShrimpParams()
    s = dummyShrimpState()
    pwm = shrimpController(p, s)
    print('Pwm', pwm)
