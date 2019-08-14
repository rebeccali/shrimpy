import numpy as np
from matplotlib import pyplot as plt

from bladeDynamics import getPropForceMoment
from shrimpClasses import ShrimpState, defaultShaftPropParams, defaultShrimpParams


def playWDiffLift():
    shrimpParams = defaultShrimpParams()
    np.set_printoptions(suppress=True)
    r_w2b_w = np.array([0, 0, 0])
    vel_w2b_w = np.array([0, 0, 0])
    euler_w2b = np.array([0, 0, 0]) * np.pi / 180.
    angvel_w2b_b = np.array([0, 0, 0])
    inflowVel = 0
    yaw_b2p = 0
    yaw_f2b = 0
    yawDot_f2b = 0
    yawDot_b2p = 0. #rotPerSec2RadiansPerSecond(500)
    state = ShrimpState.fromVecState(r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                                     yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b)
    propParams = defaultShaftPropParams()
    propParams.numBlades = 1
    (f1, m1) = getPropForceMoment(propParams, shrimpParams, state)

    print('f1', f1)
    print('m1', m1)

    def varyVel(transVel):
        vel_w2b_w = np.array([0, transVel, 0])
        # We should get differential lift regardless
        state3 = ShrimpState.fromVecState(r_w2b_w, vel_w2b_w, euler_w2b, angvel_w2b_b, inflowVel,
                                          yaw_b2p, yaw_f2b, yawDot_b2p, yawDot_f2b)
        (f3, m3) = getPropForceMoment(propParams, shrimpParams, state3)
        return (f3, m3)

    transVels = np.arange(0.1, 3, 0.1)
    forceMoments = [varyVel(v) for v in transVels]
    forces = np.array([f for (f,m) in forceMoments])
    moments = np.array([m for (f, m) in forceMoments])
    momentsOvervel = np.divide(moments.T, transVels).T
    plt.figure(figsize=(12,8))
    plt.subplot(411)
    plt.plot(transVels, forces[:,2])
    plt.ylabel('Lift')
    plt.subplot(412)
    plt.plot(transVels, -moments[:, 1])
    plt.xlabel('Translational velocity')
    plt.ylabel('- Diff Lift moment')
    plt.subplot(413)
    plt.plot(transVels, -momentsOvervel[:, 1])
    plt.xlabel('Translational velocity')
    plt.ylabel('- Diff Lift moment/V')
    plt.subplot(414)
    plt.plot(transVels, forces[:,2]/(-moments[:, 1]))
    plt.xlabel('Translational velocity')
    plt.ylabel('Lift/(- Diff Lift moment)')
    plt.show()

playWDiffLift()