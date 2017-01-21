def point_planet(plnt, t0='PyKEP.epoch(0)', N=60, units=1.0):
    from PyKEP import MU_SUN, SEC2DAY, epoch, AU
    from math import pi, sqrt
    import numpy as np


    if t0 == 'PyKEP.epoch(0)':
        t0 = epoch(0)

    # orbit period at epoch
    T = plnt.compute_period(t0) * SEC2DAY

    # points where the orbit will be plotted
    when = np.linspace(0, T, N)

    # Ephemerides Calculation for the given planet
    x = np.array([0.0] * N)
    y = np.array([0.0] * N)
    z = np.array([0.0] * N)

    for i, day in enumerate(when):
        r, v = plnt.eph(epoch(t0.mjd2000 + day))
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units

    return x,y,z


def point_kepler(r, v, t, mu, N=60, units=1):
    from PyKEP import propagate_lagrangian

    # We define the integration time ...
    dt = t / (N - 1)

    # ... and calculate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        r, v = propagate_lagrangian(r, v, dt, mu)

    tt = range(N)
    tt = [e*dt/(60*60*24) for e in tt]

    return x,y,z,tt


def point_taylor(r, v, m, u, t, mu, veff, N=60, units=1):

    from PyKEP import propagate_taylor

    # We define the integration time ...    
    dt = t / (N - 1)
    mm = []

    # ... and calcuate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        mm = mm + [m]
        r, v, m = propagate_taylor(r, v, m, u, dt, mu, veff, -10, -10)
        
    tt = range(N)
    tt = [e*dt/(60*60*24) for e in tt]

    return x,y,z,tt,mm


def point_sf_leg(leg, N=5, units=1, color='b'):
    
    from PyKEP import propagate_lagrangian, AU, DAY2SEC, G0, propagate_taylor, SEC2DAY
    import numpy as np
    from scipy.linalg import norm
    from math import exp, sqrt


    # We compute the number of segments for forward and backward propagation
    n_seg = len(leg.get_throttles())
    fwd_seg = (n_seg + 1) // 2
    back_seg = n_seg // 2

    #we see if the propulsion is nuclear or solar powered
    sp = leg.solar_powered;

    # We extract information on the spacecraft
    sc = leg.get_spacecraft()
    isp = sc.isp
    max_thrust = sc.thrust

    # And on the leg
    throttles = leg.get_throttles()
    mu = leg.get_mu()

    # Forward propagation

    # x,y,z contain the cartesian components of all points (grid+midpints)
    x_bounds = []
    y_bounds = []
    z_bounds = []
    x = []
    y = []
    z = []
    Tx = []
    Ty = []
    Tz = []
    Tmax = []
    mm = []
    tt = [] 

    state = leg.get_xi()

    # Initial conditions
    r = state.r
    v = state.v
    m = state.m
    x_bounds = x_bounds + [r[0] / units]
    y_bounds = y_bounds + [r[1] / units]
    z_bounds = z_bounds + [r[2] / units]

    segmentID = 0
    temptime = 0
    # We compute all points by propagation
    for i, t in enumerate(throttles[:fwd_seg]):
        dt = (t.end.mjd - t.start.mjd) * DAY2SEC
        alpha = min(norm(t.value), 1.0)
        # Keplerian propagation and dV application
        if leg.high_fidelity is False:
            if(sp):
                dSun = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                max_thrust = sc.get_thrust_electricSolar(dSun)

            dV = [max_thrust / m * dt * dumb for dumb in t.value]
            xi,yi,zi,ti = point_kepler(r, v, dt / 2, mu, N=N, units=units)
            r, v = propagate_lagrangian(r, v, dt / 2, mu)
            #x[-2 * i - 2] = r[0] / units
            #y[-2 * i - 2] = r[1] / units
            #z[-2 * i - 2] = r[2] / units
            x = x + xi
            y = y + yi
            z = z + zi
            tt = tt + [e + temptime for e in ti]
            temptime = tt[-1]
            mm = mm + [m]*len(ti)
            x_bounds = x_bounds + [r[0] / units]
            y_bounds = y_bounds + [r[1] / units]
            z_bounds = z_bounds + [r[2] / units]
            # the thrust is constant on the segment (1st half-segment)
            Tx = Tx + [t.value[0]]*len(ti)
            Ty = Ty + [t.value[1]]*len(ti)
            Tz = Tz + [t.value[2]]*len(ti)
            Tmax = Tmax + [max_thrust]*len(ti)
            # v= v+dV
            v = [a + b for a, b in zip(v, dV)]
            xi,yi,zi,ti = point_kepler(r, v, dt / 2, mu, N=N, units=units)
            r, v = propagate_lagrangian(r, v, dt / 2, mu)
            #x[-2 * i - 3] = r[0] / units
            #y[-2 * i - 3] = r[1] / units
            #z[-2 * i - 3] = r[2] / units
            m *= exp(-norm(dV) / isp / G0)
            x = x + xi
            y = y + yi
            z = z + zi
            tt = tt + [e + temptime for e in ti]
            temptime = tt[-1]
            mm = mm+ [m]*len(ti)
            x_bounds = x_bounds + [r[0] / units]
            y_bounds = y_bounds + [r[1] / units]
            z_bounds = z_bounds + [r[2] / units]
            # the thrust is constant on the segment (2nd half-segment)
            Tx = Tx + [t.value[0]]*len(ti)
            Ty = Ty + [t.value[1]]*len(ti)
            Tz = Tz + [t.value[2]]*len(ti)
            Tmax = Tmax + [max_thrust]*len(ti)
        else:
            if(sp):
                dSun = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                max_thrust = sc.get_thrust_electricSolar(dSun)

            u = [max_thrust * dumb for dumb in t.value]
            xi,yi,zi,ti,mi = point_taylor(r, v, m, u, dt , mu, isp * G0, N=N, units=units)
            x = x + xi
            y = y + yi
            z = z + zi
            tt = tt + [e + temptime for e in ti]
            temptime = tt[-1]
            mm = mm + mi
            r, v, m = propagate_taylor(r, v, m, u, dt , mu, isp * G0, -12, -12)
            x_bounds = x_bounds + [r[0] / units]
            y_bounds = y_bounds + [r[1] / units]
            z_bounds = z_bounds + [r[2] / units]
            # the thrust is constant of the whole segment
            Tx = Tx + [t.value[0]]*len(ti)
            Ty = Ty + [t.value[1]]*len(ti)
            Tz = Tz + [t.value[2]]*len(ti)
            Tmax = Tmax + [max_thrust]*len(ti)


    state = leg.get_xf()

    # Final conditions
    r = state.r
    v = state.v
    m = state.m
    # used to reverse the order caused by the back-propagation
    xt = []
    yt = []
    zt = []
    mmt = []
    ttt = []
    Txt = []
    Tyt = []
    Tzt = []
    Tmaxt = []
    x_bounds = x_bounds + [r[0] / units] #gonna be the wrong side over...
    y_bounds = y_bounds + [r[1] / units]
    z_bounds = z_bounds + [r[2] / units]

    temptime = 0
    for i, t in enumerate(throttles[-1:-back_seg - 1:-1]):
        dt = (t.end.mjd - t.start.mjd) * DAY2SEC
        alpha = min(norm(t.value), 1.0)
        if leg.high_fidelity is False:
            if(sp):
                dSun = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                max_thrust = sc.get_thrust_electricSolar(dSun)

            dV = [max_thrust / m * dt * dumb for dumb in t.value]
            xi,yi,zi,ti = point_kepler(r, v, -dt / 2, mu, N=N, units=units)
            r, v = propagate_lagrangian(r, v, -dt / 2, mu)
            #x[-2 * i - 2] = r[0] / units
            #y[-2 * i - 2] = r[1] / units
            #z[-2 * i - 2] = r[2] / units
            xt = xt + xi
            yt = yt + yi
            zt = zt + zi
            ttt = ttt + [e + temptime for e in ti]
            temptime = ttt[-1]
            mmt = mmt + [m]*len(ti)
            x_bounds = x_bounds + [r[0] / units]
            y_bounds = y_bounds + [r[1] / units]
            z_bounds = z_bounds + [r[2] / units]
            # the thrust is constant on the segment (1st half-segment)
            Txt = Txt + [t.value[0]]*len(ti)
            Tyt = Tyt + [t.value[1]]*len(ti)
            Tzt = Tzt + [t.value[2]]*len(ti)
            Tmaxt = Tmaxt + [max_thrust]*len(ti)
            # v= v+dV
            v = [a - b for a, b in zip(v, dV)]
            xi,yi,zi,ti = point_kepler(r, v, -dt / 2, mu, N=N, units=units)
            r, v = propagate_lagrangian(r, v, -dt / 2, mu)
            #x[-2 * i - 3] = r[0] / units
            #y[-2 * i - 3] = r[1] / units
            #z[-2 * i - 3] = r[2] / units
            m *= exp(norm(dV) / isp / G0)
            xt = xt + xi
            yt = yt + yi
            zt = zt + zi
            ttt = ttt + [e + temptime for e in ti]
            temptime = ttt[-1]
            mmt = mmt + [m]*len(ti)
            x_bounds = x_bounds + [r[0] / units]
            y_bounds = y_bounds + [r[1] / units]
            z_bounds = z_bounds + [r[2] / units]
            # the thrust is constant on the segment (2nd half-segment)
            Txt = Txt + [t.value[0]]*len(ti)
            Tyt = Tyt + [t.value[1]]*len(ti)
            Tzt = Tzt + [t.value[2]]*len(ti)
            Tmaxt = Tmaxt + [max_thrust]*len(ti)
        else:
            if(sp):
                dSun = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                max_thrust = sc.get_thrust_electricSolar(dSun)

            u = [max_thrust * dumb for dumb in t.value]
            xi,yi,zi,ti,mi = point_taylor(r, v, m, u, -dt , mu, isp * G0, N=N, units=units)
            xt = xt + xi
            yt = yt + yi
            zt = zt + zi
            ttt = ttt + [e + temptime for e in ti]
            temptime = ttt[-1]
            mmt = mmt + mi
            r, v, m = propagate_taylor(r, v, m, u, -dt , mu, isp * G0, -12, -12)
            x_bounds = x_bounds + [r[0] / units]
            y_bounds = y_bounds + [r[1] / units]
            z_bounds = z_bounds + [r[2] / units]
            # the thrust is constant on the whole segment
            Txt = Txt + [t.value[0]]*len(ti)
            Tyt = Tyt + [t.value[1]]*len(ti)
            Tzt = Tzt + [t.value[2]]*len(ti)
            Tmaxt = Tmaxt + [max_thrust]*len(ti)

    Tpx = Tx
    Tpy = Ty
    Tpz = Tz
    x = x + xt[::-1]
    y = y + yt[::-1]
    z = z + zt[::-1]
    Tx = Tx + Txt[::-1]
    Ty = Ty + Tyt[::-1]
    Tz = Tz + Tzt[::-1]
    Tmax = Tmax + Tmaxt[::-1]
    mm = mm + mmt[::-1]
    temptime = tt[-1] - ttt[-1]
    tt = tt + [e + temptime for e in ttt[::-1]]
    return x,y,z,mm,tt,Tx,Ty,Tz,Tmax,x_bounds,y_bounds,z_bounds