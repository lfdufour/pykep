from PyGMO.problem import base as base_problem
from PyKEP.core import epoch, fb_con, EARTH_VELOCITY, AU, MU_SUN
#from PyKEP.planet import jpl_lp
from PyKEP import planet
from PyKEP.sims_flanagan._sims_flanagan import leg, spacecraft, sc_state
from PyKEP.trajopt.motor import getMotor
from scipy.linalg import norm

class mga_return_lt_nep(base_problem):

    """
    This class is a PyGMO (http://esa.github.io/pygmo/) problem representing a low-thrust
    interplanetary trajectory modelled as a Multiple Gravity Assist trajectory with sims_flanagan legs

    - Yam, C.H., di Lorenzo, D., and Izzo, D.,   Low-Thrust Trajectory Design as a Constrained Global Optimization Problem,  Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering, 225(11), pp.1243-1251, 2011.

    The decision vector (chromosome) is::

      [t0] + [T1, mf1, Vxi1, Vyi1, Vzi1, Vxf1, Vyf1, Vzf1] + [T2, mf2, Vxi2, Vyi2, Vzi2, Vxf2, Vyf2, Vzf2] + ... + [throttles1] + [throttles2] + ...

    .. note::

      The resulting problem is non linearly constrained. The resulting trajectory is not time-bounded.
    """

    def __init__(self,
                 seq1=['earth','mars','ceres'],
                 seq2=['ceres','mars','earth'],
                 n_seg=[15] * 4,
                 t0=[epoch(11000), epoch(13000)],
                 tos = [200,1000],# time of stay on planet
                 tof=[[40, 900]]*4,
                 vinf_dep1 = 0.001,
                 vinf_arr1 = 0.001,
                 vinf_dep2 = 0.001,
                 vinf_arr2 = 5,
                 dm = 400.0,
                 mass_end=700.0,
                 Tmax=0.4,
                 Isp=3800.0,
                 spacecrafts = None,
                 fb_rel_vel=6,
                 multi_objective=False,
                 high_fidelity=True,
		         solar_powered = True):
      

        """
        prob = mga_lt_nep(seq = [jpl_lp('earth'),jpl_lp('venus'),jpl_lp('earth')], n_seg = [10]*2,
        t0 = [epoch(0),epoch(1000)], tof = [[200,500],[200,500]], Vinf_dep=2.5, Vinf_arr=2.0, mass=4000.0, Tmax=1.0, Isp=2000.0,
        multi_objective = False, fb_rel_vel = 6, high_fidelity=False)

        - seq: list of PyKEP.planet defining the encounter sequence for the trajectoty (including the initial planet)
        - n_seg: list of integers containing the number of segments to be used for each leg (len(n_seg) = len(seq)-1)
        - t0: list of PyKEP epochs defining the launch window
        - tof: minimum and maximum time of each leg (days)
        - vinf_dep: maximum launch hyperbolic velocity allowed (in km/sec)
        - vinf_arr: maximum arrival hyperbolic velocity allowed (in km/sec)
        - mass: spacecraft starting mass
        - Tmax: maximum thrust
        - Isp: engine specific impulse
        - fb_rel_vel = determines the bounds on the maximum allowed relative velocity at all fly-bys (in km/sec)
        - multi-objective: when True defines the problem as a multi-objective problem, returning total DV and time of flight
        - high_fidelity = makes the trajectory computations slower, but actually dynamically feasible.
        """

        # We instanciate the sequences
        if seq1[-1] != seq2[0]:
            print 'Error : last planet of inbound and first planet of outbound must be the same'
        else:
            seq1 = self.make_sequence(seq1)
            seq2 = self.make_sequence(seq2)


        # 1) We compute the problem dimensions .... and call the base problem constructor
        self.__n_legs = len(seq1) + len(seq2) - 2
        n_fb = self.__n_legs - 2 # number of fylbys
        # 1a) The decision vector length
        dim = 1 + 1 + 1 + self.__n_legs * 8 + sum(n_seg) * 3
        # 1b) The total number of constraints (mass + mismatch + fly-by + boundary + throttles
        c_dim = 1 + self.__n_legs * 7 + n_fb * 2 + 4 + sum(n_seg)
        # 1c) The number of inequality constraints (boundary + fly-by angle + throttles)
        c_ineq_dim = 4 + n_fb + sum(n_seg)
        # 1d) the number of objectives
        f_dim = multi_objective + 1
        # First we call the constructor for the base PyGMO problem
        # As our problem is n dimensional, box-bounded (may be multi-objective), we write
        # (dim, integer dim, number of obj, number of con, number of inequality con, tolerance on con violation)
        super(_mga_return_lt_nep, self).__init__(dim, 0, f_dim, c_dim, c_ineq_dim, 1e-4)

        # 2) We then define some class data members
        # public:
        self.seq1 = seq1
        self.seq2 = seq2
        # private:
        self.__n_seg = n_seg
        self.__vinf_dep1 = vinf_dep1 * 1000
        self.__vinf_arr1 = vinf_arr1 * 1000
        self.__vinf_dep2 = vinf_dep2 * 1000
        self.__vinf_arr2 = vinf_arr2 * 1000
        self.__mf = mass_end;
        self.__dm = dm
        #Spacecaft(s)
        if spacecrafts == None:
            sc = spacecraft(mass_end,Tmax,Isp)
            print "!!Spacecraft not defined (motors)!!"
            #sc = spacecraft(mass_end, Tmax, Isp, 0.0000331,-0.0077, 9000, 2500, 29400, 300) #T6 x 2
            #sc = spacecraft(mass_end, Tmax, Isp, 0.0000331,-0.0077, 9000, 2500, 25000, 300) #T6 x 2
            self.__spacecrafts = [sc, sc];
        elif len(spacecrafts)==1:
            print "Copying spacecraft twice"
            self.__spacecrafts = [spacecrafts, spacecrafts];
        elif len(spacecrafts) == 2:
            self.__spacecrafts = spacecrafts;
        else:
            print "!!Spacecraft definition does ot have the right size!!"
            #sc = spacecraft(mass_end, Tmax, Isp, 0.0000331,-0.0077, 9000, 2500, 29400, 300) #T6 x 2
            #sc = spacecraft(mass_end, Tmax, Isp, 0.0000331,-0.0077, 9000, 2500, 25000, 300) #T6 x 2
        # T6 Qinetic x 2
        #self.__sc = spacecraft(mass_end, Tmax, Isp, 0.0000331,-0.0077, 9000, 2500, 25000, 300)
        #self.__sc = spacecraft(mass_end, Tmax, Isp, 0.00003333,-0.0072, 10000, 2000, 29400, 300)
        #self.__sc = spacecraft(mass_end, Tmax, Isp,0.0000,0.2, 9000, 2500, 10000, 300)
        self.__leg = leg()
        self.__leg.set_mu(MU_SUN)
        self.__leg.set_spacecraft(self.__spacecrafts[0])
        self.__leg.high_fidelity = high_fidelity
        self.__leg.solar_powered = solar_powered
        fb_rel_vel *= 1000
        # 3) We compute the bounds
        lb = [t0[0].mjd2000] + [tos[0]] + [1000] + [0, mass_end, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel, -fb_rel_vel] * self.__n_legs + [-1, -1, -1] * sum(self.__n_seg)
        ub = [t0[1].mjd2000] + [tos[1]] + [5000] + [1, 5000, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel, fb_rel_vel] * self.__n_legs + [1, 1, 1] * sum(self.__n_seg)
        # 3a ... and account for the bounds on the vinfs......
        # 3a1 for departure 1
        lb[5:8] = [-self.__vinf_dep1] * 3
        ub[5:8] = [self.__vinf_dep1] * 3
        # 3a2 for arrival 1
        lb[5 + (len(seq1) - 2)*8 :8 + (len(seq1) -2)*8] = [-self.__vinf_arr1] * 3
        ub[5 + (len(seq1) - 2)*8 :8 + (len(seq1) -2)*8] = [self.__vinf_arr1] * 3
        # 3a3 for departure 2
        lb[5 + (len(seq1) -2 + 1)*8 :8 + (len(seq1) -2 + 1)*8] = [-self.__vinf_dep2] * 3
        ub[5 + (len(seq1) -2 + 1)*8 :8 + (len(seq1) -2 + 1)*8] = [self.__vinf_dep2] * 3
        # 3a4 for arrival 2
        lb[-sum(self.__n_seg) * 3 - 3:-sum(self.__n_seg) * 3] = [-self.__vinf_arr2] * 3
        ub[-sum(self.__n_seg) * 3 - 3:-sum(self.__n_seg) * 3] = [self.__vinf_arr2] * 3
        # 3b... and for the time of flight
        lb[3:1 + 8 * self.__n_legs:8] = [el[0] for el in tof]
        ub[3:1 + 8 * self.__n_legs:8] = [el[1] for el in tof]
        # 4) And we set the bounds
        self.set_bounds(lb, ub)

    # used during init to create sequences
    def make_sequence(self, planets):
        seq = []
        for e in planets:
            if e == 'earth':
                seq.append(planet.jpl_lp('earth'))
            elif e == 'mars':
                seq.append(planet.jpl_lp('mars'))
            elif e == 'venus':
                seq.append(planet.jpl_lp('venus'))
            elif e == 'ceres':
                seq.append(planet.mpcorb('00001    3.34  0.12 K167V 224.09531   72.81471   80.31427 '\
                '  10.59170  0.0757051  0.21400472   2.7681342  0 MPO384741  6634 113 1801-2016'\
                ' 0.60 M-v 30h MPCLINUX   0000          Ceres              20160723'))
            else:
                print 'Planet not known'
        return seq

    # Objective function
    def _objfun_impl(self, x):
        if self.f_dimension == 1:
            #return (-x[4 + (self.__n_legs - 1) * 8],)
            return (x[2],)
        else:
            return (x[2], sum(x[3:1 + 8 * self.__n_legs:8]))

    # Constraints function
    def _compute_constraints_impl(self, x):

        seq1 = self.seq1
        seq2 = self.seq2
        # 1 - We decode the chromosome extracting the time of flights
        T = list([0] * (self.__n_legs))
        for i in range(self.__n_legs):
            T[i] = x[3 + i * 8]

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 2))
        r_P = list([None] * (self.__n_legs + 2))
        v_P = list([None] * (self.__n_legs + 2))

        for i, planet in enumerate(self.seq1):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq1[i].eph(t_P[i])
        for i, planet in enumerate(self.seq2):
            t_P[i+len(seq1)] = epoch(x[0] + sum(T[0:i+len(seq1)-1]) + x[1])
            r_P[i+len(seq1)], v_P[i+len(seq1)] = self.seq2[i].eph(t_P[i+len(seq1)])

        # 3 - We iterate through legs to compute mismatches and throttles constraints
        ceq = list()
        cineq = list()
        m0 = x[2]
        #inbound
        self.__leg.set_spacecraft(self.__spacecrafts[0])
        for i in range(len(self.seq1)-1):
            # First Leg
            v = [a + b for a, b in zip(v_P[i], x[(5 + i * 8):(8 + i * 8)])]
            x0 = sc_state(r_P[i], v, m0)
            v = [a + b for a, b in zip(v_P[i + 1], x[(8 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 1], v, x[4 + i * 8])
            throttles = x[(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i], x0, throttles, t_P[i + 1], xe)
            # update mass!
            m0 = x[4 + 8 * i]
            ceq.extend(self.__leg.mismatch_constraints())
            cineq.extend(self.__leg.throttles_constraints())
        #outbound
        self.__leg.set_spacecraft(self.__spacecrafts[1])
        m0 = m0 - self.__dm # mass lost around ceres
        for i in range(len(self.seq1)-1,len(self.seq1) + len(self.seq2)-2):
            # First Leg
            v = [a + b for a, b in zip(v_P[i+1], x[(5 + i * 8):(8 + i * 8)])]
            x0 = sc_state(r_P[i+1], v, m0)
            v = [a + b for a, b in zip(v_P[i + 2], x[(8 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 2], v, x[4 + i * 8])
            throttles = x[(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i+1], x0, throttles, t_P[i + 2], xe)
            # update mass!
            m0 = x[4 + 8 * i]
            ceq.extend(self.__leg.mismatch_constraints())
            cineq.extend(self.__leg.throttles_constraints())

        # Adding the boundary constraints
        # departure1
        v_dep_con1 = (x[5] ** 2 + x[6] ** 2 + x[7] ** 2 - self.__vinf_dep1 ** 2) / (EARTH_VELOCITY ** 2)
        # arrival1
        v_arr_con1 = (x[8 + (len(seq1) -2) * 8] ** 2 + x[9 + (len(seq1) -2) * 8] ** 2 + x[10 + (len(seq1) -2) * 8] ** 2 - self.__vinf_arr1 ** 2) / (EARTH_VELOCITY ** 2)
         # departure2
        v_dep_con2 = (x[5 + (len(seq1) -1) * 8] ** 2 + x[6 + (len(seq1) -1) * 8] ** 2 + x[7 + (len(seq1) -1) * 8] ** 2 - self.__vinf_dep2 ** 2) / (EARTH_VELOCITY ** 2)
        # arrival2
        v_arr_con2 = (x[8 + (len(seq1) -1 + len(seq2)-2) * 8] ** 2 + x[9 + (len(seq1) -1 + len(seq2)-2)] ** 2 + x[10 * (len(seq1) -1 + len(seq2)-2)] ** 2 - self.__vinf_arr2 ** 2) / (EARTH_VELOCITY ** 2)
        cineq.append(v_dep_con1 * 100)
        cineq.append(v_arr_con1 * 100)
        cineq.append(v_dep_con2 * 100)
        cineq.append(v_arr_con2 * 100)


        # We add the fly-by constraints
        for i in range(len(seq1)-2):
            DV_eq, alpha_ineq = fb_con(x[8 + i * 8:11 + i * 8], x[13 + i * 8:16 + i * 8], self.seq1[i + 1])
            ceq.append(DV_eq / (EARTH_VELOCITY ** 2))
            cineq.append(alpha_ineq)
        for i in range(len(seq2)-2):
            DV_eq, alpha_ineq = fb_con(x[8 + (i+len(seq1)-1) * 8:11 + (i+len(seq1)-1) * 8], x[13 + (i+len(seq1)-1) * 8:16 + (i+len(seq1)-1) * 8], self.seq2[i + 1])
            ceq.append(DV_eq / (EARTH_VELOCITY ** 2))
            cineq.append(alpha_ineq)


        #Adding the mass constraint
        mass_con = self.__mf - m0
        ceq.append(mass_con/self.__mf)

        # Making the mismatches non dimensional
        for i in range(self.__n_legs):
            ceq[0 + i * 7] /= AU
            ceq[1 + i * 7] /= AU
            ceq[2 + i * 7] /= AU
            ceq[3 + i * 7] /= EARTH_VELOCITY
            ceq[4 + i * 7] /= EARTH_VELOCITY
            ceq[5 + i * 7] /= EARTH_VELOCITY
            ceq[6 + i * 7] /= self.__mf


        # We assemble the constraint vector
        retval = list()
        retval.extend(ceq)
        retval.extend(cineq)

        return retval

    # And this helps visualizing the trajectory
    def plot(self, x, ax=None):
        """
        ax = prob.plot(x, ax=None)

        - x: encoded trajectory
        - ax: matplotlib axis where to plot. If None figure and axis will be created
        - [out] ax: matplotlib axis where to plot

        Plots the trajectory represented by a decision vector x on the 3d axis ax

        Example::

          ax = prob.plot(x)
        """

        print("CAUTION : does not work with solar powered!!")

        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        from PyKEP import epoch, AU
        from PyKEP.sims_flanagan import sc_state
        from PyKEP.orbit_plots import plot_planet, plot_sf_leg

        seq1 = self.seq1
        seq2 = self.seq2

        # Creating the axis if necessary
        if ax is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axis = fig.gca(projection='3d')
        else:
            axis = ax

        # Plotting the Sun ........
        axis.scatter([0], [0], [0], color='y')

        # Plotting the legs .......
        # 1 - We decode the chromosome extracting the time of flights
        T = list([0] * (self.__n_legs))
        for i in range(self.__n_legs):
            T[i] = x[3 + i * 8]

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 2))
        r_P = list([None] * (self.__n_legs + 2))
        v_P = list([None] * (self.__n_legs + 2))

        for i, planet in enumerate(self.seq1):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq1[i].eph(t_P[i])
        for i, planet in enumerate(self.seq2):
            t_P[i+len(seq1)] = epoch(x[0] + sum(T[0:i+len(seq1)-1]) + x[1])
            r_P[i+len(seq1)], v_P[i+len(seq1)] = self.seq2[i].eph(t_P[i+len(seq1)])

        # 3 - We iterate through legs to compute mismatches and throttles constraints
        ceq = list()
        cineq = list()
        m0 = x[2]
        #inbound
        self.__leg.set_spacecraft(self.__spacecrafts[0])
        for i in range(len(self.seq1)-1):
            # First Leg
            v = [a + b for a, b in zip(v_P[i], x[(5 + i * 8):(8 + i * 8)])]
            x0 = sc_state(r_P[i], v, m0)
            v = [a + b for a, b in zip(v_P[i + 1], x[(8 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 1], v, x[4 + i * 8])
            throttles = x[(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i], x0, throttles, t_P[i + 1], xe)
            # update mass!
            m0 = x[4 + 8 * i]
            plot_sf_leg(self.__leg, units=AU, N=10, ax=axis)
        #outbound
        self.__leg.set_spacecraft(self.__spacecrafts[1])
        m0 = m0-self.__dm # mass lost around ceres
        for i in range(len(self.seq1)-1,len(self.seq1) + len(self.seq2)-2):
            # First Leg
            v = [a + b for a, b in zip(v_P[i+1], x[(5 + i * 8):(8 + i * 8)])]
            x0 = sc_state(r_P[i+1], v, m0)
            v = [a + b for a, b in zip(v_P[i + 2], x[(8 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 2], v, x[4 + i * 8])
            throttles = x[(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i+1], x0, throttles, t_P[i + 2], xe)
            # update mass!
            m0 = x[4 + 8 * i]
            plot_sf_leg(self.__leg, units=AU, N=10, ax=axis)
            
        # Plotting planets
        for i, planet in enumerate(self.seq1):
            plot_planet(planet, t_P[i], units=AU, legend=True, color=(0.7, 0.7, 1), ax = axis)
        for i, planet in enumerate(self.seq2):
            plot_planet(planet, t_P[i+len(seq1)], units=AU, legend=True, color=(0.7, 0.7, 1), ax = axis)

        plt.show()
        return axis

    def high_fidelity(self, boolean):
        """
        prob.high_fidelity(status)

        - status: either True or False (True sets high fidelity on)

        Sets the trajectory high fidelity mode

        Example::

          prob.high_fidelity(True)
        """
        # We avoid here that objfun and constraint are kept that have been evaluated wrt a different fidelity
        self.reset_caches()
        # We set the propagation fidelity
        self.__leg.high_fidelity = boolean
       
    # def double_segments(self, x):
        """
        x_doubled = prob.double_segments(x)

        - x: compatible trajectory as encoded by an mga_1dsm mga_lt_nep

        Returns the decision vector encoding a low trust trajectory having double the number of segments with respect to x
        and a 'similar' throttle history. In case high fidelity is True, and x is a feasible trajectory, the returned decision vector
        also encodes a feasible trajectory that can be further optimized

        Example::

          prob = traj.mga_lt_nep(nseg=[[10],[20]])
          pop = population(prob,1)
          .......OPTIMIZE.......
          x = prob.double_segments(pop.champion.x)
          prob = traj.mga_lt_nep(nseg=[[20],[40]])
          pop = population(prob)
          pop.push_back(x)
          .......OPTIMIZE AGAIN......
        """
    #    y = list()
    #    y.extend(x[:-sum(self.__n_seg) * 3])
    #    for i in range(sum(self.__n_seg)):
    #        y.extend(x[-(sum(self.__n_seg) - i) * 3:-(sum(self.__n_seg) - 1 - i) * 3] * 2)
    #    y.extend(x[-3:] * 2)
    #    return y

    def point(self, x, filename):

        from PyKEP import epoch, AU
        from PyKEP.sims_flanagan import sc_state
        from PyKEP.orbit_plots import plot_planet, plot_sf_leg
        from PyKEP.orbit_plots import point_sf_leg, point_kepler, point_taylor, point_planet
        from scipy.linalg import norm
        from math import sqrt
        import csv
        import pprint
        from PyKEP.trajopt.motor import getObjectMotor
        from PyKEP.trajopt.spacecraft import getObjectSpacecraft

        XYZ = [];
        xx = []
        yy = []
        zz = []
        tt = []
        mm = []
        x_bounds = []
        y_bounds = []
        z_bounds = []

        seq1 = self.seq1
        seq2 = self.seq2

        data = {}

        # Plotting the legs .......
        # 1 - We decode the chromosome extracting the time of flights
        T = list([0] * (self.__n_legs))
        for i in range(self.__n_legs):
            T[i] = x[3 + i * 8]

        # 2 - We compute the epochs and ephemerides of the planetary encounters
        t_P = list([None] * (self.__n_legs + 2))
        r_P = list([None] * (self.__n_legs + 2))
        v_P = list([None] * (self.__n_legs + 2))

        for i, planet in enumerate(self.seq1):
            t_P[i] = epoch(x[0] + sum(T[0:i]))
            r_P[i], v_P[i] = self.seq1[i].eph(t_P[i])
        for i, planet in enumerate(self.seq2):
            t_P[i+len(seq1)] = epoch(x[0] + sum(T[0:i+len(seq1)-1]) + x[1])
            r_P[i+len(seq1)], v_P[i+len(seq1)] = self.seq2[i].eph(t_P[i+len(seq1)])
            
        # 3 - We iterate through legs to compute mismatches and throttles constraints
        ceq = list()
        cineq = list()
        m0 = x[2]
        data['legs'] = []
        massesAlongTheWay = [m0];
        self.__leg.set_spacecraft(self.__spacecrafts[0]);
        sc = self.__spacecrafts[0];
        for i in range(len(self.seq1)-1):
            # First Leg
            v = [a + b for a, b in zip(v_P[i], x[(5 + i * 8):(8 + i * 8)])]
            x0 = sc_state(r_P[i], v, m0)
            v = [a + b for a, b in zip(v_P[i + 1], x[(8 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 1], v, x[4 + i * 8])
            throttles = x[(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i], x0, throttles, t_P[i + 1], xe)
            # update mass!
            minit = m0
            m0 = x[4 + 8 * i]
            massesAlongTheWay = massesAlongTheWay + [m0]
            xi,yi,zi,mi,ti,Txi,Tyi,Tzi,Tmaxi,x_boundsi,y_boundsi,z_boundsi = point_sf_leg(self.__leg, units=AU, N=10)
            xx = xx + xi
            yy = yy + yi
            zz = zz + zi
            ti = [e+t_P[i].jd for e  in ti]
            #for k in range(len(Txi)):
            #    print norm([Tyi[k],Txi[k],Tzi[k]])
            data['legs']=data['legs']+[{'x':[e for e in xi],'y':[e for e in yi],'z':[e for e in zi],
            'm':[e for e in mi], 't':[e for e in ti],
            'P_tot':[sc.get_totalPower(norm(e)) for e in zip(xi,yi,zi)],
            'P_engine':[sc.get_powerThrust(norm(e[0:3])*e[3]) for e in zip(Txi,Tyi,Tzi,Tmaxi)],
            'Tx':[e for e in Txi],'Ty':[e for e in Tyi],'Tz':[e for e in Tzi],'Tmax':[e for e in Tmaxi],
            'mi':minit,'mf':m0,'planet1':{'name':seq1[i].name,'date':t_P[i].jd},
            'planet2':{'name':seq1[i+1].name,'date':t_P[i+1].jd}}]   
            XYZ = XYZ + [xi]+[yi]+[zi]
            x_bounds = x_bounds + x_boundsi
            y_bounds = y_bounds + y_boundsi
            z_bounds = z_bounds + z_boundsi
        #outbound
        m0 = m0-self.__dm # mass lost around ceres
        massesAlongTheWay = massesAlongTheWay + [m0];
        self.__leg.set_spacecraft(self.__spacecrafts[1]);
        sc = self.__spacecrafts[1];
        for i in range(len(self.seq1)-1,len(self.seq1) + len(self.seq2)-2):
            # First Leg
            v = [a + b for a, b in zip(v_P[i+1], x[(5 + i * 8):(8 + i * 8)])]
            x0 = sc_state(r_P[i+1], v, m0)
            v = [a + b for a, b in zip(v_P[i + 2], x[(8 + i * 8):(11 + i * 8)])]
            xe = sc_state(r_P[i + 2], v, x[4 + i * 8])
            throttles = x[(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i])):(3 + 8 * self.__n_legs + 3 * sum(self.__n_seg[:i]) + 3 * self.__n_seg[i])]
            self.__leg.set(t_P[i+1], x0, throttles, t_P[i + 2], xe)
            # update mass!
            minit = m0
            m0 = x[4 + 8 * i]
            massesAlongTheWay = massesAlongTheWay + [m0]
            xi,yi,zi,mi,ti,Txi,Tyi,Tzi,Tmaxi,x_boundsi,y_boundsi,z_boundsi = point_sf_leg(self.__leg, units=AU, N=10)
            xx = xx + xi
            yy = yy + yi
            zz = zz + zi
            ti = [e+t_P[i+1].jd for e  in ti]
            #for k in range(len(Txi)):
            #    print norm([Tyi[k],Txi[k],Tzi[k]])
            data['legs']=data['legs']+[{'x':[e for e in xi],'y':[e for e in yi],'z':[e for e in zi],
            'm':[e for e in mi], 't':[e for e in ti],
            'P_tot':[sc.get_totalPower(norm(e)) for e in zip(xi,yi,zi)],
            'P_engine':[sc.get_powerThrust(norm(e[0:3])*e[3]) for e in zip(Txi,Tyi,Tzi,Tmaxi)],
            'Tx':[e for e in Txi],'Ty':[e for e in Tyi],'Tz':[e for e in Tzi],'Tmax':[e for e in Tmaxi],
            'mi':minit,'mf':m0,'planet1':{'name':seq2[i-len(seq1)-2].name,'date':t_P[i+1].jd},
            'planet2':{'name':seq2[i+1-len(seq1)-2].name,'date':t_P[i+2].jd}}]
            XYZ = XYZ + [xi]+[yi]+[zi]
            x_bounds = x_bounds + x_boundsi
            y_bounds = y_bounds + y_boundsi
            z_bounds = z_bounds + z_boundsi
            
        XYZplanet = []
        data['planetOrbits'] = [];
        data['planetPosition'] = []
        listAlreadySaved = [];
        for i, planet in enumerate(self.seq1):
            xp,yp,zp = point_planet(planet, t_P[i], units=AU)
            #Save orbit if it has not been done already
            if (seq1[i].name not in listAlreadySaved):
                listAlreadySaved.append(seq1[i].name)
                data['planetOrbits'].append({'name':seq1[i].name,'x':[e for e in xp],'y':[e for e in yp],'z':[e for e in zp]})
            data['planetPosition'].append({'name':seq1[i].name, 'x':xp[0], 'y':yp[0],'z': zp[0],
                'dateS':pprint.pformat(t_P[i])[0:20], 'date':t_P[i].jd,'scmass': massesAlongTheWay[i]})
        for i, planet in enumerate(self.seq2):
            xp,yp, zp = point_planet(planet, t_P[i+len(seq1)], units=AU)
            #Save orbit if it has not been done already
            if (seq2[i].name not in listAlreadySaved):
                listAlreadySaved.append(seq2[i].name)
                data['planetOrbits'].append({'name':seq2[i].name,'x':[e for e in xp],'y':[e for e in yp],'z':[e for e in zp]})
            data['planetPosition'].append({'name':seq2[i].name, 'x':xp[0], 'y':yp[0], 'z': zp[0],
                'dateS':pprint.pformat(t_P[i+len(seq1)])[0:20],'date':t_P[i+len(seq1)].jd,'scmass': massesAlongTheWay[i+3]})


        #Calculate vinf velocities
        v_dep_1 = sqrt(x[5] ** 2 + x[6] ** 2 + x[7] ** 2);
        v_arr_1 = sqrt(x[8 + (len(seq1) -2) * 8] ** 2 + x[9 + (len(seq1) -2) * 8] ** 2 + x[10 + (len(seq1) -2) * 8] ** 2);
        v_dep_2 = sqrt(x[5 + (len(seq1) -1) * 8] ** 2 + x[6 + (len(seq1) -1) * 8] ** 2 + x[7 + (len(seq1) -1) * 8] ** 2);
        v_arr_2 = sqrt(x[8 + (len(seq1) -1 + len(seq2)-2) * 8] ** 2 + x[9 + (len(seq1) -1 + len(seq2)-2)] ** 2 + x[10 * (len(seq1) -1 + len(seq2)-2)] ** 2);

        #Logging all trajectory parameters
        data['traj'] = {}
        data['traj']['high_fidelity'] = self.__leg.high_fidelity
        data['traj']['solar_powered'] = self.__leg.solar_powered
        data['traj']['outbound'] = [[seq1[i].name, t_P[i].jd] for i in range(len(seq1))] 
        data['traj']['inbound'] = [[seq2[i].name, t_P[i+len(seq1)].jd] for i in range(len(seq2))] 
        data['traj']['vinf'] = {'v_dep1':v_dep_1,'v_dep2':v_dep_2,'v_arr1':v_arr_1,'v_arr2':v_arr_2};
        data['traj']['spacecrafts'] = [getObjectSpacecraft(self.__spacecrafts[0]),
                    getObjectSpacecraft(self.__spacecrafts[1])]
        data['traj']['dm'] = self.__dm;
        data['traj']['mi'] = x[2];
        data['traj']['mf'] = self.__mf;
        data['traj']['decision_vector']= x
        data['traj']['dt'] = x[1]
        data['traj']['t_launch'] = t_P[0].jd
        data['traj']['t_launchS'] = pprint.pformat(t_P[0])[0:20]
        data['traj']['t_return'] = t_P[-1].jd
        data['traj']['t_returnS'] = pprint.pformat(t_P[-1])[0:20]
        data['traj']['duration'] = t_P[-1].jd-t_P[0].jd

        import json
        with open(filename, 'w') as outfile:
            json.dump(data, outfile)
        # File structure:
        # Number of legs, length points, number of planets to plot, length points
        # Legs trajectories
        # Planets trajectories
        # Planets names
        #with open('test1.csv', 'w') as csvfile:
        #    writer = csv.writer(csvfile)
        #    #writer.writerow([self.__n_legs, len(XYZ[0]), len(list_planet),len(XYZplanet[0])])
        #    for i in range(len(XYZ)):
        #        writer.writerow(XYZ[i])
        #    for i in range(len(XYZplanet)):
        #        writer.writerow(XYZplanet[i])
        #    #writer.writerow(list_planet)
         
        return xx,yy,zz,x_bounds,y_bounds,z_bounds, data
