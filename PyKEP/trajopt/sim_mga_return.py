def test():
	from PyGMO import problem, algorithm, island, archipelago, topology, migration
	from _mga_return_lt_nep import mga_return_lt_nep
	import glob
	from PyKEP.sims_flanagan import spacecraft
	l = list()
	algo = algorithm.scipy_slsqp(max_iter=2000, acc=1e-4)
	algo2= algorithm.mbh(algo, 5, 0.07)
	algo2.screen_output = True
	for i in range(40):
		#define problem and optimization
		prob = problemCeSaM()

		"""algo_slsqp = algorithm.scipy_slsqp(max_iter=500, acc=1e-5)
		# algo = algorithm.snopt(major_iter=500, opt_tol=1e-3, feas_tol=1e-9)
		algo2_slsqp = algorithm.mbh(algo_slsqp, 5, 0.1)
		algo2_slsqp.screen_output = True
		archi = archipelago(algo2_slsqp,prob, 4, 1, topology = topology.fully_connected())
		print("Running Monotonic Basin Hopping ....")
		for j in range(1):
			archi.evolve(2)
			print min([isl.population.champion.f[0] for isl in archi])
		archi.join();
		tmp = [isl for isl in archi];
		tmp.sort(key = lambda x: x.population.champion.f[0]);
		x_so = tmp[0].population.champion.x;

		print("Is the solution found a feasible trajectory? " +
			str(prob.feasibility_x(x_so)))

		l.append(tmp[0].population.champion);

		if(prob.feasibility_x(x_so)):
			#find new file name
			filesjson = glob.glob("Web/*.json")
			ids = [-1]
			for filejson in filesjson:
				if filejson[4:-5].isdigit():
					ids = ids + [int(float(filejson[4:-5]))]
			newid = max(ids)+1
			prob.point(x_so,'Web/' + str(newid) + '.json');"""


		isl = island(algo2, prob, 1)
		print("Running Monotonic Basin Hopping ....")
		isl.evolve(1)
		isl.join()
		print("Is the solution found a feasible trajectory? " +
			str(prob.feasibility_x(isl.population.champion.x)))
		l.append(isl.population.champion)
		print (prob._compute_constraints_impl(isl.population.champion.x) )
		if(prob.feasibility_x(isl.population.champion.x)):
			#find new file name
			filesjson = glob.glob("Web/*.json")
			ids = []
			for filejson in filesjson:
				if filejson[4:-5].isdigit():
					ids = ids + [int(float(filejson[4:-5]))]
			newid = max(ids)+1
			prob.point(isl.population.champion.x,'Web/' + str(newid) + '.json');

	return l,prob

from _mga_return_lt_nep import mga_return_lt_nep
from PyKEP.sims_flanagan import spacecraft
from motor_factory import getMotor

class problemCeSaM(mga_return_lt_nep):
	#result solved in DB/1
	def __init__(self):
		AT, BT, maxP, minP, Isp = getMotor('custom 1')
		#AT, BT, maxP, minP, Isp = getMotor('T6 QinetiQ x2')
		sc = spacecraft(700,0.2,Isp*0+3000,AT,BT,maxP,minP,21000,300); #21000 = 60m^2
		super(problemCeSaM, self).__init__(spacecrafts = [sc,sc], vinf_dep1 = 4, solar_powered = True,
			tof=[[40,800],[700,2000],[700,2000],[40,800]])

