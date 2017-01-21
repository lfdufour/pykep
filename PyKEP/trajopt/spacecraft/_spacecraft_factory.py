#from PyKEP.sims_flanagan._sims_flanagan import spacecraft

def getObjectSpacecraft(sc):
	from PyKEP.sims_flanagan import spacecraft
	minN = sc.AT*sc.minP + sc.BT;
	maxN = sc.AT*sc.maxP + sc.BT;
	return {'minP':sc.minP, 'maxP':sc.maxP, 'minN':minN, 'maxN':maxN, 'Isp':sc.isp, 'P1AU': sc.P1AU, 'Pmargin': sc.Pmargin}
