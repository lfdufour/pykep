# motor_factory:

def getMotor(name):
	# returns (AT, BT, maxP, minP, Isp)
	import json
	
	with open('motorDB.json') as data_file:
		data = json.load(data_file)
	if name in data.keys():
		AT = data[name]['AT']
		BT = data[name]['BT']
		maxP = data[name]['maxP']
		minP = data[name]['minP']
		Isp = data[name]['Isp']
		return AT, BT, maxP, minP, Isp
	else:
		print "This motor is not defined in motorDB.json"

def findName(maxN, minN, maxP, minP):
	from PyKEP.sims_flanagan import spacecraft
	import json
	with open('motorDB.json') as data_file:
		data = json.load(data_file)
	for key in data.keys():
		ATt, BTt, maxPt, minPt, Ispt = getMotor(key)
		minNt = ATt*minPt + BTt;
		maxNt = ATt*maxPt + BTt;
		if (minNt == minN) & (maxNt == maxN) & (maxPt == maxP) & (minPt == minP):
			return key
	return 'NaM'


def getObjectMotor(sc):
	from PyKEP.sims_flanagan import spacecraft
	minN = sc.AT*sc.minP + sc.BT;
	maxN = sc.AT*sc.maxP + sc.BT;
	return {'minP':sc.minP, 'maxP':sc.maxP, 'minN':minN, 'maxN':maxN, 'Isp':sc.isp}

def addMotor(name, minP, minN, maxP, maxN, Isp, IspMin, IspMax):
	# add motor to DB
	import json
	with open('motorDB.json') as data_file:
		data = json.load(data_file)
	if name in data.keys():
		print("The motor already exist in the database:")
		print("Old = " + name + ": minP = " + str(data[name]['minP']) +
		 	" - minN = " + str(data[name]['AT']*data[name]['minP']+data[name]['BT']) +
		 	" - maxP = " + str(data[name]['maxP']) + 
		 	" - maxN = " + str(data[name]['AT']*data[name]['maxP']+data[name]['BT']) + 
		 	" - Isp = " + str(data[name]['Isp']))
		print("New = " + name + ": minP = " + str(minP) +
		 	" - minN = " + str(minN) +
		 	" - maxP = " + str(maxP) + 
		 	" - maxN = " + str(maxN) + 
		 	" - Isp = " + str(Isp))
		while True:
			response = raw_input("Do you want to replace the existing motor? (y/n) ")
			if(response == 'n'):	
				return	
			elif(response == 'y'):
				break
	AT = (maxN - minN)/(maxP-minP);
	BT = minN- AT*minP;
	data[name] = {'AT':AT, 'BT':BT, 'minP':minP, 'maxP':maxP, 'Isp':Isp, 'IspMin':IspMin, 'IspMax': IspMax}
	with open('motorDB.json', 'w') as outfile:
		json.dump(data, outfile)


def createDB():

	import json
	with open('motorDB.json','w+') as data_file:
		dataInit = {'filename':'motorDB.json'};
		json.dump(dataInit, data_file);
	data_file.close();
	#NSTAR
	addMotor(name='NSTAR', minP=500.0, minN=19e-3, maxP=2300.0, maxN=92e-3, Isp=2400.0, IspMin=1900.0, IspMax =3100.0)
	#NEXT
	addMotor(name='NEXT', minP=1100.0, minN=50e-3 ,maxP=6100.0, maxN=237e-3, Isp= 3000.0, IspMin = 2210.0, IspMax = 4100.0)
	#NEXIS
	addMotor(name='NEXIS' , minP=15000.0 , minN=0.42 , maxP=25000.0 , maxN=0.475 , Isp=6800.0 ,IspMin=6000.0, IspMax=7500.0)
	#XIPS
	addMotor(name='XIPS' , minP=2200.0 , minN=79e-3 , maxP=4500.0 , maxN=168e-3 , Isp=3450.0 ,IspMin=3400.0, IspMax=3500.0)
 	#RIT 10 EVO
 	addMotor(name='RIT 10EVO' , minP=145.0 , minN=5e-3 , maxP=760.0 , maxN=25e-3 , Isp=2650.0 ,IspMin=1900.0, IspMax=3200.0)
	#RIT 2X
	addMotor(name='RIT 2X' , minP=2185.0 , minN=80e-3 , maxP=5785.0 , maxN=200e-3 , Isp=3850.0 ,IspMin=3400.0, IspMax=4300.0)
	#T6 QinetiQ
	addMotor(name='T6 QinetiQ' , minP= 2500.0, minN=75e-3 , maxP= 4500.0, maxN=145e-3 , Isp=3905.0 ,IspMin=3710.0, IspMax=4120.0)
	#T6 QinetiQ x2
	addMotor(name='T6 QinetiQ x2' , minP= 2500.0, minN=75e-3 , maxP= 4500.0*2, maxN=145e-3*2 , Isp=3905.0 ,IspMin=3710.0, IspMax=4120.0)
	#fakel SPT-140
	addMotor(name='SPT-140' , minP=1500.0 , minN=85e-3 , maxP= 5000.0, maxN=300e-3 , Isp=1850.0 ,IspMin= 1750.0, IspMax=1950.0)	
	#custom 1 - based on 2x QinetiQ T6
	addMotor(name='custom 1' , minP=1500 , minN= 40e-3, maxP= 4500*2, maxN= 145e-3*2, Isp= 3000,IspMin=2500, IspMax=3500)	
	# template addMotor(name= , minP= , minN= , maxP= , maxN= , Isp= ,IspMin=, IspMax=)	
