This fork implements:
  - the possibility to have a more realistic solar powered electrical propulsion, where the maximum thrust delivered by the engines depends on the distance of the spacecraft to the sun.
  - a realistic modelling of the ion propulsion engine, with the thrust and the Isp depending on the power delivered to the engine and its duty cycle. A list of the more common ion propulsion engine and their models can be found under pykep/trajopt/motor/ motor_factory.py
  - in PyKEP/trajopt, the mga_return_lt_nep implementes the problem (to be optimized with PyGMO) of a sample return mission, using low thrust propulsion. It also uses the newly define get_points.py methods to generate a structured file containing all the optimized trajectory parameters (position and mass as a function of time, engine operating conditions, power consumption and production, planets position, details on legs, etc...) 
  - those data file can be used with the postProcessingToolbox.py which uses plotly in Jupyter Notebook to display iteractif plots on specific trajectories, and also on analysis on a set of possible trajectories.
  

PyKEP
=====

PyKEP is a scientific library providing basic tools to perform research in interplanetary trajectory design or, more generally, in mission analysis. Algorithmic efficiency is a main focus of the library, which is written in C++ and exposed to Python using the boost::python library. At the library core is the implementation of an efficient solver for the multiple revolutions Lambert’s problem, objects representing the Sims-Flanagan low-thrust model, efficient Keplerian propagators, Taylor-integrators and more ...

Check the official documentation at http://esa.github.io/pykep/
