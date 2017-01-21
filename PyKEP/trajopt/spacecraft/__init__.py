from PyKEP import __extensions__

if (__extensions__['pygmo']):
    from PyKEP.trajopt.spacecraft._spacecraft_factory import *
