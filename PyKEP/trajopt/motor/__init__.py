from PyKEP import __extensions__

if (__extensions__['pygmo']):
    from PyKEP.trajopt.motor._motor_factory import *
