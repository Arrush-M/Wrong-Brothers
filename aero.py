import math
# ALL UNITS IN SI

dx = 0.00001
L = 0.2 # wing length
plane = 0 # 0: suzanne, 1: super dart, 2: alkonost, 3: chinese glider
density = 1.054 # from altitude and temp
viscosity = 0.00001852 # from temp

def f(x: float): # top outline of aerofoil
    return math.e * x * math.exp(-20.0*x) / 1.25

def derivative(x):
    return (f(x+dx) - f(x))/dx

def arc_length():
    t = 0.0 # parameter equal to x, not time
    arc = 0.0
    while t < L:
        arc += math.sqrt(1 + (derivative(t)*derivative(t)))*dx
        t += dx
    return arc

def lift(plane, airspeed, aspect): # Bernoulli effect lift
    lift_P = density*airspeed*airspeed*((arc_length()/L)*(arc_length()/L) - 1)
    lift_F = 0.0
    if plane == 0:
        lift_F = lift_P*L*(L*aspect)/1.6 # L*aspect = wingspan, denominators DUMMY for now
    elif plane == 1:
        lift_F = lift_P*L*(L*aspect)/2.0
    elif plane == 2:
        lift_F = lift_P*L*(L*aspect)/1.9
    elif plane == 3:
        lift_F = lift_P*L*(L*aspect)/1.4
    else:
        raise ValueError("plane must be in {0,1,2,3}")
    return lift_F

def skin_friction(plane, airspeed, aspect):
    Re = density*airspeed*L/viscosity
    c_f = 1.328/math.sqrt(Re)
    fric_f = 0.0
    if plane == 0:
        fric_f = 0.5*c_f*density*airspeed*airspeed*L*(L*aspect)/1.6
    elif plane == 1:
        fric_f = 0.5*c_f*density*airspeed*airspeed*L*(L*aspect)/2.0
    elif plane == 2:
        fric_f = 0.5*c_f*density*airspeed*airspeed*L*(L*aspect)/1.9
    elif plane == 3:
        fric_f = 0.5*c_f*density*airspeed*airspeed*L*(L*aspect)/1.4
    else:
        raise ValueError("plane must be in {0,1,2,3}")
    return fric_f

def vsquared_drag(airspeed, aspect, wing_angle):
    area_cs_wing = L*aspect*0.00018 # paper thickness DUMMY
    area_cs_fuselage = L*math.sin(wing_angle)*0.001 # fuselage thickness DUMMY
    return 0.5*airspeed*airspeed*density*(area_cs_fuselage*0.04 + area_cs_wing*0.09) #0.5Cd*rho*A*v^2 for each part
