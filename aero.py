import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Physical Constants (All in SI Units)
DENSITY = 1.18
GRAVITY = 9.78
MASS = 0.0125         
CHORD = 0.20            
WING_AREA = 0.04      
INERTIA = MASS * (CHORD**2) / 12

# Helper Functions for Numerical Safety
def safe_sigmoid(x, k=10.0):
    """
    Sigmoid function with overflow protection.
    Returns 1 / (1 + exp(-k*x))
    """
    val = -k * x
    val = np.clip(val, -50.0, 50.0)
    return 1.0 / (1.0 + np.exp(val))

def smooth_blend(v_curr, v_start, v_end):
    """
    Returns 1.0 when v < v_start, 0.0 when v > v_end, smooth transition.
    """
    if v_curr <= v_start: return 1.0
    if v_curr >= v_end: return 0.0
    t = (v_curr - v_start) / (v_end - v_start)
    return 1.0 - (t * t * (3.0 - 2.0 * t)) # Cubic ease-in-out

class PaperPlane:
    def __init__(self, name, cm0, cd0, cm_alpha, cm_q, aspect_ratio):
        self.name = name
        self.cm0_base = cm0
        self.cd0 = cd0
        self.cm_alpha = cm_alpha
        self.cm_q = cm_q
        self.AR = aspect_ratio
        # Induced drag factor K = 1 / (pi * e * AR)
        self.K = 1.0 / (3.14159 * 0.8 * self.AR)

    def get_forces(self, v, alpha, theta, q_rate, aspect):

        washout = smooth_blend(v, 4.0, 12.0) * 0.8 + 0.2

        trim_scale = smooth_blend(v, 3.0, 8.0)
        
        # AoA-dependent lift coefficient
        cl_base = 1.8 * np.sin(2.0 * alpha)
        
        # Washout reduces lift here (induced, not fixed)
        CL = cl_base * washout
        
        cd_induced = self.K * (CL**2)
        cd_stall = 1.8 * (np.sin(alpha)**2)
        
        CD = self.cd0 + cd_induced + cd_stall

        is_slow = safe_sigmoid(4.0 - v, k=2.0) # 1 if v<4
        is_nose_up = safe_sigmoid(theta - 0.3, k=10.0) # 1 if theta > 0.3 rad
        stall_torque = -0.15 * is_slow * is_nose_up
        
        current_cm0 = self.cm0_base * trim_scale

        damp_denom = max(2.0 * v, 0.5) 
        damping = self.cm_q * (CHORD * q_rate / damp_denom)
        
        Cm = current_cm0 + (self.cm_alpha * alpha) + damping + stall_torque
        
        return CL, CD, Cm

def simulate_flight(plane, v0, theta0, max_time=10.0):
    
    def derivatives(t, state):
        x, y, vx, vy, theta, omega = state

        v = np.hypot(vx, vy)
        
        # Stop if on ground (or slightly underground due to step)
        if y < 0:
            return [0.0]*6

        gamma = np.arctan2(vy, vx)
        # Normalized Angle of Attack
        alpha = (theta - gamma + np.pi) % (2*np.pi) - np.pi
        

        CL, CD, Cm = plane.get_forces(v, alpha, theta, omega, plane.AR)
        
        # Dynamic Pressure * wing area
        Q = 0.5 * DENSITY * v**2 * WING_AREA / plane.AR # True wing area, since WING_AREA is placeholder
        
        # Normalize vx, vy for forces
        vn_x, vn_y = vx/v, vy/v
        
        Lift = Q * CL
        Drag = Q * CD * plane.AR # Wing area vs cross-sectional area scaling
        
        Fx = -Drag * vn_x - Lift * vn_y
        Fy = -Drag * vn_y + Lift * vn_x - MASS * GRAVITY
        
        Moment = Q * CHORD * Cm

        ax = Fx / MASS
        ay = Fy / MASS
        alpha_acc = Moment / INERTIA
        
        return [vx, vy, ax, ay, omega, alpha_acc]

    def hit_ground(t, state):
        return state[1]
    hit_ground.terminal = True
    hit_ground.direction = -1

    state0 = [0.0, 1.5, v0 * np.cos(theta0), v0 * np.sin(theta0), theta0, 0.0]
    
    try:
        sol = solve_ivp(
            derivatives, (0, max_time), state0,
            method='Radau', # Implicit method, stable for stiff aerodynamics
            events=hit_ground,
            rtol=1e-5, atol=1e-6,
            max_step=0.05  
        )
        return sol.y[0], sol.y[1]
    except Exception as e:
        print(f"Solver failed: {e}")
        return [0], [0]

# --- Setup Planes ---
def suzanne(aspect):
    return PaperPlane("Suzanne", cm0=0.03, cd0=0.03*(aspect/3.5 + 3/aspect), cm_alpha=-0.2, cm_q=-3.0, aspect_ratio=aspect)
    # cd0 scaled with this function to reflect wing drag and gap between fuselage layers

def alkonost(aspect):
    return PaperPlane("Alkonost", cm0=0.025, cd0=0.03*(aspect/3 + 3/aspect), cm_alpha=-0.3, cm_q=-4.0, aspect_ratio=aspect)

def super_dart(aspect):
    return PaperPlane("Super Dart", cm0=0.02, cd0=0.02*(aspect/4 + 4/aspect), cm_alpha=-0.3, cm_q=-4.0, aspect_ratio=aspect)

def chinese_glider(aspect):
    return PaperPlane("Chinese Glider", cm0=0.03, cd0=0.03*(aspect/3 + 3/aspect), cm_alpha=-0.2, cm_q=-3.0, aspect_ratio=aspect)

# --- Visualization ---
plt.figure(figsize=(12, 7))

ranges_suzanne = []
for twentyasp in range(15, 101):
    plane = suzanne(twentyasp/20)
    x, y = simulate_flight(plane, v0=20.0, theta0=0.64, max_time=50)
    ranges_suzanne.append([twentyasp/20, x[-1]])

ranges_alkonost = []
for twentyasp in range(15, 101):
    plane = alkonost(twentyasp/20)
    x, y = simulate_flight(plane, v0=20.0, theta0=0.64, max_time=50)
    ranges_alkonost.append([twentyasp/20, x[-1]])

ranges_super_dart = []
for twentyasp in range(15, 101):
    plane = super_dart(twentyasp/20)
    x, y = simulate_flight(plane, v0=20.0, theta0=0.64, max_time=50)
    ranges_super_dart.append([twentyasp/20, x[-1]])

ranges_chinese_glider = []
for twentyasp in range(15, 101):
    plane = chinese_glider(twentyasp/20)
    x, y = simulate_flight(plane, v0=20.0, theta0=0.64, max_time=50)
    ranges_chinese_glider.append([twentyasp/20, x[-1]])

plt.plot([r[0] for r in ranges_suzanne], [r[1] for r in ranges_suzanne], 'b-', linewidth=2)
plt.plot([r[0] for r in ranges_alkonost], [r[1] for r in ranges_alkonost], 'r-', linewidth=2)
plt.plot([r[0] for r in ranges_super_dart], [r[1] for r in ranges_super_dart], 'g-', linewidth=2)
plt.plot([r[0] for r in ranges_chinese_glider], [r[1] for r in ranges_chinese_glider], 'm-', linewidth=2)

plt.title("Range vs Aspect Ratio")
plt.xlabel("Aspect Ratio")
plt.ylabel("Range (m)")
plt.axhline(0, color='black', linewidth=1)
plt.grid(True, linestyle='--')
plt.xlim(0, 6)
plt.legend(['Suzanne', 'Alkonost', 'Super Dart', 'Chinese Glider'])
plt.show()