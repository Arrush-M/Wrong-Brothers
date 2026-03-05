import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math

# --- Physical Constants (SI Units) ---
# Realistic parameters for a high-performance folded paper plane
DENSITY = 1.225
GRAVITY = 9.81
MASS = 0.0125            # 12.5 grams
CHORD = 0.10            # 10 cm mean chord
WING_AREA = 0.017       # Projected area
INERTIA = MASS * (CHORD**2) / 10.0 # Estimate

# --- Helper Functions for Numerical Safety ---
def safe_sigmoid(x, k=10.0):
    """
    Sigmoid function with overflow protection.
    Returns 1 / (1 + exp(-k*x))
    """
    # Clamp input to avoid exp overflow
    # np.exp(700) is approx limit. 
    # k*x should be within [-50, 50] for valid float range of 0-1
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

    def get_forces(self, v, alpha, theta, q_rate):
        # 1. --- Aeroelastic Deformation ---
        # As velocity increases, wings flatten/twist (washout).
        # This reduces the Lift Slope and Pitch Trim.
        # Washout: 1.0 at low speed, approaches 0.2 at high speed.
        
        
        # Lift Scale: Prevents loop by reducing lift at high speed (15m/s)
        # Transition from 1.0 to low value between 4m/s and 10m/s
        washout = smooth_blend(v, 4.0, 12.0) * 0.8 + 0.2
        
        # Trim Scale: Prevents continuous nose-up at high speed
        # Transition trim to 0 as speed exceeds 5 m/s
        trim_scale = smooth_blend(v, 3.0, 8.0)
        
        # 2. --- Lift Coefficient (CL) ---
        # Continuous function to avoid solver discontinuities.
        # Approximation of linear lift + stall.
        # Low aspect ratio lift slope ~ 3.0. Max CL ~ 1.0.
        # sin(2*alpha) is a good geometric approximation for plates.
        cl_base = 1.8 * np.sin(2.0 * alpha)
        
        # Apply washout
        CL = cl_base * washout
        
        # 3. --- Drag Coefficient (CD) ---
        # Parasitic + Induced
        cd_induced = self.K * (CL**2)
        # Separation/Stall drag: dominates at high alpha
        # sin(alpha)^2 fits flow separation well
        cd_stall = 1.8 * (np.sin(alpha)**2)
        
        CD = self.cd0 + cd_induced + cd_stall
        
        # 4. --- Pitching Moment (Cm) ---
        # Static Stability + Damping + Trim
        
        # Stall Recovery Nudge:
        # If speed is low (<4 m/s) and nose is high (>20 deg),
        # simulate center of pressure shift forcing nose down.
        # Use safe sigmoid for smooth transition.
        is_slow = safe_sigmoid(4.0 - v, k=2.0) # 1 if v<4
        is_nose_up = safe_sigmoid(theta - 0.3, k=10.0) # 1 if theta > 0.3 rad
        stall_torque = -0.15 * is_slow * is_nose_up
        
        current_cm0 = self.cm0_base * trim_scale
        
        # Damping (opposes q)
        # Limit denominator to avoid div/0
        damp_denom = max(2.0 * v, 0.5) 
        damping = self.cm_q * (CHORD * q_rate / damp_denom)
        
        Cm = current_cm0 + (self.cm_alpha * alpha) + damping + stall_torque
        
        return CL, CD, Cm

def simulate_flight(plane, v0, theta0, max_time=10.0):
    
    def derivatives(t, state):
        x, y, vx, vy, theta, omega = state
        
        # Compute velocity magnitude safely
        v = np.hypot(vx, vy)
        
        # Stop if on ground (or slightly underground due to step)
        if y < 0:
            return [0.0]*6
            
        # Stop if stopped
        if v < 0.1:
            return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        gamma = np.arctan2(vy, vx)
        # Normalized Angle of Attack
        alpha = (theta - gamma + np.pi) % (2*np.pi) - np.pi
        
        # Aerodynamics
        CL, CD, Cm = plane.get_forces(v, alpha, theta, omega)
        
        # Dynamic Pressure
        Q = 0.5 * DENSITY * v**2 * WING_AREA
        
        # Forces (Inertial Frame)
        # Lift acts perpendicular to Velocity (-vy, vx)
        # Drag acts opposite to Velocity (-vx, -vy)
        # Normalize direction safely
        vn_x, vn_y = vx/v, vy/v
        
        Lift = Q * CL
        Drag = Q * CD
        
        Fx = -Drag * vn_x - Lift * vn_y
        Fy = -Drag * vn_y + Lift * vn_x - MASS * GRAVITY
        
        Moment = Q * CHORD * Cm
        
        # Accelerations
        ax = Fx / MASS
        ay = Fy / MASS
        alpha_acc = Moment / INERTIA
        
        return [vx, vy, ax, ay, omega, alpha_acc]

    # Ground event check
    def hit_ground(t, state):
        return state[1]
    hit_ground.terminal = True
    hit_ground.direction = -1

    # Solver
    state0 = [0.0, 1.5, v0 * np.cos(theta0), v0 * np.sin(theta0), theta0, 0.0]
    
    try:
        sol = solve_ivp(
            derivatives, (0, max_time), state0,
            method='Radau', # Implicit method, stable for stiff aerodynamics
            events=hit_ground,
            rtol=1e-5, atol=1e-6,
            max_step=0.05   # Limit step size to catch turns
        )
        return sol.y[0], sol.y[1]
    except Exception as e:
        print(f"Solver failed: {e}")
        return [0], [0]

# --- Setup Planes ---
def suzanne(aspect):
    return PaperPlane("Suzanne", cm0=0.03, cd0=0.02, cm_alpha=-0.2, cm_q=-3.0, aspect_ratio=aspect)

# Alkonost: More stable, faster recovery.
def alkonost(aspect):
    return PaperPlane("Alkonost", cm0=0.02, cd0=0.025, cm_alpha=-0.3, cm_q=-4.0, aspect_ratio=aspect)

# --- Visualization ---
plt.figure(figsize=(12, 7))

# 1. Hard Throws (15 m/s)
# Should climb, slow down ballistically (no loop), then glide.
x1, y1 = simulate_flight(suzanne(1.6), v0=25.0, theta0=np.pi/6, max_time=25)
plt.plot(x1, y1, 'b-', linewidth=2, label="Suzanne Hard (25 m/s)")

# 2. Soft Throws (Optimization Check)
# Should float gently.
x3, y3 = simulate_flight(suzanne(1.6), v0=5.0, theta0=np.pi/6, max_time=15)
plt.plot(x3, y3, 'g--', linewidth=2, label="Suzanne Soft (5 m/s)")

plt.title("Paper Aeroplane Trajectory")
plt.xlabel("Distance (m)")
plt.ylabel("Height (m)")
plt.axhline(0, color='black', linewidth=1)
plt.grid(True, linestyle='--')
plt.axis('equal')
plt.legend()
plt.ylim(bottom=-0.5)
plt.show()