import math

def brake_power(V_v, inclination, m_v, A_v):
    # Convert velocity to m/s
    V = V_v / 3.6

    # Convert inclination to radians
    inclination = inclination * math.pi / 180

    # Constants
    C_r = 0.01  # Rolling resistance coeff.
    C_d = 0.5   # Aerodynamic drag coeff.
    rho_air = 1.225  # Air density [kg/m^3]
    g = 9.81  # Gravity [m/s^2]

    # Rolling resistance force
    F_rr = C_r * m_v * g * math.cos(inclination)

    # Grading resistance force
    F_g = m_v * g * math.sin(inclination)

    # Aerodynamic drag force
    F_d = (1/2) * C_d * rho_air * (V**2) * A_v

    # Total resistance force
    F_T = F_rr + F_g + F_d

    # Brake power [kW]
    P_b = (F_T * V) / 1000

    print('P_b:', P_b)

# Example usage:
V_v = 60  # velocity in km/h
inclination = 5  # inclination in degrees
m_v = 1500  # vehicle mass in kg
A_v = 2.5  # vehicle frontal area in m^2

brake_power(V_v, inclination, m_v, A_v)
