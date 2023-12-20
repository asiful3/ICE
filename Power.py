import numpy as np

def power(V, angle):
    mass = 24000  # total mass of bus
    A = 2.55 * 2  # frontal area of bus
    cr = 0.01  # rolling coefficient
    cd = 0.5  # drag coefficient
    g = 9.81  # gravity acceleration
    Qc = 24  # cooling load
    COP = 14.5  # COP of system
    AC = Qc / COP * 1000  # A/C participation
    V = V * 1000 / 3600  # km/h to m/s
    rho = 1.225  # density of air
    Fa = mass * 0.4 * 0  # acceleration force
    Fr = mass * g * np.cos(np.deg2rad(angle)) * cr  # rolling force
    Fd = 0.5 * cd * A * V**2 * rho  # drag force
    Fg = mass * g * np.sin(np.deg2rad(angle))  # gravity force
    P_b = (V * (Fr + Fd + Fg + Fa) + AC) / 0.9  # brake power
    return Fr, Fd, Fg, Fa, P_b

# Example usage
V = 30  # speed of vehicle in km/h
angle = 5  # road angle in degrees
Fr, Fd, Fg, Fa, P_b = power(V, angle)
print(f"Fr: {Fr}, Fd: {Fd}, Fg: {Fg}, Fa: {Fa}, P_b: {P_b}")
