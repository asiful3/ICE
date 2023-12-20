import numpy as np

def geometry_theta(theta, B, L, r, R, Ach, Ap, N, degrees=0):
    """
    Calculate various geometric parameters for an internal combustion engine.

    Inputs:
    theta: crank angle [radians];
    B:     Bore [length];
    L:     Stroke [length];
    r:     Compression ratio [unitless];
    R:     Geometry ratio = l/a [unitless];
    Ach:   Area of the cylinder head [length^2];
    Ap:    Area of the piston crown [length^2];
    N:     Rotational crank speed [revs per min (RPM)];
    degrees: Optional input, default value = 0 (theta given in radians).
             Set degrees=1 if theta is specified in degrees.

    Ouputs:
    V:     instantaneous cylinder volume [length^3];
    dV:    instantaneous rate of change of cylinder volume, dV/d(theta) [length^3/radian];
    A:     instantaneous area [length^2];
    dA:    instantaneous rate of change of cylinder area, dA/d(theta) [length^3/radian];
    Up:    instantaneous piston velocity [length/second];
    """
    # if needed, convert crank angle to radians
    T_rad = np.radians(theta) if degrees else theta

    # calculate displacement volume
    Vd = np.pi * B**2 * L / 4

    phi = np.sqrt(R**2 - np.sin(T_rad)**2)
    V = (Vd / (r - 1)) * (1 + 0.5 * (r - 1) * (R + 1 - np.cos(T_rad) - np.sqrt(R**2 - np.sin(T_rad)**2)))
    dV = 0.5 * Vd * np.sin(T_rad) * (1 + np.cos(T_rad) / np.sqrt(R**2 - np.sin(T_rad)**2))
    A = Ach + Ap + (np.pi * B * L / 2) * (R + 1 - np.cos(T_rad) - np.sqrt(R**2 - np.sin(T_rad)**2))
    dA = (np.pi * B * L / 2) * np.sin(T_rad) * (1 + np.cos(T_rad) / np.sqrt(R**2 - np.sin(T_rad)**2))
    Up = np.pi * L * (N / 60) * np.sin(T_rad) * (1 + np.cos(T_rad) / np.sqrt(R**2 - np.sin(T_rad)**2))

    return V, dV, A, dA, Up

# Example usage:
theta_value = 30  # Replace with your specific value for theta
B_value = 0.1  # Replace with your specific value for B
L_value = 0.15  # Replace with your specific value for L
r_value = 10  # Replace with your specific value for r
R_value = 5  # Replace with your specific value for R
Ach_value = 0.2  # Replace with your specific value for Ach
Ap_value = 0.3  # Replace with your specific value for Ap
N_value = 2000  # Replace with your specific value for N

result = geometry_theta(theta_value, B_value, L_value, r_value, R_value, Ach_value, Ap_value, N_value)
print(result)
