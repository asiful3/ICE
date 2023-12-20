import math

def f_haaland_turbulent(epsD):
    # Compute the fully turbulent Darcy friction factor from the Haaland equation
    # 1/sqrt(f) = -1.8 log_10 [ ( (eps/D)/3.7 )^1.11 + 6.9/Re ]
    # Neglect the Re contribution (i.e. Re -> inf) for fully turbulent flow

    out = -1.8 * math.log10((epsD / 3.7) ** 1.11)
    out = 1.0 / out
    out = out ** 2

    return out

# Example usage:
epsD_value = 0.01  # Replace with your specific value for epsD
result = f_haaland_turbulent(epsD_value)
print(result)
