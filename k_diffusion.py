def k_diffusion(T):
    # Thermal conductivity
    a = 1.52e-4
    b = 4.42e-5
    c = 8e-9
    out = a + b * T + c * T**2  # [W/m.K]
    return out
