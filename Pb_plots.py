import numpy as np
import matplotlib.pyplot as plt

def Power(V_v, inclination):
    # Your Power function implementation goes here
    # Replace this with your actual implementation
    Fr = Fd = Fg = Fa = P_b = 0
    return Fr, Fd, Fg, Fa, P_b

angle = np.linspace(0, 20, 5)
brake_power = np.zeros((150, len(angle)))

for i in range(1, 151):
    for j, a in enumerate(angle):
        Fr, Fd, Fg, Fa, P_b = Power(i, a)
        brake_power[i - 1, j] = P_b

brake_power = brake_power / 1000
velo_vect = np.linspace(1, 150, 150)
_, _, _, _, P_b = Power(30, 5)
velo_vect_2 = np.linspace(1, 150, 2)
P_b_vect = np.array([P_b, P_b]) / 1000

plt.figure()
plt.plot(velo_vect, brake_power[:, 0], '-r', velo_vect, brake_power[:, 1], '-b', velo_vect, brake_power[:, 2], '-g',
         velo_vect, brake_power[:, 3], 'y', velo_vect, brake_power[:, 4], 'm', velo_vect_2, P_b_vect, '-c')
plt.title('Power vs Speed')
plt.xlabel('Speed [km/h]')
plt.ylabel('Power [kW]')
plt.legend(['0 deg', '5 deg ', '10 deg', '15 deg', '20 deg', 'ref Power'])
plt.show()
