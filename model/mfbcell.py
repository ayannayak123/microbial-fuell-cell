import numpy as np
import matplotlib.pyplot as plt

# Constants
F = 96485  # Faraday constant, C/mol
n = 4      # electrons per mol of substrate oxidized (simplified)

# Parameters
time = np.linspace(0, 100, 1000)  # hours
dt = time[1] - time[0]

# Initial conditions
COD0 = 1000  # mg/L
X0 = 100     # mg/L biomass
S = COD0
X = X0

# Model parameters
mu_max = 0.4       # 1/hr, max specific growth rate
Ks = 100           # mg/L, half-saturation constant
Y = 0.5            # mg biomass / mg COD
b = 0.02           # 1/hr, decay coefficient
qe = 0.01          # mol e- per mg COD removed
R_int = 10         # ohm, internal resistance
R_ext = 100        # ohm, external resistance

# Store values
S_arr, X_arr, I_arr, V_arr = [], [], [], []

for t in time:
    mu = mu_max * S / (Ks + S)  # Monod equation
    dX = (mu * X - b * X) * dt
    dS = -(1 / Y) * mu * X * dt
    
    X += dX
    S += dS
    S = max(S, 0)

    # Electrons released (mol/s)
    re = -dS * qe / dt  # mol e-/L/hr
    re_Cs = re * 1e-3 * 3600  # mol e-/L/s

    # Current: I = n * F * re
    I = n * F * re_Cs  # A

    # Voltage from Ohm's Law: V = I * R_ext
    V = I * R_ext

    # Store
    S_arr.append(S)
    X_arr.append(X)
    I_arr.append(I)
    V_arr.append(V)

# Plotting
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(time, S_arr, label='COD (mg/L)', color='blue')
plt.xlabel('Time (h)')
plt.ylabel('COD')
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(time, X_arr, label='Biomass (mg/L)', color='green')
plt.xlabel('Time (h)')
plt.ylabel('Biomass')
plt.grid(True)

plt.subplot(2, 2, 3)
plt.plot(time, I_arr, label='Current (A)', color='orange')
plt.xlabel('Time (h)')
plt.ylabel('Current')
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(time, V_arr, label='Voltage (V)', color='red')
plt.xlabel('Time (h)')
plt.ylabel('Voltage')
plt.grid(True)

plt.tight_layout()
plt.suptitle('Realistic Microbial Fuel Cell Model', fontsize=16, y=1.02)
plt.show()
