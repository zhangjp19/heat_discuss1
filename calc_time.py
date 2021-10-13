import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
import h5py
import sys

def calc_time(P_CRIT, P_A_INIT, P_B_INIT, P_B_FIN, T_A_INIT, T_B_INIT, V_A, V_B, q2, t_step):
    s_B2 = PropsSI("S", "P", P_B_INIT, "T", T_B_INIT, "REFPROP::CO2")
    d_B2_INIT = PropsSI("D", "P", P_B_INIT, "T", T_B_INIT, "REFPROP::CO2")
    d_B2_CRIT = PropsSI("D", "P", P_CRIT, "S", s_B2, "REFPROP::CO2")
    t_part1 = V_B * (d_B2_INIT - d_B2_CRIT) / q2

    t = 0
    t_list = []
    P_B2_list = []
    P_A_list = []
    T_B2_list = []
    T_A_list = []

    P_B2 = P_B_INIT * 1
    P_A = P_A_INIT * 1
    d_B2 = d_B2_INIT * 1
    while d_B2 > d_B2_CRIT:
        t_list.append(t)
        P_B2_list.append(P_B2)
        P_A_list.append(P_A)
        T_B2_list.append(PropsSI("T", "P", P_B2, "S", s_B2, "REFPROP::CO2"))
        T_A_list.append(T_A_INIT)
        delta_d_B2 = -q2 * t_step / V_B
        d_B2 += delta_d_B2
        P_B2 = PropsSI("P", "D", d_B2, "S", s_B2, "REFPROP::CO2")
        t += t_step

    t_part2 = 0
    s_A = PropsSI("S", "P", P_A_INIT, "T", T_A_INIT, "REFPROP::HELIUM")
    d_A = PropsSI("D", "P", P_A_INIT, "T", T_A_INIT, "REFPROP::HELIUM")
    h_A = PropsSI("H", "S", s_A, "D", d_A, "REFPROP::HELIUM")
    V_B1 = 0
    d_B1 = d_A * 1
    u_B1 = PropsSI("U", "D", d_B1, "P", P_CRIT, "REFPROP::HELIUM")
    P_A = P_A_INIT * 1
    while (P_A > P_CRIT and V_B1 < V_B):
        delta_V_B1 = q2 / d_B2_CRIT * t_step
        V_B1 += delta_V_B1
        u_B1_partial_d_B1 = (PropsSI("U", "D", d_B1 * 1.01, "P", P_CRIT, "REFPROP::HELIUM") - PropsSI("U", "D", d_B1 * 0.99, "P", P_CRIT, "REFPROP::HELIUM")) / (d_B1 * 0.02)
        delta_d_B1 = (-(u_B1 - h_A) * d_B1 * delta_V_B1 - P_CRIT * delta_V_B1) / ((u_B1 - h_A) * V_B1 + d_B1 * V_B1 * u_B1_partial_d_B1)
        delta_d_A = -(V_B1 * delta_d_B1 + d_B1 * delta_V_B1) / V_A
        d_B1 += delta_d_B1
        d_A += delta_d_A
        h_A = PropsSI("H", "S", s_A, "D", d_A, "REFPROP::HELIUM")
        P_A = PropsSI("P", "S", s_A, "D", d_A, "REFPROP::HELIUM")
        u_B1 = PropsSI("U", "D", d_B1, "P", P_CRIT, "REFPROP::HELIUM")
        t_part2 += t_step
        t_list.append(t)
        P_A_list.append(P_A)
        P_B2_list.append(P_CRIT)
        T_A_list.append(PropsSI("T", "P", P_A, "S", s_A, "REFPROP::HELIUM"))
        T_B2_list.append(PropsSI("T", "P", P_B2, "S", s_B2, "REFPROP::CO2"))
        t += t_step
    
    print("first step t:", t_part1)
    print("second step t:", t_part2)
    print("total t:", t_part1 + t_part2)
    print("V_B1:", V_B1)

    data = {}
    data["t"] = t_list
    data["P_A"] = P_A_list
    data["P_B2"] = P_B2_list
    data["T_A"] = T_A_list
    data["T_B2"] = T_B2_list
    
    return data

P_A_INIT = 50000000
P_B_INIT = 20000000
P_B_FIN = 7500000
T_A_INIT = 293.15
T_B_INIT = 293.15
V_A = 0.032
V_B = 0.015
q2 = 0.1

P_CRIT = float(sys.argv[1]) * 1000000
t_step = 0.1

filename = sys.argv[2]

data = calc_time(P_CRIT, P_A_INIT, P_B_INIT, P_B_FIN, T_A_INIT, T_B_INIT, V_A, V_B, q2, t_step)
print("final T_A:", data["T_A"][-1])

fig, ax = plt.subplots(2,1)

ax[0].plot(data["t"], data["T_A"], label="T_A")
ax[0].plot(data["t"], data["T_B2"], label="T_B2")
ax[1].plot(data["t"], data["P_A"], label="P_A")
ax[1].plot(data["t"], data["P_B2"], label="P_B2")
ax[0].legend()
ax[0].set_xlabel("Time(s)")
ax[0].set_ylabel("Temperature(K)")
ax[1].set_xlabel("Time(s)")
ax[1].set_ylabel("Pressure(Pa)")
ax[1].legend()
ax[0].grid(True)
ax[1].grid(True)
fig.tight_layout()
plt.show()

fig.savefig(filename + ".png", format="png")

with h5py.File(filename + ".h5", "w") as out_file:
    out_file["t"] = data["t"]
    out_file["P_A"] = data["P_A"]
    out_file["P_B2"] = data["P_B2"]
    out_file["T_A"] = data["T_A"]
    out_file["T_B2"] = data["T_B2"]

