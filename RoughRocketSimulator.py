import matplotlib.pyplot as plt
import numpy as np

g = 9.81
R = 287
dt = 0.01

m = 55.3
m_dry = 46.3 #N2 thrust is negligable
thrust = 3000
isp = 228
mdot = thrust / g / isp

T = 0
h = 0
v = 0

TList = []

hList = []
vList = []
machList = []

aList = []
dragList = []
pressureList = []

massList = []

def pressure(h):
    if h < 11000:
        return 1.225 * (1 + -0.0065 * h / 288.15)**(-g / R / -0.0065 - 1)
    elif h < 20000:
        return 0.3636 * np.exp(-g / R / 216.66 * (h - 11000))
    else:
        return 0.0879 * (1 + 0.001 * (h-20000) / 216.66)**(-g / R / 0.001 - 1)



while h >= 0:
    TList.append(T)
    hList.append(h)
    vList.append(v)

    Fdrag = - np.sign(v) * 0.5 * pressure(h) * v * v * 0.049 * np.pi * 0.12 * 0.12
    athrust = 0
    if m > m_dry:
        athrust = thrust / m
        m = m - mdot * dt

    aList.append(-g + Fdrag / m + athrust)
    machList.append(v / np.sqrt(1.4 * R * 288))
    dragList.append(Fdrag)
    pressureList.append(pressure(h))

    h = h + v * dt
    v = v + (-g + Fdrag / m + athrust) * dt

    T = T + dt

fig, axs = plt.subplots(3,2)
axs[0,0].plot(TList, hList)
axs[0,0].set(xlabel='height')
axs[1,0].plot(TList, vList)
axs[1,0].set(xlabel='velocity')
axs[2,0].plot(TList, machList)
axs[2,0].set(xlabel='mach')
axs[0,1].plot(TList, aList)
axs[0,1].set(xlabel='acceleration')
axs[1,1].plot(TList, dragList)
axs[1,1].set(xlabel='drag')
axs[2,1].plot(TList, pressureList)
axs[2,1].set(xlabel='pressure')
plt.show()