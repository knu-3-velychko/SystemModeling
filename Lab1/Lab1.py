# %%
import numpy as np
# %%
import matplotlib.pyplot as plt
# %%
with open('f19.txt') as file:
    inputData = np.array([float(val) for val in file.read().split()])

T = 5
dt = 0.01
time = np.arange(0, T+dt, dt)
plt.grid(True)
plt.rc('figure', figsize=(30.0, 10.0))
plt.plot(time, inputData)
# %%
n = time.shape[0]
frequency = np.zeros(n)
for pointID in range(n):
    sinFrequency = 0
    cosFrequency = 0
    for signal in range(n):
        sinFrequency += inputData[signal] * \
            np.sin(2.*np.pi*pointID*signal/float(n))
        cosFrequency += inputData[signal] * \
            np.cos(2.*np.pi*pointID*signal/float(n))
    sinFrequency /= float(n)
    cosFrequency /= float(n)
    frequency[pointID] = np.sqrt(sinFrequency**2+cosFrequency**2)
plt.grid(True)
plt.plot(time, frequency)

# %%
biggestValue = []
for i in range(3, n // 2):
    if np.max(frequency[i-3:i+3]) == frequency[i]:
        biggestValue.append(i)
        print(frequency[i])
mainFrequency = biggestValue[0]/T
#%%
b = np.array([np.sum(inputData * time ** 3), np.sum(inputData * time ** 2), np.sum(inputData * time),
              np.sum(inputData * np.sin(2. * np.pi * mainFrequency * time)), np.sum(inputData)])