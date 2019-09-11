# To add a new cell, type '#%%'
# To add a new markdown cell, type '#%% [markdown]'
#%%
print("It works");
#%%
import matplotlib.pyplot as plt
#%%
import numpy as np
#%%
with open('f19.txt') as file:
    inputData = np.array([float(val) for val in file.read().split()]) 


time = np.arange(0,5.01, 0.01)
plt.grid(True)
plt.rc('figure', figsize=(30.0, 10.0))
plt.plot(time,inputData);
#%%
