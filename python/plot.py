import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df1=pd.read_csv("observation_data/observation_data.csv",header=None,sep=" ")
df2=pd.read_csv("true_data/true_data.csv",header=None,sep=" ")

y1=df1.iloc[:100,19]
y2=df2.iloc[:100,19]

x=np.arange(0,25,0.25)

plt.figure(figsize=(12,8))
plt.plot(x,y1,label="observation")
plt.plot(x,y2,label="true")

plt.legend()
plt.show()