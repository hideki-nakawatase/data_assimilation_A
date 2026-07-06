import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')

df1 = pd.read_csv(
    "EnKF_data/PO_data.csv", header=None, sep=" "
)
df2 = pd.read_csv("KF_data/KF_data_delete.csv", header=None, sep=" ")
df3 = pd.read_csv(
    "EnKF_data/SRF_data.csv", header=None, sep=" "
)

PO = df1.iloc[0, ::-1]
x = np.arange(len(PO))
KF = df2.iloc[0, ::-1]
SRF = df3.iloc[0,::-1]

plt.figure(figsize=(12, 8))
plt.plot(x, PO,color='orange', label="EnKF PO")
plt.plot(x, SRF,color='green', label="EnKF SRF")
plt.plot(x, KF,color='blue', label="KF")
plt.title("EnKF PO method")
plt.xlabel("observation data num")
plt.ylabel("rms")
plt.ylim(top=6)
plt.legend()
plt.savefig("EnKF_data/PO_rms_time_avg.jpg")
plt.show()
