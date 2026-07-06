import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df1 = pd.read_csv(
    "EnKF_data/PO_data.csv", header=None, sep=" "
)
df2 = pd.read_csv("KF_data/KF_data_delete.csv", header=None, sep=" ")

variational_method = df1.iloc[0, ::-1]
x = np.arange(len(variational_method))
KF = df2.iloc[0, ::-1]

plt.figure(figsize=(12, 8))
plt.plot(x, variational_method,color='orange', label="EnKF PO")
plt.plot(x, KF,color='blue', label="KF")
plt.title("EnKF PO method")
plt.xlabel("observation data num")
plt.ylabel("rms")
plt.ylim(top=6)
plt.legend()
plt.savefig("EnKF_data/PO_rms_time_avg.jpg")
plt.show()
