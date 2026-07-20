import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.use("Agg")
df = pd.read_csv("EnKF_data/localized_transform_data.csv", header=None, sep=" ")
df1 = pd.read_csv("EnKF_data/PO_data.csv", header=None, sep=" ")
# df2 = pd.read_csv("KF_data/KF_data_delete.csv", header=None, sep=" ")
# df3 = pd.read_csv("EnKF_data/SRF_data.csv", header=None, sep=" ")
df4 = pd.read_csv("EnKF_data/transform_data.csv", header=None, sep=" ")

PO = df1.iloc[0, ::-1]
x = np.arange(len(PO))
# KF = df2.iloc[0, ::-1]
# SRF = df3.iloc[0, ::-1]
transform = df4.iloc[0, ::-1]
localozed_transform = df.iloc[:, ::-1]

plt.figure(figsize=(12, 8))

for i, row in localozed_transform.iterrows():
    plt.plot(x, row, label=f"L={(i+1)*0.5*4.0}")

# plt.plot(x, PO, color="orange", label="EnKF PO")
# plt.plot(x, SRF, color="green", label="EnKF SRF")
# plt.plot(x, KF, color="blue", label="KF")
plt.plot(x, transform, color="red", label="ETKF")
plt.title("EnKF")
plt.xlabel("observation data num")
plt.ylabel("rms")
plt.ylim(top=6)
plt.legend()
plt.savefig("EnKF_data/rms_time_avg.jpg")
plt.show()
