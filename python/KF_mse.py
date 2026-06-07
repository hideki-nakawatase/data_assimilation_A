import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("KF_data/KF_rms_delta.csv", header=None, sep=",")
plt.figure(figsize=(12, 8))
sum = 0
x = df.iloc[:, 0]
rms = df.iloc[:, 1]
plt.plot(x, rms, label="kf_true_rms")
plt.xlabel("time_step")
plt.ylabel("rms")
plt.savefig("KF_data/KF_rms_delta.jpg")
plt.show()
