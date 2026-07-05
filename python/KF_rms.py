import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

matplotlib.use('Agg')

df = pd.read_csv("EnKF_data/PO_rms.csv", header=None, sep=",")
plt.figure(figsize=(12, 8))
sum = 0
x = df.iloc[:, 0]
rms = df.iloc[:, 1]
avg=rms.mean()
print(avg)
plt.plot(x, rms, label="kf_true_rms")
plt.xlabel("time_step")
plt.ylabel("rms")
plt.savefig("EnKF_data/PO_mse.jpg")
plt.show()
