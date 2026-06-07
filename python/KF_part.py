import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df1 = pd.read_csv("observation_data/observation_data.csv", header=None, sep=" ")
df2 = pd.read_csv("true_data/true_data.csv", header=None, sep=" ")
df3 = pd.read_csv("KF_data/KF_data.csv", header=None, sep=" ")
df4 = pd.read_csv("KF_data/no_KF_data.csv", header=None, sep=" ")

y1 = df1.iloc[:1460, 19]
y2 = df2.iloc[:1460, 19]
y3 = df3.iloc[:1460, 19]
y4 = df4.iloc[:1460, 19]


x1 = np.arange(0, 75, 0.25)
x2 = np.arange(75, 150, 0.25)
x3 = np.arange(150, 225, 0.25)
x4 = np.arange(225, 300, 0.25)
x5 = np.arange(300, 365, 0.25)


fig = plt.figure(num=4, figsize=(14, 8))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax1.plot(x1, y2[:300], label="true")
ax1.plot(x1, y3[:300], label="KF")
ax1.plot(x1, y4[:300], label="no_KF")
ax1.set_title("~75days")

ax2.plot(x2, y2[300:600])
ax2.plot(x2, y3[300:600])
ax2.plot(x2, y4[300:600])
ax2.set_title("75~150days")

ax3.plot(x3, y2[600:900])
ax3.plot(x3, y3[600:900])
ax3.plot(x3, y4[600:900])
ax3.set_title("150~225days")

ax4.plot(x4, y2[900:1200])
ax4.plot(x4, y3[900:1200])
ax4.plot(x4, y4[900:1200])
ax4.set_title("225~300days")

fig.suptitle("x_20")
fig.legend()
fig.show()
fig.savefig("KF_data/KF_data_with_no_KF.png")


fig2 = plt.figure(figsize=(14, 8))
ax1 = fig2.subplots()
ax1.plot(x1, y2[:300], label="true")
ax1.plot(x1, y3[:300], label="KF")
ax1.plot(x1, y4[:300], label="no_KF")

fig2.suptitle("x_20")
fig2.legend()
fig2.show()
fig2.savefig("KF_data/KF_data_75days_with_no_KF.png")

fig3 = plt.figure(figsize=(14, 8))
ax = fig3.subplots()
ax.plot(x5, y2[1200:1460], label="true")
ax.plot(x5, y3[1200:1460], label="KF")
ax.plot(x5, y4[1200:1460], label="no_KF")

fig3.suptitle("x_20")
fig3.legend()
fig3.show()
fig3.savefig("KF_data/KF_data_last_65days_with_no_KF.png")
