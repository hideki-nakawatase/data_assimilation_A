import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("plob2_data\lyapunov_10_points.csv",header=None,sep=' ')

fig = plt.figure(num=3, figsize=(10, 12))
ax1 = fig.add_subplot(111)

dt = 0.005
point = 40
time_steps = 99999
x=range(time_steps)
y=[]
for i in range(10):
  y.append(df.iloc[:time_steps, i])
y_avg=df.iloc[:time_steps,10]

for y_i in y:
  ax1.plot(x,y_i,ls='-.',linewidth=0.5)

ax1.plot(x,y_avg,linewidth=1.0)
fig.suptitle("Lorenz 96: lyapunov progression")

plt.xlabel("Time Step")

plt.savefig("plob2_data\lyapunov_avg.jpg")
plt.show()
