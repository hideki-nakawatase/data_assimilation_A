import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("-A/lyapunov.csv",header=None)

fig = plt.figure(num=3, figsize=(10, 12))
ax1 = fig.add_subplot(111)

dt = 0.01
point = 40
time_steps = 400
x=range(time_steps)
y=df.iloc[:time_steps, 0]

ax1.plot(x,y,color='black',linewidth=0.5,label='Node 20')

fig.suptitle("Lorenz 96: lyapunov progression")

plt.xlabel("Time Step")

plt.savefig("lyapunov.jpg")
plt.show()
