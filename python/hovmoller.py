import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("plob1.csv", header=None)

fig = plt.figure(num=3, figsize=(10, 12))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# 最初の500ステップ、全40地点を表示
dt = 0.01
point = 40
first_steps = 1000
second_step = 5000
time_steps = 100000
offset = 1  # グラフを上にずらす幅（データの振幅に合わせて調整）

for i in range(0, first_steps, 50):
    ax1.plot(df.iloc[i, :point] * 3 + (i * offset), color="black", linewidth=0.5)

for i in range(first_steps, second_step, 100):
    ax2.plot(df.iloc[i, :point] + (i / 10 * offset), color="black", linewidth=0.5)

for i in range(second_step, time_steps, 5000):
    ax3.plot(
        df.iloc[i, :point] * 1 / 3 + (i / 1000 * offset + 100 * offset),
        color="black",
        linewidth=0.5,
    )

fig.suptitle("Lorenz 96: Spatiotemporal Propagation (Offset Plot)")
ax1.set_title("t=0 ~ t=10")
ax2.set_title("t=10 ~ t=50")
ax3.set_title("t=50 ~ t=1000")
ax1.get_yticks([i * offset for i in range(0, int(first_steps * dt), 50 * dt)])
ax2.get_yticks(
    [i * offset for i in range(int(first_steps * dt), second_step * dt, 100 * dt)]
)
ax3.get_yticks(
    [i * offset for i in range(second_step * dt, time_steps * dt, 5000 * dt)]
)
plt.ylabel("Time Step")
plt.xlabel("Grid Point (with Offset)")
plt.xticks([i * offset for i in range(0, 40, 5)], [f"x{i}" for i in range(0, 40, 5)])
plt.grid(axis="x", linestyle="--", alpha=0.7)

plt.savefig("hovmoller_dialog_f_6.jpg")
plt.tight_layout()
plt.show()
