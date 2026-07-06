import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df1 = pd.read_csv(
    "3D_variational_method_data/random_delete_rms_avg.csv", header=None, sep=" "
)
df2 = pd.read_csv("KF_data/KF_data_delete.csv", header=None, sep=" ")

variational_method = df1.iloc[0, ::-1]
x = np.arange(len(variational_method))
KF = df2.iloc[0, ::-1]

plt.figure(figsize=(12, 8))
plt.plot(x, variational_method,color='orange', label="3D_variation_method")
plt.plot(x, KF,color='blue', label="KF")
plt.title("KF, 3D variational method with insufficient observation")
plt.xlabel("observation data num")
plt.ylabel("rms")
plt.ylim(top=6)
plt.legend()
plt.savefig("3D_variational_method_data/3D_variation_KF_delete.jpg")
plt.show()
