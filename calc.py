import pandas as pd

df = pd.read_csv("KF_data/no_KF_data_rms.csv")

rms_vec = df.iloc[60 * 4 :, 1]
sum = 0
for rms in rms_vec:
    sum += rms
avg = sum / len(rms_vec)
print(avg)
