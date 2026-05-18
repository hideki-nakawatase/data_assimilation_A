import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde, norm

df=pd.read_csv("random/random.csv",header=None)
y=df.iloc[:,0]
y=sorted(y)

kde=gaussian_kde(y)

x_grid=np.linspace(-4,4,1000)

y_kde=kde(x_grid)

y_theoretical = norm.pdf(x_grid, loc=0, scale=1)

plt.plot(x_grid,y_kde)
plt.plot(x_grid,y_theoretical)
plt.show()