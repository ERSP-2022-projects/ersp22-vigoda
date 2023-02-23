# dont edit this here, just copy it, i used jupyter notebook

import numpy as np
import matplotlib.pyplot as plt

xvals = [500, 1000, 1500, 2000]
data = [
    [350, 390, 550, 270, 230, 390, 470, 310, 350, 430],
    [350, 470, 270, 310, 510, 310, 350, 350, 430, 310],
    [390, 390, 270, 350, 350, 270, 430, 350, 390, 470],
    [390, 430, 270, 270, 230, 230, 270, 470, 510, 310]
]

plt.boxplot(data, labels=xvals)
plt.title("trial 2")
plt.xlabel("# of chars")
plt.ylabel("# of generations")
palette = ['r', 'g', 'b', 'y'] # 4 colors, add as needed

for i in range(len(xvals)):
    plt.scatter(i+1, np.mean(data[i]), color='k') # show mean as black dot
    for y in data[i]:
        x = i + 1;
        x += np.random.normal(0,0.04) # optional: add a bit of horizontal scatter
        plt.scatter(x, y, alpha=0.4, color=palette[i])

plt.show()