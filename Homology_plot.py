import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('homology.csv')
df = df.pivot(index='sequence', columns='score', values='homology').plot(kind='bar' , figsize=(15,8))


plt.xlabel('align test')
plt.ylabel('Homology')
plt.title('homology on different score')
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()
