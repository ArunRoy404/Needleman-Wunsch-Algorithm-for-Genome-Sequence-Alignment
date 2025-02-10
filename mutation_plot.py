import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('nucleotide.csv')
df = df.pivot(index='sequence', columns='score', values='nucleotide_mutation').plot(kind='bar' , figsize=(15,8))


plt.xlabel('align test')
plt.ylabel('Nucleotide')
plt.title('Mutated Nucleotide on different score')
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()
