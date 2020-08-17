import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates
# Open .bed file
content = []
with open("merged_simple.bed")as f:
    for line in f:
        content.append(line.strip().split())
# Set first row as the header
df = pd.DataFrame(content)
df.columns = df.iloc[0]
df = df[1:]
# Remove "RepeatUnit" duplications
df = df.loc[:,~df.columns.duplicated()]
# Create "ID" column
df["Loci"] = df[["#Chromosome", "start", "stop"]].apply(lambda x: '  '.join(x), axis=1)

# Index "AnchoredIrrCount"
df.set_index(['#Chromosome','start','stop','Loci','RepeatUnit', 'AnchoredIrrCount_fa','AnchoredIrrCount_mo','AnchoredIrrCount_p1', 'AnchoredIrrCount_s1'],inplace=True)
# Remove unused columns
df = df.drop(columns=['AnchorCopies_fa','AnchorCopies_mo','AnchorCopies_p1','AnchorCopies_s1'])
# Change strings to floats
df = df.astype(float)
# Remove rows where any "IrrPairCount" != 0
df = df[~(df != 0).any(axis=1)]
# Change "AnchoredIrrCount" from index to columns
df.reset_index(level=['AnchoredIrrCount_fa','AnchoredIrrCount_mo','AnchoredIrrCount_p1','AnchoredIrrCount_s1'], inplace=True)
df = df.astype(float)
# Remove unused columns
d3 = df.drop(columns=['IrrPairCount_fa','IrrPairCount_mo','IrrPairCount_p1','IrrPairCount_s1'])
# Remove all values that are between the str.fa and str.mo size for s1
d1 = df.drop(df[(df.AnchoredIrrCount_p1 <= df['AnchoredIrrCount_mo']) | (df.AnchoredIrrCount_p1 <= df['AnchoredIrrCount_fa'])].index)
d2 = df.drop(df[(df.AnchoredIrrCount_p1 >= df['AnchoredIrrCount_mo']) | (df.AnchoredIrrCount_p1 >= df['AnchoredIrrCount_fa'])].index)
d3 = pd.concat([d1,d2])
# Remove all values that are between the str.fa and str.mo size for p1
d4 = df.drop(df[(df.AnchoredIrrCount_s1 <= df['AnchoredIrrCount_mo']) | (df.AnchoredIrrCount_s1 <= df['AnchoredIrrCount_fa'])].index)
d5 = df.drop(df[(df.AnchoredIrrCount_s1 >= df['AnchoredIrrCount_mo']) | (df.AnchoredIrrCount_s1 >= df['AnchoredIrrCount_fa'])].index)
d6 = pd.concat([d4,d5])
# Move 'ID' and 'RepeatUnit' back to columns
d3.reset_index(level= ['#Chromosome','start','stop','Loci', 'RepeatUnit'], inplace=True)
d6.reset_index(level= ['start','stop','#Chromosome','Loci', 'RepeatUnit'], inplace=True)
d3 = pd.concat([d3,d6]).drop_duplicates()
########
d3 = d3.iloc[range(0,50)]
#######
# Set index
d3.index = range(len(d3['Loci']))
d3.reset_index(level=0, inplace=True)
# Scale index
Dataframe_max = max(list(df.max(axis = 0)))
max_index = (max(d3['index']))
scale_value = np.ceil(Dataframe_max / max_index)
# Add "ID" to y axis
d3['index'] = d3['index'] * scale_value
my_yticks =  list(d3.Loci)
y = d3['index']
plt.yticks(y, my_yticks)
# Create parallel plot
ax = parallel_coordinates(d3, 'RepeatUnit',cols = ['index','AnchoredIrrCount_fa','AnchoredIrrCount_mo','AnchoredIrrCount_p1','AnchoredIrrCount_s1'],\
                          colormap=plt.get_cmap("Set3"))
ax.grid(False)
# Add second yaxis
plt.subplots_adjust(wspace=0)
ax2= ax.twinx()
y2 = d3['index']
ax2.plot(y2, color = 'white', alpha=0)
ax2.grid(True)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height * 1])
ax2.set_position([box.x0, box.y0, box.width * 0.9, box.height *1])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), prop={'size': 9})
#plt.title('p1 STR Genotypes', size = 30, y = 1.04)
plt.title('Short Repeat Genotypes', size = 30, y = 1.04)
plt.ylabel('Number of Anchored IRR Reads', size = 15, x = 1, rotation = 270)
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
figure = plt.gcf() # get current figure
figure.set_size_inches(20, 15)
# when saving, specify the DPI
plt.savefig("short_repeat.png", dpi = 100)
plt.show()
#print (d3)
