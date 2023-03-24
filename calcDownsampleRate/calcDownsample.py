#!/usr/bin/python3

import pandas as pd
import sys

#print(str(sys.argv))
table = sys.argv[1]
ex = sys.argv[2]
s = sys.argv[3]
name = sys.argv[4]
#name = "test.txt"
#table = "Idxstats.txt"
#ex = 'chrHsv1_s17'
#s = 'all'
df = pd.read_table(table)
#df2 = pd.read_table(table)
exclude = ex.split(',')
samples = s.split(',')
df = df[~df.contigt.isin(exclude)]
df = df[df.contigt!='*']
if s!="all":
    df = df[df['sample'].isin(samples)]
sampleList = list(set(df['sample']))
map = []
for sample in sampleList:
    map = map + [df.loc[df['sample'] == sample, 'mapped'].sum()]
    #print(df.loc[df['sample'] == sample, 'mapped'].sum())
print(map)
min = min(map)
#df = df2[~df2.contigt.isin(exclude)]
#df = df[df.contigt!='*']
factor = min/map
data = {'sample': sampleList,
        'factor': factor}
dat = pd.DataFrame(data)
print(dat)
#excludes = df2[df2.contigt.isin(exclude)]
#excludes['factor'] = 1
#all = pd.concat([df, excludes],axis=0)

dat.to_csv(name, header=True, index=None, sep='\t')