import pandas as pd
import numpy as np
import re
from repDNA.nac import Kmer

# kmer = Kmer(k=2)
# print(kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']))
# kmer = Kmer(k=2, upto=True)
# print( kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']))

# 整理hairpin.fa的数据集，提取出Homo sapiens的基因信息
data = pd.read_csv("hairpin.fa", header=None)
df = np.array(data)
# print(df)
df = df.flatten()
# print('---------------')
# print(df)
arr = ''.join(df)
arr = arr.split('>')[1:]
# print(arr)
b = []
for item in arr:
    a = item.split(' ')
    b.append(a)
dfData = pd.DataFrame(b)
# print(dfData[2])
dfData[2] = dfData[2] + ' ' + dfData[3]
dfData[3] = dfData[4]
dfData[4] = dfData[5]
dfData = dfData.loc[:, [0, 1, 2, 3, 4]]
# print(dfData)
hs = dfData[dfData[2].isin(['Homo sapiens'])]
# print(hs)
# print(dfData)  #0列是name，1列是Derives_from
o = []
q = []
for m in hs[4]:
    n = m[0:9]
    o.append(n)
    p = m[9:]
    q.append(p)
hs[4] = o
hs[5] = q
# print(hs[4])
t = []
for q in hs[5]:
    t.append(q.replace('U', 'T'))
hs[6] = t
# print(hs[1])
# hs = hs.drop(columns=[4])
hss = hs.copy()
hss = hss.drop(columns=[4])
# Short	Symbol	Accession	Family	Species	length	Sequence	CDNA
new_col = ['Symbol', 'Accession', 'species', 'Family', 'sequences', 'CDNA']
hss.columns = new_col
hss['length'] = hss['sequences'].str.len()
g = []
for i in hss['Symbol']:
    i = i[:3]
    g.append(i)
hss['short'] = g
hss = hss[['short','Symbol', 'Accession', 'Family','species', 'length', 'sequences', 'CDNA']]
hss.to_csv("hsahairpin.csv", index=False, header=True)

# 整理mature.fa的数据集，提取出Homo sapiens的基因信息
# data1 = pd.read_csv("mature.fa",header=None)
with open("mature.fa", "r") as f:
    data1 = f.read().splitlines()
# print(type(data1))
de = np.array(data1)
de = de.flatten()
arr1 = ' '.join(de)
arr1 = arr1.split('>')[1:]
c = []
for item1 in arr1:
    a1 = item1.split(' ')
    c.append(a1)
dfData1 = pd.DataFrame(c)
# print(dfData1.loc[:,[5]])
dfData1[2] = dfData1[2] + ' ' + dfData1[3]
dfData1[3] = dfData1[4]
dfData1[4] = dfData1[5]
# dfData1 = dfData1.loc[:,[0,1,2]]
hs1 = dfData1[dfData1[2].isin(['Homo sapiens'])]
hs2 = dfData1[dfData1[2].isin(['Homo sapiens'])]
# print(hs1)
# hs1.to_csv("cleanmature1.csv",index=False,header=True)
hs1 = hs1.drop(columns=[5, 6])
hs2 = hs2.drop(columns=[5, 6])
# hs2.dropna(inplace = True)
# print(hs1)
g = []
for i in hs2[0]:
    i = i[:3]
    g.append(i)
hs2[3] = g

new_col = ['Symbol', 'Accession', 'species', 'short', 'sequences']
hs2.columns = new_col
hs2['length'] = hs2['sequences'].str.len()
# 添加一列
f = []
for q in hs2['sequences']:
    f.append(q.replace('U', 'T'))
hs2['CDNA'] = f
hs2 = hs2[['short','Symbol', 'Accession', 'species','length', 'sequences','CDNA']]
hs2.to_csv("hsamature.csv", index=False, header=True)

# 整理hsagff3.fa数据集，利用ID联系hairpin.fa和mature.fa基因信息
data2 = pd.read_csv("hsa.gff3.txt", header=None)
dh = np.array(data2)
# print(data2)
dh = dh.flatten()
arr2 = ''.join(dh)
arr2 = arr2.split('chr')[1:]
d = []
for item2 in arr2:
    a2 = re.split('[\t;]', item2)
    d.append(a2)
dfData2 = pd.DataFrame(d)
# print(dfData2)
# print(dfData2[9])
e = []  # 截取出ID中的有效字段，删除无效列(‘.’)
for i in dfData2[8]:
    i = i[3:]
    e.append(i)
dfData2[8] = e
e = []  # 截取出Alias中的有效字段，删除无效列(‘.’)
for i in dfData2[9]:
    i = i[6:]
    e.append(i)
dfData2[9] = e
dfData2 = dfData2.drop(columns=[1, 5, 7], axis=1)
dfData2 = dfData2.T.reset_index(drop=True).T
dfData2.to_csv("cleanhsagff3.csv", index=False, header=True)

# 提取miRNA_primary_transcript
dfData2_1 = dfData2[dfData2[1].isin(['miRNA_primary_transcript'])]
dfData2_1.to_csv("hairpin_hsagff3_1.csv",index=False,header=True)

dfData3 = pd.DataFrame()
# dfData3 = dfData2_1.copy()
# print(hs[1])
# 将ID 为 hairpin的筛选出来,hs为hsahairpin.csv，若用hs[1]，则输出的id顺序是hairpin.fa文件里面的顺序
# 用dfData2_1[5],则输出的id顺序是按照hsa.gff3.txt
# for i in hs[1]:
# dfData3 = dfData3.append(dfData2.loc[dfData2[5]==i],ignore_index=True)
for i in dfData2_1[5]:
    dfData3 = dfData3.append(dfData2_1.loc[dfData2_1[5] == i], ignore_index=True)
# print(dfData3)
# cnt = dfData2.loc[dfData2[8].value(i)]
# hairpin_ID_hsagff3 = dfData2[0:cnt]
dfData3.to_csv("hairpin_ID_hsagff3_1.csv",index=False,header=True)
# maoshengdong >>
g = []
for i in dfData3[0]:
    a = "chr", i
    a = ''.join(a)
    g.append(a)
chr = pd.DataFrame(g)
dfData3[0] = chr
# print(dfData3[2])
dfData3 = dfData3.drop(columns=[1])

# dfData3.dropna(inplace = True)
# 输出列数
# print(dfData3.shape[1])
# 重新排列列
dfData3 = dfData3.T.reset_index(drop=True).T
dfData3 = dfData3.drop(columns=[5, 7])
dfData3 = dfData3.T.reset_index(drop=True).T
# 截取出ID中的有效字段，删除无效列(‘.’)
g = []
for i in dfData3[5]:
    i = i[5:]
    g.append(i)
dfData3[5] = g
new_col = ['Chromosome', 'Start', 'End', 'Strand', 'Accession', 'Name']
dfData3.columns = new_col
dfData3 = dfData3[['Accession', 'Name', 'Chromosome', 'Strand', 'Start', 'End']]
dfData3.to_csv("hairpin_hsagff3.csv", index=False, header=True)
# maoshengdong <<

dfData4 = pd.DataFrame()
# 提取miRNA
dfData2_2 = dfData2[dfData2[1].isin(['miRNA'])]
dfData2_2.to_csv("mature_hsagff3.csv", index=False, header=True)
# 若用hs[1]，则输出的id顺序是hairpin.fa文件里面的顺序
# 用dfData2_1[5],则输出的id顺序是按照hsa.gff3.txt
# for k in hs1[1]:
#     dfData4 = dfData4.append(dfData2.loc[dfData2[5]==k],ignore_index=True)
for k in dfData2_2[5]:
    dfData4 = dfData4.append(dfData2_2.loc[dfData2_2[5] == k], ignore_index=True)
# mature_ID_hsagff3 = dfData2[0:cnt]
# 删除Alias=这一列
# print(dfData4)
dfData4 = dfData4.drop(columns=[1, 5])
dfData4 = dfData4.T.reset_index(drop=True).T
dfData4.to_csv("mature_ID_hsagff3_1.csv", index=False, header=True)
# 截取出ID中的有效字段，删除无效列(‘.’)
g = []
for i in dfData4[5]:
    i = i[5:]
    g.append(i)
dfData4[5] = g
g = []
for i in dfData4[6]:
    i = i[13:]
    g.append(i)
dfData4[6] = g
# 添加chr
g = []
for i in dfData4[0]:
    a = "chr", i
    a = ''.join(a)
    g.append(a)
chr_p = pd.DataFrame(g)
dfData4[0] = chr_p
# 插入新的一列，列名为Arm
# print(dfData4)
g = []
for i in dfData4[5]:
    i = i[-2:]
    g.append(i)
dfData4[7] = g
# 添加列名，并更换列的顺序
new_col = ['Chromosome', 'Start', 'End', 'Strand', 'Accession', 'Name', 'Derives_from','Arm']
dfData4.columns = new_col
dfData4 = dfData4[['Derives_from', 'Accession', 'Name','Arm', 'Chromosome', 'Strand', 'Start', 'End']]
dfData4.to_csv("mature_ID_hsagff3.csv", index=False, header=True)

f = []
for q in hs1[4]:
    f.append(q.replace('U', 'T'))
hs1[5] = f
# print(hs1[5])
hs1.to_csv("cleanmature.csv", index=False, header=True)

# hs1指cleanmature.csv,hss指hsahairpin.csv,hs2指hsamature.csv
# 提取hsamature.csv中没有3p/5p的序列 >>
hs2_1 = hs2.copy()
sp = ['3p', '5p']
symbols = '|'.join(sp)
hs2_1 = hs2_1[~(hs2_1['Symbol'].str.contains(symbols))]
hs2_1.to_csv("hs2_1.csv", index=False, header=True)
# <<

# >> 从dfData4中提取包含有hs2_1的Accession的信息
dfData4_1 = pd.DataFrame()
for i in hs2_1['Accession']:
    dfData4_1 = dfData4_1.append(dfData4.loc[dfData4['Accession'] == i], ignore_index=True)
dfData4_1.to_csv("dfData4_1.csv", index=True, header=True)
# <<
# >> 从hsahairpin.csv中提取dfData4_1中Derives_from的信息
hss_1 = pd.DataFrame()
for i in dfData4_1['Derives_from']:
    hss_1 = hss_1.append(hss.loc[hss['Accession'] == i], ignore_index=True)
hss_1.to_csv("hss_1.csv", index=True, header=True)
# <<
# >> 从hsamature.csv中提取出包含dfData4_1.csv的Accession的信息
hs2_2 = pd.DataFrame()
for i in dfData4_1['Accession']:
    hs2_2 = hs2_2.append(hs2.loc[hs2['Accession'] == i], ignore_index=True)
# hs2_2.to_csv("hs2_2.csv",index=True,header=True)
# <<
# >> 定位hs2_2.csv中的sequences在hss_1.csv的sequences的始终坐标，
e = []
for i, j in zip(hs2_2['sequences'], hss_1['sequences']):
    e.append(re.search(i, j).start())
hs2_2['indexOfHairpin'] = e
# <<
# >> 根据坐标来判断mature的Arm,右为3p，左为5p
hs2_2['Arm'] = np.where(hs2_2['indexOfHairpin'] > hss_1['length']//2, '3p','5p')

hs2_2.to_csv("hs2_2.csv", index=True, header=True)
# <<
# >> 把hs2_2['Arm']放入mature_ID_hsagff3.csv

for i,j in zip(hs2_2['Accession'],hs2_2['Arm']):
    # print(j)
    dfData4.loc[i == dfData4['Accession'], 'Arm'] = j

dfData4.to_csv("mature_ID_hsagff3.csv", index=False, header=True)
# <<


hs1 = hs1.reset_index(drop=True)
hs1.dropna(inplace=True)
ss = 'U'
hs1['isis'] = hss['sequences'].str.find(ss, 0)

hs1.to_csv("cleanmature.csv", index=True, header=True)
# g = []
# for i in dfData4[8]:
#     i=i[13:]
#     g.append(i)
# dfData4[8] = g
# dfData4.to_csv("mature_ID_hsagff3.csv",index=False,header=True)

# hs1指cleanmature.csv，dfData4指mature_ID_hsagff3.csv
# print(hs1[1])
h = []
dfData5 = pd.DataFrame()
#     dfData5 = dfData4.append(dfData4.loc[dfData4[8] == i], ignore_index=True)
# print(dfData5)
# print(hs1)
# print(hs)
# for i in dfData4[11]:
#     h.append(i)
# print(h)
# for j in dfData4[11]:
#     for i in hs[1]:
#         if j == i:
#             print(j)
#         break
#
# hs指cleanhairpin.csv,hs1指cleanmature.csv,dfData4指mature_ID_hsagff3.csv
# print(hs[1])
for j in dfData4['Derives_from']:
    dfData5 = dfData5.append(hs.loc[hs[1] == j], ignore_index=True)
    # dfData5 = dfData5.drop_duplicates(subset=1)
dfData5.to_csv("matureIDhsagff3_to_cleanhairpin.csv", index=False, header=True)
# print(dfData4[11])


