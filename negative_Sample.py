from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import random

# 读取FASTA文件
records = list(SeqIO.parse("Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa", "fasta"))

# 将序列字母全部转换为大写
for record in records:
    record.seq = record.seq.upper()

# 将修改后的记录写回FASTA文件
with open("mirBASE_uppercase.fasta", "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")

# 打开CSV文件并读取数据

with open('hsahairpin1.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    # 生成Fasta格式的数据
    records = []
    for row in reader:
        # row = row.seq
        # row = row.tomutable()
        # 生成SeqRecord对象
        # Seq(row[5])
        record = SeqIO.SeqRecord(
            seq=Seq(row[5]),
            id=row[0],
            description=' '.join(row[1:5])
        )
        records.append(record)

# 将SeqRecord对象列表写入Fasta文件
with open('hsahairpin.fasta', 'w') as outfile:
    SeqIO.write(records, outfile, 'fasta')

# 读取hsahairpin.fasta文件
from Bio.SeqRecord import SeqRecord

b_records = SeqIO.to_dict(SeqIO.parse("hsahairpin.fasta", "fasta"))

target_records = []
# 遍历mirBASE_uppercase.fasta文件
for a_record in SeqIO.parse("mirBASE_uppercase.fasta", "fasta"):
    # 在a_record中查找b_record的序列
    for b_record in b_records.values():
        target_seq = b_record.seq
        if target_seq in a_record.seq:
            # 定位到目标序列
            start = a_record.seq.find(target_seq)
            end = start + len(target_seq)
            # 截取目标序列前后10个字母
            flank_len = 100
            flank_start = max(0, start-flank_len)
            flank_end = min(len(a_record.seq), end+flank_len)
            target_seq_with_flank = a_record.seq[flank_start:flank_end]

            # 创建新的SeqRecord对象
            target_record = a_record[start:end]
            target_record.id = b_record.id
            target_record.description = ""
            target_record.seq = Seq(target_seq_with_flank)
            target_records.append(target_record)
            # 输出目标序列信息
            print("ID:", target_record.id)
            print("Target sequence with flank:", target_seq_with_flank)
            # new_record = SeqRecord(target_records, id=f"{b_record.id}",
            #                        description=f"Negative sample ")
            SeqIO.write(target_records, "from_hsa_mirBase.fasta", "fasta")

# 对hsahairpin.fasta的文件序列生成负样本（每条序列生成10个负样本），使用random.shuffle的方法
# 打开输入和输出FASTA文件
input_file = "hsahairpin.fasta"
output_file = "sample_hsahairpin.fasta"
with open(output_file, "w") as outfile:
    # 读取输入FASTA文件中的每个序列
    for record in SeqIO.parse(input_file, "fasta"):
        # 生成10个负样本序列
        for i in range(10):
            # 随机重排碱基
            new_seq = record.seq.tomutable()
            random.shuffle(new_seq)
            # 创建新的SeqRecord对象并写入输出FASTA文件
            new_record = SeqRecord(new_seq.toseq(), id=f"{record.id}_neg{i+1}", description=f"Negative sample random.shuffle {i+1}")
            SeqIO.write(new_record, outfile, "fasta")


# 对from_hsa_mirBase.fasta的文件序列生成负样本（每条序列生成10个负样本），使用random.shuffle的方法
# 打开输入和输出FASTA文件
input_file = "from_hsa_mirBase.fasta"
output_file = "sample_from_hsa_mirBase_random_shuffle.fasta"
with open(output_file, "w") as outfile:
    # 读取输入FASTA文件中的每个序列
    for record in SeqIO.parse(input_file, "fasta"):
        # 生成10个负样本序列
        for i in range(10):
            # 随机重排碱基
            new_seq = record.seq.tomutable()
            random.shuffle(new_seq)
            # 创建新的SeqRecord对象并写入输出FASTA文件
            new_record = SeqRecord(new_seq.toseq(), id=f"{record.id}_neg{i+1}", description=f"Negative sample random.shuffle {i+1}")
            SeqIO.write(new_record, outfile, "fasta")

# 对hsahairpin.fasta的文件序列生成负样本，使用reverse_complement的方法
# 打开输入和输出FASTA文件
input_file = "hsahairpin.fasta"
output_file = "sample_hsahairpin_reverse_complement.fasta"
with open(output_file, "w") as outfile:
    # 读取输入FASTA文件中的每个序列
    for i, record in enumerate(SeqIO.parse(input_file, "fasta"), 1):
        # 反向互补序列
        new_seq = record.seq.tomutable()
        new_seq.reverse_complement()
        # 创建新的SeqRecord对象并写入输出FASTA文件
        new_record = SeqRecord(Seq(new_seq), id=f"{record.id}", description=f"Negative sample by reverse_complement {i}")
        SeqIO.write(new_record, outfile, "fasta")

# 对from_hsa_mirBase.fasta的文件序列生成负样本，使用reverse_complement的方法
# 打开输入和输出FASTA文件
input_file = "from_hsa_mirBase.fasta"
output_file = "sample_from_hsa_mirBase_reverse_complement.fasta"
with open(output_file, "w") as outfile:
    # 读取输入FASTA文件中的每个序列
    for i, record in enumerate(SeqIO.parse(input_file, "fasta"), 1):
        # 反向互补序列
        new_seq = record.seq.tomutable()
        new_seq.reverse_complement()
        # 创建新的SeqRecord对象并写入输出FASTA文件
        new_record = SeqRecord(Seq(new_seq), id=f"{record.id}", description=f"Negative sample by reverse_complement {i}")
        SeqIO.write(new_record, outfile, "fasta")