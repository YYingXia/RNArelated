import csv
import os
import pandas as pd

def dot_bracket_to_bpseq(sequence, structure):
    stack = []
    bpseq = []

    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i + 1)  # 1-based index
        elif char == ')':
            j = stack.pop()
            bpseq.append((j, i + 1))

    # 将结果整理为BPSEQ格式
    bpseq_output = []
    pairs = {i: j for i, j in bpseq}
    pairs.update({j: i for i, j in bpseq})

    for i, nucleotide in enumerate(sequence):
        paired_with = pairs.get(i + 1, 0)  # 1-based index
        bpseq_output.append((i + 1, nucleotide, paired_with))

    return bpseq_output

def save_bpseq_to_csv(bpseq, csv_file):
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        writer.writerows(bpseq)  # 写入BPSEQ数据

def bpseq2fasta(file):
    df = pd.read_csv(file, names=['id', 'base', 'pair'], sep=' ')
    seq = ''.join(df['base'].tolist())
    return seq
    

def rfamid2seqid():
    files = os.listdir('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/fastaFiles')
    for file in files:
        if file.endswith('fasta'):
            with open('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/fastaFiles/' + file, 'r') as f:
                text = f.readlines()
            seqid = ''.join(text[0].strip().split('_')[1:])
            sequence = text[1].strip()

            with open('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/fasta2/{}.fasta'.format(seqid), 'w') as f:
                f.writelines(text)

            filenum = file.strip('.fasta').split('_')[-1]
            filebpseq = file.replace(filenum, str(int(filenum)//2 + 1)).replace('fasta', 'bpseq')
            sequence2 = bpseq2fasta('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/bpseqFiles/' + filebpseq)

            if sequence != sequence2:
                print(file, filebpseq)
                print(sequence)
                print(sequence2)
            else:
                with open('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/bpseqFiles/' + filebpseq, 'r') as f:
                    text = f.readlines()
                with open('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/bpseq2/{}.bpseq'.format(seqid), 'w') as f:
                    f.writelines(text)

# rfam seed 
def rfamseed_process():
    bpRNA1m_list = []
    for file in os.listdir('/data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/bpseq2'):
        bpRNA1m_list.append(file.strip('.bpseq'))

    bpRNAnew_list = []
    for file in os.listdir('/data3/xiaying/model/RNAsearch/dataset/bpRNAnew'):
        bpRNAnew_list.append(file.strip('.bpseq'))  


    with open('/data3/xiaying/model/RNAsearch/dataset/rfam/Rfam.seed','r', encoding='utf-8', errors='ignore') as f:
        text = "".join(f.readlines())
    seq_texts = text.split('//\n# STOCKHOLM 1.0')
    for seq_text in seq_texts:
        family_name = seq_text.split('#=GF AC')[1].split('#=GF ID')[0].strip('\n').strip(' ')
        # ss = seq_text.split('#=GC SS_cons')[1].split('#=GC RF')[0].strip('\n').strip(' ')
        seed_seqs = seq_text.split('#=GF SQ')[1].split('#=GC SS_cons')[0].strip('\n').strip(' ').split('\n')[1:]
        saved_pth = '/data3/xiaying/model/RNAsearch/dataset/cmset/rfam/bpseq/' + family_name
        for item in seed_seqs:
            if len(item)>1:
                seqid = item.split(' ')[0]
                if seqid in bpRNA1m_list:
                    if not os.path.exists(saved_pth):
                        os.mkdir(saved_pth)
                    os.system('cp /data3/xiaying/model/RNAsearch/dataset/bpRNA-1m/rfam/bpseq2/{}.bpseq {}'.format(seqid, saved_pth))
                elif seqid in bpRNAnew_list:
                    if not os.path.exists(saved_pth):
                        os.mkdir(saved_pth)
                    os.system('cp /data3/xiaying/model/RNAsearch/dataset/bpRNAnew/{}.bpseq {}'.format(seqid, saved_pth))
        

def r2dtrfam():
    files = os.listdir('/data3/xiaying/model/R2DT/data/rfam')
    for file in files:
        if file.startswith('RF'):
            with open('/data3/xiaying/model/R2DT/data/rfam/{}/{}-traveler.fasta'.format(file, file), 'r') as f:
                text = f.readlines()
            bpseq = dot_bracket_to_bpseq(text[1].strip(), text[2].strip())
            save_bpseq_to_csv(bpseq, '/data3/xiaying/model/RNAsearch/dataset/cmset/rfam/bpseq/{}.bpseq'.format(file))
            a=0

def rfamCM():
    with open('/data3/xiaying/model/RNAsearch/dataset/rfam/Rfam.cm', 'r') as f:
        text = ''.join(f.readlines())
    cms = text.split('INFERNAL1/a')[1:]
    cm_dict = {}
    for cm in cms:
        family_id = cm.split('ACC')[1].split('DESC')[0].strip()
        cm_dict[family_id] = 'INFERNAL1/a' + cm

    for file in os.listdir('/data3/xiaying/model/RNAsearch/dataset/cmset/rfam/bpseq'):
        if file.endswith('bpseq'):
            family_id = file.split('.')[0]
            if family_id in cm_dict.keys():
                with open('/data3/xiaying/model/RNAsearch/dataset/cmset/rfam/cms/{}.cm'.format(family_id), 'w') as f:
                    f.writelines(cm_dict[family_id])

def getfasta():
    directory = '/data3/xiaying/model/RNAsearch/dataset/cmset/rfam'
    for file in os.listdir(directory + '/bpseq'):
        if not file.endswith('.bpseq'):
            continue
        seqid = file.strip('.bpseq')
        with open(directory + '/bpseq/' + file, 'r') as f:
            text = f.readlines()
        for i in range(len(text)):
            if text[i].startswith('1 '):
                break
        text2 = text[i:]
        # with open(directory + '/bpseq2/' + file, 'w') as f:
            # f.writelines(text2)

        sequence = ''  
        for line in text2:
            base = line.split(' ')[1]
            sequence += base
        
        with open(directory + '/fasta/{}.fasta'.format(seqid), 'w') as f:
            f.write('>{}\n{}\n'.format(seqid, sequence))
                

# getfasta()
