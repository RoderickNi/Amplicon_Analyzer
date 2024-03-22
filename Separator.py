import os
import sys


'''
USAGE:

python Separator.py -fq <input: fq> -Pp <input: primers_info> -OD <output: Result>

arg_lst:

-fq
-Pp
-OD


'''

#序列反向互补函数
def RC_seq(seq):
    '''
    Input:
    string(DNAseq)
    
    =============================================================================================
    Return:
    string(DNAseq)
    '''
    compl = str.maketrans('ACGTacgt', 'TGCAtgca')                       # Create a replacement mapping table
    sequence_RC = seq.translate(compl)[::-1]
    return sequence_RC

#引物信息读取
def read_primers(primers):
    '''
    Input:
    The path of primers file in which the column should be splited by '\t'
    For example:
    
    PrimerPair-1    ACAGTGTGCCATAAAGTGCGTCTGGT    ACAGTGTCCCAATAAGATTGCAACAAGAACT    50    PopName1
    ...
    
    =============================================================================================
    Return:
    The list of primers info
    e.g.:
    
    [
    ('PopName1','ACAGTGTGCCATAAAGTGCGTCTGGT','AGTTCTTGTTGCAATCTTATTGGGACACTGT',50,5,5),
    ('PopName2','TCAGTGTGCCATAAAGTGCGTCTGGT','AGTTCTTGTTGCAATCTTATTGGGACACTGA',50,5,5)
    ...
    ]
    
    '''
    pri = open(primers,'r').read().split('\n')
    while '' in pri:
        pri.remove('')
    lst=[]
    for line in pri:
        pop=str(line.split('\t')[0])
        F1=str(line.split('\t')[1])
        R1=RC_seq(str(line.split('\t')[2]))
        sample_size=int(line.split('\t')[3])
        F_Trim=int(line.split('\t')[4])
        R_Trim=int(line.split('\t')[5])
        lst.append((pop,F1,R1,sample_size,F_Trim,R_Trim))
    return lst



if __name__ in '__main__':

    
    # 序列分选
    for i in range(1,len(sys.argv)):
        if r'-OD' in sys.argv[i]:
            OutputDir=sys.argv[i+1]
            os.mkdir(OutputDir)  
        elif r'-fq' in sys.argv[i]:
            fq=sys.argv[i+1]

    
    for i in range(1,len(sys.argv)):
        if r'-Pp' in sys.argv[i]:
            print('Loading primers ... ')
            primers=sys.argv[i+1]
            primer_lst=read_primers(primers)
            print('Separating Reads ...')
            for primer_pair in primer_lst:
                print(''.join([str(x) for x in primer_pair[0]]))
                Pop=primer_pair[0]
                F1=primer_pair[1]
                R1=primer_pair[2]
                R2=RC_seq(F1)
                F2=RC_seq(R1)
                write_in_1=os.path.join(OutputDir,Pop+'.fq')
                write_in_2=os.path.join(OutputDir,Pop+'_RC.fq')
                os.system(f'''
                touch {write_in_1}
                touch {write_in_2}
                ''')
                os.system(f'''
                echo Getting reads: head-{F1} tail-{R1} from {fq} 
                cat {fq} | grep -i -B 1 -A 2 ^{F1}.*.{R1}$ | sed 's/^--$//g' | sed '/^$/d' >> {write_in_1} 
                echo Getting reads: head-{F2} tail-{R2} from {fq} 
                cat {fq} | grep -i -B 1 -A 2 ^{F2}.*.{R2}$ | sed 's/^--$//g' | sed '/^$/d' >> {write_in_2} 
                seqtk seq -r {write_in_2} >> {write_in_1}
                rm {write_in_2}
                ''')
    # 标签序列剔除
    for i in range(1,len(sys.argv)):
        if r'-Pp' in sys.argv[i]:
            print('Loading Trim info ... ')
            primers=sys.argv[i+1]
            primer_lst=read_primers(primers)
            print('Trimming Reads ...')
            for Trim_pair in primer_lst:
                F_Trim = Trim_pair[4]
                R_Trim = Trim_pair[5]
                for file in os.listdir(OutputDir):
                    if file.endswith('.fq') and file.replace(r'.fq','') == Trim_pair[0]:
                        print(f'''Trimming {file} with 5'-{str(F_Trim)}bp ... 3'-{str(R_Trim)}bp''')
                        fastq=open(os.path.join(OutputDir,file),'r').read().split('\n')
                        while '' in fastq:
                            fastq.remove('')
                        OTPT=open(os.path.join(OutputDir,file).replace('.fq','.fastq'),'a',encoding='utf-8')
                        
                        for i in range(0,len(fastq),4):
                            name=fastq[i]
                            seq=fastq[i+1][F_Trim:]
                            seq=seq[:-R_Trim]
                            three=fastq[i+2]
                            quality=fastq[i+3][F_Trim:]
                            quality=quality[:-R_Trim]
                            print(name+'\n'+seq+'\n'+three+'\n'+quality,file = OTPT)
                        OTPT.close()
                        os.system(f'''
                        rm {os.path.join(OutputDir,file)}
                        ''')
    

# python Separator.py -fq ./B_Overlapping/out.extendedFrags.fastq -Pp ./primers.tab -OD C_PopReads