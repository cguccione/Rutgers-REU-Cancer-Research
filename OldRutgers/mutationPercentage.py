file=open('TRF068722.xls', 'r')
#file=open('TRF046732.xls', 'r')
lines=file.readlines()
n, comp_purity = (lines[8]).split(":")
n, path_purity = (lines[7]).split(":")

'''
**Splits up the code into lists so that it can be filtered
easily in the future. **
    -   Right now the data that has
    applified, loss or fusion is not included in the lists
    because they don't follow the normal format.
    -   Count starts at 1 instead of 0 so it matches up with
    the xls file numbering
'''
count =1
gene_list=[]
allele_list=[]
LOHGICpath_list=[]
LOHGICcomp_list=[]
test=0
for line in lines:
    #This section just contains the mutations
    if test == 0:
        none,TRF=line.split(':')
        TRF=TRF.strip('\n')
        test=1
    if count >=17:
        test=0
        for word in line.split('\t'):
            test+=1
        if test >= 12:
            chr, start, end, gene, aminoAcid, mutation, allele, depth, strand, CN, LOHGICpath, LOHGICcomp = line.split('\t')
            gene_list.append(gene)
            allele_list.append(allele)
            LOHGICpath_list.append(LOHGICpath)
            LOHGICcomp_list.append(LOHGICcomp)
    count=count+1

'''
This section tests to see if there is any TP53 alleles
and if so it takes the average of thier allele freq.
'''
TP53_allele_list=[]
TPcount=0
for x in gene_list:
    if x == 'TP53':
        TP53_allele_list.append(allele_list[TPcount])
    TPcount+=1

TP53_avg=0
total_TP53=0
for y in TP53_allele_list:
    TP53_avg+=float(y)
    total_TP53+=1
if total_TP53 != 0:
    TP53_avg=TP53_avg/total_TP53

'''
Option1= Cmut=1, LOH not present
Option2= Cmut=2, LOH present
Option3= Cmut=1, LOH present
'''

#Comp Purity
Option1_compPurity=[]
Option2_compPurity=[]
Option3_compPurity=[]

#Path Purity
Option1_pathPurity=[]
Option2_pathPurity=[]
Option3_pathPurity=[]

#Number of cells tested with, currently set at 100 but can be modifed
cells=10

def CCF_calc (VAF, comp_purity, path_purity):

    comp_purity=float(comp_purity)*0.01
    path_purity=float(path_purity)*0.01
    VAF=float(VAF)*0.01

    #Number of cancer cells per 10 cells
    total_cancer_cells_comp=comp_purity*cells
    total_cancer_cells_path=path_purity*cells

    #Alleles with mutation for Options 1 & 2
    alles_w_mut_o12=(cells*2)*VAF

    #Option1
    CCF_comp = alles_w_mut_o12/total_cancer_cells_comp
    CCF_path = alles_w_mut_o12/total_cancer_cells_path

    #Assuming that mutated allels are only on
    #'cancer' alleles
    Option1_compPurity.append(round((CCF_comp*100),2))
    Option1_pathPurity.append(round((CCF_path*100),2))

    #Option2
    CCF_comp = alles_w_mut_o12/(2*total_cancer_cells_comp)
    CCF_path = alles_w_mut_o12/(2*total_cancer_cells_path)

    #Assuming if a cell has cancer then both alllels are mutated
    Option2_compPurity.append(round((CCF_comp*100),2))
    Option2_pathPurity.append(round((CCF_path*100),2))

    #Option3
    alles_w_mut_o3_comp=((cells*2)-total_cancer_cells_comp)*VAF
    alles_w_mut_o3_path=((cells*2)-total_cancer_cells_path)*VAF

    CCF_comp = alles_w_mut_o3_comp/total_cancer_cells_comp
    CCF_path = alles_w_mut_o3_path/total_cancer_cells_path

    Option3_compPurity.append(round((CCF_comp*100),2))
    Option3_pathPurity.append(round((CCF_path*100),2))

    #These are also all 'worse case' in terms of number
    #of cancer cells. They all assume that even if two
    #alleles are present, unless they are both cancerous
    #alleles, the cancer allele will move into a new cell
    return

for VAF in allele_list:
    CCF_calc(VAF, comp_purity, path_purity)

def CCF_calc_TP53 (VAF, TP53_purity):

    TP53_purity=float(TP53_purity)*0.01
    VAF=float(VAF)*0.01

    #Number of cancer cells per 10 cells
    total_cancer_cells_TP53=TP53_purity*cells

    #Alleles with mutation for Options 1 & 2
    alles_w_mut_o12=(cells*2)*VAF

    #Option1
    CCF_TP53 = alles_w_mut_o12/total_cancer_cells_TP53
    #Assuming that mutated allels are only on
    #'cancer' alleles
    Option1_TP53.append(round((CCF_TP53*100),2))

    #Option2
    CCF_TP53 = alles_w_mut_o12/(2*total_cancer_cells_TP53)
    Option2_TP53.append(round((CCF_TP53*100),2))

    #Option3
    alles_w_mut_o3_TP53=((cells*2)-total_cancer_cells_TP53)*VAF
    CCF_TP53 = alles_w_mut_o3_TP53/total_cancer_cells_TP53
    Option3_TP53.append(round((CCF_TP53*100),2))

    #These are also all 'worse case' in terms of number
    #of cancer cells. They all assume that even if two
    #alleles are present, unless they are both cancerous
    #alleles, the cancer allele will move into a new cell
    return

#TP53 Purity
if TP53_avg !=0:
    Option1_TP53=[]
    Option2_TP53=[]
    Option3_TP53=[]
    for VAF in allele_list:
        CCF_calc_TP53(VAF, TP53_avg)
    #Print TP53 Purity
    print(Option1_TP53)
    print(Option2_TP53)
    print(Option3_TP53)

temp='CCF_'+ TRF+'.xls'
new_file=open(temp, 'w+')

new_file.write("Protein \t Option1: Comp \t Option2: Comp \t Option3: Comp \t Option1: Path \t Option2: Path \t Option3: Path \t Option1: TP53 \t Option2: TP53 \t Option3: TP53")

l=len(Option1_compPurity)
x=0

while x < l:
    new_file.write("\n")
    a=str(Option1_compPurity[x])
    b=str(Option2_compPurity[x])
    c=str(Option3_compPurity[x])
    d=str(Option1_pathPurity[x])
    e=str(Option2_pathPurity[x])
    f=str(Option3_pathPurity[x])
    new_file.write(gene_list[x] + "\t" + a + "\t" + b + "\t" + c + "\t")
    new_file.write(d+ "\t" + e+ "\t" + f + "\t")
    if TP53_avg != 0:
        g=str(Option1_TP53[x])
        h=str(Option2_TP53[x])
        i=str(Option3_TP53[x])
        new_file.write(g+ "\t" + h+ "\t" + i)
    x=x+1

#Print Comp Purity
print(Option1_compPurity)
print(Option2_compPurity)
print(Option3_compPurity)

#Print Path Purity
print(Option1_pathPurity)
print(Option2_pathPurity)
print(Option3_pathPurity)

#Test cased, same one used in the notebook example
#CCF_calc(39.54, 28, 1)
#CCF_calc(23, 41, 1)
