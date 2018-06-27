import math
from statsmodels.stats.proportion import proportion_confint as binofit

#file=open('TRF068722.xls', 'r')
#file=open('TRF046732.xls', 'r')
file=open('TRF351303.xls', 'r')
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
depth_list=[]
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
            depth_list.append(depth)
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

#Total Number of cells
cells=10

def CCF_calc (VAF, depth, comp_purity, path_purity):

    comp_purity=float(comp_purity)*0.01
    P_comp=0.05
    path_purity=float(path_purity)*0.01
    P_path=0.05
    VAF=float(VAF)*0.01
    depth=float(depth)

    #Calculates the correct bounds on VAF using the bionomial parameter
    X = binofit((depth*VAF), depth, 0.01, method='binom_test')
    X=str(X)
    X=X.strip('()')
    Xup, Xlow=X.split(',')
    Xup=float(Xup)
    Xlow=float(Xlow)
    V=((Xup-Xlow)/2)


    #Number of cancer cells per 10 cells
    total_cancer_cells_comp=comp_purity*cells
    T_comp=(P_comp/comp_purity)* total_cancer_cells_comp
    total_cancer_cells_path=path_purity*cells
    T_path=(P_path/path_purity)* total_cancer_cells_path

    #Alleles with mutation for Options 1 & 2
    alles_w_mut_o12=float((cells*2)*VAF)
    A_o12=(V/VAF)* alles_w_mut_o12

    #Option1
    CCF_comp = alles_w_mut_o12/total_cancer_cells_comp
    C_comp=(math.sqrt(math.pow((T_comp/total_cancer_cells_comp),2) + math.pow((A_o12/alles_w_mut_o12),2)))*CCF_comp
    CCF_path = alles_w_mut_o12/total_cancer_cells_path
    C_path=(math.sqrt(math.pow((T_path/total_cancer_cells_path),2) + math.pow((A_o12/alles_w_mut_o12),2)))*CCF_path

    #Assuming that mutated allels are only on
    #'cancer' alleles
    Option1_compPurity.append(str(round((CCF_comp*100),2)) + '±' + str(round((C_comp*100),2)))
    Option1_pathPurity.append(str(round((CCF_path*100),2)) + '±' + str(round((C_path*100),2)))

    #Option2
    CCF_comp = alles_w_mut_o12/(2*total_cancer_cells_comp)
    y_comp=2*total_cancer_cells_comp
    Y_comp=2*T_comp
    C_comp=(math.sqrt(math.pow((Y_comp/y_comp),2) + math.pow((A_o12/alles_w_mut_o12),2)))*CCF_comp
    CCF_path = alles_w_mut_o12/(2*total_cancer_cells_path)
    y_path=2*total_cancer_cells_path
    Y_path=2*T_path
    C_path=(math.sqrt(math.pow((Y_path/y_path),2) + math.pow((A_o12/alles_w_mut_o12),2)))*CCF_path

    #Assuming if a cell has cancer then both allels are mutated
    Option2_compPurity.append(str(round((CCF_comp*100),2)) + '±' + str(round((C_comp*100),2)))
    Option2_pathPurity.append(str(round((CCF_path*100),2)) + '±' + str(round((C_path*100),2)))

    #Option3

    #Calculate the total number of alleles with mutation
    x_comp=(cells*2)-total_cancer_cells_comp
    X_comp=T_comp
    alles_w_mut_o3_comp=x_comp*VAF
    A_o3_comp=(math.sqrt(math.pow((X_comp/x_comp),2) + math.pow((V/VAF),2)))*alles_w_mut_o3_comp

    x_path=(cells*2)-total_cancer_cells_path
    X_path=T_path
    alles_w_mut_o3_path=x_path*VAF
    A_o3_path=(math.sqrt(math.pow((X_path/x_path),2) + math.pow((V/VAF),2)))*alles_w_mut_o3_path

    #Calculate the CCF
    CCF_comp = alles_w_mut_o3_comp/total_cancer_cells_comp
    C_comp=(math.sqrt(math.pow((T_comp/total_cancer_cells_comp),2) + math.pow((A_o3_comp/alles_w_mut_o3_comp),2)))*CCF_comp
    CCF_path = alles_w_mut_o3_path/total_cancer_cells_path
    C_path=(math.sqrt(math.pow((T_path/total_cancer_cells_path),2) + math.pow((A_o3_path/alles_w_mut_o3_path),2)))*CCF_path

    Option3_compPurity.append(str(round((CCF_comp*100),2)) + '±' + str(round((C_comp*100),2)))
    Option3_pathPurity.append(str(round((CCF_path*100),2)) + '±' + str(round((C_path*100),2)))

    #These are also all 'worse case' in terms of number
    #of cancer cells. They all assume that even if two
    #alleles are present, unless they are both cancerous
    #alleles, the cancer allele will move into a new cell
    return

count =0;
while count < len(allele_list):
    CCF_calc(allele_list[count], depth_list[count], comp_purity, path_purity)
    count=count+1

def CCF_calc_TP53 (VAF, depth, TP53_purity):

    TP53_purity=float(TP53_purity)*0.01
    P=0.05
    VAF=float(VAF)*0.01
    depth=float(depth)

    #Calculates the correct bounds on VAF using the bionomial parameter
    X = binofit((depth*VAF), depth, 0.01, method='binom_test')
    X=str(X)
    X=X.strip('()')
    Xup, Xlow=X.split(',')
    Xup=float(Xup)
    Xlow=float(Xlow)
    V=((Xup-Xlow)/2)

    #Number of cancer cells per 10 cells
    total_cancer_cells_TP53=TP53_purity*cells
    T=(P/TP53_purity) * total_cancer_cells_TP53

    #Alleles with mutation for Options 1 & 2
    alles_w_mut_o12=float((cells*2)*VAF)
    A_o12=(V/VAF)*alles_w_mut_o12

    #Option1
    CCF_TP53 = alles_w_mut_o12/total_cancer_cells_TP53
    C=(math.sqrt(math.pow((T/total_cancer_cells_TP53),2) + math.pow((A_o12/alles_w_mut_o12),2)))*CCF_TP53
    #Assuming that mutated allels are only on
    #'cancer' alleles
    Option1_TP53.append(str(round((CCF_TP53*100),2)) + '±' + str(round((C*100),2)))

    #Option2
    CCF_TP53 = alles_w_mut_o12/(2*total_cancer_cells_TP53)
    y=2*total_cancer_cells_TP53
    Y=2*T
    C=(math.sqrt(math.pow((Y/y),2) + math.pow((A_o12/alles_w_mut_o12),2)))*CCF_TP53
    Option2_TP53.append(str(round((CCF_TP53*100),2)) + '±' + str(round((C*100),2)))

    #Option3
    #Calculate the total number of alleles with mutation
    x=(cells*2)-total_cancer_cells_TP53
    X=T
    alles_w_mut_o3=x*VAF
    A_o3=(math.sqrt(math.pow((X/x),2) + math.pow((V/VAF),2)))*alles_w_mut_o3

    CCF_TP53 = alles_w_mut_o3/total_cancer_cells_TP53
    C=(math.sqrt(math.pow((T/total_cancer_cells_TP53),2) + math.pow((A_o3/alles_w_mut_o3),2)))*CCF_TP53
    Option3_TP53.append(str(round((CCF_TP53*100),2)) + '±' + str(round((C*100),2)))

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

    count =0;
    while count < len(allele_list):
        CCF_calc_TP53(allele_list[count], depth_list[count], TP53_avg)
        count=count+1

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
