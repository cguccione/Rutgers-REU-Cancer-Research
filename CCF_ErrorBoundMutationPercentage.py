#Caitlin Guccione
#Summer 2018
#Rutgers Cancer Institute of New Jersey
#DIMACS REU

import math
#This program finds the parameter estimates and confidence intervals for binomial data
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

#Total Number of cells: doesn't matter is divided out later
cells=10

def choose_option (LOHGIC):
    #This function determines which option is correct for the
    #given LOGHIC output and returns option #, path/comp

    #a is the information we want form LOHGIC and
    #b is the part we don't need right now
    a,b=LOHGIC.split(";")
    return_test=False
    #This assumes that if logic outputs multiple options, in order
    #for the mutation to chose an option it must be great than 70% confident
    if ',' in a:
        #c is the confidence level
        a,b=a.split(',')
        b,c=a.split('(')
        c,b=c.split(')')
        if (float(c) < 0.7):
            return_Test=True
            return 0
    #sg is somatic or germline
    #b,c,d is information not needed
    #cmut is the Cmut number
    sg,a=a.split('=')
    cmut,b=a.split('(')
    #Option 1:Somatic CNmut=1
    #Option 2:Somatic CNmut=2 *So LOH occurred here but LOHGIC doesn't say that
    #Option 3: Somatic LOH CNmut=1
    #Option4= Germline Cmut=1
    #Option5= Germline Cmut=2, *So LOH occurred here but LOHGIC doesn't say that
    #Option6= Germline LOH Cmut=1
    #Option0= Inconclusive
    if return_test == False:
        if sg == "Germline LOH CNmut ":
            return 6
        elif sg == "Somatic LOH CNmut ":
            return 3
        elif sg == "Germline CNmut ":
            if int(cmut) == 1:
                    return 4
            else:
                    return 5
        elif sg == "Somatic CNmut ":
            if int(cmut) == 1:
                    return 1
            else:
                    return 2
    else:
        return 0

def CCF_calc (VAF, depth, purity, option):
    '''
    These are also all 'worse case' in terms of number
    of cancer cells. They all assume that even if two
    alleles are present, unless they are both cancerous
    alleles, the cancer allele will move into a new cell
    Option1= Somatic Cmut=1, LOH not present
    Option2= Somatic Cmut=2, LOH present
    Option3= Somatic Cmut=1, LOH present
    Option4= Germline Cmut=1, LOH not present
    Option5= Germline Cmut=2, LOH present
    Option6= Germline Cmut=1, LOH present
    Option0= Inconclusive
    '''

    purity=float(purity)*0.01
    #P=0.05
    P=purity*0.1
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

    option=int(option)
    if option == 1 or option == 2 or option == 4 or option == 5:
        Y=2
    else:
        Y=1
    if option ==1 or option == 3 or option == 4 or option == 6:
        cnmut=1
    else:
        cnmut=2

    #Simplifying variable names for following equation
    p=purity
    f=VAF
    F=V

    #Notes in notebook and LaTex on logic behind this
    if option < 4:
        #Somatic Mutations
        a=1-p
        A=(P/p)*a
        b=2*a
        B=2*A
        c=Y*p
        C=P*Y
        d=b+c
        D=math.sqrt((B)**2 + (C)**2)
        e=f*d
        E=math.sqrt((F/f)**2 + (D/d)**2) * e
        h=cnmut * p
        H = P * cnmut
        g = e/h
        G= math.sqrt((E/e)**2 + (H/h)** 2) * g

        purity_prec=str(round((g*100),2)) + '±' + str(round((G*100),2))

    else:
        print("f", f)
        print("p", p)
        print('cnmut', cnmut)
        print("Y", Y)
        #Germline mutations
        a=1-p
        A=(P/p)*a
        b=2*a
        B=2*A
        c=Y*p
        C=Y*P
        d=b+c
        D=math.sqrt((B)**2 + (C)**2)
        e=f*d
        E=math.sqrt((F/f)**2 + (D/d)**2) * e
        g=1-p
        G=(P/p)*g
        h=e-g
        H=math.sqrt((E)**2 + (G)**2)
        i=cnmut *p
        I=cnmut * P
        j=h/i
        J=math.sqrt((H/h)**2 + (I/i)**2) * j

        purity_prec=str(round((j*100),2)) + '±' + str(round((J*100),2))

    print(option)
    if option == 0:
        purity_prec = 0

    return purity_prec, option

def compare_purity (VAF, depth, lpath, lcomp):
    '''
    #Path purity
    option=choose_option(lpath)
    path_purity_prec, path_option = CCF_calc(VAF, depth, path_purity, option)
    return path_purity_prec, path_option, path_purity, 'Path'
    '''

    #Comp purity
    option=choose_option(lcomp)
    comp_purity_prec, comp_option = CCF_calc(VAF, depth, comp_purity, option)
    return comp_purity_prec, comp_option, comp_purity, 'Comp'


temp='CCF_'+ TRF+'.xls'
new_file=open(temp, 'w+')

new_file.write("Protein \t VAF \t Option # \t Purity \t Mutation Percentage \t LOHGIC Path Purity \t LOHGIC Comp Purity")

count =0;
while count < len(allele_list):
    purity_perc, option, purity, purity_type=compare_purity(allele_list[count], depth_list[count], LOHGICpath_list[count], LOHGICcomp_list[count])
    a=(LOHGICcomp_list[count]).strip('\n')
    purity_type=purity_type + ' - ' + purity
    purity_type=purity_type.strip('\n')
    if purity_perc == 0:
        new_file.write("\n")
        new_file.write(gene_list[count] + "\t" + str(allele_list[count]) + "\t NA \t NA \t NA \t" + str(LOHGICpath_list[count]) + "\t" + a)
    else:
        new_file.write("\n")
        new_file.write(gene_list[count] + "\t" + str(allele_list[count]) + "\t" + str(option) + "\t" + purity_type + "\t"+ str(purity_perc) + "\t")
    count=count+1
