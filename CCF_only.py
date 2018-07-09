#Caitlin Guccione
#Summer 2018
#Rutgers Cancer Institute of New Jersey
#DIMACS REU

import math
#This program finds the parameter estimates and confidence intervals for binomial data
from statsmodels.stats.proportion import proportion_confint as binofit

'''
****************************************NO ERROR CALC
*Same as CCF_error_only but calculations don't output error

Input:
VAF: Allele Frequency for particular mutation
depth: Depth for particular mutation
Y: Total number of alleles per cell
cmut: Total number of cancerous alleles per cell
purity: entire sample purity
sg: Somatic or Germline input either: 's' or 'g'

Output:
purity_perC: The Cancer Cell frequency for the particular mutation

Example:
CCF_calc (34.11, 601, 65.9, 2, 1, 's') -> 103.52
which is CCF 103.52
'''

def CCF_calc (VAF, depth, purity, Y, cnmut, sg):
    '''
    They all assume that the cancer will follow the sample
    pattern when it comes to number of cancerous alleles per
    cancer cell.

    Also as general note, all capital letters are error bounds
    and all lowecase letters are the values with no error
    '''

    purity=float(purity)*0.01
    P=purity*0.1
    VAF=float(VAF)*0.01
    depth=float(depth)
    Y=int(Y)
    cnmut=int(cnmut)

    #Simplifying variable names for following equation
    p=purity
    f=VAF

    #Notes in notebook and LaTex on logic behind this
    if sg == 's':
        #Somatic Mutations
        a=1-p
        b=2*a
        c=Y*p
        d=b+c
        e=f*d
        h=cnmut * p
        g = e/h

        purity_prec=str(round((g*100),2))

    elif sg == 'g':
        #Germline mutations
        a=1-p
        b=2*a
        c=Y*p
        d=b+c
        e=f*d
        g=1-p
        h=e-g
        i=cnmut *p
        j=h/i

        purity_prec=str(round((j*100),2))

    else:
        #Tells the users if the input was incorrect
        purity_prec = 'Incorrect Input'

    return purity_prec
