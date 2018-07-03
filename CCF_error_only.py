#Caitlin Guccione
#Summer 2018
#Rutgers Cancer Institute of New Jersey
#DIMACS REU

import math
#This program finds the parameter estimates and confidence intervals for binomial data
from statsmodels.stats.proportion import proportion_confint as binofit

'''Input:
VAF: Allele Frequency for particular mutation
depth: Depth for particular mutation
Y: Total number of alleles per cell
cmut: Total number of cancerous alleles per cell
purity: entire sample purity
sg: Somatic or Germline input either: 's' or 'g'

Output:
purity_perC: The Cancer Cell frequency for the particular mutation
It is currently written in the form CCF; error

Example:
CCF_calc (34.11, 601, 65.9, 2, 1, 's') -> 103.52;19.92
which is CCF 103.52±19.92
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

    #Calculates the correct bounds on VAF using the bionomial parameter
    X = binofit((depth*VAF), depth, 0.01, method='binom_test')
    X=str(X)
    X=X.strip('()')
    Xup, Xlow=X.split(',')
    Xup=float(Xup)
    Xlow=float(Xlow)
    V=((Xup-Xlow)/2)


    #Simplifying variable names for following equation
    p=purity
    f=VAF
    F=V

    #Notes in notebook and LaTex on logic behind this
    if sg == 's':
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

        purity_prec=str(round((g*100),2)) + ';' + str(round((G*100),2))

    elif sg == 'g':
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

        purity_prec=str(round((j*100),2)) + ';' + str(round((J*100),2))

    else:
        #Tells the users if the input was incorrect
        purity_prec = 'Incorrect Input'

    return purity_prec
