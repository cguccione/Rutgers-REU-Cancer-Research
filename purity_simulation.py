#Caitlin Guccione
#Summer 2018
#Rutgers Cancer Institute of New Jersey
#DIMACS REU

'''
This file creates multiple fake mutations and gives them an allele frequency,
depth and purity. The allele frequency (VAF) and depth are pushed into one
file and all the correct purity's are placed into another file so that they
may be referenced later.
'''

import math
import random
from scipy.stats import binom

#Only change this number if you want an entire new set of data
set_count = 3
#Creats the anwser key file
temp='SIM_DATA_ANS_'+ str(set_count) +'.xls'
ans_file=open(temp, 'w+')
ans_file.write("File Number \t" + "Purity" + "\n")

def mut_calc(mut_o1, mut_o2, mut_o3, mut_o4, mut_o5, mut_o6, mut_o7, mut_o8, purity):
    '''Calculates and write the VAF and depth to the sample data files '''
    purity=round((purity* 0.01), 2)
    while mut_o1 > 0:
        Y= 2
        cnmut=1
        #Calculates the allele frequency
        f= purity/2
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o1 = mut_o1 - 1
    while mut_o2 > 0:
        Y= 2
        cnmut=2
        #Calculates the allele frequency
        f= purity/(2-purity)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o2 = mut_o2 - 1
    while mut_o3 > 0:
        Y= 1
        cnmut=1
        #Calculates the allele frequency
        f= purity
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o3 = mut_o3 - 1
    while mut_o4 > 0:
        Y= 2
        cnmut=1
        #Calculates the allele frequency
        f= 1/2
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o4 = mut_o4 - 1
    while mut_o5 > 0:
        Y= 2
        cnmut=2
        #Calculates the allele frequency
        f= 1/(2-purity)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o5 = mut_o5 - 1
    while mut_o6 > 0:
        Y= 1
        cnmut=1
        #Calculates the allele frequency
        f= (1 + purity)/2
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o6 = mut_o6 - 1
    while mut_o7 > 0:
        #Selects a random number of alleles from 3 to 8
        Y= random.randint(3,8)
        cnmut=Y
        #Calculates the allele frequency
        f= (cnmut * purity)/((2*(1-purity)) + (Y * purity))
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o7 = mut_o7 - 1
    while mut_o8 > 0:
        #Selects a random number of alleles from 3 to 8
        Y= random.randint(3,8)
        cnmut=Y
        #Calculates the allele frequency
        f= ((1-purity) + (cnmut * purity))/((2 * (1 - purity)) + (Y * purity))
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds noies to the allele frequency to make it more realistic
        VAF=round(((binom.rvs(depth, f))/depth),2)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) +  "\n")
        mut_o8 = mut_o8 - 1

#Changes the amount of samples you want per anwser file
samples=10
samples_count=1

while samples_count <= samples:

    #Selects a random purity from 10-90
    purity=random.randint(10,90)
    #Adds the correct purity to the anwser file
    ans_file.write(str(set_count) + "." + str(samples_count) + "\t" + str(purity) + "\n")

    #Creates the sample data files
    temp='SIM_DATA_'+ str(set_count) + "." + str(samples_count) +'.xls'
    sim_file=open(temp, 'w+')
    sim_file.write("Allele Freq. \t" + "Depth \n")

    #Change the number of mutations per sample here:
    mut_count=random.randint(5,20)

    #Divides the mutations into somatic/germline mutations
    mut_germline = int(round(mut_count/3))
    mut_somatic = mut_count - mut_germline

    #Divides the somatic mutations into various CNmuts/LOH
    #Option 1:Somatic CNmut=1
    mut_o1= int(math.ceil(mut_somatic/2))
    #Option 2:Somatic LOH CNmut=2
    mut_o2= int(math.floor((mut_somatic-mut_o1)/3))
    #Option 3: Somatic LOH CNmut=1
    mut_o3= int(math.ceil((mut_somatic-mut_o1))/2)
    #Option7: Somatic Cmut=?, LOH
    mut_o7= int(math.floor((mut_somatic-mut_o1)- mut_o2)/2)

    #Divides the germline mutations into various CNmuts/LOH
    #Option4: Germline Cmut=1
    mut_o4= int(math.ceil(mut_germline/2))
    #Option5: Germline LOHN Cmut=2
    mut_o5= int(math.floor((mut_germline-mut_o4)/3))
    #Option6: Germline LOH Cmut=1
    mut_o6= int(math.ceil((mut_germline-mut_o4)- mut_o5)/2)
    #Option8: Germline Cmut=?, LOH
    mut_o8= int(math.floor((mut_germline-mut_o4)- mut_o5)/2)

    mut_calc(mut_o1, mut_o2, mut_o3, mut_o4, mut_o5, mut_o6, mut_o7, mut_o8, purity)

    samples_count = samples_count+1
