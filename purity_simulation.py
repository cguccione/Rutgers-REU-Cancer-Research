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

#Only change this number if you want an entire new set of data
set_count = 2
#Creats the anwser key file
temp='SIM_DATA_ANS_'+ str(set_count) +'.xls'
ans_file=open(temp, 'w+')
ans_file.write("File Number \t" + "Purity" + "\n")

def mut_calc(mut_o1, mut_o2, mut_o3, mut_o4, mut_o5, mut_o6, purity):
    '''Calculates and write the VAF and depth to the sample data files '''
    purity=round((purity* 0.01), 2)
    while mut_o1 > 0:
        #Calculates the allele frequency
        f= purity/2
        #Adds noies to the allele frequency to make it more realistic
        VAF= round(random.normalvariate(f, 0.035), 2)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o1 = mut_o1 - 1
    while mut_o2 > 0:
        #Calculates the allele frequency
        f= purity/(2-purity)
        #Adds noies to the allele frequency to make it more realistic
        VAF= round(random.normalvariate(f, 0.035), 2)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o2 = mut_o2 - 1
    while mut_o3 > 0:
        #Calculates the allele frequency
        f= purity
        #Adds noies to the allele frequency to make it more realistic
        VAF= round(random.normalvariate(f, 0.035), 2)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o3 = mut_o3 - 1
    while mut_o4 > 0:
        #Calculates the allele frequency
        f= 1/2
        #Adds noies to the allele frequency to make it more realistic
        VAF= round(random.normalvariate(f, 0.035), 2)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o4 = mut_o4 - 1
    while mut_o5 > 0:
        #Calculates the allele frequency
        f= 1/(2-purity)
        #Adds noies to the allele frequency to make it more realistic
        VAF= round(random.normalvariate(f, 0.035), 2)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o5 = mut_o5 - 1
    while mut_o6 > 0:
        #Calculates the allele frequency
        f= (1 + purity)/2
        #Adds noies to the allele frequency to make it more realistic
        VAF= round(random.normalvariate(f, 0.035), 2)
        #The depth is just a random number from 300-1000
        depth=random.randint(300, 1000)
        #Adds the data to the file
        sim_file.write(str(VAF*100) + "\t" + str(depth) + "\n")
        mut_o6 = mut_o6 - 1

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
    sim_file.write("Allele Freq. \t" + "Depth" + "\n")

    #Change the number of mutations per sample here:
    mut_count=random.randint(5,12)

    #Divides the mutations into somatic/germline mutations
    mut_germline = int(round(mut_count/3))
    mut_somatic = mut_count - mut_germline

    #Divides the somatic mutations into various CNmuts/LOH
    #Option 1:Somatic CNmut=1
    mut_o1= int(math.ceil(mut_somatic/2))
    #mut_calc(mut_o1, 1, purity)
    #Option 2:Somatic LOH CNmut=2
    mut_o2= int(math.ceil((mut_somatic-mut_o1)/2))
    #mut_calc(mut_o2, 2, purity)
    #Option 3: Somatic LOH CNmut=1
    mut_o3= int(math.floor((mut_somatic-mut_o1))/2)
    #mut_calc(mut_o3, 3, purity)

    #Divides the germline mutations into various CNmuts/LOH
    #Option4: Germline Cmut=1
    mut_o4= int(math.ceil(mut_germline/2))
    #mut_calc(mut_o4, 4, purity)
    #Option5: Germline LOHN Cmut=2
    mut_o5= int(math.ceil((mut_germline-mut_o4)/2))
    #mut_calc(mut_o5, 5, purity)
    #Option6: Germline LOH Cmut=1
    mut_o6= int(math.floor((mut_germline-mut_o4))/2)
    mut_calc(mut_o1, mut_o2, mut_o3, mut_o4, mut_o5, mut_o6, purity)
    samples_count = samples_count+1
