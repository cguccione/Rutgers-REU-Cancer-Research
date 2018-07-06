#Caitlin Guccione
#Summer 2018
#Rutgers Cancer Institute of New Jersey
#DIMACS REU

'''
This file creates multiple fake mutations and gives them an allele frequency,
depth and purity. The allele frequency (VAF) and depth are pushed into one
file and all the correct purity's are placed into another file so that they
may be referenced later.

**This version gives every mutation an equal amount of being chosen and adds in
a few subclonal mutaions
'''

import math
import random
from scipy.stats import binom
import os.path

#Only change this number if you want an entire new set of data
set_count = 1

#Creates a new folder for all the files in the set to go into
current_directory = os.getcwd()
file_name="PuritySim_Set" + str(set_count)
final_directory = os.path.join(current_directory, file_name)
if not os.path.exists(final_directory):
   os.makedirs(final_directory)

#Creats the anwser key file
temp=final_directory + '\SIM_DATA_ANS_'+ str(set_count) +'.xls'
ans_file=open(temp, 'w+')
ans_file.write("File Number \t" + "Purity" + "\n")

def mut_calc(mut_o1, mut_o2, mut_o3, mut_o4, mut_o5, mut_o6, mut_o7, mut_o8, purity):
    '''Calculates and write the VAF and depth to the sample data files '''
    purity=round((purity* 0.01), 2)
    while mut_o1 > 0:
        Y= 2
        cnmut=1

        #Randomly chooses between clonal and subcloanl mutaions only for this
        #Option1. There is currently a bit less than 1/4 chance that this mutation turns
        #subclonal. This can be adjusted by changing the random int and while loop below
        sub=random.randint(1,8)
        if sub == 1 and purity > 0.2:
            purity= purity - (random.randint(5,10) * 0.01)
        if sub == 2 :
            purity = purity /2

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
    temp=final_directory + '\SIM_DATA_'+ str(set_count) + "." + str(samples_count) +'.xls'
    sim_file=open(temp, 'w+')
    sim_file.write("Allele Freq. \t" + "Depth \n")

    #Change the number of mutations per sample here:
    mut_count=random.randint(5,20)

    #The type of mutation is randomly selected
    #Option 1:Somatic CNmut=1
    #   If option 1 is choosen there is a 1/5 chance it will become subclonal
    #Option 2:Somatic LOH CNmut=2
    #Option 3: Somatic LOH CNmut=1
    #Option4: Germline Cmut=1
    #Option5: Germline LOHN Cmut=2
    #Option6: Germline LOH Cmut=1
    #Option7: Somatic Cmut=?, LOH
    #Option8: Germline Cmut=?, LOH
    mut_o1 =0
    mut_o2 =0
    mut_o3 =0
    mut_o4 =0
    mut_o5 =0
    mut_o6 =0
    mut_o7 =0
    mut_o8 =0

    count=0
    while count < mut_count:
        mutation=random.randint(1,8)
        if mutation == 1:
            mut_o1 = mut_o1 + 1
        elif mutation == 2:
            mut_o2 = mut_o2 + 1
        elif mutation == 3:
            mut_o3 = mut_o3 + 1
        elif mutation == 4:
            mut_o4 = mut_o4 + 1
        elif mutation == 5:
            mut_o5 = mut_o5 + 1
        elif mutation == 6:
            mut_o6 = mut_o6 + 1
        elif mutation == 7:
            mut_o7 = mut_o7 + 1
        else:
            mut_o8 = mut_o8 + 1
        count = count +1
        
    mut_calc(mut_o1, mut_o2, mut_o3, mut_o4, mut_o5, mut_o6, mut_o7, mut_o8, purity)

    samples_count = samples_count+1
