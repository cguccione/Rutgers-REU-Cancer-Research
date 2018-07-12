#Caitlin Guccione
#Summer 2018
#Rutgers Cancer Institute of New Jersey
#DIMACS REU

'''
**This version outputs the data as a legit matrix 
This program monitors mutations and thier progression over time

Input:
3 TRF files from the same paitent at diffrent estimates

Output:
A matrix with the changes across muatations reflected that
can then be created into a tree

*This program uses hamming distance, or just simply if the
mutation is present or absent in order to track mutations.
In the future things such as the type of mutation and depth
may be taken into account.
'''

import numpy as np

#Opens the correct file and reads it into a variable for future usage
t1_file=open('TRF351303.xls', 'r')
t2_file=open('F2TRF351303.xls', 'r')
t3_file=open('F3TRF351303.xls', 'r')
t1_lines=t1_file.readlines()
t2_lines=t2_file.readlines()
t3_lines=t3_file.readlines()

def split_data(lines, germline):
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
    geneA_list=[]
    for line in lines:
        #This loop gets the TRF value from the first line
        if count == 1:
            none,TRF=line.split(':')
            TRF=TRF.strip('\n')
            TRF=TRF.strip('\t')
        #This section just contains the mutations
        elif count >=17:
            #This for loop makes sure the line is long enough to work, for example
            #things like applifications are cut out of the data
            test=0
            for word in line.split('\t'):
                test+=1
            #Now the clean mutations are filtered and used in the rest of the program
            if test >= 12:
                chro, start, end, gene, aminoAcid, mutation, allele, depth, strand, CN, LOHGICpath, LOHGICcomp = line.split('\t')
                geneA= str(gene) + '_' + str(aminoAcid)
                '''
                This finds all the germline mutations and sets them
                into a diffrent section since they will be present in
                every sample. We assume that all germline mutations will
                be found in the first time sample.

                *We also assume comp purity is correct and use LOHGIC's
                output from the COMP  purity
                '''
                trash, LOHGIC = LOHGICcomp.split(';')
                germ, somatic, gLOH, sLOH, gcmut, scmut = LOHGIC.split(',')
                if float(germ) > 0.5:
                        if int(germline) == 1:
                            germline_list.append(geneA)
                geneA_list.append(geneA)
        count=count+1
    return TRF, geneA_list

def compare_mutations(list1, list2):
    temp1=list1.copy()
    temp2=list2.copy()
    '''
    Finds diffrences in mutations
    Input: Two sets of mutations
    Outputs: A list of mutations appear in only one list and the length of that list
    '''
    count1= len(temp1)
    #Loops through the first list
    while count1 > 0:
        count1=count1-1
        count2= len(temp2)
        #Loops through the second list
        while count2 > 0:
            count2= count2 -1
            #Checks to see if any elements in the first list also appear in the second
            if temp1[count1] == temp2[count2]:
                #If a mutation appears in both lists it is removed from both lists
                bad=temp1[count1]
                temp1.remove(bad)
                temp2.remove(bad)
                break
    #The remaning elements found in only one list are combined
    final_list= temp1 + temp2
    count=len(final_list)
    return final_list, count

#Creates an empty list where all germline gene will be placed
germline_list=[]

#Places all gene and amino acid names into a list
t1_TRF, t1_geneA_list = split_data(t1_lines, 1)
t2_TRF, t2_geneA_list = split_data(t2_lines, 0)
t3_TRF, t3_geneA_list = split_data(t3_lines, 0)

'''
For matrix input:
x = Time 1
y = Time 2
z = Time 3
g = Germline
'''

#xg: Time 1 vs. Germline
xg_list, xg = compare_mutations(t1_geneA_list, germline_list)
#yg: Time 2 vs. Germline
yg_list, yg = compare_mutations(t2_geneA_list, germline_list)
#zg: Time 3 vs. Germline
zg_list, zg = compare_mutations(t3_geneA_list, germline_list)
#xy: Time 1 vs. Time 2
xy_list, xy = compare_mutations( t1_geneA_list, t2_geneA_list)
#xz: Time 1 vs. Time 3
xz_list, xz = compare_mutations( t1_geneA_list, t3_geneA_list)
#yz: Time 2 vs. Time 3
yz_list, yz = compare_mutations( t2_geneA_list, t3_geneA_list)

#Places the calculated diffrences into the matrix
matrix_L1=str(0) + ' ' + str(xg) + ' ' + str(yg) + ' ' + str(zg) + '; '
matrix_L2= str(xg) + ' ' + str(0) + ' ' + str(xy) + ' ' + str(xz) + '; '
matrix_L3= str(yg) + ' ' + str(xy) + ' ' + str(0) + ' ' + str(yz) + '; '
matrix_L4= str(zg) + ' ' + str(xz) + ' ' + str(yz) + ' ' + str(0)
matrix_input = matrix_L1 + matrix_L2 + matrix_L3 + matrix_L4
matrix_output = np.matrix(matrix_input)

#Creates combination of all three file names and writes the matrix to the file
t_combo= str(t1_TRF) + '_' + str(t2_TRF) + '_' + str(t3_TRF) + '.txt'
matrix_file=open(t_combo, 'w+')
matrix_file.write(str(matrix_output))
