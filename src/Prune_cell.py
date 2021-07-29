
import numpy as np
#import pandas as pd 

thredshold_propotion = 0.6
pruned_cell_label = 0
dataset_dir = 'GSE114459' 

archetype_propotion = [] 

f = open('/output_dir/archetypal-explicit-function.txt', 'r')
line = f.readline()
line = f.readline()
while line:
    line = line.split(' ')
    line[-1] = line[-1].strip('\n')
    line = [float(each) for each in line[1:]]
    archetype_propotion.append(line)
    line = f.readline()
f.close()

assigned_archetype = [] 
assigned_cluster = []

def read_cluster():
    f = open('/output_dir/assigned-subpopulation.txt' , 'r')
    cluster_name = [] 
    line = f.readline()
    line = f.readline()
    while line:
        line = line.split(' ')
        line = line[1].strip('\n').strip('"').strip('c')
        cluster_name.append(int(line))
        line = f.readline()
    f.close()
    return cluster_name 

assigned_cluster = read_cluster()

new_assigned_cluster = []
was_pruned_count = 0
for i in range(len(archetype_propotion)):
    max_ = np.max(archetype_propotion[i])
    tmp = archetype_propotion[i].index(max_)

    if max_ < thredshold_propotion:
        new_assigned_cluster.append(pruned_cell_label)
        was_pruned_count += 1
    else:
        if dataset_dir =='PBMC' and assigned_cluster[i]==9:
            new_assigned_cluster.append(pruned_cell_label)
        else:
            new_assigned_cluster.append(assigned_cluster[i])

print("total pruned out: " + str(was_pruned_count))

f = open('/output_dir/pruned-assigned-subpopulation.txt', 'w')
for i in range(len(new_assigned_cluster)):
    f.write('\"' + str(i+1) + '\" ' +  str(new_assigned_cluster[i]) + '\n')
f.close()

