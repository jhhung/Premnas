import csv
import itertools
import math
import pickle
import sys
import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None

global threshold
global con_threshold

class Drug:
    
    def __init__(self, name, dose, time, inst, effect):
        self.name = name
        self.dose = dose
        self.time = time
        self.id = '{}_{:.3f}_{}'.format(self.name, self.dose, self.time)
        self.inst = inst
        self.effect = np.array(effect)

    def __eq__(self, other):
        return self.dose == other.dose and self.time == other.time

    def __lt__(self, other):
        return (self.dose > other.dose)

    def __gt__(self, other):
        return (self.dose < other.dose)

    def __str__(self):
        return '{} {}uM {}h'.format(self.name, self.dose, self.time)


def append_df_effect(df, drug_id, result_effect):
    n_killed_all = len([x for x in result_effect if x <= threshold])
    df.loc[drug_id] = result_effect + [n_killed_all]
    return df


def read_metadata(filename):
    df_metadata = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
    grouped_plate = df_metadata.groupby(['det_plate'])[['pert_iname','pert_dose','pert_dose_unit','pert_time','pert_time_unit']]
    list_plate_df = [grouped_plate.get_group(x) for x in grouped_plate.groups]
    return list_plate_df


def add_consistency_info(df, DICT_DRUG):
    for cluster_id in df.columns[:-1]:
        colname = 'con_{}'.format(cluster_id)
        df[colname] = 0
        for drug_id in df.index:
            if(df.loc[drug_id, cluster_id] <= threshold):
                drug = DICT_DRUG[drug_id]
                keys = [x for x in DICT_DRUG.keys() if DICT_DRUG[x].name == drug.name and DICT_DRUG[x].dose >= drug.dose and DICT_DRUG[x].time >= drug.time]
                if(df.loc[keys,cluster_id].quantile(q=0.5) <= con_threshold):
                    df.loc[drug_id, colname] = 1
    return df


def calEffectConsistency(df):
    mid_col = int(len(df.columns)/2)
    df_eff = df[df.columns[:mid_col+1]]
    df_con = df[df.columns[mid_col+1:]]
    for i, j in zip(df_eff.columns.to_list()[:-1], df_con.columns.to_list()):
        tmp = df_eff[i] * df_con[j]
        df_eff[i] = tmp.values
    df_eff['kill_all_count'] = df_eff.apply(lambda row: len([x for x in row if x <= threshold]), axis=1)
    return df_eff


def update_df_effect(df, removed_clusters = []):
    #print('\nremoved_clusters: ', removed_clusters)
    df = df.drop(columns = removed_clusters)
    df['kill_all_count'] = df.apply(lambda row: len([x for x in row if x <= threshold]), axis=1)
    return df

def choose_least_dose(drug_ids):
    if len(drug_ids) <= 1 :
        return drug_ids
    else:
        drug_list = [DICT_DRUG[x] for x in drug_ids]
        least_dose = max(drug_list).dose
        return [x.id for x in drug_list if x.dose == least_dose]

def select_candidate_drugs(df, _max):
    if(_max <= 0): return []
    #print('\n1st step: kill same amount clusters ( %d clusters) : ' %  _max)
    same_ability_drugs = df[df['kill_all_count'].values == _max].index.to_list()
    #choose the one with minimum dose
    min_dose_drugs = choose_least_dose(same_ability_drugs)
    return min_dose_drugs


def find_drug(df, solution=[], LIST_SOLUTION=[]):
    # All cell types were killed
    if (len(df.columns) == 1) :
        LIST_SOLUTION.append(sorted(solution))
    # If no perturbation could induce complete apoptosis of any subpopulation
    elif (df.min().min() > threshold):
        i_col = 0
        while (df[df.columns[i_col]].min() > -0.5 and df.columns[i_col] != df.columns[-1]):
            #print('\nNo candidates for %s' % df.columns[i_col])
            i_col += 1
        col = df.columns[i_col]
        if(col != df.columns[-1]):
            candidates = [x for x in df.index if df.loc[x,col] == df[col].min()]
            for drug in candidates:
                solution.append(drug)
                df_tmp = update_df_effect(df, df.columns[:i_col+1])
                find_drug(df_tmp, solution, LIST_SOLUTION)
                solution.pop()
        df = update_df_effect(df, df.columns[:-1])
        find_drug(df, solution, LIST_SOLUTION)
    else :
        candidates = select_candidate_drugs(df, df['kill_all_count'].values.max())
        for drug in candidates:
            solution.append(drug)
            killed_clusters = [x for x in df.columns if df.loc[drug,x] <= threshold]
            df_tmp = update_df_effect(df, killed_clusters)
            find_drug(df_tmp, solution, LIST_SOLUTION)      
            solution.pop()



if __name__ == '__main__':
    if(len(sys.argv)<4):
        print("Usage: python3 .py $meta_data.txt $compostion.csv $threshold $consistency_threshold")
        exit(0)

    print('--------Preprocessing--------')
#####
    # with open(sys.argv[1], 'rb') as fh:
    #     df_effect = pickle.load(fh)
    # with open(sys.argv[2], 'rb') as fh:
    #     DICT_DRUG = pickle.load(fh)
#####

    threshold = float(sys.argv[3])
    con_threshold = float(sys.argv[4])
    print('threshold: {}, con_threshold: {} '.format(threshold, con_threshold))
    
    print('Reading the metadata file and the decomposition results...')
    # read metadata
    list_df_metadata = read_metadata(sys.argv[1])

    # read decomposition results
    if( '.csv' in sys.argv[2]):
        df_comp = pd.read_csv(sys.argv[2], index_col=0)
    else:
        df_comp = pd.read_csv(sys.argv[2], index_col=0, sep='\t')
    df_comp = df_comp.reindex(sorted(df_comp.columns), axis=1)
    df_percent = pd.DataFrame(data=0, index = df_comp.index, columns = sorted(df_comp.columns[:-3].to_list()))
    n_clusters = len(df_comp.columns) - 3 
    
    col_names = df_comp.columns[:-3].to_list()+ ['kill_all_count']

    df_effect = pd.DataFrame(columns=col_names)
    DICT_DRUG_PRE = {}

    for df_metadata in list_df_metadata:
        df_result = pd.merge(df_metadata, df_percent, how='left', left_index=True, right_index=True, validate='one_to_one')
        #print(df_result.shape)
        
        ctrl_inst_ids = df_result[df_result['pert_dose']<0].index
        pert_inst_ids = df_result[df_result['pert_dose']>0].index
        df_result = df_result[df_result['pert_dose']>0]

        # average ctrl results
        df_ctrl_composition = df_comp.loc[ctrl_inst_ids,df_comp.columns[:n_clusters]]
        ctrl_composition = df_ctrl_composition.mean()
        ctrl_composition = ctrl_composition/np.sum(ctrl_composition.values)

        #print('after dropping control vehicle: ',df_result.shape)
        
        # fill the result dataframe
        for perturbation in df_result.index:
            drug_composition = df_comp.loc[perturbation,df_comp.columns[:n_clusters]]
            change_list = (drug_composition-ctrl_composition)/(drug_composition+ctrl_composition)

            # fill nan (occurs when 0/0) with 2 
            change_list = [ 2 if math.isnan(x) else x for x in change_list]
            df_result.loc[perturbation,df_result.columns[-n_clusters:]] = change_list

            drug = Drug(df_metadata.loc[perturbation,'pert_iname'],\
                        df_metadata.loc[perturbation,'pert_dose'],\
                        df_metadata.loc[perturbation,'pert_time'],\
                        perturbation,\
                        change_list)
            drug_id = '{}_{:.3f}_{}'.format(\
                df_metadata.loc[perturbation,'pert_iname'],\
                df_metadata.loc[perturbation,'pert_dose'],\
                df_metadata.loc[perturbation,'pert_time'])

            if drug_id not in DICT_DRUG_PRE:
                DICT_DRUG_PRE[drug_id] = [drug]
            else:
                DICT_DRUG_PRE[drug_id].append(drug)

    #print('DICT_DRUG_PRE: ', len(DICT_DRUG_PRE))

    print('Calculating drug effects...')
    # average replicates
    DICT_DRUG = {}
    for drug in DICT_DRUG_PRE:
        # only keep those have replicates
        if (len(DICT_DRUG_PRE[drug]) <=2):
            continue
        avg_effect = np.array([0.0]*n_clusters)
        inst_list = []
        counts = np.array([0.0]*n_clusters)
        
        for x in DICT_DRUG_PRE[drug]:
            for i in range(0,len(x.effect)):
                if(x.effect[i] < 2): # The subpopulation existed in the control samples
                    avg_effect[i] += x.effect[i]
                    counts[i] += 1
        for i in range(0, len(avg_effect)):
            if counts[i] >= 3:
                avg_effect[i] /= counts[i]
            else:
                avg_effect[i] = 0
        
        x = DICT_DRUG_PRE[drug][0]
        avg_drug = Drug(x.name, x.dose, x.time, inst_list, avg_effect)
        DICT_DRUG[avg_drug.id] = avg_drug
        df_effect = append_df_effect(df_effect, avg_drug.id, avg_effect.tolist())

    df_effect = add_consistency_info(df_effect, DICT_DRUG)


    # store output
    fh = open(file = '{}_t{}_ct{}_df_effect.pickle'.format(sys.argv[2].rsplit('.',0)[0], threshold, con_threshold), mode='wb')
    pickle.dump(df_effect, fh) 
    fh.close()
    fh = open(file = '{}_t{}_ct{}_DICT_DRUG_PRE.pickle'.format(sys.argv[2].rsplit('.',0)[0], threshold, con_threshold), mode='wb')
    pickle.dump(DICT_DRUG_PRE, fh)
    fh = open(file = '{}_t{}_ct{}_DICT_DRUG.pickle'.format(sys.argv[2].rsplit('.',0)[0], threshold, con_threshold), mode='wb')
    pickle.dump(DICT_DRUG, fh)
    fh.close()


    print('--------Subpopulation Analysis--------')
    # subpopulation analysis
    df_effect = calEffectConsistency(df_effect)
    grouped = df_effect.groupby('kill_all_count')
    print(grouped.size())
    for name, group in grouped:
        print(name, group)
    

    # subpopulation analysis
    for cluster in df_effect.columns[:-1]:
        print(cluster)
        pert = df_effect[df_effect[cluster] <= threshold].index.to_list()
        print('can be killed by {} perturbations'.format(len(pert)))
        print('peturbation with best efficacy: {} - {}'.format(df_effect[cluster].idxmin(), df_effect[cluster].min()))

    

    print('--------Treatment Selection--------')
    LIST_SOLUTION = []

    # find a cocktail therapy
    find_drug(df_effect, LIST_SOLUTION=LIST_SOLUTION)
    
  
    LIST_SOLUTION = sorted(LIST_SOLUTION)
    list_result = list(LIST_SOLUTION for LIST_SOLUTION,_ in itertools.groupby(LIST_SOLUTION))

#    with open('{}_solution_list_t{}_cont{}.csv'.format(sys.argv[2].rsplit('.',0)[0], threshold, con_threshold),"w+") as f:
    with open('Treatment-selection-output.csv',"w+") as f:
        writer = csv.writer(f)
        writer.writerows(list_result)

#    print('Done! Cocktail therapy is stored in {}_solution_list_t{}_cont{}.csv'.format(sys.argv[2].rsplit('.',0)[0], threshold, con_threshold))




