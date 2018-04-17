import sys, os

from glob import glob
from functools import reduce

import pandas as pd
import numpy as np


def parse_single_run(files): 
	dic = {}
	tags = []

	for file in files: 
	    
	    tag = file.split('/')[-1].rstrip(".nodes.tsv")
	    tags.append(tag)
	    
	    with open(file, 'r') as f: 
	        
	        f.readline()
	        
	        for line in f.readlines(): 
	            
	            protein, degree, _, prize, node_type = line.rstrip().split('\t')
	            
	            if protein in dic: 
	                
	                dic[protein]["membership"].append(tag)
	            
	            else: 
	                
	                dic[protein] = {}
	                dic[protein]["membership"] = []
	                dic[protein]["log_degree"] = np.log2(int(float(degree)))
	                dic[protein]["prize"] = prize
	                dic[protein]["type"] = node_type

	out = []
	header = "\t".join(["gene", "prize", "type", "log_degree"]+tags)

	out.append(header)

	for protein in dic: 
	    row = "\t".join([protein, dic[protein]["prize"], dic[protein]["type"], str(dic[protein]["log_degree"])]+
	                    ["1" if tag in dic[protein]["membership"] else "0" for tag in tags])
	    out.append(row)

	return out


def parse_multi_runs(files, top_n=1000): 

	dfs_robustness = []
	dfs_specificity = []
	dic = {}

	for file in files: 
	    
	    tag = file.split('/')[-1].rstrip(".nodes.tsv")
	    
	    df = pd.read_csv(file, sep='\t')
	    df_robustness = df[["protein", "robustness"]]
	    df_robustness.columns = ["protein", tag]
	    dfs_robustness.append(df_robustness)
	    
	    df_specificity = df[["protein", "specificity"]]
	    df_specificity.columns = ["protein", tag]
	    dfs_specificity.append(df_specificity)
	    
	    with open(file, 'r') as f: 
	        f.readline()
	        for line in f.readlines(): 
	            protein, degree, _, prize, robustness, specificity, node_type = line.rstrip().split('\t')
	            if protein not in dic: 
	                dic[protein] = {}
	                dic[protein]["log_degree"] = np.log2(int(float(degree)))
	                dic[protein]["prize"] = prize
	                dic[protein]["type"] = node_type
	                
	                
	attribs = pd.DataFrame.from_dict(dic, orient="index")    
	    
	df_robustness_final = reduce(lambda left,right: pd.merge(left,right,on='protein', how='outer'), dfs_robustness)
	df_robustness_final.fillna(0, inplace=True)
	if df_robustness_final.shape[0] > top_n: 
		cutoff = sorted(df_robustness_final.sum(axis=1).tolist(), reverse=True)[top_n]
		df_robustness_final = df_robustness_final[df_robustness_final.sum(axis=1)>cutoff]

	df_specificity_final = reduce(lambda left,right: pd.merge(left,right,on='protein', how='outer'), dfs_specificity)
	df_specificity_final.fillna(0, inplace=True)
	if df_specificity_final.shape[0] > top_n: 
		cutoff = sorted(df_specificity_final.sum(axis=1).tolist(), reverse=True)[top_n]
		df_specificity_final = df_specificity_final[df_specificity_final.sum(axis=1)>cutoff]

	df_robustness_final = df_robustness_final.merge(attribs, left_on="protein", right_index=True, how='left')
	df_specificity_final = df_specificity_final.merge(attribs, left_on="protein", right_index=True, how='left')

	return df_robustness_final, df_specificity_final



def main(): 

	path = sys.argv[1]

	forest_paths = glob(os.path.join(path, "forest/*nodes.tsv"))
	robust_paths = glob(os.path.join(path, "robust_network/*nodes.tsv"))
	augmented_paths = glob(os.path.join(path, "augmented_forest/*nodes.tsv"))


	if len(forest_paths) > 0: 
		print("{} forest paths found.".format(len(forest_paths)))
		membership = parse_single_run(forest_paths)
		with open(os.path.join(path, "param_sweep_membership.tsv"), 'w') as f: 
		    f.write("\n".join(membership))

	if len(robust_paths) > 0: 
		print("{} robust paths found.".format(len(robust_paths)))
		robustness, specificity = parse_multi_runs(robust_paths)
		robustness.to_csv(os.path.join(path, "robust_network.robustness_attributes.tsv"), sep='\t', index=False)
		specificity.to_csv(os.path.join(path, "robust_network.specificity_attributes.tsv"), sep='\t', index=False)

	if len(augmented_paths) > 0: 
		print("{} augmented paths found.".format(len(augmented_paths)))
		robustness, specificity = parse_multi_runs(augmented_paths)
		robustness.to_csv(os.path.join(path, "augmented_forest.robustness_attributes.tsv"), sep='\t', index=False)
		specificity.to_csv(os.path.join(path, "augmented_forest.specificity_attributes.tsv"), sep='\t', index=False)





main()