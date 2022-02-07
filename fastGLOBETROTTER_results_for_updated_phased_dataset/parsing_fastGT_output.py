'''
Script name: parsing_fastGT_output.py
Author: Jiawei Zhao
Description:
    This python script parses ".main.txt" and ".boot.txt" output files for historical admixture inference from
    fastGT (https://github.com/sahwa/fastGLOBETROTTER) with a .csv file of annotation for surrogate
    populations and their corresponding regions. Ultimately, a row (for "one-date" as conclusion) or multiple
    rows (for "one-date-multiway" or "multiple-dates" as conclusion) of summary would be output to a 
    user-specified .csv file. The row(s) of summary essentially consists of five-parts:
    1. Conclusion of (a successful) fastGT inference: one-date, one-date-multiway or multiple-dates
    2. The order of admixture event, date1, date2 etc
    3. The order of principal-component, 1st, 2nd etc 
    4. Components of the admixture event by regions of surrogate populations listed at the .csv file of 
    annotation for surrogate populations specified at "--surr-pop-info", listing summed-up percentage for
    surrogate populations from each region like ADMIXTURE (https://dalexander.github.io/admixture/). The 
    percentages could be easily used for stacked barplots
    5. Inferred value of admixture date (N generations ago) at ".main.txt" and 95%CI of bootstrapping values
    at ".boot.txt" (only available for date1)
Required Non-Standard Python Libraries: pandas, numpy, scipy
Usage:
    python parsing_fastGT_output.py --input FastGT_Poland-AJ --pop Poland-AJ --output fastGT_Jews_Turks \
    --surr-pop-info 180pops_info.csv 
'''

import pandas as pd
import numpy as np
import argparse
import re
from scipy.stats import t
from scipy.stats import sem
from os.path import exists as fexists

def parsed_args(): # this function aims at parsing arguments with the library "argparse"
    ## Construct an argument parser
    argparser = argparse.ArgumentParser()
    ## Define the five args that our script needs
    argparser.add_argument("--input", help="Specifies the prefix of .main.txt and .boot.txt input files to parse",\
                           required=True)
    argparser.add_argument("--pop", help="Specifies the name of fastGLOBETROTTER-inferred population",\
                           required=True)
    argparser.add_argument("--output", help="Specifies the prefix of .csv file for parsed output", required=True)
    argparser.add_argument("--surr-pop-info", help="Specifies the prefix of .csv file listing information of \
                           surrogate populations in fastGT inference", required=True)
    ## Parse the arguments
    args = argparser.parse_args()
    input_prefix = str(args.input) 
    inferred_pop = str(args.pop)
    output_prefix = str(args.output)
    surr_pop_info = str(args.surr_pop_info)
    return (input_prefix, inferred_pop, output_prefix, surr_pop_info)

def parse_main_txt(input_f: str, inferred_pop: str, surr_pop_info_df: pd.DataFrame): 
    # returns to a dataframe consisting of inferred admixture date & proportions
    output_info_dict_list = []
    with open(input_f, "r") as f:
        lines = [line.strip() for line in f.readlines()]
        concl_line = lines[0]
        # if conclusion is either "unclear signal" or "uncertain", don't output anything
        if re.search("(unclear|uncertain)", concl_line): return ""
        # otherwise, start to parse $surr_pop_info_df
        pop2group_dict = surr_pop_info_df.set_index("pop")["pop_group"].to_dict()
        #print(pop2group_dict) # only for testing
        # parse admixture components according to conclusion
        if re.search("one-date", concl_line):
            onedate_gen = round(float(lines[4].split(" ")[0]), 2) # inferred single-event admixture date
            surr_PC1_comp1_perc = float(lines[4].split(" ")[1]) 
            surr_PC1_comp2_perc = 1 - surr_PC1_comp1_perc
            surr_PC1_comp1_pops = lines[11].split(" ")[1:] # list of surrogate populations at comp1 of PC1
            surr_PC1_comp1_pops_perc = lines[12].split(" ")[1:] # perc for each surr pop at $surr_PC1_comp1_pops
            # obtain info for comp2 of PC1
            surr_PC1_comp2_pops = lines[13].split(" ")[1:]
            surr_PC1_comp2_pops_perc = lines[14].split(" ")[1:]
            group_perc_dict = {} # use this dict to record subsequent percentage to output
            # loop through comp1 of PC1
            for pop, perc in zip(surr_PC1_comp1_pops, surr_PC1_comp1_pops_perc):
                if pop not in pop2group_dict: 
                    print(pop)
                    return "."
                perc = float(perc)
                perc = round(perc * surr_PC1_comp1_perc, 4)
                group = pop2group_dict[pop]
                group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
            # loop through comp2 of PC1
            for pop, perc in zip(surr_PC1_comp2_pops, surr_PC1_comp2_pops_perc):
                if pop not in pop2group_dict: 
                    print(pop)
                    return "."
                perc = float(perc)
                perc = round(perc * surr_PC1_comp2_perc, 4)
                group = pop2group_dict[pop]
                group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
            # fill the columns of "unpresent" groups with 0
            for el in set(list(pop2group_dict.values())):
                if not re.search("Jew", el):    
                    group_perc_dict[el] = 0 if el not in group_perc_dict else round(group_perc_dict[el], 4)
            # if we only need to output PC1, only one row would be returned
            admix_info_dict = group_perc_dict
            admix_info_dict['N_generations_ago'] = onedate_gen
            admix_info_dict['event'] = 1
            admix_info_dict['PC'] = 1
            admix_info_dict['conclusion'] = "one-date"
            output_info_dict_list.append(admix_info_dict) # add the row for PC1
            # in this case, we also need to add info for PC2. In this case, two rows would be returned
            if re.search("multiway\)$", concl_line): 
                admix_info_dict['conclusion'] = "one-date-multiway"
                output_info_dict_list = [] # reset the list to empty
                output_info_dict_list.append(admix_info_dict) # add the row for PC1
                admix_info_dict, group_perc_dict = {}, {} # reset this temp variables to empty
                surr_PC2_comp1_perc = float(lines[4].split(" ")[-3]) 
                surr_PC2_comp2_perc = 1 - surr_PC2_comp1_perc
                surr_PC2_comp1_pops = lines[17].split(" ")[1:]
                surr_PC2_comp2_pops = lines[19].split(" ")[1:]
                surr_PC2_comp1_pops_perc = lines[18].split(" ")[1:]
                surr_PC2_comp2_pops_perc = lines[20].split(" ")[1:]
                for pop, perc in zip(surr_PC2_comp1_pops, surr_PC2_comp1_pops_perc):
                    if pop not in pop2group_dict: 
                        print(pop)
                        return "."
                    perc = float(perc)
                    perc = round(perc * surr_PC2_comp1_perc, 4)
                    group = pop2group_dict[pop]
                    group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
                for pop, perc in zip(surr_PC2_comp2_pops, surr_PC2_comp2_pops_perc):
                    if pop not in pop2group_dict: 
                        print(pop)
                        return "."
                    perc = float(perc)
                    perc = round(perc * surr_PC2_comp2_perc, 4)
                    group = pop2group_dict[pop]
                    group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
                # fill the columns of "unpresent" groups with 0
                for el in set(list(pop2group_dict.values())):
                    if not re.search("Jew", el):
                        group_perc_dict[el] = 0 if el not in group_perc_dict else round(group_perc_dict[el], 4)
                admix_info_dict = group_perc_dict 
                admix_info_dict['N_generations_ago'] = onedate_gen
                admix_info_dict['event'] = 1
                admix_info_dict['PC'] = 2
                admix_info_dict['conclusion'] = "one-date-multiway"
                output_info_dict_list.append(admix_info_dict) # add the row for PC2
        elif re.search("multiple-dates", concl_line):
            twodate_gen1 = round(float(lines[8].split(" ")[0]), 2) # inferred admixture date1
            twodate_gen2 = round(float(lines[8].split(" ")[1]), 2) # inferred admixture date2
            # obtain info for comp1 of PC1 at date1
            surr_date1_comp1_perc = float(lines[24].split(" ")[0]) 
            # list of surrogate populations at comp1 of PC1 at date1
            surr_date1_comp1_pops = lines[23].split(" ")[1:] 
            # perc for each surr pop at $surr_PC1_comp1_pops at date1
            surr_date1_comp1_pops_perc = lines[24].split(" ")[1:]             
            # obtain info for comp2 of PC1 at date1
            surr_date1_comp2_perc = 1 - surr_date1_comp1_perc
            surr_date1_comp2_pops = lines[25].split(" ")[1:]
            surr_date1_comp2_pops_perc = lines[26].split(" ")[1:]
            ## obtain information for date2 admixture event
            surr_date2_comp1_perc = float(lines[30].split(" ")[0]) 
            # list of surrogate populations at comp1 of PC1 at date2
            surr_date2_comp1_pops = lines[29].split(" ")[1:] 
            # perc for each surr pop at $surr_PC1_comp1_pops at date2
            surr_date2_comp1_pops_perc = lines[30].split(" ")[1:]             
            # obtain info for comp2 of PC1 at date2
            surr_date2_comp2_perc = 1 - surr_date2_comp1_perc
            surr_date2_comp2_pops = lines[31].split(" ")[1:]
            surr_date2_comp2_pops_perc = lines[32].split(" ")[1:]
            group_perc_dict, admix_info_dict = {}, {} # use this dict to record subsequent percentage to output
            # gather information for date1 event
            for pop, perc in zip(surr_date1_comp1_pops, surr_date1_comp1_pops_perc):
                if pop not in pop2group_dict: 
                    print(pop)
                    return "."
                perc = float(perc)
                perc = round(perc * surr_date1_comp1_perc, 4)
                group = pop2group_dict[pop]
                group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc 
            for pop, perc in zip(surr_date1_comp2_pops, surr_date1_comp2_pops_perc):
                if pop not in pop2group_dict: 
                    print(pop)
                    return "."
                perc = float(perc)
                perc = round(perc * surr_date1_comp2_perc, 4)
                group = pop2group_dict[pop]
                group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
            # fill the columns of "unpresent" groups with 0
            for el in set(list(pop2group_dict.values())):
                if not re.search("Jew", el):
                    group_perc_dict[el] = 0 if el not in group_perc_dict else round(group_perc_dict[el], 4)
            # fill the row of information for admixture event1
            admix_info_dict = group_perc_dict 
            admix_info_dict['N_generations_ago'] = twodate_gen1
            admix_info_dict['event'] = 1
            admix_info_dict['PC'] = 1
            admix_info_dict['conclusion'] = "multiple-dates"
            output_info_dict_list.append(admix_info_dict) # add the row for date1
            # gather information for date2 event
            group_perc_dict, admix_info_dict = {}, {} 
            for pop, perc in zip(surr_date2_comp1_pops, surr_date2_comp1_pops_perc):
                if pop not in pop2group_dict: 
                    print(pop)
                    return "."
                perc = float(perc)
                perc = round(perc * surr_date2_comp1_perc, 4)
                group = pop2group_dict[pop]
                group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
            for pop, perc in zip(surr_date2_comp2_pops, surr_date2_comp2_pops_perc):
                if pop not in pop2group_dict: 
                    print(pop)
                    return "."
                perc = float(perc)
                perc = round(perc * surr_date2_comp2_perc, 4)
                group = pop2group_dict[pop]
                group_perc_dict[group] = group_perc_dict[group]+perc if group in group_perc_dict else perc
            # fill the columns of "unpresent" groups with 0
            for el in set(list(pop2group_dict.values())):
                if not re.search("Jew", el):
                    group_perc_dict[el] = 0 if el not in group_perc_dict else round(group_perc_dict[el], 4)
            admix_info_dict = group_perc_dict 
            admix_info_dict['N_generations_ago'] = twodate_gen2
            admix_info_dict['event'] = 2
            admix_info_dict['PC'] = 1
            admix_info_dict['conclusion'] = "multiple-dates"
            output_info_dict_list.append(admix_info_dict) # add the row for date2
        else:
            return ""
    return pd.DataFrame.from_dict(output_info_dict_list)
    
def parse_boot_txt(input_f: str, fastGT_main_info: pd.DataFrame): 
    # returns to a dataframe consisting of mean/95% CI for bootstrapping values for admixture dates
    boot_df = pd.read_csv(input_f, header=0, sep=" ")
    gen_boot_col = boot_df.columns[1]
    boot_gens = np.array(boot_df[gen_boot_col])
    CI95_scale = t.ppf((1+0.95)/2, boot_gens.shape[0]-1)
    SE_gens = sem(boot_gens)
    low_95CI = np.mean(boot_gens) - CI95_scale*SE_gens
    high_95CI = np.mean(boot_gens) + CI95_scale*SE_gens
    fastGT_main_info['boot_gen_low95_CI'] = ""
    fastGT_main_info['boot_gen_high95_CI'] = ""
    for i in range(fastGT_main_info.shape[0]):
        row_i = fastGT_main_info.iloc[i]
        if row_i["event"] == 1:
            fastGT_main_info['boot_gen_low95_CI'].iloc[i] = round(low_95CI, 2)
            fastGT_main_info['boot_gen_high95_CI'].iloc[i] = round(high_95CI, 2)
    return fastGT_main_info

if __name__ == "__main__":
    ## Obtain the input args for this script
    input_prefix, inferred_pop, output_prefix, surr_pop_info = parsed_args()
    surr_pops_info_df = pd.read_csv(surr_pop_info)
    srr_pop_df_cols = {el for el in surr_pops_info_df.columns}
    # input check control
    if "pop" not in srr_pop_df_cols or "pop_group" not in srr_pop_df_cols:
        print("Info file for surrogate populations is incorrect!")
        exit(1)
    input_main_f = f"{input_prefix}.main.txt"
    input_boot_f = f"{input_prefix}.boot.txt"
    # input check control
    if not fexists(input_main_f) and not fexists(input_boot_f):
        print("Input fastGT info file(s) do not exist!")
        exit(1)
    fastGT_main_info = parse_main_txt(input_main_f, inferred_pop, surr_pops_info_df)
    if isinstance(fastGT_main_info, str) and fastGT_main_info == "":
        print("Unclear or uncertain historical admixture!")
        exit(0)
    if isinstance(fastGT_main_info, str) and fastGT_main_info == ".":
        print("Surrogate pop is not found at '--surr-pop-info' annotation file!")
        exit(1)
    fastGT_main_info = pd.DataFrame.from_dict(fastGT_main_info) # convert a list of dicts to a dataframe
    fastGT_main_info_cols = fastGT_main_info.columns
    non_pop_cols = {'conclusion','N_generations_ago', 'event', 'PC'}
    group_names = {}
    for col in fastGT_main_info.columns:
        if col not in non_pop_cols: group_names[col] = ""
    group_names = sorted(list(group_names))
    fastGT_info = parse_boot_txt(input_boot_f, fastGT_main_info)
    fastGT_info['inferred_pop'] = inferred_pop
    new_ordered_cols = ['inferred_pop', 'conclusion','N_generations_ago', 'event', 'PC', 'boot_gen_low95_CI', \
                       'boot_gen_high95_CI'] + group_names
    fastGT_info = fastGT_info[new_ordered_cols]
    output_f = f"{output_prefix}.csv"
    if fexists(output_f):
        output_df = pd.read_csv(output_f, header=0)
        head_cols = list(output_df.columns)
        # check if the columns of $fastGT_info_rows match the already-existent output CSV file
        if head_cols == list(fastGT_info.columns):
            fastGT_info.to_csv(output_f, mode="a", header=False, index=False)
        else: # if unmatched, print the columns to the output file as well
            fastGT_info.to_csv(output_f, mode="a", header=True, index=False)
    else:
        fastGT_info.to_csv(output_f, index=False)
