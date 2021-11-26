#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:57:14 2021

@author: sujan
"""
import json
import sys
import os
import copy
import utils.shared_utils as shut

sys.path.append(os.getcwd())
def add_qc_vars(var, qc_var_dic_in, qc):
    qc_var_dic = qc_var_dic_in.copy()
    qc_var_dic.pop("QC")
    qc_var_dic["bounds"] = [0, 1]
    qc_var_dic["sourceVariableName"] = var+"_QC"
    qc_var_dic["source2sindbadUnit"] = 1
    qc_var_dic["variableUnit"] = "adimensional"
    qc_var_dic["sourceVariableUnit"] = "adimensional"
    qc_var_dic["isEnergy"] = False
    qc_var_dic["nameLong"] = "quality flag for " + qc_var_dic_in["nameLong"]
    if qc == "QC":
        qc_var_dic["sourceVariableName"] = var+"_QC"
    else:
        qc_var_dic["sourceVariableName"] = qc
    qc_var_dic["nameShort"] = var+"_QC"
    if len(qc_var_dic_in["weightVar"]) > 0:
        qc_var_dic["weightVar"]= qc_var_dic_in["weightVar"]

    return qc_var_dic

def get_cliff_dirname(fnversion, temp_reso):
    if fnversion == 'FLUXNET2015':
        opath = 'BRK15'
    elif fnversion == 'LaThuile':
        opath = 'LTL07'
    else:
        opath=''

    if temp_reso == 'daily':
        otime = 'DD'
    else:
        otime = 'HR'
    return f'fluxnetBGI2021.{opath}.{otime}'

def setup_out_dir(exp_config):
    out_main_dir = exp_config['OutPath']
    out_sub_path = get_cliff_dirname(exp_config["FLUXNET_version"], exp_config["temporal_resolution"])

    exp_config['OutPath'] = {}
    exp_config['OutPath']['main'] =  os.path.join(out_main_dir, out_sub_path)
    opath = out_sub_path.split('.')[-2]
    otime = out_sub_path.split('.')[-1]
    exp_config['OutPath']['nc_file'] = os.path.join(exp_config['OutPath']['main'], 'data')
    os.makedirs(exp_config['OutPath']['nc_file'], exist_ok=True)
    # figures
    if exp_config['diagnostic_plots']:
        exp_config['OutPath']['figs'] = os.path.join(exp_config['OutPath']['main'], 'figs_diagno')
        os.makedirs(exp_config['OutPath']['figs'], exist_ok=True)
    # configuration output
    exp_config['OutPath']['expName'] = f'fluxnetBGI2021.{opath}.{otime}'
    exp_config['OutPath']['info_config'] = os.path.join(exp_config['OutPath']['main'], 'config_info')
    exp_config['OutPath']['log'] = os.path.join(exp_config['OutPath']['main'], 'logs')
    os.makedirs(exp_config['OutPath']['log'], exist_ok=True)
    os.makedirs(exp_config['OutPath']['info_config'], exist_ok=True)
    return exp_config

def merge_var_info(base_info, sel_var_info):
    out_info = {}
    for b_field in list(base_info.keys()):
        out_info[b_field] = base_info[b_field]
        if b_field in sel_var_info:
            out_info[b_field] = sel_var_info[b_field]
    return out_info

def get_selected_list(sel, full):
    olist = full
    if isinstance(sel, list):
        if "all" in sel:
           return olist
        else:
            return sel
    else:
        if sel == "all":
            return olist
        else:
            sel_list = [x.strip() for x in sel.split(" ")]
            if "except" in sel_list:
                remove_list = [x.strip() for x in sel.split("except")[-1].strip().split(" ")]
                return(list(set(olist)-set(remove_list)))
            elif "all" in sel_list:
                return olist
            else:
                return sel_list


def add_suffix_to_sel_variables(var_list, suffix):
    if len(suffix) > 0:
        return [_l+suffix for _l in var_list]
    else:
        return var_list

def get_variable_info(_config):
    _config_out = _config.copy()
    variables_all = {}
    for data_key in _config["sel_datasets"]:
        data_info = _config['dataset'][data_key]
        variables = {}
        var_sfx = data_info["var_suffix"].strip()
        if len(var_sfx) > 0:
            var_sfx = "_" + var_sfx
        var_info = data_info["variables"]
        sel_vars = get_selected_list(data_info["sel_variables"], list(data_info["variables"].keys()))
        data_info['sel_variables'] = sel_vars
        for var_base in data_info["sel_variables"]:
            var = var_base + var_sfx
            merged_var_info = merge_var_info(data_info["base_info"], var_info[var_base])
            variables[var] = merged_var_info.copy()
            # if len(data_info["cube_data_path"].strip()) > 0:
            #     data_path = data_info["cube_data_path"].strip()
            #     variables[var]["data_path"]= data_path
            qc_vars = variables[var]["QC"].copy()

            if _config["temporal_resolution"] == "daily":
                variables[var]["nameLong"] = var_info[var_base]["nameLong"]
            else:
                variables[var]["nameLong"] = var_info[var_base]["nameLong"]

            for QC in qc_vars:
                var_field_name = var_base + "_QC" + var_sfx
                variables[var_field_name] = add_qc_vars(var_base, merged_var_info.copy(), QC)
                variables[var]['QC'] = variables[var_field_name]["sourceVariableName"] 
            if _config["temporal_resolution"] == "daily":
                field_daystats = {"DayMin": "daily minimum",
                "DayMax": "daily maximum", "DayTime": "daytime mean", "DaySum": "daily sum", "wDayMean":"weighted daily mean"}
                for fi_da, fi_prx in field_daystats.items():
                    if fi_da in var_info[var_base] and var_info[var_base][fi_da]:
                        var_field_name = var_base + "_" + fi_da + var_sfx
                        variables[var_field_name] = merged_var_info.copy()
                        if fi_da == "wDayMean":
                            if len(variables[var_field_name]["weightVar"]) == 0:
                                variables[var_field_name]["weightVar"] = var_base+"_QC"
                            variables[var_field_name]["nameLong"] = f'{fi_prx} {var_info[var_base]["nameLong"]} [using {variables[var_field_name]["weightVar"]}]'
                        else:
                            variables[var_field_name]["nameLong"] = f'{fi_prx} {var_info[var_base]["nameLong"]}'

                        # variables[var_field_name] = merged_var_info
                        for QC in qc_vars:
                            var_field_name_qc = var_base + "_QC_" + fi_da +  var_sfx
                            variables[var_field_name_qc] = add_qc_vars(var_base, variables[var_field_name], QC)
                            variables[var_field_name]['QC'] = variables[var_field_name_qc]["sourceVariableName"]
                            
            variables_all.update(variables)
        data_info["sel_variables"] = add_suffix_to_sel_variables(sel_vars, var_sfx)

        for var, var_det in variables.items():
            for unwant in "wDayMean DayMin DayMax DaySum DayTime".split():
                if unwant in var_det:
                    variables[var].pop(unwant)
        _config_out['dataset'][data_key].pop('base_info')
        _config_out["dataset"][data_key]['variables'] = variables
    return variables_all, _config_out

def get_prov_configuration(conf_file):
    with open(conf_file, 'r') as config_file:
        prov_rf=config_file.read()
    prov_conf = json.loads(prov_rf)
    return prov_conf

def check_units_consistency(full_info):
    for name in full_info['sel_datasets']:
        dataset = full_info['dataset'][name]
        print(f'dataset: {name} | variables: {dataset["sel_variables"]}')
        for var in dataset['sel_variables']:
            # print(var, dataset['variables'])
            info = dataset['variables'][var]
            src_unit = info["sourceVariableUnit"]
            tar_unit = info["variableUnit"]
            unit_scalar = info["source2sindbadUnit"]
            check_unit = info["skipUnitcheck"]
            if check_unit == False:
                if src_unit != tar_unit:
                    if float(unit_scalar) == 1:
                        sys.exit(f'dataset: {name}: {var} is in {src_unit} in source [sourceVariableUnit] and {tar_unit} in target [variableUnit], but the unit scalar [source2sindbadUnit] is {unit_scalar}. Fix inconsistency setting the units and coversion correctly or using skipUnitcheck: true')

def get_inp_config(exp_config_path):
    with open(exp_config_path, 'r') as config_file:
        exp_config_rf=config_file.read()
    exp_config = json.loads(exp_config_rf)
    return exp_config

def get_exp_configuration(exp_config_path, fn_version=None, temporal_resolution=None):
    # get input configuration
    exp_config_r = get_inp_config(exp_config_path)

    # replace version and time scale configs if passed
    if fn_version is not None:
        exp_config_r["FLUXNET_version"] = fn_version
    if temporal_resolution is not None:
        exp_config_r["temporal_resolution"] = temporal_resolution
    exp_config = copy.deepcopy(exp_config_r)

    # setup directory
    exp_config = setup_out_dir(exp_config)
    exp_config['sel_datasets'] = get_selected_list(exp_config['sel_datasets'], list(exp_config['dataset'].keys()))

    # get provider information for each of the selected datasets
    for prov in exp_config['sel_datasets']:
        prov_dict = get_prov_configuration(exp_config['dataset'][prov])
        exp_config['dataset'][prov] = prov_dict

    # get variable information for each of the selected provider
    variables, exp_config = get_variable_info(exp_config)

    # check if the units are consistent with scalars
    check_units_consistency(exp_config)

    # save configuration
    shut.save_json(os.path.join(exp_config['OutPath']['info_config'], "variable_info.json"), variables)
    shut.save_json(os.path.join(exp_config['OutPath']['info_config'], "exp_info.json"), exp_config)

    return exp_config

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Utilities: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Config generator: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)

