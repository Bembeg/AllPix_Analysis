#!/usr/bin/python3

from ROOT import TFile, TEfficiency, TGraphErrors
import numpy as np
from math import floor
from scipy.stats import sem
from os import system
import gc
import cProfile
from guppy import hpy

def ImportROOTFile(input_name):
    root_file = TFile("data/raw/" + input_name)
    n_particles = int(str(root_file.config.Get("Allpix").Get("number_of_events")))

    hit_tree = root_file.PixelCharge
    i = 0
    hit_dict = dict()
    for event in hit_tree:
        if i/n_particles * 100 % 1 == 1:
            print("Processing events:", i/n_particles*100, end="\r")
        charges_dict = dict()
        for stripHit in event.dut:
            charges_dict[stripHit.getIndex().Y()] = stripHit.getCharge()
        hit_dict[i] = charges_dict
        i += 1
    print()
    root_file.Close()
    return hit_dict

def ScanThresholdDict(hit_dict, threshold):
    cluster_list = list()
    for event in hit_dict.values():
        hits = sum(charge >= threshold for charge in event.values())
        cluster_list.append(hits)
    
    return cluster_list
        

def ScanThresholdTree(hit_tree, threshold):
    cluster_list = list()
    for event in hit_tree:
        charges_dict = dict()
        for stripHit in event.dut:
            charges_dict[stripHit.getIndex().Y()] = stripHit.getCharge()
        hits = sum(charge >= threshold for charge in charges_dict.values())
        cluster_list.append(hits)
    return cluster_list 

def RunAnalysis(input_name, output_name=""):
    if not output_name: 
        output_name = input_name.split("_")[0] + "_analysed.root"
    print("INPUT:", input_name, "\nOUTPUT:", output_name)

    # root_file = TFile("data/raw/" + input_name)

    # hit_tree = root_file.PixelCharge
    # n_strips = int(str(root_file.models.Get("atlas17_dut").Get("number_of_pixels")).split(" ")[1])
    # n_particles = int(str(root_file.config.Get("Allpix").Get("number_of_events")))

    (thr_start, thr_end, thr_step) = (0, 8, 0.1)
    thr_range = np.arange(thr_start, thr_end+thr_step, thr_step) 
    n_thr = len(thr_range)

    hit_dict = ImportROOTFile(input_name)

    write_file = TFile("data/" + output_name, "recreate") 
    write_file.cd()

    eff = TEfficiency("Efficiency", "Efficiency;Threshold [fC];Efficiency", 80, 0, 8)
    clus_graph = TGraphErrors(n_thr)
    clus_graph.SetNameTitle("Average_cluster_size", "Average_cluster_size")
    clus_graph.GetXaxis().SetTitle("Threshold [fC]")
    clus_graph.GetYaxis().SetTitle("Average cluster size") 

    for thr in thr_range:
        thrE = thr * 6242.2
        print("Threshold scanning:", round(thr, 1), end="\r")

        cluster_list = ScanThresholdDict(hit_dict, thrE)
        
        for cluster in cluster_list:
            eff.Fill(bool(cluster), thr)
        clus_graph.SetPoint(int(np.where(thr_range == thr)[0][0]), thr, np.mean(cluster_list))  
        clus_graph.SetPointError(int(np.where(thr_range == thr)[0][0]), ex=0, ey=sem(cluster_list))
    print()
    eff.Write()
    clus_graph.Write()
    # write_file.Write()
    write_file.Close()


# ImportROOTFile("0deg-lin_output.root")
# RunAnalysis("0deg-lin_output.root", "test_analysed.root")
# RunAnalysis("0deg-EF_output.root", "0deg-EF_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-WF4-EF_output.root", "0deg-WF4-EF_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
RunAnalysis("0deg-WF-EF_output.root", "0deg-WF-EF_analysed.root")
# RunAnalysis("0deg-EF-CTint_output.root", "0deg-EF-CTint_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
