#!/usr/bin/python3

from ROOT import TFile, TEfficiency, TGraphErrors, TCanvas
import numpy as np
from math import floor
from scipy.stats import sem


def ImportROOTFile(input_name):
    root_file = TFile("data/raw/" + input_name)
    n_particles = int(str(root_file.config.Get("Allpix").Get("number_of_events")))
    hit_tree = root_file.PixelCharge
    i = 0
    hit_dict = dict()
    for event in hit_tree:
        if (i+2) % 10 == 1:
            print("Processing events:", round((i+1)/n_particles*100,1),"%",  end="\r")
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


def IntegrateCharge(modules_file_name, q_low = 0, q_high = 50):
    input_file = TFile("data/raw/" + modules_file_name)
    cluster_charge = input_file.Get("DetectorHistogrammer").Get("dut").Get("cluster_charge")
    bin_low = cluster_charge.FindBin(q_low)
    bin_high = cluster_charge.FindBin(q_high)

    total_charge = sum([cluster_charge.GetBinContent(bin)*cluster_charge.GetBinCenter(bin) for bin in range(bin_low, bin_high+1)])
    
    print(modules_file_name, ": ", round(total_charge,2), "ke")


def DrawCharge(modules_file_names):
    file_paths = ["data/raw/" + file_name for file_name in modules_file_names]
    input_files = [TFile(file_path, "r") for file_path in file_paths]
    histograms = [input_file.Get("DetectorHistogrammer").Get("dut").Get("cluster_charge") for input_file in input_files]
    canvas = TCanvas("canvas", "canvas", 600, 600)
    for i in range(len(histograms)):
        histograms[i].Rebin(20)
        if i == 0:
            histograms[i].Draw()
            histograms[i].GetXaxis().SetRangeUser(0, 70)
            histograms[i].GetXaxis().SetTitleOffset(1.0)
            histograms[i].GetYaxis().SetTitleOffset(1.0)
        else:
            histograms[i].Draw("same")
        histograms[i].SetLineColor(i+1)    
    canvas.SaveAs("plot.pdf")


def RunAnalysis(input_name, output_name=""):
    if not output_name: 
        output_name = input_name.split("_")[0] + "_analysed.root"
    print("INPUT:", input_name, "\nOUTPUT:", output_name)
    
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
        print("Threshold scanning:", round(thr/thr_end*100, 1), "%", end="\r")

        cluster_list = ScanThresholdDict(hit_dict, thrE)
        
        for cluster in cluster_list:
            eff.Fill(bool(cluster), thr)      

        cluster_list = [np.nan if cluster == 0 else cluster for cluster in cluster_list]
        
        try:
            cluster_size = np.nanmean(cluster_list)
            clus_graph.SetPoint(int(np.where(thr_range == thr)[0][0]), thr, cluster_size)  
            clus_graph.SetPointError(int(np.where(thr_range == thr)[0][0]), ex=0, ey=np.nanstd(cluster_list))
        except RuntimeWarning:
            pass        
            
    print("\nDone.")
    eff.Write()
    clus_graph.Write()
    # write_file.Write()
    write_file.Close()


# ImportROOTFile("0deg-WF-EF_output.root")
# RunAnalysis("0deg-lin_output.root", "test_analysed.root")
# RunAnalysis("0deg-EF_output.root")
# RunAnalysis("0deg-WF-EF_output.root", "0deg-WF-EF_analysed.root")
# RunAnalysis("0deg-EF-CTint_output.root", "0deg-EF-CTint_analysed.root")

IntegrateCharge("0deg-EF_modules.root")
IntegrateCharge("0deg-WF-EF_modules.root")

# DrawCharge(["0deg-EF_modules.root", "0deg-WF-EF_modules.root"])