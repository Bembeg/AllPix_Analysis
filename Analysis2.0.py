#!/usr/bin/python3

from ROOT import TFile, TEfficiency, TGraphErrors, TCanvas, TF1, TString, TDirectory
import numpy as np
from math import floor, sqrt
from scipy.stats import sem


def GetHitDict(input_name):
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


def ScanThreshold(hit_dict, threshold):
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


def DrawEfficiencyCluster(file_names, output_name, ref_file=""):
    """A function to plot efficiency and cluster size.
    
    Parameters
    ----------
    file_names : list
        List with file names
    output_name : str
        File name of the produced output file
    ref_file : str
        Optional name of the file with reference results
    """

    file_paths = ["data/" + file_name for file_name in file_names]
    input_files = [TFile(file_path, "r") for file_path in file_paths]

    eff = [input_file.Get("Efficiency") for input_file in input_files]
    clus = [input_file.Get("Average_cluster_size") for input_file in input_files]
    titles = [input_file.Get("Info").Get("title") for input_file in input_files]

    canvas = TCanvas("canvas1", "canvas1", 800, 600)
    for i in range(len(eff)):
        if i == 0:
            eff[i].Draw()
        else:
            eff[i].Draw("same")
        eff[i].SetLineColor(i+1)
        eff[i].GetListOfFunctions().FindObject("Efficiency_fit").SetLineColor(i+1)
    if ref_file:
        eff_ref = TFile("data/" + ref_file).Get("efficiency_vs_threshold_time_corrected")
        eff_ref.SetLineColor(51)
        eff_ref.SetMarkerColor(51)
        eff_ref.GetListOfFunctions().FindObject("erfcFit_timing").SetLineColor(51)
        eff_ref.Draw("PEX0same")

    canvas.SaveAs("results/" + output_name + "_eff.pdf")

    canvas = TCanvas("canvas2", "canvas2", 800, 600)
    for i in range(len(clus)):
        if i == 0:
            clus[i].Draw()
        else:
            clus[i].Draw("same")
        clus[i].SetLineColor(i+1)
    if ref_file:
        clus_ref = TFile("data/" + ref_file).Get("cluster_size_vs_threshold")
        clus_ref.SetLineColor(51)
        clus_ref.Draw("same")
    canvas.SaveAs("results/" + output_name + "_clus.pdf")


def RunAnalysis(input_name, output_name=""):
    # Check output name, set by default if not passed to the function
    if not output_name: 
        output_name = input_name.split("_")[0] + "_analysed.root"
    print("INPUT:", input_name, "\nOUTPUT:", output_name)
    
    # Configuring thresholds
    (thr_start, thr_end, thr_step) = (0.3, 8, 0.1)
    thr_range = np.arange(thr_start, thr_end+thr_step, thr_step) 
    n_thr = len(thr_range)

    # Get hit dictionary
    hit_dict = GetHitDict(input_name)
    
    # Open root file to write the results into
    write_file = TFile("data/" + output_name, "recreate") 
    write_file.cd()

    # Initialize efficiency and cluster size objects
    eff = TEfficiency("Efficiency", "Efficiency;Threshold [fC];Efficiency", n_thr, thr_start, thr_end)
    clus_graph = TGraphErrors(n_thr)
    clus_graph.SetNameTitle("Average_cluster_size", "Average_cluster_size")
    clus_graph.GetXaxis().SetTitle("Threshold [fC]")
    clus_graph.GetYaxis().SetTitle("Average cluster size") 

    # Perform threshold scanning
    for thr in thr_range:
        thrE = thr * 6242.2
        print("Threshold scanning:", round(thr/thr_end*100, 1), "%", end="\r")

        # Scan a threshold and obtain a list of clusters
        cluster_list = ScanThreshold(hit_dict, thrE)
        
        # Analyze cluster list and fill TEfficiency object
        for cluster in cluster_list:
            eff.Fill(bool(cluster), thr)      

        # Replace zeros with NaN to obtain proper average later
        cluster_list = [np.nan if cluster == 0 else cluster for cluster in cluster_list]
        
        try:
            # Calculate average cluster size and add points to clus_graph
            cluster_size = np.nanmean(cluster_list)
            N = np.count_nonzero(~np.isnan(cluster_list))
            err = np.nanstd(cluster_list)/sqrt(N-1)
            clus_graph.SetPoint(int(np.where(thr_range == thr)[0][0]), thr, cluster_size)  
            clus_graph.SetPointError(int(np.where(thr_range == thr)[0][0]), ex=0, ey=err)
        except RuntimeWarning:
            pass        
            
    print("\nDone.")

    # Efficiency fit
    fit_form = "0.5*[0]*TMath::Erfc((x-[1])/(TMath::Sqrt(2)*[2])*(1-0.6*TMath::TanH([3]*(x-[1])/TMath::Sqrt(2)*[2])))"
    fit_func = TF1("Efficiency_fit", fit_form, 0, 8)
    # Set initial parameter values and limits
    fit_func.SetParameters(1, 4, 1, 1, 0.5)
    fit_func.SetParLimits(1, 0, 5)
    fit_func.SetParLimits(2, 0, 2)
    fit_func.SetParLimits(3, 0, 2)
    fit_func.SetParLimits(4, 0.5, 0.7)
    # Perform fit
    eff.Fit(fit_func, "R")
    
    # Write outputs to file
    eff.Write()
    clus_graph.Write()

    # Collect info
    source = "allpix"
    angle = input_name.split("-")[0]
    descr = input_name.strip("_output.root").strip("0deg-")
    n_events = str(int(eff.GetTotalHistogram().GetEntries() / n_thr))
    title = source + "," + angle + "," + descr + "(" + str(n_events) + "ev)"
    vt50 = str(fit_func.GetParameter(1))
    vt50_err = str(fit_func.GetParError(1))
    thr_range = str(thr_start) + ":" + str(thr_end) + ":" + str(thr_step)
    
    # Write info to Info directory
    info_dir = write_file.mkdir("Info")
    info_dir.cd()
    info_dir.WriteObject(TString(source), "source")
    info_dir.WriteObject(TString(angle), "angle")
    info_dir.WriteObject(TString(descr), "descr")
    info_dir.WriteObject(TString(n_events), "n_events")
    info_dir.WriteObject(TString(title), "title")
    info_dir.WriteObject(TString(vt50), "vt50")
    info_dir.WriteObject(TString(vt50), "vt50_err")
    info_dir.WriteObject(TString(thr_range), "thr_range")
    write_file.Close()


# RunAnalysis("0deg-EF_output.root", "test.root")
# RunAnalysis("0deg-EF_output.root")
# RunAnalysis("0deg-lin_output.root")
# RunAnalysis("0deg-WF-EF_output.root")
# RunAnalysis("0deg-EF-CTint_output.root")
# RunAnalysis("0deg-histat_output.root")

# IntegrateCharge("0deg-EF_modules.root")
# IntegrateCharge("0deg-WF-EF_modules.root")

# DrawCharge(["0deg-EF_modules.root", "0deg-WF-EF_modules.root"])

DrawEfficiencyCluster(["0deg-WF-EF_analysed.root", "ser005.root"], output_name="test", ref_file="ref-0deg-testbeam.root")
