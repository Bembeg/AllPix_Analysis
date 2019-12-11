#!/usr/bin/python3

from ROOT import TFile, TCanvas, TH1D, gStyle, TBrowser, TLegend, TMath, TF1, TGraph, Double, TLine, TLatex

def InterpolateHist(hist, x):
    """
    Interpolates between two bins of a histogram "hist", returns x,y position of left and right "point" and the interpolated y value to the corresponding "x" input.
    """
    xBin = hist.FindBin(x)
    binWidth = hist.GetBinWidth(hist.GetBin(xBin))
    # Looking for left bin
    for i in range(10):
        currBin = hist.FindBin(x-i*binWidth) 
        if hist.GetBinContent(currBin) > 0.01:
            if i == 0:
                return (0,0,0,0,hist.GetBinContent(currBin))
            xL = hist.GetBinCenter(currBin)
            yL = hist.GetBinContent(currBin)
            break 

    # Looking for right bin
    for i in range(10):
        currBin = hist.FindBin(x+i*binWidth) 
        if hist.GetBinContent(currBin) > 0.01:
            xR = hist.GetBinCenter(currBin)
            yR = hist.GetBinContent(currBin)
            break
                       
    k = (yR-yL)/(xR-xL)
    q = yR - k*xR       
    y = k*x + q
    return (xL,yL,xR,yR,y)

def PlotEfficiency (fileNames=0, legendEntries=0, refFileNames=[], refLegendEntries=[], plotName=0, legendHeader=0, refOption=0, axisTitleX="Threshold [fC]", axisTitleY="Efficiency"):
    """
    Pass a list of root files to have the efficiency plotted along with preferred legend entries for these plots (if not provided, they will be assumed from the file names). Reference root file (eg. with testbeam data) can be also passed along with the appropriate legend entry (or else assumed from the file name) and will be plotted as well.
    """
    canvas = TCanvas("c1", "c1", 500,500)
    gStyle.SetOptStat(0)            #hides stat table
    gStyle.SetOptTitle(0)           #hides title

    fitForm = "0.5*[0]*TMath::Erfc((x-[1])/(TMath::Sqrt(2)*[2])*(1-0.6*TMath::TanH([3]*(x-[1])/TMath::Sqrt(2)*[2])))"

    color = [1,2,4,6,9,9,6,4,2,1]
    lineStyle = [2,2,2,2,2,1,1,1,1,1]   
    markerStyle = [21,22,23,33,34,28,27,32,26,25]
      
    axisRangeXLow = 1
    axisRangeXHigh = 6
    axisRangeYLow = 0
    axisRangeYHigh = 1.02

    markerSize = 0.8
    textSize = 0.030
    effHist = []
    effTitle = []
    rootFile = []
    fitFunction = []
    refFile = []
    effRefPoints = []
    effRefFunc = []
      
    # Preparation and checks
    if fileNames == 0: 
        print("No file names were passed.")
        return 1
        
    if legendEntries == 0:
        print("No legend entries passed, assuming entries from file names.")
        return 1
        
    if refFileNames != 0 and refLegendEntries == 0:
        print("No reference legend entry passed, assuming from file name.")
        return 1
        
    if refFileNames != 0 and refOption == 0:
        print("No reference option selected (choose \"time\", \"center\") or \"edge\".")    
        return 1   

    if legendHeader == 0:
        print("No legend header provided.")
        return 1

    if plotName == 0:
        print("No plot title provided.")
        return 1

    if len(refFileNames) != len(refLegendEntries):
        print("Number of reference files and reference legend entries different.")
        return 1

    if len(fileNames) != len(legendEntries):
        print("Number of files and legend entries different.")
        return 1

    print ("PLOTTING", len(fileNames), "FILES ", end="")
    if refFileNames != 0: print("AND REFERENCE", end="")
    print("\n")

    # Efficiency histograms and fitting
    for i in range(len(fileNames)):
        rootFile.append(TFile("data/" + fileNames[i]))
        effTitle.append("Efficiency - " + fileNames[i].split("_")[0])
        effHist.append(rootFile[i].Get(effTitle[i]))
        
        effHist[i].SetMarkerSize(markerSize)
        effHist[i].SetMarkerColor(color[i])
        effHist[i].SetMarkerStyle(markerStyle[i])
        effHist[i].SetLineColor(color[i])
        
        if i == 0:              
            effHist[i].SetTitle("Efficiency - " + plotName)
            effHist[i].GetXaxis().SetTitle(axisTitleX)
            effHist[i].GetYaxis().SetTitle(axisTitleY)
            effHist[i].SetAxisRange(axisRangeYLow, axisRangeYHigh, "Y")
            effHist[i].SetAxisRange(axisRangeXLow, axisRangeXHigh, "X")
            effHist[i].Draw()            #() for points, ("AXIS") for just the axes
        else:                   
            effHist[i].Draw("same")

        fitFunction.append(TF1("fitFunc"+str(i), fitForm, 0, 8))
        fitFunction[i].SetParameters(1, 4, 1, 1, 0.5)
        fitFunction[i].SetParLimits(4, 0.5, 0.7)
        # fitFunction[i].FixParameter(0,1)
        fitFunction[i].SetParLimits(1, 0, 5)
        fitFunction[i].SetParLimits(2, 0, 2)
        fitFunction[i].SetParLimits(3, 0, 2)
        fitFunction[i].SetLineColor(color[i])
        fitFunction[i].SetLineStyle(lineStyle[i])
        fitFunction[i].SetMarkerSize(markerSize)
        fitFunction[i].SetMarkerStyle(markerStyle[i])
        fitFunction[i].SetMarkerColor(color[i])
        print("Fitting", fileNames[i])
        effHist[i].Fit("fitFunc"+str(i), "R")
        fitFunction[i].Draw("same")

        # Reference data plot
        if refFileNames != 0:
            if refOption == "time":
                histName = "efficiency_vs_threshold_time_corrected"
                functionName = "erfcFit_timing"
            if refOption == "center":
                histName = "efficiency_vs_threshold_strip_ctr"
                functionName = "erfcFit_ctr"
            if refOption == "edge":
                histName = "efficiency_vs_threshold_strip_edge"
                functionName = "erfcFit_edge"
            
            chisquares = []
            lines = []
            interpLines = []
            for i in range(len(refFileNames)):                  
                refFile.append(TFile ("data/" + refFileNames[i]))
                effRefFunc.append(refFile[i].Get(histName).GetListOfFunctions().FindObject(functionName))
                effRefPoints.append(refFile[i].Get(histName))
                refFile[i].Close()
        
                effRefPoints[i].SetMarkerSize(markerSize)
                effRefPoints[i].SetMarkerStyle(markerStyle[-i-1])
                effRefPoints[i].SetMarkerColor(color[-i-1])
                effRefPoints[i].SetLineColor(color[-i-1])
                effRefPoints[i].SetLineStyle(lineStyle[-i-1])
                effRefPoints[i].SetLineWidth(2)
                effRefPoints[i].Draw("Psame")
                effRefFunc[i].SetLineColor(color[-i-1])
                effRefFunc[i].SetLineStyle(lineStyle[-i-1])
                effRefFunc[i].Draw("Lsame")
                        
                # get chi2 of Data-MC
                doChi2 = 1
                if doChi2 == 1:
                    chisquares.append(0)
                    lines.append([])
                    interpLines.append([])
                    # effRefPoints[i].Print()
                    for j in range(effRefPoints[i].GetN()):
                        (x,y,ey) = (effRefPoints[i].GetX()[j],effRefPoints[i].GetY()[j],effRefPoints[i].GetEY()[j])      
                        (xL, yL, xR, yR, fy) = InterpolateHist(effHist[i], x)  
                        #fy = fitFunction[i].Eval(x)
                        # interpLines[-1].append(TLine(xL, yL, xR, yR))
                        # interpLines[-1][-1].SetLineColor(2)
                        # interpLines[-1][-1].SetLineWidth(2)
                        # interpLines[-1][-1].Draw("same")
                        chisquares[-1] += (y-fy)**2/ey**2
                        # lines[-1].append(TLine(x,y,x,fy))
                        # lines[-1][-1].SetLineColor(2)
                        # lines[-1][-1].SetLineWidth(2)
                        # lines[-1][-1].Draw("same")
                    print("chi2/NDF =", chisquares[i], "/", effRefPoints[i].GetN())
                     
    # Legend
    legendWidth = 0.9 - 0.013*max([len(max(legendEntries, key=len))+5, len(legendHeader), len(refLegendEntries)+5])
    legendHeight = 0.9 - (len(fileNames)+len(refFileNames))*0.05-0.05
    legend = TLegend(legendWidth, legendHeight, 0.9, 0.9)
    legend.SetHeader(legendHeader,"C")
    legend.SetTextSize(textSize)
    for i in range(len(fileNames)):
        if refFileNames !=0 and i <= len(refFileNames) - 1: 
            legend.AddEntry(effRefPoints[i], refLegendEntries[i], "pl")
        legend.AddEntry(fitFunction[i], legendEntries[i], "pl")
    legend.Draw("same")

    #Print Chi2/NDF
    if doChi2 == 1:     
        textHeader = TLatex(4.8, legendHeight-0.12, " #bf{#chi^{2}/ndf_{(data - sim)}}")
        textHeader.SetTextSize(textSize)
        
        textHeader.Draw("same")
        texts = []
        for i in range(len(chisquares)):
            texts.append( TLatex(4.8, legendHeight-(i+1)*0.05-0.13, "#bf{" + str(round(chisquares[i],0)) + "/" + str(effRefPoints[i].GetN()) + " = "+ str(round(chisquares[i]/effRefPoints[i].GetN(),1)) + "}" ))
            texts[i].SetTextSize(textSize)
            texts[i].SetTextColor(color[i])
            texts[i].Draw("same")

            # draw box around chi2/ndf text
            # lineLeft = TLine(5.22, legendHeight-(len(chisquares)+1)*0.05-0.09, 5.22, legendHeight-0.06)
            # lineLeft.Draw("same")
            # lineTop = TLine(5.22, legendHeight-0.06, 7.04, legendHeight-0.06)
            # lineTop.Draw("same")
            # lineBottom = TLine(5.22, legendHeight-(len(chisquares)+2)*0.05-0.04, 7.04, legendHeight-(len(chisquares)+2)*0.05-0.04)
            # lineBottom.Draw("same")
            # lineHeader = TLine(5.22, legendHeight-0.125, 7.04, legendHeight-0.125)
            # lineHeader.Draw("same")

    # Print and save
    canvas.SaveAs("results/" + plotName + "_eff.pdf")


def PlotClusterSize (fileNames=0, legendEntries=0, refFileNames=0, refLegendEntries=0, plotName=0, legendHeader=0, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size"):
    """
    Pass a list of root files to have the cluster size plotted along with preferred legend entries for these plots (if not provided, they will be assumed from the file names). Reference root file (eg. with testbeam data) can be also passed along with the appropriate legend entry (or else assumed from the file name) and will be plotted as well.
    """
    canvas = TCanvas("c1", "c1", 600,600)
    gStyle.SetOptStat(0)            #hides stat table
    gStyle.SetOptTitle(0)           #hides title

    color = [1,2,4,6,9,9,6,4,2,1]
    lineStyle = [2,2,2,2,2,1,1,1,1,1]   
    markerStyle = [21,22,23,33,34,28,27,32,26,25]
    markerSize = 0.6
    textSize = 0.030

    axisRangeXLow = 0
    axisRangeXHigh = 6
    axisRangeYLow = 1
    axisRangeYHigh = 1.7
    
    clusTitle = []
    rootFile = []
    clusRef = []
    clusHist = []
    refFile = []

    # Preparation and checks
    if fileNames == 0: 
        print("No file names were passed.")
        return 1
    
    if legendEntries == 0:
        print("No legend entries passed, assuming entries from file names.")
        legendEntries = []
        for fileName in fileNames:
            legendEntries.append(fileName.split("_")[0])
    
    if refFileNames != 0 and refLegendEntries == 0:
        print("No reference legend entry passed, assuming from file name.")
        refLegendEntries = refFileName.split(".")[0].strip("ref_")

    if legendHeader == 0:
        print("No legend header provided.")
        return 1

    if plotName == 0:
        print("No plot title provided, setting to \"Cluster_size_plot_00.pdf\".")
        plotName = "Cluster_size_plot_00.pdf"

    if len(fileNames) != len(legendEntries):
        print("Number of files and legend entries different.")
        return 1
    else: 
        nOfPlots = len(fileNames)
        print ("PLOTTING", nOfPlots, "FILES ", end="")
        if refFileNames != 0: print("AND REFERENCE", end="")
        print("\n")

    # Cluster size histograms
    for i in range(nOfPlots):
        rootFile.append(TFile("data/" + fileNames[i]))
        clusTitle.append("Cluster Size - " + fileNames[i].split("_")[0] + "_pfx")
        clusHist.append(rootFile[i].Get(clusTitle[i]))
        
        clusHist[i].SetMarkerSize(markerSize)
        clusHist[i].SetMarkerColor(color[i])
        clusHist[i].SetMarkerStyle(markerStyle[i])
        clusHist[i].SetLineColor(color[i])

        if i == 0:           
            clusHist[i].SetTitle("Cluster Size - " + plotName)
            clusHist[i].GetXaxis().SetTitle(axisTitleX)
            clusHist[i].GetYaxis().SetTitleOffset(1.4)
            clusHist[i].GetYaxis().SetTitle(axisTitleY)
            clusHist[i].SetAxisRange(axisRangeYLow, axisRangeYHigh, "Y")
            clusHist[i].SetAxisRange(axisRangeXLow, axisRangeXHigh, "X")
            clusHist[i].Draw("")
        else:              
            clusHist[i].Draw("same")

    # Reference data plot
    if refFileNames != 0:
        histName = "cluster_size_vs_threshold"
        for i in range(len(refFileNames)):
            refFile.append(TFile ("data/" + refFileNames[i]))
            clusRef.append(refFile[i].Get(histName))
            refFile[i].Close()
            clusRef[i].SetMarkerSize(markerSize)
            clusRef[i].SetMarkerStyle(markerStyle[-1-i])
            clusRef[i].SetMarkerColor(color[-1-i])
            clusRef[i].SetLineColor(color[-1-i])  
            clusRef[i].SetLineStyle(lineStyle[-1-i])
            clusRef[i].SetLineWidth(2)
            clusRef[i].Draw("samePLE")

    # Legend
    legendWidth = 0.9 - 0.013*max([len(max(legendEntries, key=len))+5, len(legendHeader), len(refLegendEntries)+5])
    legendHeight = 0.9 - (len(fileNames)+len(refFileNames))*0.05-0.05
    legend = TLegend(legendWidth, legendHeight, 0.9, 0.9)
    legend.SetHeader(legendHeader,"C")
    legend.SetTextSize(textSize)
    for i in range(len(fileNames)):
        if refFileNames !=0 and i <= len(refFileNames) - 1:
            legend.AddEntry(clusRef[i], refLegendEntries[i], "pl")
        legend.AddEntry(clusHist[i], legendEntries[i], "p")
            
    legend.Draw("same")

    # Print and save
    canvas.SaveAs("results/" + plotName + "_clus.pdf")


# ------------------------------------------------------------------------------

fileNames =  ["rot0deg-300um-noCT_analysed.root", "rot0deg-300um-CT_analysed.root", "0deg-290um-864e-CT_analysed.root"]
legendEntries = ["Sim - 300um", "AUW", "new Sim."]
refFileNames = ["ref-0deg-testbeam.root"]
refLegendEntries = ["Test beam data"]            
plotName = "diff"      
legendHeader = "Data points"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames =  ["rot0deg-300um-noCT_analysed.root"]
legendEntries = ["300um, CT 0 %"]
refFileNames = ["ref-0deg-testbeam.root"]
refLegendEntries = ["TB data, CT ?"]            
plotName = "rot0deg_crosstalk_1"      
legendHeader = "Crosstalk:"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames =  ["rot0deg-300um-noCT_analysed.root", "rot0deg-300um-CT_analysed.root", "rot0deg-300um-CThigh_analysed.root"]
legendEntries = ["Simulation", "Sim. + crosstalk", "Sim. + 2x crosstalk"]
refFileNames = ["ref-0deg-testbeam.root"]
refLegendEntries = ["Test beam data"]
plotName = "CT_over"      
legendHeader = "Crosstalk:"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames = ["0deg-290um-864e-CT_analysed.root", "y5deg-290um-864e-CT_analysed.root", "y12deg-290um-864e-CT_analysed.root"] 
legendEntries = ["0 deg - Simulation",  "5 deg - Simulation", "12 deg - Simulation"]
refFileNames =  ["ref-0deg-testbeam.root", "ref-5degy-testbeam.root", "ref-12degy-testbeam.root"]
refLegendEntries = ["0 deg - Test beam", "5 deg - Test beam", "12 deg - Test beam"]    
plotName = "chi_testing_Y2"
legendHeader = "Sensor rotation (y-axis):"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames = ["0deg-290um-864e-CT_analysed.root", "x23deg-290um-864e-CT_analysed.root"]
legendEntries = ["0 deg - Simulation",  "23 deg - Simulation"]
refFileNames =  ["ref-0deg-testbeam.root", "ref-23degx-testbeam.root"]
refLegendEntries = ["0 deg - Test beam", "23 deg - Test beam"]    
plotName = "rotXaxis"
legendHeader = "Sensor rotation (x-axis):"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames =  ["roty12deg-300um-noCT_analysed.root", "roty12deg-300um-CT_analysed.root"]
legendEntries = ["Simulation", "Sim. + crosstalk"]
refFileNames = ["ref-12degy-testbeam.root"]
refLegendEntries = ["Test beam data"]            
plotName = "roty12deg_crosstalk"      
legendHeader = "Crosstalk:"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames =  ["rotx23deg-300um-noCT_analysed.root", "rotx23deg-300um-CT_analysed.root"]
legendEntries = ["Simulation", "Sim. + crosstalk"]
refFileNames = ["ref-23degx-testbeam.root"]
refLegendEntries = ["Test beam data"]            
plotName = "rotx23deg_crosstalk"      
legendHeader = "Crosstalk:"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")


fileNames =  ["y5deg-290um-864e-CT_analysed.root"]
legendEntries = ["Simulation"]
refFileNames = ["ref-12degy-testbeam.root"]
refLegendEntries = ["Testbeam"]
plotName = "chi_testing_1"      
legendHeader = "Data points:"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency")
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, axisTitleX="Threshold [fC]", axisTitleY="Average cluster size")
