#!/usr/bin/python3

from ROOT import TFile, TCanvas, TH1D, gStyle, TBrowser, TLegend, TMath, TF1, TGraph, Double, TLine, TLatex, TPad, TGraphErrors
from datetime import datetime as date
from math import ceil, sqrt


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
                return (0,0,0,0,hist.GetBinContent(currBin), hist.GetBinError(currBin))
            xL = hist.GetBinCenter(currBin)
            yL = hist.GetBinContent(currBin)
            eyL = hist.GetBinError(currBin)
            break 

    # Looking for right bin
    for i in range(10):
        currBin = hist.FindBin(x+i*binWidth) 
        if hist.GetBinContent(currBin) > 0.01:
            xR = hist.GetBinCenter(currBin)
            yR = hist.GetBinContent(currBin)
            eyR = hist.GetBinError(currBin)
            break
                       
    k = (yR-yL)/(xR-xL)
    q = yR - k*xR       
    y = k*x + q
    ey = (eyL*(xR-x) + eyR*(x-xL))/(xR-xL)
    # print("x=", x, "xL=", xL, "xR=", xR, "eyL=", eyL, "eyR=", eyR, "ey=", ey)
    # print("y=",round(y,5),"ey=",round(ey,5))
    return (xL,yL,xR,yR,y,ey)


def PlotEfficiency (fileNames=0, legendEntries=0, refFileNames=[], refLegendEntries=[], plotName=0, legendHeader=0, refOption="time", axisTitleX="Threshold [fC]", axisTitleY="Efficiency", plotRatio=0):
    """
    Pass a list of root files to have the efficiency plotted along with preferred legend entries for these plots (if not provided, they will be assumed from the file names). Reference root file (eg. with testbeam data) can be also passed along with the appropriate legend entry (or else assumed from the file name) and will be plotted as well.
    """

    canvas = TCanvas("c1", "c1", 500,500)
    
    if plotRatio == 1: 
        mainPadYlow = 0.275
        ratioPadYhigh = 0.275
        mainPadMargin = 0.012
        ratioPadMargin = 0
        ratioPadBotMargin = 0.25
        mainToRatio = (1 - mainPadYlow + ratioPadMargin + mainPadMargin*(1-mainPadYlow+ratioPadBotMargin)) / (ratioPadYhigh)

        mainPad = TPad("mainPad", "mainPad", 0, mainPadYlow, 1, 1)
        mainPad.SetBottomMargin(mainPadMargin)
        mainPad.Draw()
         
        ratioPad = TPad("ratioPad","ratioPad", 0, 0, 1, ratioPadYhigh)
        ratioPad.SetTopMargin(ratioPadMargin)
        ratioPad.SetBottomMargin(ratioPadBotMargin)
        ratioPad.Draw()
        mainPad.cd()

    gStyle.SetOptStat(0)            #hides stat table
    gStyle.SetOptTitle(0)           #hides title

    fitForm = "0.5*[0]*TMath::Erfc((x-[1])/(TMath::Sqrt(2)*[2])*(1-0.6*TMath::TanH([3]*(x-[1])/TMath::Sqrt(2)*[2])))"

    color = [1,2,4,6,9,9,6,4,2,1]
    lineStyle = [2,2,2,2,2,1,1,1,1,1]   
    markerStyle = [21,22,23,33,34,28,27,32,26,25]
    
    axisRangeXLow = 0
    axisRangeXHigh = 6
    axisRangeYLow = 0
    axisRangeYHigh = 1.02

    doChi2 = 0

    markerSize = 0.6
    textSize = 0.030
    effHist = []
    effTitle = []
    rootFile = []
    fitFunction = []
    refFile = []
    effRefPoints = []
    effRefFunc = []
    effRatio = []
    
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

    if len(fileNames) == len(refFileNames):
        chi2Mode = "toEach"
    else:
        chi2Mode = "toOne"

    # if len(refFileNames) != 1 and plotRatio == 1:
    #     print("Plot ratio impossible with other than 1 reference file.")
    #     return 1

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
            if plotRatio == 0:
                effHist[i].GetXaxis().SetTitle(axisTitleX)
            if plotRatio == 1:
                effHist[i].GetXaxis().SetLabelSize(0)
            effHist[i].GetYaxis().SetTitle(axisTitleY)
            effHist[i].SetAxisRange(axisRangeYLow, axisRangeYHigh, "Y")
            effHist[i].SetAxisRange(axisRangeXLow, axisRangeXHigh, "X")
            effHist[i].GetXaxis().SetTitleOffset(1.1)
            effHist[i].Draw("X0")            #() for points, ("AXIS") for just the axes
        else:                   
            effHist[i].Draw("X0same")

        fitFunction.append(TF1("fitFunc"+str(i), fitForm, 0, 8))
        fitFunction[i].SetParameters(1, 4, 1, 1, 0.5)
        fitFunction[i].SetParLimits(4, 0.5, 0.7)
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
            effRefPoints[i].Draw("PE0same")
            effRefFunc[i].SetLineColor(color[-i-1])
            effRefFunc[i].SetLineStyle(lineStyle[-i-1])
            effRefFunc[i].Draw("Lsame")
                        
            # Get Chi2 of Data-MC for "toEach" mode
            if doChi2 == 1 and chi2Mode == "toEach":
                chisquares.append(0)
                lines.append([])
                interpLines.append([])
                
                if plotRatio == 1:
                    effRatio.append(TGraph(effRefPoints[i].GetN()))

                for j in range(effRefPoints[i].GetN()):
                    (x, y, ey) = (effRefPoints[i].GetX()[j], effRefPoints[i].GetY()[j], effRefPoints[i].GetEY()[j])      
                    (xL, yL, xR, yR, fy, fEy) = InterpolateHist(effHist[i], x)  
                    # fy = fitFunction[i].Eval(x)
                    # interpLines[-1].append(TLine(xL, yL, xR, yR))
                    # interpLines[-1][-1].SetLineColor(2)
                    # interpLines[-1][-1].SetLineWidth(2)
                    # interpLines[-1][-1].Draw("same")
                    chisquares[-1] += (y-fy)**2/(ey+fEy)**2
                    # lines[-1].append(TLine(x,y,x,fy))
                    # lines[-1][-1].SetLineColor(2)
                    # lines[-1][-1].SetLineWidth(2)
                    # lines[-1][-1].Draw("same")
                    if plotRatio == 1:
                        effRatio[-1].SetPoint(j, x, fy/y)
                        E = fy/y * sqrt((ey/y)**2+(fEy/fy)**2)
                        effRatio[-1].SetPointError(j, 0, E)
                # print("chi2/NDF =", chisquares[i], "/", effRefPoints[i].GetN())

    # Get Chi2 of Data-MC for "toOne" mode
    if doChi2 == 1 and chi2Mode == "toOne":       
        for i in range(len(effHist)):
            if plotRatio == 1:  effRatio.append(TGraphErrors(effRefPoints[0].GetN()))
            chisquares.append(0)
            for j in range(effRefPoints[0].GetN()):
                (x, y, ey) = (effRefPoints[0].GetX()[j], effRefPoints[0].GetY()[j], effRefPoints[0].GetEY()[j])
                (xL, yL, xR, yR, fy, fEy) = InterpolateHist(effHist[i], x) 
                chisquares[-1] += (y-fy)**2/(ey+fEy)**2
                if plotRatio == 1:
                    effRatio[-1].SetPoint(j, x, fy/y)
                    E = fy/y * sqrt((ey/y)**2+(fEy/fy)**2)
                    effRatio[-1].SetPointError(j, 0, E)

    # Legend
    if plotRatio == 1: legendWidthCoeff = 0.011
    else: legendWidthCoeff = 0.014
    legendWidth = 0.9 - legendWidthCoeff*max([len(max(legendEntries, key=len))+4, len(legendHeader)])#, len(max(refLegendEntries, key=len))+4])
    legendHeight = 0.9 - (len(fileNames)+len(refFileNames))*0.05#-0.05
    legend = TLegend(legendWidth, legendHeight, 0.9, 0.9)
    # legend.SetHeader(legendHeader,"C")
    legend.SetTextSize(textSize)
    for i in range(len(fileNames)):
        if refFileNames !=0 and i <= len(refFileNames) - 1: 
            legend.AddEntry(effRefPoints[i], refLegendEntries[i], "pl")
        legend.AddEntry(fitFunction[i], legendEntries[i], "pl")
    legend.Draw("same")

    # Print Chi2/NDF
    if doChi2 == 1 and chi2Mode == "toEach":     
        textHeader = TLatex(4.6, legendHeight-0.12, " #bf{#chi^{2}/ndf_{(data - sim)}}")
        textHeader.SetTextSize(textSize)
        textHeader.Draw("same")
        texts = []
        for i in range(len(chisquares)):
            texts.append( TLatex(4.6, legendHeight-(i+1)*0.05-0.13, "#bf{" + str(round(chisquares[i],0)) + "/" + str(effRefPoints[i].GetN()) + " = "+ str(round(chisquares[i]/effRefPoints[i].GetN(),1)) + "}" ))
            texts[i].SetTextSize(textSize)
            texts[i].SetTextColor(color[i])
            texts[i].Draw("same")
    
    elif doChi2 == 1 and chi2Mode == "toOne":
        textHeader = TLatex(4.8, legendHeight-0.12, " #bf{#chi^{2}/ndf_{(data - sim)}}")
        textHeader.SetTextSize(textSize)
        # print(chisquares)
        textHeader.Draw("same")
        texts = []
        for i in range(len(chisquares)):
            texts.append( TLatex(4.8, legendHeight-(i+1)*0.05-0.13, "#bf{       " + str(round(chisquares[i]/effRefPoints[0].GetN(),1)) + "}" ))
            texts[i].SetTextSize(textSize)
            texts[i].SetTextColor(color[i])
            texts[i].Draw("same")

    # Plot ratio graphs
    if plotRatio == 1:
        ratioPad.cd()
        # Axis center 1 for division, 0 for subtraction
        axisCenter = 1
        axisShift = ceil(10*max([max(abs(gr.GetHistogram().GetMinimum())-axisCenter, gr.GetHistogram().GetMaximum()-axisCenter) for gr in effRatio]))/10      
        
        for i in range(len(effRatio)):
            effRatio[i].SetLineColor(color[i])
            effRatio[i].SetMarkerStyle(20)
            effRatio[i].SetMarkerColor(color[i])
            effRatio[i].SetMarkerSize(markerSize*0.7)
            effRatio[i].SetMarkerStyle(markerStyle[i])
            effRatio[i].SetLineWidth(2)
            if i == 0:
                effRatio[i].GetXaxis().SetTitle(axisTitleX)
                effRatio[i].GetXaxis().SetTitleSize(textSize*mainToRatio)
                effRatio[i].GetXaxis().SetTitleOffset(1.2)
                effRatio[i].GetYaxis().SetRangeUser(axisCenter-axisShift, axisCenter+axisShift)
                effRatio[i].GetYaxis().SetTitle("Sim / Data")
                effRatio[i].GetYaxis().SetTitleSize(textSize*mainToRatio)
                effRatio[i].GetYaxis().SetTitleOffset(0.52)
                effRatio[i].GetYaxis().SetLabelSize(textSize*mainToRatio)
                effRatio[i].GetYaxis().SetNdivisions(6,5,1)
            
                effRatio[i].GetXaxis().SetLimits(axisRangeXLow, axisRangeXHigh + 0.04)
                effRatio[i].GetXaxis().SetLabelSize(textSize*mainToRatio)
                effRatio[i].GetXaxis().SetTickSize(0.08)
                effRatio[i].GetXaxis().SetLabelOffset(0.02)
                effRatio[i].Draw("ALPE")
            else:
                effRatio[i].Draw("sameLPE")
        oneLine = TLine(axisRangeXLow+0.1, axisCenter, axisRangeXHigh+0.04, axisCenter)
        oneLine.SetLineColor(13)
        oneLine.Draw("same")   


    # Print and save
    canvas.SaveAs("results/" + plotName + "_eff.pdf")
    logFile = open("log_fits.txt", "a")
    logFile.write("\nPLOT:\t\t\t\t" + plotName + "_eff.pdf")
    logFile.write("\nDATE:\t\t\t\t" + str(date.now()))
    for i in range(len(fitFunction)):    
        logFile.write("\n\tFCT:\t\t\t" + fileNames[i].split("_")[0])
        logFile.write("\n\t\tChi2:\t\t" + str(round(fitFunction[i].GetChisquare(),2)))
        logFile.write("\n\t\tNDF:\t\t" + str(fitFunction[i].GetNDF()))
        logFile.write("\n\t\tChi2/NDF:\t" + str(round(fitFunction[i].GetChisquare()/fitFunction[i].GetNDF(),2)))
        for j in range(fitFunction[i].GetNpar()):
            logFile.write("\n\t\tP" + str(j) + ":\t\t\t" + str(round(fitFunction[i].GetParameter(j),3)) + " +- " + str(round(fitFunction[i].GetParError(j),3)))
    logFile.write("\n-----\n")
    logFile.close()


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
        return 1

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
plotName = "toOne_div"      
legendHeader = "Data points"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, plotRatio=1)
# PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader)


fileNames =  ["0deg-athena_analysed.root"]
legendEntries = ["Athena"]
refFileNames = []
refLegendEntries = []            
plotName = "Athena"      
legendHeader = "Data points:"                        
PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, plotRatio=0)
PlotClusterSize(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader)


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
plotName = "toEach_div"
legendHeader = "Sensor rotation (y-axis):"                        
# PlotEfficiency(fileNames, legendEntries, refFileNames, refLegendEntries, plotName, legendHeader, refOption="time", plotRatio=1)
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

