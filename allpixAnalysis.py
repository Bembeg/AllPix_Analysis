#!/usr/bin/python3

from ROOT import TFile, TH1D, TH2D, TMath, TF1, gStyle
import numpy as np
import sys
import os

def RunAnalysis(inputName, crosstalkSide, crosstalkBack):  
    # Check if input file exists
    if not os.path.exists(inputName):
        print("Input file doesn't exist.")
        return 1
    readFile = TFile(inputName, "read")

    # Check if output file already exists
    if os.path.exists("analysed.root"):
        if input("Output file already exists, overwrite? [Y/n] ") == "n":
            return 1
    writeFile = TFile("analysed.root", "recreate") 

    # Fetch PixelCharge objects, simulation parameters and close the input file
    hitTree = readFile.PixelCharge
    nOfStrips = int(str(readFile.models.Get("atlas17_dut").Get("number_of_pixels")).split(" ")[0])
    nOfEvents = int(str(readFile.config.Get("Allpix").Get("number_of_events")))

    # Set threshold scan parameters
    (thrStartFC, thrEndFC, thrStepFC) = (0.0, 8.0, 0.2)

    # Prepare histograms for efficiency and cluster size
    effHist = TH1D("Efficiency", "Efficiency", 200, thrStartFC, thrEndFC)
    clusHist = TH2D("Cluster Size", "Cluster Size", 200, thrStartFC, thrEndFC, 200, 0, 10)


    # Perform threshold scan 
    for thr in np.arange(thrStartFC, thrEndFC, thrStepFC):
        # Convert threshold from fC to e
        thrE = thr * 6242.2
        print "Current threshold [fC]:", round(thr, 1), "/", thrEndFC,

        for event in hitTree:  # Iterating over all events
            stripCharge = np.zeros(nOfStrips)
            for stripHit in event.dut:    # Iterating over all strips hit in that event
                stripCharge[stripHit.getIndex().X()] = stripHit.getCharge()     #Fill array with charges 

            stripChargeAdj = np.copy(stripCharge)
            # Check if cross talk enabled
            if crosstalkSide > 0 or crosstalkBack > 0:   
                stripMask = np.nonzero(stripCharge)    # Get strip indices with non-zero charge
                for stripIndex in stripMask:
                    # Calculate amount of charge to share by cross talk
                    chargeSide = stripCharge[stripIndex] * crosstalkSide
                    chargeBack = stripCharge[stripIndex] * crosstalkBack

                    # Transfer charges by cross talk
                    try:
                        stripChargeAdj[stripIndex] -= 2 * chargeSide + chargeBack
                        stripChargeAdj[stripIndex - 1] += chargeSide
                        stripChargeAdj[stripIndex + 1] += chargeSide
                    except IndexError:  # Handle list index going out of range 
                        print("INDEXERROR")
                        pass

            # Get cluster size as a number of strips with charge above threshold
            cluster = len(np.where(stripChargeAdj > thrE)[0])

            # If at least 1 hit, fill efficiency and cluster size histograms
            if cluster > 0:
                effHist.Fill(thr)
                clusHist.Fill(thr, cluster)

    effHist.Scale(1 / nOfEvents)
    clusHist = clusHist.ProfileX()
    print("Analysis done.                                         \n")

    # Aesthetic changes to the plots
    gStyle.SetOptStat("n")
    gStyle.SetOptFit(1)  
    gStyle.SetOptTitle(0)
    effHist.GetXaxis().SetTitle("Threshold [fC]")
    effHist.GetYaxis().SetTitle("Efficiency")
    effHist.SetMarkerStyle(20)
    effHist.SetMarkerColor(4)
    clusHist.GetXaxis().SetTitle("Threshold [fC]")
    clusHist.GetYaxis().SetTitle("Average cluster size")
    clusHist.SetMarkerStyle(20)
    clusHist.SetMarkerColor(4)

    # Fit efficiency with a skewed complementary error function
    fitForm = "0.5*[0]*TMath::Erfc((x-[1])/(TMath::Sqrt(2)*[2])*(1-0.6*TMath::TanH([3]*(x-[1])/TMath::Sqrt(2)*[2])))"
    fitFunction = TF1("fitFunc", fitForm, 0, 8)
    fitFunction.SetParameters(1, 4, 1, 1, 0.5)
    fitFunction.SetParLimits(4, 0.5, 0.7)
    fitFunction.SetParLimits(1, 0, 5)
    fitFunction.SetParLimits(2, 0, 2)
    fitFunction.SetParLimits(3, 0, 2)
    fitFunction.SetParNames("Max Efficiency","Median Charge","Width (sigma)","Skew")
    effHist.Fit("fitFunc", "quietR")

    # Write the histograms to a file and close.
    writeFile.Write()
    writeFile.Close()
    readFile.Close()

args = sys.argv[1:]
# Check validity of passed arguments
if len(args) == 0 or len(args) > 5 or len(args)%2 != 1:
    print("Invalid arguments.\nUsage: python3 allpixAnalysis.py INPUT_FILE [--side crosstalkSide] [--back crosstalkBack]")
else:
    # By default no cross talk
    crosstalkSide = 0.0
    crosstalkBack = 0.0

    # Get input name and cross talk values from passed args
    inputName = args[0]
    for i in range(len(args)):
        if args[i] == "--side":
            # Check that cross talk value is a number
            try:
                crosstalkSide = float(args[i+1])
            except ValueError:
                print("Cross talk value not a number, reverting to 0.")
                crosstalkSide = 0.0

            # Check that cross talk value is non-negative
            if crosstalkSide < 0 or crosstalkSide > 1:
                print("Invalid cross talk value (negative or >1), reverting to 0.")
                crosstalkSide = 0.0

        if args[i] == "--back":
            # Check that cross talk value is a number
            try:
                crosstalkBack = float(args[i+1])
            except ValueError:
                print("Cross talk value not a number, reverting to 0.")
                crosstalkBack = 0.0

            # Check that cross talk value is non-negative
            if crosstalkBack < 0 or crosstalkBack > 1:
                print("Invalid cross talk value (negative or >1), reverting to 0.")
                crosstalkBack = 0.0

    RunAnalysis(inputName, crosstalkSide, crosstalkBack) 
