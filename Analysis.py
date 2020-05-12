#!/usr/bin/python3

from ROOT import TFile, TH1D, TH2D, TMath
import numpy as np
from datetime import datetime as date
 

def RunAnalysis(inputName, outputName="", source="", CT_StS=0.0, CT_StBP=0.0):
    if outputName == 0: 
        outputName = inputName.split("_")[0] + "_analysed.root"
    print("INPUT:", inputName, "\nOUTPUT:", outputName)
    rootFile = TFile("data/raw/" + inputName)
    writeFile = TFile("data/" + outputName, "recreate") 
    
    (thrStartFC, thrEndFC, thrStepFC) = (0, 8, 0.1)
    effTitle = "Efficiency - " + outputName.split("_")[0]
    effHist = TH1D(effTitle, effTitle, 200, thrStartFC, thrEndFC)
    clusTitle = "Cluster Size - " + outputName.split("_")[0]
    clusHist = TH2D(clusTitle, clusTitle, 200, thrStartFC, thrEndFC, 1600, 0, 10)

    if source == "athena":
        hitTree = rootFile.Get("SCT_RDOAnalysis").Get("SCT_RDOAna")
        nOfStrips = 1280 
        nOfParts = 0
        for event in hitTree:
            if len(list(event.charge)) > 0:
                nOfParts += 1 
    elif source == "allpix":
        hitTree = rootFile.PixelCharge
        nOfStrips = int(str(rootFile.models.Get("atlas17_dut").Get("number_of_pixels")).split(" ")[0])
        nOfParts = int(str(rootFile.config.Get("Allpix").Get("number_of_events")))
    else:
        print("Unknown source.")
        return 1

    print("CONFIG: Source:", source, ",  Events:", nOfParts, ",  CT_StS:", CT_StS, ",  CT_StBP:", CT_StBP)
        
    for thr in np.arange(thrStartFC, thrEndFC, thrStepFC):
        thrE = thr * 6242.2
        print("Current threshold [fC]:", round(thr, 1), end="\r")
        for event in hitTree:  # iterating over all events
            stripCharge = np.zeros(nOfStrips)
            if source == "allpix":
                for stripHit in event.dut:    # iterating over all strips hit in that event of Allpix simulation
                    stripCharge[stripHit.getIndex().X()] = stripHit.getCharge()
            elif source == "athena":
                nonEmptyStrips = list(event.strip_sdo)
                nonEmptyStripsCharge = list(event.charge)
                for i in range(len(nonEmptyStrips)):    # iterating over all strips hit in that event of Athena simulation
                    stripCharge[nonEmptyStrips[i]] = nonEmptyStripsCharge[i]
                
            stripChargeAdj = np.copy(stripCharge)
            if CT_StBP > 0 and CT_StS > 0:
                chargeMask = np.nonzero(stripCharge)
                for stripIndex in chargeMask:
                    chargeStoS = stripCharge[stripIndex] * CT_StS
                    chargeStoBP = stripCharge[stripIndex] * CT_StBP
                    try:
                        stripChargeAdj[stripIndex] -= 2 * chargeStoS + chargeStoBP
                        stripChargeAdj[stripIndex - 1] += chargeStoS
                        stripChargeAdj[stripIndex + 1] += chargeStoS
                    except IndexError:
                        pass

            cluster = len(np.where(stripChargeAdj > thrE)[0])
            if cluster > 0:
                effHist.Fill(thr)
                clusHist.Fill(thr, cluster)

    effHist.Scale(1 / nOfParts)
    clusHist = clusHist.ProfileX()
    print("Analysis done.                                         \n")
    rootFile.Close()
    writeFile.Write()
    writeFile.Close()

    logFile = open("log_analysis.txt", "a")
    # logFile.writelines("\nOUTPUT:\t" + outputName + "\nINPUT:\t" + inputName + "\nDATE:\t" + str(date.now())+ "\nCONFIG:\t" + "thickness=" + str(thickness) + "  events=" + str(nOfParts) + "  scaleXfactor=" + str(scaleXFactor) + "  CT_StS=" + str(CT_StS) +  "  CT_StBP=" + str(CT_StBP) + "  thrStep=" + str(thrStepFC))
    logFile.writelines(("\nOUTPUT:\t" + outputName + "\nINPUT:\t" + inputName + "\nDATE:\t" + str(date.now())+ "\nCONFIG:\t" + "  events=" + str(nOfParts) + "  CT_StS=" + str(CT_StS) +  "  CT_StBP=" + str(CT_StBP) + "  thrStep=" + str(thrStepFC)))
    logFile.write("\n-----\n")
    logFile.close()


# ------------------------------------------------------------------------
# CT l
# Cc = 24, Cb = 0.25+0.2, Ci = 0.8
# CT_StS = 1,58%, CT_StBP = 1,78%

# CT h
# Cc = 20, Cb = 0.45, Ci = 0.8
# CT_StS = 1,88%, CT_StBP = 2,11%

#CT final
# Cc = 25, Cb = 0.25, Ci = 0.8
# CT_StS = 1,53%, CT_StBP = 0,96%

# RunAnalysis("0deg-280um-864e-TCAD_output.root", "0deg-280um-864e-TCAD-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
# RunAnalysis("0deg-280um-864e-TCAD_output.root", "0deg-280um-864e-TCAD_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
RunAnalysis("y5deg-280um-864e-TCAD_output.root", "y5deg-280um-864e-TCAD-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("y12deg-280um-864e-TCAD_output.root", "y12deg-280um-864e-TCAD-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("x23deg-280um-864e-TCAD_output.root", "x23deg-280um-864e-TCAD-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)

RunAnalysis("0deg-25um-864e_output.root", "0deg-25um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("0deg-50um-864e_output.root", "0deg-50um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("0deg-100um-864e_output.root", "0deg-100um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("0deg-150um-864e_output.root", "0deg-150um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("0deg-200um-864e_output.root", "0deg-200um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
RunAnalysis("0deg-250um-864e_output.root", "0deg-250um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)

# RunAnalysis("0deg-270um-864e_output.root", "0deg-270um-864e_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-280um-864e_output.root", "0deg-280um-864e_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-290um-864e_output.root", "0deg-290um-864e_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-864e_output.root", "0deg-300um-864e_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-310um-864e_output.root", "0deg-310um-864e_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)

# RunAnalysis("0deg-270um-864e_output.root", "0deg-270um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
# RunAnalysis("0deg-280um-864e_output.root", "0deg-280um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
# RunAnalysis("0deg-290um-864e_output.root", "0deg-290um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
# RunAnalysis("0deg-300um-864e_output.root", "0deg-300um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)
# RunAnalysis("0deg-310um-864e_output.root", "0deg-310um-864e-CT_analysed.root", source="allpix", CT_StS=0.0153, CT_StBP=0.0096)

# RunAnalysis("0deg-280um-athena-cut50_output.root", "0deg-280um-athena-cut50um_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-280um-athena-cut15_output.root", "0deg-280um-athena-cut15um_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)

# RunAnalysis("0deg-25um-athena_output.root", "0deg-25um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-50um-athena_output.root", "0deg-50um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-100um-athena_output.root", "0deg-100um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-150um-athena_output.root", "0deg-150um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-200um-athena_output.root", "0deg-200um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-250um-athena_output.root", "0deg-250um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-270um-athena_output.root", "0deg-270um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)

# RunAnalysis("0deg-290um-athena_output.root", "0deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-athena_output.root", "0deg-300um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-310um-athena_output.root", "0deg-310um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
