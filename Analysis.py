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

RunAnalysis("0deg-290um-864e_output.root", "0deg-290um-864e-CT_analysed.root", source="allpix", CT_StS=0.0188, CT_StBP=0.0211)
RunAnalysis("y5deg-290um-864e_output.root", "y5deg-290um-864e-CT_analysed.root", source="allpix", CT_StS=0.0188, CT_StBP=0.0211)
RunAnalysis("y12deg-290um-864e_output.root", "y12deg-290um-864e-CT_analysed.root", source="allpix", CT_StS=0.0188, CT_StBP=0.0211)
RunAnalysis("x23deg-290um-864e_output.root", "x23deg-290um-864e-CT_analysed.root", source="allpix", CT_StS=0.0188, CT_StBP=0.0211)

# RunAnalysis("rot0deg-300um_output.root","rot0deg-300um-noCT_analysed.root", CT_StS=0.00, CT_StBP=0.00)
# RunAnalysis("rot0deg-300um_output.root","rot0deg-300um-CT_analysed.root", CT_StS=0.042, CT_StBP=0.012)

# RunAnalysis("rotx23deg-300um_output.root","rotx23deg-300um-noCT_analysed.root", CT_StS=0.00, CT_StBP=0.00)
# RunAnalysis("rotx23deg-300um_output.root","rotx23deg-300um-CT_analysed.root", CT_StS=0.042, CT_StBP=0.012)

# RunAnalysis("roty5deg-300um_output.root","roty5deg-300um-noCT_analysed.root", CT_StS=0.00, CT_StBP=0.00)
# RunAnalysis("roty5deg-300um_output.root","roty5deg-300um-CT_analysed.root", CT_StS=0.042, CT_StBP=0.012)

# RunAnalysis("roty12deg-300um_output.root","roty12deg-300um-noCT_analysed.root", CT_StS=0.00, CT_StBP=0.00)
# RunAnalysis("roty12deg-300um_output.root","roty12deg-300um-CT_analysed.root", CT_StS=0.42, CT_StBP=0.012)

# RunAnalysis("roty5deg-270um_output.root")
# RunAnalysis("roty12deg-270um_output.root")
# RunAnalysis("roty15deg-270um_output.root")
# RunAnalysis("roty18deg-270um_output.root")
# RunAnalysis("rotx23deg-270um_output.root")

# RunAnalysis("rot0deg-300um_output.root", "rot0deg-LIV_analysed.root", CT_StS=0.042, CT_StBP=0.01)
# RunAnalysis("rot0deg-PEN_output.root", CT_StS=0.042, CT_StBP=0.01)
# RunAnalysis("rot0deg-EMV_output.root", CT_StS=0.042, CT_StBP=0.01)
# RunAnalysis("rot0deg-EMZ_output.root", CT_StS=0.0159, CT_StBP=0.01)

# RunAnalysis("0deg-290um-0e_output.root","0deg-290um-0e-athCT_analysed.root", source="allpix", CT_StS=0.07, CT_StBP=0.1)
# RunAnalysis("y5deg-290um-0e_output.root","y5deg-290um-0e-athCT_analysed.root", source="allpix", CT_StS=0.07, CT_StBP=0.1)
# RunAnalysis("y12deg-290um-0e_output.root","y12deg-290um-0e-athCT_analysed.root", source="allpix", CT_StS=0.07, CT_StBP=0.1)
# RunAnalysis("x23deg-290um-0e_output.root","x23deg-290um-0e-athCT_analysed.root", source="allpix", CT_StS=0.07, CT_StBP=0.1)


# RunAnalysis("0deg-athena-fullCT_output.root", "0deg-athena-fullCT_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-athena-noCT_output.root", "0deg-athena-noCT_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)

# RunAnalysis("0deg-290um-athena_output.root", "0deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("a5deg-290um-athena_output.root", "a5deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("b5deg-290um-athena_output.root", "b5deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("a12deg-290um-athena_output.root", "a12deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("b12deg-290um-athena_output.root", "b12deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("a23deg-290um-athena_output.root", "a23deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("b23deg-290um-athena_output.root", "b23deg-290um-athena_analysed.root", source="athena", CT_StS=0.0, CT_StBP=0.0)

# RunAnalysis("0deg-300um-step0.1_output.root", "0deg-300um-step0.1_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-step5_output.root", "0deg-300um-step5_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-FTFPBERT_output.root", "0deg-300um-FTFPBERT_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-QGSPBERT_output.root", "0deg-300um-QGSPBERT_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-FTFPBERTPEN_output.root", "0deg-300um-FTFPBERTPEN_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-FTFPBERTEMV_output.root", "0deg-300um-FTFPBERTEMV_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-FTFPBERTEMZ_output.root", "0deg-300um-FTFPBERTEMZ_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-noPAI_output.root", "0deg-300um-noPAI_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-chPerStep10_output.root", "0deg-300um-chPerStep10_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um-chPerStep100_output.root", "0deg-300um-chPerStep100_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-270um_output.root", "0deg-270um_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-280um_output.root", "0deg-280um_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-290um_output.root", "0deg-290um_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-300um_output.root", "0deg-300um_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)
# RunAnalysis("0deg-310um_output.root", "0deg-310um_analysed.root", source="allpix", CT_StS=0.0, CT_StBP=0.0)

# RunAnalysis("0deg-290um_output.root", "0deg-290um-ctL_analysed.root", source="allpix", CT_StS=0.0158, CT_StBP=0.0178)
# RunAnalysis("0deg-290um_output.root", "0deg-290um-ctH_analysed.root", source="allpix", CT_StS=0.0188, CT_StBP=0.0211)
