#!/usr/bin/python3

from ROOT import TFile, TH1D, TH2D, TMath
import numpy as np


def RunAnalysis(inputName, outputName="", CT_StS=0.0, CT_StBP=0.0):
    # path = "/home/b/pCloudDrive/Work/MgrThesis/Prog/AllPix/tb/output/"
    rootFile = TFile("data/raw/" + inputName)
    if outputName == 0: outputName = inputName.split("_")[0] + "_analysed.root"
    print("INPUT:", inputName, "\nOUTPUT:", outputName)

    nOfParts = int(str(rootFile.config.Get("Allpix").Get("number_of_events")))
    thickness = rootFile.models.Get("atlas17_dut").Get("sensor_thickness")  # itk_strip2_dut or atlas17_dut
    scaleXFactor = 1  # otherwise int(str(thickness).strip("um"))
    nOfStrips = int(str(rootFile.models.Get("atlas17_dut").Get("number_of_pixels")).split(" ")[0])  # itk_strip2_dut or atlas17_dut
    # orientation = rootFile.detectors.Get("dut").Get("orientation")
    print("CONFIG: Sensor:", thickness, ",  Events:", nOfParts, ", scaleXFactor:", scaleXFactor, ", nOfStrips:",
          nOfStrips, ", CT_StS:", CT_StS, ", CT_StBP:", CT_StBP)

    writeFile = TFile("data/" + outputName, "recreate")
    hitTree = rootFile.PixelCharge

    (thrStartFC, thrEndFC, thrStepFC) = (0, 8, 0.2)
    effTitle = "Efficiency - " + outputName.split("_")[0]
    effHist = TH1D(effTitle, effTitle, 200, thrStartFC / scaleXFactor, thrEndFC / scaleXFactor)
    clusTitle = "Cluster Size - " + outputName.split("_")[0]
    clusHist = TH2D(clusTitle, clusTitle, 200, thrStartFC / scaleXFactor, thrEndFC / scaleXFactor, 1600, 0, 10)

    for thr in np.arange(thrStartFC, thrEndFC, thrStepFC):
        thrE = thr * 6242.2
        print("Current threshold [fC]:", round(thr, 1), end="\r")
        for event in hitTree:  # iterating over all events
            stripCharge = np.zeros(nOfStrips)
            for stripHit in event.dut:  # iterating over all strips hit in that event
                stripCharge[stripHit.getIndex().X()] = stripHit.getCharge()

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
                effHist.Fill(thr / scaleXFactor)
                clusHist.Fill(thr / scaleXFactor, cluster)

    effHist.Scale(1 / nOfParts)
    clusHist = clusHist.ProfileX()
    print("Analysis done.                           \n")
    rootFile.Close()
    writeFile.Write()
    writeFile.Close()


# ------------------------------------------------------------------------

# Cc = 24, Cb = 0.25+0.2, Ci = 0.8
# CT_StS = 2,48%, CT_StBP = 1,78%

# Cc = 24, Cb = 0.25, Ci = 0.8
# CT_Sts = 2.09%, CT_StBP = 1%

# Cc = 20, Cb = 0.45, Ci = 0.8
# CT_StS = 2,94%, CT_StBP = 2,11%

#RunAnalysis("0deg-290um-864e_output.root", "0deg-290um-864e-CT_analysed.root", CT_StS=0.0248, CT_StBP=0.018)
#RunAnalysis("y5deg-290um-864e_output.root", "y5deg-290um-864e-CT_analysed.root", CT_StS=0.0248, CT_StBP=0.018)
#RunAnalysis("y12deg-290um-864e_output.root", "y12deg-290um-864e-CT_analysed.root", CT_StS=0.0248, CT_StBP=0.018)
#RunAnalysis("x23deg-290um-864e_output.root", "x23deg-290um-864e-CT_analysed.root", CT_StS=0.0248, CT_StBP=0.018)

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

RunAnalysis("0deg-290um-864e_output.root","test.root", CT_StS=0.0, CT_StBP=0.0)
