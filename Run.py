from os import system as cmd

def ModifyGeom(angle):
    if angle[0] == "x":
        orientation = angle.strip("x") + " 0 0"
    elif angle[0] == "y":
        orientation = "0 " + angle.strip("y") + " 0"
    elif angle[0] == "z":
        orientation = "0 0 " + angle.strip("z")
    else:
        orientation = "0 0 0"

    for line in geomCont:
        if "orientation" in line:
            geomCont[geomCont.index(line)] = "orientation = " + orientation + "\n"
    writeFile = open(geomPath + "geom.conf","w")
    writeFile.writelines(geomCont)        
    writeFile.close()  


def ModifyConf(noise, nOfEvents):
    for line in configCont:
        if "electronics_noise" in line:
            configCont[configCont.index(line)] = "electronics_noise = " + noise + "\n"
        if "number_of_events" in line: 
            configCont[configCont.index(line)] = "number_of_events = " + nOfEvents + "\n"
    writeFile = open(configPath + "cfg.conf","w")
    writeFile.writelines(configCont)        
    writeFile.close()
    

def ModifyModel(thickness):
    for line in modelCont:
        if "sensor_thickness" in line:
            modelCont[modelCont.index(line)] = "sensor_thickness = " + thickness + "\n"
    writeFile = open(modelPath + "atlas17.conf","w")
    writeFile.writelines(modelCont)        
    writeFile.close()     


def RunSimulation():
    allpixPath = "/afs/cern.ch/user/r/rprivara/Allpix-" + allpixVers
    cmd(allpixPath + "/bin/allpix -c cfg.conf")
    print("Running simulation.", end = " ")
    


def MoveOutput(angle,noise,thickness):
    fileName = outputPath + angle + "-" + thickness + "-" + noise + "_output.root"
    cmd("cp " + outputPath + "output.root " + fileName)
    print("Moving output.")

#-------------------------------------------------------------------------------------------
allpixVers = "1.3"
configPath = geomPath = "/afs/cern.ch/user/r/rprivara/tb/"
modelPath = "/afs/cern.ch/user/r/rprivara/Allpix-" + allpixVers + "/models/"
outputPath = "/afs/cern.ch/user/r/rprivara/tb/output/"

#configPath = geomPath = "/home/b/pCloudDrive/Work/MgrThesis/Prog/AllPix/testing/"
#modelPath = "/home/b/pCloudDrive/Work/MgrThesis/Prog/AllPix/testing/"

# open, read and close config files
configFile = open(configPath + "cfg_def.conf", "r")
configCont = configFile.readlines()
configFile.close()
modelFile = open(modelPath + "atlas17_def.conf", "r")
modelCont = modelFile.readlines()
modelFile.close()
geomFile = open(geomPath + "geom_def.conf", "r")
geomCont = geomFile.readlines()
geomFile.close()

angles = ["0deg", "y5deg", "y12deg", "x23deg"]#, "x23deg"]#, "y10deg", "z15deg"]
noises = ["864e"]#, "700e", "900e"]
thicknesses = ["290um"]#, "300um", "305um", "310um", "315um"]
nOfEvents = "50000"

for angle in angles:
    for noise in noises:
        for thickness in thicknesses:
            print("noise:",noise,"thickness:",thickness, "angle:", angle)
            ModifyGeom(angle)
            ModifyConf(noise, nOfEvents)
            ModifyModel(thickness)
            RunSimulation()
            MoveOutput(angle, noise, thickness) 
