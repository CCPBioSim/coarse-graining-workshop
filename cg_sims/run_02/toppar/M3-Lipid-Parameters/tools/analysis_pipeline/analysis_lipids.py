#!/usr/bin/env python3

# WARNING this is a simple analysis script for analysis of the auto-built lipid test simulations. 
# This script is not meant to be generalized but provided as an example.   

import math
import os, sys
import do_order_module_itp as do_order_module
import area_lipid_g5
import subprocess
import shlex
import re
import statistics
import multiprocessing
import numpy as np
from numpy import linspace,exp
from numpy.random import randn
import qwrap
import diffusion

#fpFS is the complete path of the directory that holds fatslim, please make sure
#to change this to the directory where you have fatslim installed
fpFS = "/g/g90/helgi/.venvs/devbending_specific/bin/fatslim"
# File with all lipid .itp's - used in do_order to get all @BONDS info
lipidItpFileM3 = "../top/martini_v3.0_lipids_all_alpha_v1.itp" # Martini 3 file
lipidItpFileM2 = "../top/martini2/martini_v2.0_lipids_all_201506.itp" # Martini 2 file

'''
 Simple analysis pipeline used for analysis of small Martini 2 and Martini 3 lipid bilayer simulations.
 Code adapted from larger analysis_pipeline described in:
   Blumer et al. Simulations of Asymmetric Membranes Illustrate Cooperative Leaflet Coupling and Lipid Adaptability.
   Front Cell Dev Biol 2020, DOI: 10.3389/fcell.2020.00575
   git: https://github.com/llnl-hmc-clinic/Bilayer-Simulation-Pipeline
'''
# Functions for handling xvg files
def parse_xvg(fname, sel_columns='all'):
    """Parses XVG file legends and data"""

    _ignored = set(('legend', 'view'))
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')

    metadata = {}
    num_data = []

    metadata['labels'] = {}
    metadata['labels']['series'] = []

    ff_path = os.path.abspath(fname)
    if not os.path.isfile(ff_path):
        raise IOError('File not readable: {0}'.format(ff_path))

    with open(ff_path, 'r') as fhandle:
        for line in fhandle:
            line = line.strip()
            if line.startswith('@'):
                tokens = shlex.split(line[1:])
                if tokens[0] in _ignored:
                    continue
                elif tokens[0] == 'TYPE':
                    if tokens[1] != 'xy':
                        raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
                elif _re_series.match(tokens[0]):
                    metadata['labels']['series'].append(tokens[-1])
                elif _re_xyaxis.match(tokens[0]):
                    metadata['labels'][tokens[0]] = tokens[-1]
                elif len(tokens) == 2:
                    metadata[tokens[0]] = tokens[1]
                else:
                    print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
            elif line[0].isdigit():
                num_data.append(map(float, line.split()))

    num_data = list(zip(*num_data))

    if not metadata['labels']['series']:
        for series in range(len(num_data) - 1):
            metadata['labels']['series'].append('')

    # Column selection if asked
    if sel_columns != 'all':
        sel_columns = map(int, sel_columns)
        x_axis = num_data[0]
        num_data = [x_axis] + [num_data[col] for col in sel_columns]
        metadata['labels']['series'] = [metadata['labels']['series'][col - 1] for col in sel_columns]

    return metadata, num_data

def running_average(data, metadata, window=10):
    """
    Performs a running average calculation over all series in data.
    Assumes the first series is the x-axis.
    Appends the series and a new label to the original data and label arrays.
    """

    weights = np.repeat(1.0, window)/window
    s_labels = metadata['labels']['series']
    for n_series, series in enumerate(data[1:]):
        series_rav = np.convolve(series, weights, 'valid')
        s_labels.append('{0} (Av)'.format(s_labels[n_series]))
        data.append(series_rav)
    return metadata, data


def avg_and_bse(num_data):
    ''' given a list of measurements, find their mean and find their block standard
        error '''
    dataMean = sum(num_data)/len(num_data)
    blockLengths = np.linspace(1, int(len(num_data)/5), num=int(len(num_data)/5))
    bse = []
    for n in blockLengths:
        numBlocks = len(num_data)/n
        blockAverages = []
        for i in range(int(numBlocks)):
            avgi = sum(num_data[int(i*n):int((i+1)*n)])/n
            blockAverages.append(avgi)
        # now we have means for all our blocks: compute standard error
        if len(blockAverages) > 1:
            blockse = statistics.stdev(blockAverages)/math.sqrt(numBlocks)
        else:
            blockse = 0
        bse.append(blockse)
    bseError = max(bse)
    return dataMean, bseError

# Adapted form Kasper Busk Pedersen (calc_2Dc.py, calc_DB.py, calc_DHH.py)
def calc_DHH(x,dens):
    x_neg_ndx = np.where(x < 0)[0]   #Split the density in two, to find both maxima in the electron density
    x_pos_ndx = np.where(x > 0)[0]
    x_neg = x[x_neg_ndx]
    x_pos = x[x_pos_ndx]

    x_neg_max = x_neg[np.argmax(dens[x_neg_ndx])]
    x_pos_max = x_pos[np.argmax(dens[x_pos_ndx])]
  
    return abs(x_pos_max - x_neg_max)

def calc_2Dc(x,dens):
    h = dens - (max(dens) / 2 )
    indexes = np.where(h > 0)[0]
    return abs(x[indexes[-1]] - x[indexes[0]])

def calc_DB(x,dens):
    d_z = np.abs(x[-1]-x[0])         # repeat distance
    bulkdens = np.mean(dens[0:10])
    norm_dens = np.divide(dens,bulkdens, out=np.zeros_like(dens), where=bulkdens!=0)
    h   = d_z - np.trapz(norm_dens, x)
    return h

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_density(ext_func, filenames, n_blocks):
    data = [np.loadtxt("{}_{}.xvg".format(filenames, i), comments=["#", "@"]) for i in range(n_blocks)]
    val = []
    for d in data:
        x = d[:,0]
        dens = d[:,1]

        smooth_dens = smooth(dens,5) #Smoothing needed due to noisy density in CG - window of 5 found by visual inspection, should be very robust
        x_cut = x[5:-6]              #Cut the first and last 5 datapoints to avoid artefacts from smoothing
        smooth_dens_cut = smooth_dens[5:-6]

        val.append(ext_func(x_cut,smooth_dens_cut))

    val_final_mean = np.mean(val)
    val_final_std = np.std(val,ddof=1) #Bessels correction
    val_final_sem = val_final_std / np.sqrt(n_blocks)

    return [val_final_mean, val_final_std, val_final_sem]


# input a list of lipid names in the same order as they appear in description.txt,
# the total time of the simulation (in ps) and the time skip between frames that are saved
# (again in ps)
def analyze(dirList, resultsFP):
    resultFile = open(resultsFP,"w+")
    resultArray = []
    resultArrayLegend = []
    beginAnalysis = 200000 # start analysis at 200 ns in ps
    precision = 0
    dt = 0
    nsteps = 0
    temperature = 310
    timeBetweenFrames = 1 # 1 ns
    finalSimulationTime = 1000000 # end analysis at 1us in ps

    # for each run, find the asymmetry value, the 20 fs tpr, gro, xtc, edr, mdout, trr  and index file
    for runDir in dirList:
        groFP = runDir + "confout.gro"
        tprFP = runDir + "topol.tpr"
        xtcFP = runDir + "traj_comp.xtc"
        edrFP = runDir + "ener.edr"
        indFP = runDir + "index.ndx"
        topFP = runDir + "system.top"
        mdoFP = "mdout.mdp"
        trrFP = runDir + "traj.trr"
        # Assuming first two letters for simulation folders are M2 or M3 for Martini 2 or Martini 3
        currentSim = os.path.basename(os.path.dirname(runDir)) # the provided inputs are always xxx/ so we unwrapp one time
        if currentSim[0:2] == "M2":
            lipidItpFile = lipidItpFileM2
        elif currentSim[0:2] == "M3":
            lipidItpFile = lipidItpFileM3
        else:
            print(f"ERROR - Martini version {currentSim[0:2]} not supported only M2 and M3 from sim {runDir}")
            sys.exit(-1)
        print('Running 1, dir {} result {} using .itp {} for order parameter bonds'.format(runDir, resultsFP, lipidItpFile))

        # read values out of topology file
        lipidNameList = []
        leaflet1 = []
        leaflet2 = []
        topFile = open(topFP, "r")
        # WARNING will pickup anything as lipid in system.top that has 3 or 4 characters
        p3 = re.compile('[A-Z]{2}[A-Z\d] ')
        p4 = re.compile('[A-Z]{3}[A-Z\d] ')
        for line in topFile.readlines():
            if bool(p3.match(line)) or bool(p4.match(line)):
                if line[0:4] in lipidNameList:
                    i = lipidNameList.index(line[0:4])
                    if len(leaflet2) < i:
                        while len(leaflet2) < i:
                            leaflet2.append(0)
                    leaflet2.append(int(line[4:]))
                else:
                    lipidNameList.append(line[0:4])
                    leaflet1.append(int(line[4:]))
        while len(leaflet2) < len(leaflet1):
            leaflet2.append(0)
        topFile.close()
        print(f" L1 names {lipidNameList} lipids {leaflet1} total {sum(leaflet1)}")
        print(f" L2 names {lipidNameList} lipids {leaflet2} total {sum(leaflet2)}")

        # Found Lipids per leaflet
        for i in range(len(lipidNameList)):
            resultFile.write(f"Lipid counts {lipidNameList[i]} lower {leaflet1[i]} upper {leaflet2[i]}\n")
            resultArray.append([leaflet1[i], leaflet2[i]])
            resultArrayLegend.append(f"Lipid counts {lipidNameList[i]} [lower,upper]")

        # Leaflet asymmetry
        asym = abs(sum(leaflet2) - sum(leaflet1))
        resultFile.write(f"Asymmetry {asym}\n")
        resultArray.append(asym)
        resultArrayLegend.append("asymmetry")
        topFile.close()
        # now calculate property values

        # find area per lipid and bilayer thickness
        # create index file
        if True:
            subprocess.call(["mkdir", runDir + "/analysis"])
            # make joined file - orginal shoudl work for all lipids 
            headgroupsFP = runDir + "/analysis/headgroups.ndx"
            if currentSim[0:2] == "M2" and any(x[2:4] == "CE" for x in lipidNameList): # In Martini 2 Ceramide (CE) lipids have no headgroup so use linker beads if CE lipids in sim
                command = "(echo del 0-100; echo a AM1 AM2 PO4 NC3 NH3 CNO ROH GL0 PS1 PS2 COH; echo name 0 headgroups; echo q) | gmx make_ndx -f " + groFP + " -o " + headgroupsFP
            else:  # all other case use full headgroup list - will pick what is needed - NOTE PIPs would only pickup PO4 which shoudl be enough
                command = "(echo del 0-100; echo a PO4 NC3 NH3 CNO ROH GL0 PS1 PS2 COH; echo name 0 headgroups; echo q) | gmx make_ndx -f " + groFP + " -o " + headgroupsFP
            subprocess.call(command, shell=True)
            # make PO4 index (warning will not work for CHOL, CER, DAG etc)
            headgroupsPO4 = runDir + "/analysis/index_PO4.ndx"
            command = "(echo del 0-100; echo a PO4; echo name 0 PO4; echo q) | gmx make_ndx -f " + groFP + " -o " + headgroupsPO4
            subprocess.call(command, shell=True)
            # make 2Dc index (warning will not work for CHOL)
            headgroups2Dc = runDir + "/analysis/index_2Dc.ndx"
            command = "(echo del 0-100; echo 'a C* | a D* | a T*'; echo name 0 DC2; echo q) | gmx make_ndx -f " + groFP + " -o " + headgroups2Dc
            subprocess.call(command, shell=True)
            # make membrane index - assumes simulation only has W NA CL and membrane 
            indexMem = runDir + "/analysis/index_mem.ndx"
            command = "(echo keep 0; echo 'a W | a NA | a CL'; echo !1; echo name 2 membrane; echo 'a W'; echo name 3 W; echo q) | gmx make_ndx -f " + groFP + " -o " + indexMem
            subprocess.call(command, shell=True)


        print('Running 2, asym {} finalSimulationTime {}'.format(asym, finalSimulationTime))

        # Do analysis on last 800ns (skip first 200ns)
        if False:
            command = f'(echo 2; echo 0) | gmx trjconv -f {xtcFP} -o {runDir}/traj_centered.xtc -n {indFP} -center -pbc mol -s {tprFP} -b {beginAnalysis}'
            subprocess.call(command, shell=True)
            xtcFP = runDir + '/traj_centered.xtc'

        # Make the trajectory whole for CoM computation and unwrapping
        if True:
            command = f"echo \"System System\" | gmx trjconv -pbc whole -center -s {tprFP} -f {xtcFP} -o {runDir}/traj_whole.xtc -b {beginAnalysis} -e {finalSimulationTime}"
            subprocess.call(command, shell=True)
            xtcFP = runDir + '/traj_whole.xtc'

        print('Running 3, xtcFP {}'.format(xtcFP))

        if True:
            # call FatSlim commands
            subprocess.call(["python3", fpFS, "membranes", "-c", groFP, "-n", headgroupsFP, "--output-index", runDir + "/analysis/bilayer_leaflet.ndx"])
            subprocess.call(["python3", fpFS, "thickness", "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(beginAnalysis), "-e", str(finalSimulationTime), "--plot-thickness", runDir + "/analysis/thickness.xvg", "--nthreads", "2"])
            subprocess.call(["python3", fpFS, "apl",       "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(beginAnalysis), "-e", str(finalSimulationTime), "--plot-apl", runDir + "/analysis/apl.xvg", "--plot-area", runDir + "/analysis/area.xvg", "--nthreads", "2"])
            subprocess.call(["python3", fpFS, "apl",       "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(beginAnalysis), "-e", str(finalSimulationTime), "--plot-apl", runDir + "/analysis/apl_type.xvg", "--nthreads", "2", "--apl-by-type"])
            # new FatSlim commands to match Kasper calls
            subprocess.call(["python3", fpFS, "apl",       "-c", groFP, "-n", headgroupsPO4,"-t", xtcFP, "-b", str(beginAnalysis), "-e", str(finalSimulationTime), "--hg-group", "PO4", "--plot-area", runDir + "/analysis/apl_PO4.xvg", "--nthreads", "2", "--apl-cutoff", "10.0"])
            subprocess.call(["python3", fpFS, "thickness", "-c", groFP, "-n", headgroupsPO4,"-t", xtcFP, "-b", str(beginAnalysis), "-e", str(finalSimulationTime), "--hg-group", "PO4", "--plot-thickness", runDir + "/analysis/thickness_PO4.xvg", "--nthreads", "2", "--thickness-cutoff", "8.0"])

        if True:
            # parse the output xvg
            metadata, num_data = parse_xvg(runDir + "/analysis/thickness.xvg")
            thickness, thicknesserror = avg_and_bse(num_data[1])
            resultFile.write("Bilayer Thickness = %s with a standard error of %s\n" % (thickness, thicknesserror))
            resultArray.append([thickness, thicknesserror])
            resultArrayLegend.append("FS thickness [avg,SE]")

            metadata, num_data = parse_xvg(runDir + "/analysis/apl.xvg")
            apl_down, aplerror_down = avg_and_bse(num_data[2])
            apl_up, aplerror_up = avg_and_bse(num_data[3])
            resultFile.write("Area per Lipid (Lower) = %s with a standard error of %s\n" % (apl_down, aplerror_down))
            resultFile.write("Area per Lipid (Upper) = %s with a standard error of %s\n" % (apl_up, aplerror_up))
            resultArray.append([apl_down, aplerror_down, apl_up, aplerror_up])
            resultArrayLegend.append("FS avg APL [lower,SE,upper,SE]")

            metadata, num_data = parse_xvg(runDir + "/analysis/apl_type.xvg")
            apl_L1, aplerror_L1 = avg_and_bse(num_data[1])
            apl_L2, aplerror_L2 = avg_and_bse(num_data[2])
            apl_L3, aplerror_L3 = avg_and_bse(num_data[3])
            resultFile.write("APL L1 = %s with a standard error of %s\n" % (apl_L1, aplerror_L1))
            resultArray.append([apl_L1, aplerror_L1])
            resultArrayLegend.append("FS L1 APL [avg,SE]")
            resultFile.write("APL L2 = %s with a standard error of %s\n" % (apl_L2, aplerror_L2))
            resultArray.append([apl_L2, aplerror_L2])
            resultArrayLegend.append("FS L2 APL [avg,SE]")
            if len(leaflet2) > 2:
                apl_L3, aplerror_L3 = avg_and_bse(num_data[3])
                resultFile.write("APL L3+ = %s with a standard error of %s\n" % (apl_L3, aplerror_L3))
                resultArray.append([apl_L3, aplerror_L3])
                resultArrayLegend.append("FS L3 APL [avg,SE]")

            # Get props form new FatSlim commands to match Kasper calls
            try:
                metadata, num_data = parse_xvg(runDir + "/analysis/thickness_PO4.xvg")
                thickness, thicknesserror = avg_and_bse(num_data[1])
            except Exception as error:
                print(' Error in reading - thickness_PO4.xvg, this can happen e.g. if no lipids with PO4')
                thickness = -1
                thicknesserror = -1
            resultFile.write("Bilayer Thickness PO4 = %s with a standard error of %s\n" % (thickness, thicknesserror))
            resultArray.append([thickness, thicknesserror])
            resultArrayLegend.append("FS PO4 thickness [avg,SE]")

            try:
                metadata, num_data = parse_xvg(runDir + "/analysis/apl_PO4.xvg")
                apl, aplerror = avg_and_bse(num_data[1])
            except Exception as error:
                print(' Error in reading - apl_PO4.xvg this can happen e.g. if no lipids with PO4')
                apl = -1 
                aplerror = -1
            resultFile.write("Area per Lipid PO4 = %s with a standard error of %s\n" % (apl, aplerror))
            resultArray.append([apl, aplerror])
            resultArrayLegend.append("FS PO4 avg APL [avg,SE]")

            print('Running 4, fatslim done, apl_down {}'.format(apl_down))

        if True:
            # Get thickness values from density profiles - adopted from Kasper Busk Pedersen 
            # DHH (peak-peak) thickness, DB (Luzzati) thickness, and 2Dc (hydrocarbon) thickness
            n_blocks=8 # 200-1000 ns normal analysis use 8 100 ns windows
            block_length=100000 # 100 ns in ps
            electron_data = "../tools/analysis_pipeline/martini_electrons.dat" # fix path
            # DHH (peak-peak) thickness 
            for i in range(n_blocks): 
                command = f"echo \"membrane System\" | gmx density -sl 200 -dens electron -ei {electron_data} -ng 1 -n {indexMem} -f {xtcFP} -s {tprFP} -center -symm -relative -o {runDir}/analysis/DHH_dens_{i}.xvg -b {beginAnalysis + (block_length * i)} -e {beginAnalysis + (block_length * (i+1))}"
                subprocess.call(command, shell=True) 
            # DB (Luzzati) thickness
            for i in range(n_blocks):
                command = f"echo \"membrane W\" | gmx density -sl 200 -dens number -ng 1 -n {indexMem} -f {xtcFP} -s {tprFP} -center -symm -relative -o {runDir}/analysis/W_dens_{i}.xvg -b {beginAnalysis + (block_length * i)} -e {beginAnalysis + (block_length * (i+1))}"
                #print(f"HII debug run W {i} c: {command}")
                subprocess.call(command, shell=True) 
            # 2Dc (hydrocarbon) thickness
            for i in range(n_blocks):
                command = f"echo \"DC2 DC2\" | gmx density -sl 200 -dens number -ng 1 -n {headgroups2Dc} -f {xtcFP} -s {tprFP} -center -symm -relative -o {runDir}/analysis/2Dc_dens_{i}.xvg -b {beginAnalysis + (block_length * i)} -e {beginAnalysis + (block_length * (i+1))}"
                #print(f"HII debug run C {i} c: {command}")
                subprocess.call(command, shell=True) 

            # Read data for all 
            [mean, sd, se] = get_density(calc_DHH, runDir + "/analysis/DHH_dens", n_blocks) 
            resultFile.write("DHH thickness = %s with a standard error of %s\n" % (mean, se))
            resultArray.append([mean, se])
            resultArrayLegend.append("DHH thickness [avg,SE]")
            [mean, sd, se] = get_density(calc_DB, runDir + "/analysis/W_dens", n_blocks)
            resultFile.write("DB thickness = %s with a standard error of %s\n" % (mean, se))
            resultArray.append([mean, se])
            resultArrayLegend.append("DB thickness [avg,SE]")
            [mean, sd, se] = get_density(calc_2Dc, runDir + "/analysis/2Dc_dens", n_blocks)
            resultFile.write("2Dc thickness = %s with a standard error of %s\n" % (mean, se))
            resultArray.append([mean, se])
            resultArrayLegend.append("2Dc thickness [avg,SE]")


        if True:
            # find area compressibility and APL from box size - note here only doing sum(leaflet1) as assumign simulation is NOT asymetric
            area = area_lipid_g5.calculate_area_compressibility(edrFP, runDir + "/analysis", beginAnalysis, finalSimulationTime, sum(leaflet1), temperature)
            resultFile.write("Area compressibility box = {} with a standard error of {} and SD of {}\n".format(area[3], area[5], area[4]))
            resultArray.append([area[3], area[5], area[4]])
            resultArrayLegend.append("Box Ka [avg,SE,SD]")
            resultFile.write("Area per lipid box = {} with a standard error of {} and SD of {}\n".format(area[0], area[2], area[1]))
            resultArray.append([area[0], area[2], area[1]])
            resultArrayLegend.append("Box APL [avg,SE,SD]")
            print('Running 5, area box {} lipid in a leaflet {}'.format(area, sum(leaflet1)))

        if True:
            # create index files for order param
            command = f"(echo del 1; echo q) | gmx make_ndx -f {tprFP} -n {runDir}/analysis/bilayer_leaflet_0000.ndx -o {runDir}/analysis/up.ndx"
            subprocess.call(command, shell=True)
            command = f"(echo del 0; echo q) | gmx make_ndx -f {tprFP} -n {runDir}/analysis/bilayer_leaflet_0000.ndx -o {runDir}/analysis/down.ndx"
            subprocess.call(command, shell=True)
            print('Running 6, index file created')

        up = -1
        down = -1
        numfailures = 0
        if True:
            # find order parameter for each lipid type
            # NOTE, these can fail for a number or reasons - still need to track down more
            #  - If any lipids flip-flops between the leaflets this one fails as the lipids count is hardcoded per leaflet - as what FS determined
            for i in range(len(lipidNameList)):
                if lipidNameList[i] == 'CHOL':
                    continue
                if leaflet2[i] > 0:
                    try:
                        up = do_order_module.computeOrderParameter(runDir + "/analysis/", xtcFP, tprFP, runDir + "/analysis/up.ndx", beginAnalysis, finalSimulationTime, 5, 0, 0, 1, leaflet2[i], lipidNameList[i], lipidItpFile, f"{lipidNameList[i]}_up_order.dat")
                        resultFile.write("Average order parameter (Upper) for {} = {}\n".format(lipidNameList[i], up))
                        resultArray.append(up)
                        resultArrayLegend.append(f"OrderP up {lipidNameList[i]} [avg]")
                    except Exception as error:
                        numfailures += 1
                        resultFile.write("Average order parameter (Upper) for {} = {}\n".format(lipidNameList[i], "n/a"))
                        resultArray.append(-1)
                        resultArrayLegend.append(f"OrderP up {lipidNameList[i]} [avg]")
                        #raise error
                if leaflet1[i] > 0:
                    try:
                        down = do_order_module.computeOrderParameter(runDir + "/analysis/", xtcFP, tprFP, runDir + "/analysis/down.ndx", beginAnalysis, finalSimulationTime, 5, 0, 0, 1, leaflet1[i], lipidNameList[i], lipidItpFile, f"{lipidNameList[i]}_down_order.dat")
                        resultFile.write("Average order parameter (Lower) for {} = {}\n".format(lipidNameList[i], down))
                        resultArray.append(down)
                        resultArrayLegend.append(f"OrderP down {lipidNameList[i]} [avg]")
                    except Exception as error:
                        numfailures += 1
                        resultFile.write("Average order parameter (Lower) for {} = {}\n".format(lipidNameList[i], "n/a"))
                        resultArray.append(-1)
                        resultArrayLegend.append(f"OrderP down {lipidNameList[i]} [avg]")
                        #raise error

            print('Running 7, order params - failures {} up {} down {}'.format(numfailures, up, down))

        # Diffusion analysis - added by Balazs Fabian <fbalazsf@gmail.com>
        numfailures = 0
        if True:
            # Compute CoM and unwrap - WARNING xtcFP needs to be -pbc whole and -center 
            qwrap.qwrap(runDir+"topol.tpr", xtcFP, runDir+"analysis/com.gro", runDir+"analysis/unwrapped.xtc", ' '.join(lipidNameList))

            # Create groups of [lipidname, leaflet]
            headGroups = []
            for i in range(len(lipidNameList)):
                lipid = lipidNameList[i]
                if leaflet1[i] > 0:
                    headGroups.append([lipid, 1])

                if leaflet2[i] > 0:
                    headGroups.append([lipid, 0])

            # Use the GLS estimator (https://aip.scitation.org/doi/full/10.1063/5.0008312)
            glsIntervalLower = 1
            glsIntervalUpper = 40 # was 50 for 1000 ns - changed to 40 for 800 ns 
            for group in headGroups:
                lipid = group[0]
                if group[1] == 1:
                    leaflet_name = "up"
                elif group[1] == 0:
                    leaflet_name = "down"
                else:
                    print(f"ERROR leaflet_group={group[1]} not supported")
                try:
                    # group: lipid name and leaflet (1 or 0)
                    print(f"  diffusion.diffusion({runDir}analysis/com.gro, {runDir}analysis/unwrapped.xtc, {runDir}analysis/gls-{group[0]}-{group[1]}, {glsIntervalLower}, {glsIntervalUpper}, {group[0]}, {group[1]})")
                    diffusion.diffusion(runDir+"analysis/com.gro", runDir+"analysis/unwrapped.xtc", f"{runDir}analysis/gls-{group[0]}-{group[1]}", glsIntervalLower, glsIntervalUpper, group[0], group[1])
                    # load and write results
                    data = np.loadtxt(f"{runDir}analysis/gls-{group[0]}-{group[1]}.dat",skiprows=11,max_rows=(glsIntervalUpper-glsIntervalLower))
                    data = data[9][1:3]
                    coeff = data[0] * 1000          # nm^2/ns
                    err   = np.sqrt(data[1]) * 1000 # nm^2/ns
                    resultFile.write("Lateral Diffusion Coef (leaflet %s) for %s = %s with a standard error of %s [nm^2/ns]\n" % (leaflet_name, lipid, coeff, err))
                    resultArrayLegend.append(f"Diffusion {leaflet_name} {lipid} [coeff,err]")
                    resultArray.append([coeff, err])
                except Exception as error:
                    numfailures += 1
                    resultFile.write("Lateral Diffusion Coef (leaflet %s) for %s = n/a\n" % (leaflet_name, lipid))
                    resultArrayLegend.append(f"Diffusion {leaflet_name} {lipid} [coeff,err]")
                    resultArray.append(-1)
                    #raise error

            print(f'Running 8, Diffusion analysis done - failures {numfailures}')  # BF

    resultFile.close()
    np.savez(resultsFP[0:-4], legend=np.array(resultArrayLegend), data=np.array(resultArray))

analyze([sys.argv[1]], f"{sys.argv[1]}/analysisResults.txt")
