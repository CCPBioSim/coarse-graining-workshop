#!/usr/bin/env python3

"""make-itps-from-table-v02.py reads in a .csv table with info on how to make lipids - creates all lipdis and can setup and run simulations, use make-itps-from-table-v02.py -h for description"""

__author__  = "Helgi I. Ingolfsson"
__status__  = "Development"
__version__ = "M3.l01"
__email__   = "ingolfsson1@llnl.gov"

import subprocess as subp
import sys
import math
import random
import os
import csv
from timeit import default_timer as timer
import multiprocessing as mp

# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self):
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

# Description
desc = """
    TBD
Use:
       """

# Options
options = [
"""
Options:""",
("-table",    Option(str,    1,        "CG-Martini3-lipids-naming.csv", "Name of lipid definition .csv input file")),
("-ofolder",  Option(str,    1,        "lipid_itps", "Name of output folder for all single lipids .itp files generated")),
("-makeItps", Option(str,    1,        "no", "Make all lipid .itp's  yes/no (def. no)")),
("-setupSims",Option(str,    1,        "no", "Setup simulations  yes/no (def. no)")),
("-setupCustomSims",Option(str,    1,  "no", "Setup custom simulations  yes/no (def. no)")),
("-runSims",  Option(str,    1,        "no", "Run simulations in current dir yes/no (def. no)")),
("-analyseSims", Option(str, 1,        "no", "Analyse simulations in current dir yes/no (def. no)")),
("-sRunC",    Option(str,    1,        "2", "Numer of sims to run")),
("-sRunId",   Option(str,    1,        "0", "Simulation number to start with")),
("-baseLipid",Option(str,    1,        "pure", "Base (80%) lipid, options are: {pure (def.), POPC, DLPE, SSM-CHOL}")),
("-CBT",      Option(str,    1,        "yes", "Use CBT dihedrals potentials yes/no (def. yes)")),
          ]

# Options
itpGenerator = "./lipid-itp-generator-Martini3-01.py "
singleLipiditpFolder = "lipid_itps"

# Need to source GROMACX for script to work
# e.g. source ~/.local/gromacs-2020.04/bin/GMXRC

# Setup simulations for all lipids in table (with x or y in auto column)
# Note - these are of one simulation type (pure, mixed 20% with one of x3 referances), with and without CBT

# For base lipids we want x4 implmented options "pure", "POPC", "DIPE" or "DSSM-CHOL"
#baseLipid = "pure" # tested and works all referance lipid no base
#baseLipid = "POPC" # 20% lipid x in 80% POPC - HII testing this now
#baseLipid = "DLPE" # 20% lipid x in 80% DLPE (DLiPE in Martini2 called DIPE) - HII testing this now
#baseLipid = "DSSM-CHOL" # 20% lipid x in 50% DSSM and 30% Cholesterol - HII still missign

gitDir = "../.."  # This is relative from one up from cDir (in current sim dir)
resources = f"{gitDir}/tools/resources"
# Read in name of .csv file with lipid input definitons
lipidTableFile = ""


''' tableFile file format for read in .txt or .csv file
cData[0]  -> auto == x (else other), now using -1 as notes and category seperators
cData[1]  -> main category
cData[2]  -> sub category name - full
cData[3]  -> sub category acrynom
cData[4]  -> for letter name - main name
cData[5]  -> tail chain b
cData[6]  -> tail chain a
cData[7]  -> tail chain c (only use in TAG lipids)
cData[8]  -> tail name
cData[9]  -> common name
cData[10]  -> description
cData[11]  -> keywords
cData[12] -> parameterization
cData[13] -> refs
cData[14] -> created
cData[15] -> authors
cData[16] -> modified
cData[17] -> ref area per lipid
cData[18] -> Warnings / Notes
cData[19] -> not used
cData[20] -> -alhead
cData[21] -> -allink
cData[22] -> -altail
cData[23] -> -alcharge
cData[24] -> not used
'''

def makeitps(tableData, options):
    cMasteritpFileName = "temp_inital"
    if os.path.isdir(singleLipiditpFolder) != True:
        os.mkdir(singleLipiditpFolder)

    for cline in tableData:
        #cline = cline.replace('\ufeff', '') # remove BOM/utf-8-sig tag from Excel
        #cData = cline.split(',') # WARNING will get funky if there are any , in the Excel field
        cData = cline # .csv reader already seprated and parased
        initNote = cData.pop(0)

        #for i in range(len(cData)):
        #    print(f"cData {i} -> {cData[i]}")
        #if True:  # End here if only to print commands
        #    continue

        if cData[0]=='-1':
            # if -1, then this is a new master .itp file, read in name and comment
            print(f"MASTER -> {cData}")

            cMasteritpFileName = cData[1]
            cMasteritpFile = open(cMasteritpFileName,"w")
            print(f"Creating itp lipid collection {cMasteritpFile.name}")
            print(';', file=cMasteritpFile)
            print('; ' + cData[2], file=cMasteritpFile)
            print('; ' + cData[3], file=cMasteritpFile)
            print(';', file=cMasteritpFile)
            print("; Please cite:\n; " + cData[13].replace('\\n','\n;  ') + "\n;", file=cMasteritpFile)
            print('; Last updated:  ' + cData[16] + ";", file=cMasteritpFile)
            print('; Authors: ' + cData[15], file=cMasteritpFile)
            print('\n\n', file=cMasteritpFile)
            cMasteritpFile.close()

        elif cData[0]=='x' or cData[0]=='y':
            # Use this line to autogenerate lipid
            print(f"MAKEitp -> {cData}")

            # Now get all the "coorect data" for the lipid
            current = itpGenerator
            itpFileName =  f"{singleLipiditpFolder}/martini_v3.0_{cData[4]}_01.itp"
            current += "-alname " + cData[4] + " -alhead \"" + cData[20] + "\" -allink \"" + cData[21] + "\" -altail \"" + str(cData[22]).rstrip() + "\" -o " + itpFileName + " -name \""+ cData[9] +"\" "
            current += "-desc \"" + cData[10] + "\" -keyw \"" + cData[11] + "\" -parm \"" + cData[12] + "\" -refs \"" + cData[13] + "\" -crea \"" + cData[14] + "\" -auth \"" + cData[15] + "\"  "
            current += "-modi \"" + cData[16] + "\" -area \"" + cData[17] + "\" -warn \"" + cData[18] + "\" "
            print(f"  CALL -> {current}")

            # Make lipid .itp using lipid-martini-itp-vXXX
            subp.call(current, shell=True, cwd=cDir)

            # Append current lipid to active master itp file
            subp.call(f"cat {itpFileName} >> {cMasteritpFileName}", shell=True, cwd=cDir)

        else:
            #Found line not with -1,x,y so currently not used for autogenerated lipids
            print(f"NOT USED: {cData}")
# end makeitps

def printSim(sim):
    csim, csimDir, cData = sim
    print(f" -> print {csim}")
    return csim

def eqSim(sim):
    csim, csimDir, options, cData = sim
    baseLipid = options["-baseLipid"].value
    print(f" -> setup {csim} - started - full path {csimDir} - with baselipid {baseLipid}", flush = True)
    start = timer()

    os.mkdir(csimDir)

    subp.call(f"cp {gitDir}/tools/resources/martini_v2.x_new-rf-prod.mdp .", shell=True, cwd=csimDir)
    use_CBT = options["-CBT"].value == "yes"
    if use_CBT:
        subp.call(f'sed -ie "s/x_define_x/define                   = -DM3_CBT/g" martini_v2.x_new-rf-prod.mdp', shell=True, cwd=csimDir)
    else:
        subp.call(f'sed -ie "s/x_define_x/;define                   =/g" martini_v2.x_new-rf-prod.mdp', shell=True, cwd=csimDir)

    cbtsim = csim.replace("noCBT","CBT")
    if (not use_CBT) and os.path.isdir(cbtsim):
        # check if CBT true is done for this sim
        # use doen sim and just redo last topol
        subp.call(f'cp ../{cbtsim}/system.top . ', shell=True, cwd=csimDir)
        subp.call(f'cp ../{cbtsim}/lipids-water-eq4.gro . ', shell=True, cwd=csimDir)
        subp.call(f'cp ../{cbtsim}/index.ndx . ', shell=True, cwd=csimDir)
    else:
        # do regular setup

        # copy in resources
        subp.call(f"cp {gitDir}/tools/resources/system-cg-default-M3.top system.top", shell=True, cwd=csimDir)
        subp.call('sed -ie "s/xPATHx/\..\/..\/top/g" system.top', shell=True, cwd=csimDir)  # @WARNING gitDir hardcoded some complication with / and \/ in f statments

        if len(cData) < 10:
            # this is used in custom build mode
            string_insane = f"python2 {gitDir}/tools/insane_M3_lipids.py -x 8 -y 8 -z 10 -pbc square -sol W -solr 0.5 -o lipids-water.gro -l {cData[0]}:{cData[2]} -l {cData[1]}:{cData[3]} -salt 0.15 -asym 0 -rand 0.095 2>&1 | tee -a system.top"
        else:
            # this is normal build mode
            headgroup_autobuild_insane = ["C", "E", "G", "S", "PS1", "PS2", "P", "O", "COH"]
            lipidautoscaffold = ""
            # warning baseLipids POPC and DLPE needs to be in insane, and if referance one of those insane definition will be overwritten and therefore should be the same
            # Also all lipids not in auto scaffold e.g. PI PIPs etc need to be defined in insane alreaddy
            if all(elem in headgroup_autobuild_insane for elem in cData[20].split()):
                lipidautoscaffold = f" -alname {cData[4]} -alhead \"{cData[20]}\" -allink \"{cData[21]}\" -altail \"{str(cData[22]).rstrip().replace('c','C').replace('d','D').replace('F','D').replace('t','T')}\" -alcharge \"{cData[23]}\""
            # construct simulation
            if baseLipid == "pure":
                string_insane = f"python2 {gitDir}/tools/insane_M3_lipids.py -x 8 -y 8 -z 10 -pbc square -sol W -solr 0.5 -o lipids-water.gro -l {cData[4]} -salt 0.15 -asym 0 -rand 0.095 {lipidautoscaffold} 2>&1 | tee -a system.top"
            elif baseLipid == "POPC" or baseLipid == "DLPE":
                if baseLipid == cData[4]: # These help with the analysis and reading of each leaflet - simulations actually not needed as also part of pure
                    string_insane = f"python2 {gitDir}/tools/insane_M3_lipids.py -x 8 -y 8 -z 10 -pbc square -sol W -solr 0.5 -o lipids-water.gro -l {baseLipid}:100 -salt 0.15 -asym 0 -rand 0.095 {lipidautoscaffold} 2>&1 | tee -a system.top"
                else:
                    string_insane = f"python2 {gitDir}/tools/insane_M3_lipids.py -x 8 -y 8 -z 10 -pbc square -sol W -solr 0.5 -o lipids-water.gro -l {baseLipid}:80 -l {cData[4]}:20 -salt 0.15 -asym 0 -rand 0.095 {lipidautoscaffold} 2>&1 | tee -a system.top"
            elif baseLipid == "SSM-CHOL":
                # warning SSM and CHOL needs to be in insane
                string_insane = f"python2 {gitDir}/tools/insane_M3_lipids.py -x 8 -y 8 -z 10 -pbc square -sol W -solr 0.5 -o lipids-water.gro -l SSM:50 -l CHOL:30 -l {cData[4]}:20 -salt 0.15 -asym 0 -rand 0.095 {lipidautoscaffold} 2>&1 | tee -a system.top"
            else:
                print(f" -> setup {csim} - ERROR - baseLipid {baseLipid} not supported", flush = True)
                sys.exit(-1)

        # @TODO other base examples -l ${cBase}:51 -l ${clipid}:13 -o lipids-water.gro 2>&1 | tee -a system.top
        print(f" -> setup {csim} - insane string: {string_insane}", flush = True)
        subp.call(string_insane, shell=True, cwd=csimDir)

        # Run energy minimization (steepest descent, with protein constreins)
        subp.call(f'gmx grompp -c lipids-water.gro -f {resources}/martini_v2.x_new-rf-em.mdp -p system.top -o topol.tpr -maxwarn 5', shell=True, cwd=csimDir)
        # @NOTE: some gmx versions are compiled without thread-MPI and do not support setting the number of threads "-nt"
        subp.call('gmx mdrun -nt 1 -v -c lipids-water-em.gro', shell=True, cwd=csimDir)

        # Return here if just doing .itp testing
        #return csim

        # Make index file
        string_index = '''
        echo "del 1 - 200" > index-selection.txt
        echo "a W | a NA | a CL" >> index-selection.txt
        echo "name 1 Solvent" >> index-selection.txt
        echo '!1' >> index-selection.txt
        echo "name 2 Rest" >> index-selection.txt
        echo "q" >> index-selection.txt
        gmx make_ndx -f lipids-water-em.gro -o index.ndx < index-selection.txt
        '''
        subp.call(string_index, shell=True, cwd=csimDir)

        # Run EQ1 -nt 8 -v
        subp.call(f'gmx grompp -c lipids-water-em.gro -r lipids-water-em.gro -f {resources}/martini_v2.x_new-rf-eq1.mdp -p system.top -o topol.tpr -n index.ndx -maxwarn 5', shell=True, cwd=csimDir)
        subp.call('gmx mdrun -nt 1 -c lipids-water-eq1.gro >> mdrun-eq1.log 2>&1', shell=True, cwd=csimDir)

        # Run EQ2  -nt 8 -v
        subp.call(f'gmx grompp -c lipids-water-eq1.gro -r lipids-water-eq1.gro -f {resources}/martini_v2.x_new-rf-eq2.mdp -p system.top -o topol.tpr -n index.ndx -maxwarn 5', shell=True, cwd=csimDir)
        subp.call('gmx mdrun -nt 1 -c lipids-water-eq2.gro >> mdrun-eq2.log 2>&1', shell=True, cwd=csimDir)

        # Run EQ3  -nt 8 -v
        subp.call(f'gmx grompp -c lipids-water-eq2.gro -r lipids-water-eq2.gro -f {resources}/martini_v2.x_new-rf-eq3.mdp -p system.top -o topol.tpr -n index.ndx -maxwarn 5', shell=True, cwd=csimDir)
        subp.call('gmx mdrun -nt 1 -c lipids-water-eq3.gro >> mdrun-eq3.log 2>&1', shell=True, cwd=csimDir)

        # Run EQ4  -nt 8  -nt 4
        subp.call(f'gmx grompp -c lipids-water-eq3.gro -r lipids-water-eq3.gro -f {resources}/martini_v2.x_new-rf-eq4.mdp -p system.top -o topol.tpr -n index.ndx -maxwarn 5', shell=True, cwd=csimDir)
        subp.call('gmx mdrun -nt 1 -c lipids-water-eq4.gro >> mdrun-eq4.log 2>&1', shell=True, cwd=csimDir)

    # Set up production
    subp.call(f'gmx grompp -c lipids-water-eq4.gro -r lipids-water-eq4.gro -f martini_v2.x_new-rf-prod.mdp -p system.top -o topol.tpr -n index.ndx -maxwarn 5', shell=True, cwd=csimDir)
    #subp.call('gmx mdrun -nt 8 -v -maxh 0.02', shell=True, cwd=csimDir)
    # 6316.905 with -nt 8 on my laptop ~4h to finish 1us

    # Cleanup temp files and eq restart files
    subp.call('rm \#*', shell=True, cwd=csimDir)
    subp.call('rm dd_dump_err*.pdb', shell=True, cwd=csimDir)
    subp.call('rm step*.pdb', shell=True, cwd=csimDir)
    subp.call('rm state.cpt state_prev.cpt traj_comp.xtc traj.trr md.log ener.edr', shell=True, cwd=csimDir) # clean last eq run

    # Finished setup
    if os.path.isfile(os.path.join(csimDir, 'lipids-water-eq4.gro')) == True:
        subp.call('echo "setup done" > FLAG_SETUP_DONE', shell=True, cwd=csimDir)
    else:
        subp.call('echo "setup error" > FLAG_SETUP_ERROR', shell=True, cwd=csimDir)
    print(f" -> setup {csim} - ended in {timer() - start} s", flush = True)  # Setup takes ~44 sec on my laptop per sim
    return csim

def setupSims(cDir, lipidTableFile, tableData, options):
    baseLipid = options["-baseLipid"].value
    print(f"Setup simulations in {cDir} for lipid in {lipidTableFile} with baselipid {baseLipid}", flush = True)

    setsimArray = []
    for cline in tableData:
        cData = cline # .csv reader already seprated and parased
        initNote = cData.pop(0)

        if cData[0]!='x':  # only setup sims marked with x
            continue

        # Make name of current sim
        string_CBT = "noCBT"
        use_CBT = options["-CBT"].value == "yes"
        if use_CBT:
            string_CBT = "CBT"
        csim = f"M3_{string_CBT}_{baseLipid}_{cData[4]}"
        csimDir = os.path.join(cDir, csim)

        #os.chdir(cDir) # reset to inital dir
        if os.path.isdir(csimDir) == True:
            # csim alreaddy exist
            if os.path.isfile(os.path.join(csim, "FLAG_SETUP_DONE")) == True:
                print(f"-setup {csim} - sim found in {cDir}, sim is done skip setup", flush = True)
                continue
            else:
                print(f"-setup {csim} - ERROR - dir found in {cDir} but not done, skipped here", flush = True)
                continue
        setsimArray.append([csim, csimDir, options, cData]) # add to list as we want to setup this sim

    #eqSim(setsimArray[0]) # run first one for testing

    #for sims in setsimArray:
    #    print(f" -> HII run {sims[0]}", flush = True)
    #print(f" -> HII cpu count {mp.cpu_count()}", flush = True)
    thredcount = 36 # 5 # mp.cpu_count()
    pool = mp.Pool(min(thredcount, len(setsimArray))) # number of workers
    print(f"-Setup {len(setsimArray)} simulations using {thredcount} CPUs", flush = True)
    #results = pool.map(printSim, setsimArray, chunksize=1)
    results = pool.map(eqSim, setsimArray, chunksize=1)
    print(f"-All done for {results}", flush = True)
    pool.close()
    pool.join()

# end setupSims

def setupCustomSims(cDir, options):
    baseLipid = options["-baseLipid"].value
    print(f"Setup custom simulations in {cDir}", flush = True)

    # Get CBT selection
    string_CBT = "noCBT"
    use_CBT = options["-CBT"].value == "yes"
    if use_CBT:
        string_CBT = "CBT"

    setsimArray = []
    for cLipids in [["POPC", "POPE"], ["POPC", "POPG"], ["POPC", "POPS"], ["POPE", "POPG"], ["POPE", "POPS"], ["POPG", "POPS"]]:
        for cRatio in [[2, 1], [1, 1], [1, 2]]:
            csim = f"M3_{string_CBT}_{cLipids[0]}_{cRatio[0]}_{cLipids[1]}_{cRatio[1]}"
            csimDir = os.path.join(cDir, csim)

            if os.path.isdir(csimDir) == True:
                # csim alreaddy exist
                if os.path.isfile(os.path.join(csim, "FLAG_SETUP_DONE")) == True:
                    print(f"-setup {csim} - sim found in {cDir}, sim is done skip setup", flush = True)
                    continue
                else:
                    print(f"-setup {csim} - ERROR - dir found in {cDir} but not done, skipped here", flush = True)
                    continue
            setsimArray.append([csim, csimDir, options, cLipids+cRatio])

    #eqSim(setsimArray[0]) # run first one for testing
    #for sims in setsimArray:
    #    print(f" -> HII run {sims[0]}", flush = True)
    #    eqSim(sims) # run first one for testing

    thredcount = 36 # 5 # mp.cpu_count()
    pool = mp.Pool(min(thredcount, len(setsimArray))) # number of workers
    print(f"-Setup {len(setsimArray)} simulations using {thredcount} CPUs", flush = True)
    #results = pool.map(printSim, setsimArray, chunksize=1)
    results = pool.map(eqSim, setsimArray, chunksize=1)
    print(f"-All done for {results}", flush = True)
    pool.close()
    pool.join()

# end setupCustomSims

def runWorker(sim):
    csim, csimDir, cData = sim
    print(f" -> Running {csim}, full dir {csimDir}", flush = True)
    subp.call('gmx mdrun -nt 2 -cpi >> mdrun.log 2>&1', shell=True, cwd=csimDir)
    #subp.call('gmx mdrun -nt 1 -cpi -maxh 0.02 >> mdrun.log 2>&1', shell=True, cwd=csimDir)
    subp.call('rm \#*', shell=True, cwd=csimDir) # cleanup backup file from eq or previus rounds
    subp.call('echo "run done" > FLAG_RUN_DONE', shell=True, cwd=csimDir) # mark sim as done
    return csim

def runSims(cDir, runner_id, runner_count, options):
    print(f"Run unfinished simulations in {cDir}, runner {runner_id} of {runner_count}", flush = True)

    setsimArray = []
    simsDirs = 0
    simsCount = 0
    for cfile in os.listdir(cDir):
        csimDir = os.path.join(cDir,cfile)
        if cfile.startswith("M3_") and os.path.isdir(csimDir):
            # is sim folder check status
            simsDirs += 1
            if os.path.isfile(os.path.join(csimDir, "FLAG_SETUP_DONE")) == False:
                print(f"-run {cfile} - WARNING simulation setup not finished run will skip", flush = True)
            elif os.path.isfile(os.path.join(csimDir, "FLAG_RUN_DONE")) == True:
                print(f"-run {cfile} - alreaddy finished", flush = True)
            else:
                # add to run list
                simsCount += 1
                if (simsCount % runner_count) == runner_id:
                    setsimArray.append([cfile, csimDir, []])

    thredcount = 16 # mp.cpu_count()
    print(f"-Run {len(setsimArray)} simulations using {thredcount} CPUs of total {simsCount} - found {simsDirs-simsCount} simulations that are not readdy or done", flush = True)
    if len(setsimArray) == 0:
        print(f"-All done no simulation found to run", flush = True)
    else:
        pool = mp.Pool(min(thredcount, len(setsimArray))) # number of workers
        #results = pool.map(printSim, setsimArray, chunksize=1)
        results = pool.map(runWorker, setsimArray, chunksize=1)
        print(f"-All done running {len(setsimArray)} simulations: {results}", flush = True)
        pool.close()
        pool.join()
# end runSims

def analysisWorker(sim):
    csim, csimDir, cData = sim
    print(f" -> Analysing {csim}, full dir {csimDir}", flush = True)
    # Note part of analysis uses x2 CPUs
    subp.call(f'python3 /p/lustre1/helgi/martini3/tools/analysis_pipeline/analysis_lipids.py {csimDir}/ >> {csimDir}/analysisLog.log 2>&1', shell=True) # Change to relative path and make varaibale
    #print(f'python3 /p/lustre1/helgi/martini3/tools/analysis_pipeline/analysis_lipids.py {csimDir}')
    subp.call('rm \#*', shell=True, cwd=csimDir) # cleanup backup file from eq or previus rounds
    subp.call('echo "run done" > FLAG_ANALYSES_DONE', shell=True, cwd=csimDir) # mark sim as done
    return csim

def analyseSims(cDir, options):
    print(f"Analyse all simulations that have not alreaddy been analysed - Note delete flag if you want to rerurn", flush = True)

    setsimArray = []
    simsDirs = 0
    simsCount = 0
    for cfile in os.listdir(cDir):
        csimDir = os.path.join(cDir,cfile)
        if cfile.startswith("M3_") and os.path.isdir(csimDir):
            # is sim folder check status
            simsDirs += 1
            if os.path.isfile(os.path.join(csimDir, "FLAG_SETUP_DONE")) == False:
                print(f"-run {cfile} - WARNING simulation setup not finished will skip analysis", flush = True)
            elif os.path.isfile(os.path.join(csimDir, "FLAG_RUN_DONE")) == False:
                print(f"-run {cfile} - WARNING simulation run not finished will skip analysis", flush = True)
            elif os.path.isfile(os.path.join(csimDir, "FLAG_ANALYSES_DONE")) == True:
                print(f"-run {cfile} - analysis alreaddy finished", flush = True)
            else:
                # add to run list
                simsCount += 1
                setsimArray.append([cfile, csimDir, []])

    thredcount = 16 # mp.cpu_count()
    print(f"-Do analysis for {len(setsimArray)} simulations using {thredcount} threads", flush = True)
    if len(setsimArray) == 0:
        print(f"-All done no simulation found to analyse", flush = True)
    else:
        if thredcount == 1:
            for sim in setsimArray:
                analysisWorker(sim)
            print(f"-All done analysing {len(setsimArray)} simulations: {setsimArray}", flush = True)
        else:
            pool = mp.Pool(min(thredcount, len(setsimArray))) # number of workers
            #results = pool.map(printSim, setsimArray, chunksize=1)
            results = pool.map(analysisWorker, setsimArray, chunksize=1)
            print(f"-All done analysing {len(setsimArray)} simulations: {results}", flush = True)
            pool.close()
            pool.join()
# end analyseSims

if __name__ == "__main__":
    # Get arguments
    args = sys.argv[1:]

    # Print help
    if '-h' in args or '--help' in args:
        print("\n", __file__)
        print(desc)
        for thing in options:
            print(type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing)
        print()
        sys.exit()

    # Convert the option list to a dictionary, discarding all comments
    options = dict([i for i in options if not type(i) == str])
    # Process the command line
    while args:
        ar = args.pop(0)
        options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])


    # Read in name of .csv file with lipid input definitons
    lipidTableFile = options["-table"].value
    print(f"Reading lipid input definitions from {lipidTableFile}", flush = True)
    cDir = os.getcwd()
    tableFile = open(f"{lipidTableFile}","r")
    tableData = csv.reader(tableFile, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)

    singleLipiditpFolder = options["-ofolder"].value
    use_CBT = options["-CBT"].value == "yes"
    baseLipid = options["-baseLipid"].value
    print(f" running itp's in {lipidTableFile}, output itp folder {singleLipiditpFolder}, CBT {use_CBT}, baseLipid {baseLipid}", flush = True)

    if options["-makeItps"].value == "yes":
        # Make itp files for all lipids in table (with x or y in auto column)
        print(f"Make all lipid itp's in {lipidTableFile}", flush = True)
        makeitps(tableData, options)

    if options["-setupSims"].value == "yes":
        # Run simulation for all lipids in current directory (Note not in table but all readdy "same type" simulations)
        print(f"Setup simulations for lipids in {lipidTableFile} using base {baseLipid}", flush = True)
        setupSims(cDir, lipidTableFile, tableData, options)

    if options["-runSims"].value == "yes":
        runner_id = int(options["-sRunId"].value)
        runner_count = int(options["-sRunC"].value)
        print(f"Run simulations in current folder starting with M3 running {runner_count} sims starting with {runner_id}", flush = True)
        runSims(cDir, runner_id, runner_count, options)

    if options["-analyseSims"].value == "yes":
        # Run analysis script for all simulation in folder
        print(f"Analyse all simulations in current folder starting with M3", flush = True)
        analyseSims(cDir, options)

    if options["-setupCustomSims"].value == "yes":
        print(f"Setup custom simulations", flush = True)
        setupCustomSims(cDir, options)
