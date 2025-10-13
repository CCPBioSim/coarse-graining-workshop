#!/usr/bin/env python3

"""lipid-itp-generator-Martini3 creates a customized Martini lipid topologies for Martini 3 lipids, use lipid-itp-generator-Martini3.py -h for description"""

__author__  = "Helgi I. Ingolfsson, Kasper B. Pedersen, and Tsjerk A. Wassenaar"
__status__  = "Development"
__version__ = "M3.l01"
__email__   = "ingolfsson1@llnl.gov"

import sys,math

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
This script creates a customized Martini 3 lipid topology based on the head, linker and
tail specification strings provided. The tails, linkers and simple heagroup parameters follow
the main optimized (v2.0) Martini 3 lipid definitions as described in:
  K.B. Pedersen et al., The Martini 3 Lipidome: Expanded and Refined Parameters Improve Lipid Phase Behavior*, 
  ACS Central Science, 2025. doi: 10.1021/acscentsci.5c00755
The improved Martini 3 parameters for PI and PIPs headgroups are as specified in:
  L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo.
  Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the Martini 3
  Coarse-Grain Force Field, JCTC, 2021. doi:10.1021/acs.jctc.1c00615
All parameters are based on the Martini3 force field:
  P.C.T. Souza et al. Martini 3: a general purpose force field for coarse-grained molecular dynamics,
  Nat. Methods; 2021. doi: 10.1038/s41592-021-01098-3
This scripts is based on the Martini 2 lipid topology builder script lipid-martini-itp-v06.py
described in:
  T.A. Wassenaar, H.I. Ingolfsson, R.A. Bockmann, D.P. Tieleman, S.J. Marrink. Computational
  lipidomics with insane: a versatile tool for generating custom membranes for molecular simulations.
  JCTC, 150410125128004, 2015. doi:10.1021/acs.jctc.5b00209

For increased clarity and easy of modififcation the majority of the lipdi paramters are defined as  
variables that are defined in martini_v3.0.0_ffbonded_v2.itp which needs to be imported in .top file. 
  #include "xPATHx/martini_v3.0.0_ffbonded_v2.itp"  

WARNING:
  This script can generate topologies for numerous lipids many of which are untested and/or unrealistic,
  please use with discretion. 


The lipid descriptions supported are as follows:

Heads (-alhead):
  Please provide a list of lipid head beads. The left most bead will be on top, they are
  connected in a sequence from left to right and the right most bead is connected to the
  first bead in the linker.

  head bead types supported:
    C = NC3 = Choline         - bead Q1,  charge +1
    E = NH3 = Ethanolamine    - bead Q4p, charge +1
    G = GL0 = Glycerol        - bead P4r, charge  0
    S = CNO = Serine          - bead P6,  charge  0
    P = PO4 = Phosphate       - bead Q5,  charge -1
    O = PO4 = Phosphate       - bead Q5,  charge -2
  COH = COH = DAG capp        - bead TN6, charge  0
   PI = phosphatidylinositol  - x5 beads, charge -1
   P1 = PIP_1(3)              - x6 beads, charge -3
   P2 = PIP_2(3,4)            - x7 beads, charge -4
   P3 = PIP_3(3,4,5)          - x8 beads, charge -5
   P4 = PIP_1(4)              - x6 beads, charge -3
   P5 = PIP_1(5)              - x6 beads, charge -3
   P6 = PIP_2(4,5)            - x7 beads, charge -4
   P7 = PIP_2(3,4)            - x7 beads, charge -4

  Examples of lipid heads:
    "C P" -> 'NC3 PO4' - PC - PhosphatidylCholine
    "E P" -> 'NH3 PO4' - PE - PhosphatidylEthanolamine
    "G P" -> 'GL0 PO4' - PG - PhosphatidylGlycerol
    "S P" -> 'CNO PO4' - PS - PhosphatidylSerine
    "P"   -> 'PO4 ---' - PA - Phosphatidic acid (charge -1, use "O" for charge -2)
    "O"   -> 'PO4 ---' - PA - Phosphatidic acid, not protonated (charge -2)
    "COH" -> 'COH ---' - DG - Capping bead for top of diacylglycerols (DAG)

Linkers (-allink):
  Currently lipids with Glycerol linker (ester, ether and plasmalogens) and sphingosine 
  backbone are supported.

  linker beads types supported:
    G = Glycerol backbone    - normally G G -> GL1 GL2 with bead type SN4a
    L = Plasmalogen glycerol - normally G L -> GL1 PL2 with bead types SN4a and SN4as 
    E = Ether glycerol       - normally E E -> ET1 ET2 with bead type SN3a
    A = Sphingosine backbone - normally A1 A2 -> OH1 AM2 with bead type SP1 and SP2

  Examples of lipid linkers:
    "G G" -> 'GL1 GL2' - A glycerol linker
    "A A" -> 'OH1 AM2' - A sphingosine backbone

Tails (-altail):
  One lipid tail definition should be provided for each linker, separated with a space;
  extra spaces are ignored. Each tail can have an arbitrary number of tail beads.

  tail bead types supported:
    C = straight chain, corresponding roughly to a linear combination of x4 CH2 groups 
        (CH2-CH2-CH2-CH2). Represented with a C1 bead.
    c = short straight chain, corresponding roughly to a linear combination of x2 CH2
        groups (CH2-CH2). Represented with a SC1 bead.
    D = chain with double bond. Corresponding roughly to a linear combination of x2 CH2 
        and x2 CH groups (CH2-CH=CH-CH2), where the double bound is a cis bond. 
        Represented with a C4h bead.
    F = Chain with more than one double bond (normally 1.5). Represented with a C5h bead.
    T = Same as D except with a trans double bond. Represented with a C4h bead.
    t = Same as T except small. Represented with a SC4h bead.

  Examples of tails:
    Saturated tails:
    "cC     cC    " - C8:0  - diOctanoyl - DT 
    "CC     CC    " - C10:0 - diDecanoyl - DJ
    "cCCC   cCCC  " - C16:0 - diPalmitic acid - DP  
    "CCCC   CCCC  " - C18:0 - diStearoyl - DS
    "CCCCCC CCCCCC" - C26:0 - diHexacosanoyl - DC 
    Unsaturated tails:
    "CDC    CDC   " - C14:1(9c) - diMyristoleoyl - DR
    "CDCC   CDCC  " - C18:1(9c) - diOleoyl - DO 
    "CDDC   CDDC  " - C18:2(9c,12c) - diLinoleoyl - DL 
    "cFFDC  cFFDC " - C20:4(5c,8c,11c,14c) diArachidonoyl - DA 
    Mixed tails:
    "CDCC   cCCC  " - C16:0/C18:1(9c) - PO
    "CDDC   cCCC  " - C16:0/C18:2(9c,12c) - PL
    "DFFDD  cCCC  " - C16:0/C22:6(4c,7c,10c,13c,16c,19c) - PD
    Trans tails:
    "CTCC   CTCC  " - C18:1(9t) - diElaidoyl
    "tCCC   cCCC  " - palmytoyl sphingomyeline tail (AM1 contains the start of the tail)

    NOTE: the first tail (tail A) is connected to linker 1 closer to head (this is sn-2 for GLY linker lipids), which is reverse order
    compared to how regular lipid names are written. The second tail is tail B (for GLY linker lipids this is sn-1)

Use e.g.:
  ./lipid-itp-generator-Martini3-01.py -alhead 'C P' -allink 'G G' -altail "CDCC cCCC" -alname POPC -o POPC-lipid.itp
  ./lipid-itp-generator-Martini3-01.py -alhead 'E P' -allink 'G G' -altail "CDCC cCCC" -alname POPE -o POPE-lipid.itp
  ./lipid-itp-generator-Martini3-01.py -alhead 'P6' -allink 'G G' -altail "cFFDC CCCC" -alname SAP6 -o SA-PIP2-4-5-lipid.itp
"""

# Options
options = [
"""
Options:""",
("-o",       Option(str,    1,        "Martini-lipid.itp", "Output speciffic Martini lipid topology")),
("-alname",  Option(str,    1,        "POPC", "Four letter lipid name")),
("-alhead",  Option(str,    1,        "C P", "Lipid heads, see description")),
("-allink",  Option(str,    1,        "G G", "Lipid linkers, see description")),
("-altail",  Option(str,    1,        "CDCC cCCC", "Lipid tails, see description")),
("-name",    Option(str,    1,        "POPC", "A common name of the lipid, only use in comments")),
("-desc",    Option(str,    1,        "This is a ...", "A general description of what the FF is / represents, only use in comments")),
("-keyw",    Option(str,    1,        "", "List of keywords, only use in comments")),
("-parm",    Option(str,    1,        "Was modeled on ...", "How the FF was parameterized, only use in comments")),
("-refs",    Option(str,    1,        "", "List of references for the FF, only use in comments")),
("-crea",    Option(str,    1,        "", "FF created on, only use in comments")),
("-auth",    Option(str,    1,        "", "FF author, only use in comments")),
("-modi",    Option(str,    1,        "", "List of modifications to the FF, only use in comments")),
("-area",    Option(str,    1,        "", "Reference area per lipid, only use in comments")),
("-warn",    Option(str,    1,        "", "Warning(s)/Note(s) for the FF, only use in comments"))
          ]

# Define supported lipid head beads
# Lists all supported head bead types. One letter name mapped to type, atom name and charge
headMapp = {
    # beadtype, beadname, charge, ffbonded name
    "C":  ['Q1',   'NC3', '1.0',  'NC3'],  # NC3 = Choline
    "E":  ['Q4p',  'NH3', '1.0',  'NH3'],  # NH3 = Ethanolamine
    "G":  ['P4r',  'GL0', '0.0',  'GL0'],  # GL0 = Glycerol (used for PG lipid head group)
    "S":  ['P6',   'CNO', '0.0',  'CNO'],  # CNO = Serine
  "PS1":  ['SP2q', 'PS1', '-0.6', 'PS1'],  # PS1 is bead one of two for an alt. in dev. PS, represents a COO group
  "PS2":  ['Q4p',  'PS2', '0.6',  'PS2'],  # PS2 is bead two of two for an alt. in dev. PS, represents a NH3 group
    "P":  ['Q5',   'PO4', '-1.0', 'PO4'],  # PO4 = Phosphate
    "O":  ['Q5',   'PO4', '-2.0', 'PO4'],  # PO4 = Phosphate (one bond x2 charges can be used e.g. when making unprotonated PA lipids)
  "COH":  ['TN6',  'COH', '0.0',  'COH'],  # COH = Capping bead for top of diacylglycerols and ceramides
    }

linkerMapp = {
    # beadtype, beadname, charge, ffbonded name
    "G":   ['SN4a', 'GL?',  '0', 'GL' ],  # G = glycerol, normally GL1 and GL2
    "L":   ['SN4ah','PL?',  '0', 'ET' ],  # L = glycerol with alkenyl-acyl normally GL1 and PL2 for plasmalogen lipids 
    "E":   ['SN3a', 'ET?',  '0', 'ET' ],  # E = glycerol linked to ether lipid tails normally ET1 and ET2 for ether lipids 
    "A1":  ['SP1',  'OH1',  '0', 'OH1'],  # A = Shingosine A1 is mapped to OH1 (only supported in A1 A2 pair)
    "A2":  ['SP2',  'AM2',  '0', 'AM2'],  # A = Shingosine A2 is mapped to AM2 (only supported in A1 A2 pair)
    }

tailMapp = {
    # beadtype, beadname, charge, ffbonded bondname and angle/dih name
    "C":  ['C1',   'C??', '0',  'C1',  'C1'],  # C = straight chain
    "c":  ['SC1',  'C??', '0', 'SC1',  'C1'],  # c = short straight chain
    "D":  ['C4h',  'D??', '0',  'C4',  'C4'],  # D = chain with double bond
    "F":  ['C5h',  'D??', '0',  'C4',  'C4'],  # F = chain with more than one double bond (normally 1.5)
    "T":  ['C4h',  'T??', '0',  'C4',  'C1'],  # T = chain with trans double bond - normal sized bead - only used in shingosine top bead (SM and CER)
    "t":  ['SC4h', 'T??', '0', 'SC4',  'C1'],  # t = chain with trans double bond - small sized bead - only used in shingosine top bead (SM and CER)
    }

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

# Get ouput .itp file name
itpFileName  = options["-o"].value

# Get lipid description
lipidHead  = options["-alhead"].value
lipidLinker  = options["-allink"].value
lipidTail  = options["-altail"].value
if lipidLinker==None or lipidLinker==None or lipidTail==None:
    print("You have to provide a header, linker and tail lipid description, if one should be missing provide an empty string", file=sys.stderr)
    sys.exit()
lipidName = options["-alname"].value

lipidCommonName = options["-name"].value
lipidDesc = options["-desc"].value
lipidParm = options["-parm"].value
if lipidCommonName==None or lipidDesc==None or lipidParm==None:
    print("You have to provide a common name, description and list how the FF was parameterized.", file=sys.stderr)
    sys.exit()
lCharge = 0  # Update when adding charged beads

progString = "The Martini lipid itp generator version " + __version__ + "  Args are: -alname %s -alhead '%s' -allink '%s' -altail '%s'" % (lipidName, lipidHead, lipidLinker, lipidTail)
print(f"{progString} to file {itpFileName}")

headsArray = lipidHead.split()
linkersArray = lipidLinker.split()
linkersIndex = []
tailsArray = lipidTail.split()
if len(tailsArray)>len(linkersArray):
    print("A linker definition has to be provided for each tail", file=sys.stderr)
    sys.exit()


drawingArray = [[]]
drawingHeadsize = 0
bondsArray = []
anglesArray = []
beadArray = []
dihedralsArray = []
constraintsArray = []
exclusionsArray = []

# Set linker type - use for some bond and angle names
if (linkersArray[0] == "G"): # first linker is GLY
    if (linkersArray[1] == "L"): # second linker is L so this is an plasmalogen linker
        # Note G L - used for plasmalogens
        linkerType = "plasm"
        lipidType = "def" # _def is default parameters (not specifid to lipid type) 
    else: # assumign normal G G
        linkerType = "glyc"
        lipidType = "def" # _def is default parameters (not specifid to lipid type) 
elif (linkersArray[0] == "A1"): # SM for sphingomyelin and ceramide
    linkerType = "sm"
    lipidType = "sm"
elif (linkersArray[0] == "E"): # Ether lipid linker
    # Note e.g. E E - used for ether lipids 
    linkerType = "ether"
    lipidType = "def" # _def is default parameters (not specifid to lipid type) 
else: # A linker
    print(f"Linker type '{linkersArray[0]}' not supported", file=sys.stderr)
    sys.exit()


# If speciall head insert now all beads, bonds, angles, dihedrals, constreints etc
index = 1
if len(headsArray)>0 and headsArray[0]=='PI': # Add PI head
    # PI headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName, 'C1', 1, 0,    ''])
    beadArray.append([2,'SP4r', 1, lipidName, 'C2', 2, 0,    ''])
    beadArray.append([3,'SP4r', 1, lipidName, 'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName, 'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5, 'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    index += 5
    lCharge += -1.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.3720', '30000', ''])
    bondsArray.append([1, 3, 1, '0.3696', '30000', ''])
    bondsArray.append([2, 3, 1, '0.4044', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.33 ', ' 5000', 'C1 PO4'])
    bondsArray.append([5, 6, 1, '0.368', '2250', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 7, 1, '0.518', ' 600', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,  '100.0', '15.0', 'C2 C1 PO4'])
    anglesArray.append([1, 5, 6, 10,  '100.0', ' 4.0', 'C1 PO4 GL1'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -148.0,  95.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 6,  1,   180.0,   2.0,  2,  'C3 C2 C1 PO4'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.3720', 'C1 C2   Expanded 20% to account for SASA'])
    constraintsArray.append([1, 3, '0.3696', 'C1 C3'])
    constraintsArray.append([2, 3, '0.4044', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P1': # Add PIP_1(3) head
    # PIP_1(3) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2,'SP4r', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3, 'SP1', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,   'D', 1, lipidName,  'P3', 6, -2.0,  ''])
    index += 6
    lCharge += -3.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.305', '30000', ''])
    bondsArray.append([1, 3, 1, '0.302', '30000', ''])
    bondsArray.append([2, 3, 1, '0.324', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 3, 1, '0.303',  '10000', 'P3 C3'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.33 ',  '6000', 'C1 PO4'])
    bondsArray.append([5, 7, 1, '0.36 ',  '1800', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 8, 1, '0.515',  '1000', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,  '106.0',  '35.0', 'C2 C1 PO4'])
    anglesArray.append([1, 5, 7, 10,  '110.0',  '10.0', 'C1 PO4 GL1'])
    anglesArray.append([6, 3, 1, 10,  ' 99.0',  '82.5', 'P3 C3 C1'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -142.0,  125.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 7,  1,   180.0,   3.0,   2,  'C3 C2 C1 PO4'])
    dihedralsArray.append([6, 3, 1, 2,  2,  -167.0,  110.0, '',  'P3 C3 C1 C2'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.305', 'C1 C2'])
    constraintsArray.append([1, 3, '0.302', 'C1 C3'])
    constraintsArray.append([2, 3, '0.324', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P2': # Add PIP_2(3,4) head
    # PIP_2(3,4) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2,'SP4r', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3,'SN4a', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,   'D', 1, lipidName,  'P3', 6, -1.5,  ''])
    beadArray.append([7,   'D', 1, lipidName,  'P4', 7, -1.5,  ''])
    index += 7
    lCharge += -4.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.3075', '30000', ''])
    bondsArray.append([1, 3, 1, '0.2845', '30000', ''])
    bondsArray.append([2, 3, 1, '0.3030', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 3, 1, '0.340',   '9000', 'P3 C3'])
    bondsArray.append([6, 1, 1, '0.435',   '1750', 'P3 C1'])
    bondsArray.append([7, 3, 1, '0.330',   '9000', 'P4 C3'])
    bondsArray.append([6, 7, 1, '0.525',   ' 800', 'P3 P4'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.34 ',  '5000', 'C1 PO4'])
    bondsArray.append([5, 8, 1, '0.363',  '1750', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 9, 1, '0.520',  '1000', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,  '112.0',  '70.0', 'C2 C1 PO4'])
    anglesArray.append([1, 5, 8, 10,  ' 90.0',  ' 6.0', 'C1 PO4 GL1'])
    anglesArray.append([6, 3, 1, 10,  ' 88.0',  '55.0', 'P3 C3 C1'])
    anglesArray.append([7, 3, 2, 10,  ' 83.0',  '67.5', 'P4 C3 C2'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -140.0,  120.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 8,  1,   180.0,    3.5,  2,  'C3 C2 C1 PO4'])
    dihedralsArray.append([6, 3, 1, 2,  2,  -166.0,   95.0, '',  'P3 C3 C1 C2'])
    dihedralsArray.append([7, 3, 2, 1,  2,  -170.0,  120.0, '',  'P4 C3 C2 C1'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.3075', 'C1 C2'])
    constraintsArray.append([1, 3, '0.2845', 'C1 C3'])
    constraintsArray.append([2, 3, '0.3030', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6 7'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P3': # Add PIP_3(3,4,5) head
    # PIP_3(3,4,5) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2, 'SP1', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3,'SN4a', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,  'Q5', 1, lipidName,  'P3', 6, -1.3,  ''])
    beadArray.append([7,   'D', 1, lipidName,  'P4', 7, -1.4,  ''])
    beadArray.append([8,  'Q5', 1, lipidName,  'P5', 8, -1.3,  ''])
    index += 8
    lCharge += -5.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.290', '30000', ''])
    bondsArray.append([1, 3, 1, '0.284', '30000', ''])
    bondsArray.append([2, 3, 1, '0.307', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 3, 1, '0.306',   '7500', 'P3 C3'])
    bondsArray.append([6, 1, 1, '0.505',   '1250', 'P3 C3'])
    bondsArray.append([7, 3, 1, '0.325',   '5500', 'P4 C3'])
    bondsArray.append([8, 2, 1, '0.395',   '3500', 'P5 C2'])
    bondsArray.append([8, 1, 1, '0.500',   '1250', 'P5 C1'])
    bondsArray.append([7, 8, 1, '0.505',    '800', 'P4 P5'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.341',  '5000', 'C1 PO4'])
    bondsArray.append([5, 9, 1, '0.36 ',  '2000', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5,10, 1, '0.51 ',  ' 900', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,  '105.0',  '17.5', 'C2 C1 PO4'])
    anglesArray.append([1, 5, 9, 10,  ' 90.0',  ' 3.5', 'C1 PO4 GL1'])
    anglesArray.append([6, 3, 1, 10,  '110.0',  '10.0', 'P3 C3 C1'])
    anglesArray.append([7, 3, 2, 10,  ' 95.0',  ' 2.0', 'P4 C3 C2'])
    anglesArray.append([8, 2, 3, 10,  '110.0',  '20.0', 'P5 C2 C3'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -180.0,    0.5, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 9,  1,   180.0,    1.5,  2,  'C2 C1 PO4 GL1'])
    dihedralsArray.append([6, 3, 1, 2,  2,   130.0,   10.0, '',  'P3 C3 C1 C2'])
    dihedralsArray.append([7, 3, 2, 1,  2,   130.0,   15.0, '',  'P4 C3 C2 C1'])
    dihedralsArray.append([8, 2, 3, 1,  2,   115.0,   17.5, '',  'P5 C2 C3 C1'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.290', 'C1 C2'])
    constraintsArray.append([1, 3, '0.284', 'C1 C3'])
    constraintsArray.append([2, 3, '0.307', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6 7 8'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P4': # Add PIP_1(4) head
    # PIP_1(4) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2,'SP4r', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3, 'SP1', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,   'D', 1, lipidName,  'P4', 6, -2.0,  ''])
    index += 6
    lCharge += -3.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.3037', '30000', ''])
    bondsArray.append([1, 3, 1, '0.2945', '30000', ''])
    bondsArray.append([2, 3, 1, '0.3370', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 3, 1, '0.325', '5250', 'P4 C3'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.335', '5500', 'C1 PO4'])
    bondsArray.append([5, 7, 1, '0.36 ', '2250', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 8, 1, '0.52 ', ' 800', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,   '98.0',  '27.5', 'C2  C1 PO4'])
    anglesArray.append([1, 5, 7, 10,   '90.0',  ' 3.0', 'C1 PO4 GL1'])
    anglesArray.append([6, 3, 2, 10,   '92.0',  '55.0', 'P4 C3 C2'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -145.0,  130.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 7,  1,   180.0,   2.0,   2,  'C3 C2 C1 PO4'])
    dihedralsArray.append([6, 3, 2, 1,  2,  -170.0,  110.0, '',  'P4 C3 C2 C1'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.3037', 'C1 C2'])
    constraintsArray.append([1, 3, '0.2945', 'C1 C3'])
    constraintsArray.append([2, 3, '0.3370', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P5': # Add PIP_1(5) head
    # PIP_1(5) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2, 'SP1', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3,'SP4r', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,   'D', 1, lipidName,  'P5', 6, -2.0,  ''])
    index += 6
    lCharge += -3.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.297', '30000', ''])
    bondsArray.append([1, 3, 1, '0.311', '30000', ''])
    bondsArray.append([2, 3, 1, '0.324', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 2, 1, '0.335', '10000', 'P5 C2'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.337', '5000', 'C1 PO4'])
    bondsArray.append([5, 7, 1, '0.354', '1750', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 8, 1, '0.52 ', ' 600', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,   '89.0',  '19.0', 'C2  C1 PO4'])
    anglesArray.append([1, 5, 7, 10,   '90.0',  ' 2.5', 'C1 PO4 GL1'])
    anglesArray.append([6, 2, 3, 10,   '82.5', '110.0', 'P5  C2 C3'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -150.0,  115.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 7,  1,   180.0,    2.0,  2,  'C3 C2 C1 PO4'])
    dihedralsArray.append([6, 2, 3, 1,  2,   165.0,  100.0, '',  'P3 C2 C3 C1'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.297', 'C1 C2'])
    constraintsArray.append([1, 3, '0.311', 'C1 C3'])
    constraintsArray.append([2, 3, '0.324', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P6': # Add PIP_2(4,5) head
    # PIP_2(4,5) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2, 'SP1', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3, 'SP1', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,   'D', 1, lipidName,  'P4', 6, -1.5,  ''])
    beadArray.append([7,   'D', 1, lipidName,  'P5', 7, -1.5,  ''])
    index += 7
    lCharge += -4.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.2770', '30000', ''])
    bondsArray.append([1, 3, 1, '0.2865', '30000', ''])
    bondsArray.append([2, 3, 1, '0.3360', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 3, 1, '0.325',  '4250', 'P4 C3'])
    bondsArray.append([7, 2, 1, '0.337',  '2500', 'P5 C2'])
    bondsArray.append([6, 7, 1, '0.525',  '1000', 'P4 P5'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.34 ',  '5000', 'C1 PO4'])
    bondsArray.append([5, 8, 1, '0.36 ',  '1700', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 9, 1, '0.507',  '1000', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,  '100.0',  '50.0', 'C2 C1 PO4'])
    anglesArray.append([1, 5, 8, 10,  ' 90.0',  ' 2.0', 'C1 PO4 GL1'])
    anglesArray.append([6, 3, 2, 10,  '101.0',  ' 7.0', 'P4 C3 C2'])
    anglesArray.append([7, 2, 3, 10,  ' 97.5',  '12.5', 'P5 C2 C3'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -140.0,  110.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 8,  1,   180.0,    0.5,  2,  'C3 C2 C1 PO4'])
    dihedralsArray.append([6, 3, 2, 1,  2,   177.0,   50.0, '',  'P3 C3 C1 C2'])
    dihedralsArray.append([7, 2, 3, 1,  2,   130.0,    9.0, '',  'P5 C2 C3 C1'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.2770', 'C1 C2'])
    constraintsArray.append([1, 3, '0.2865', 'C1 C3'])
    constraintsArray.append([2, 3, '0.3360', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6 7'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='P7': # Add PIP_2(3,5) head
    # PIP_2(3,5) headgroup based on Martini3 PI parameters by
    # L. Borges-Araujo, P.C.T. Souza, F. Fernandes and M.N. Melo,
    # Improved Parameterization of Phosphatidylinositide Lipid Headgroups for the
    # Martini 3 Coarse-Grain Force Field, J. Chem. Theory Comput., 2021,
    # doi:10.1021/acs.jctc.1c00615

    beadArray.append([1, 'SP1', 1, lipidName,  'C1', 1, 0,    ''])
    beadArray.append([2, 'SP1', 1, lipidName,  'C2', 2, 0,    ''])
    beadArray.append([3, 'SP1', 1, lipidName,  'C3', 3, 0,    ''])
    beadArray.append([4, 'TC4', 1, lipidName,  'C4', 4, 0, 0, 'Massless virtual bead'])
    beadArray.append([5,  'Q5', 1, lipidName, 'PO4', 5, -1.0,  ''])
    beadArray.append([6,   'D', 1, lipidName,  'P3', 6, -1.5,  ''])
    beadArray.append([7,   'D', 1, lipidName,  'P5', 7, -1.5,  ''])
    index += 7
    lCharge += -4.0 # Keep track of overall lipid charge

    bondsArray.append([-2, '#ifdef FLEXIBLE'])
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, 1, '0.2910', '30000', ''])
    bondsArray.append([1, 3, 1, '0.3050', '30000', ''])
    bondsArray.append([2, 3, 1, '0.3150', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([-1, 'Headgroup phosphates'])
    bondsArray.append([6, 3, 1, '0.325',  '5000', 'P3 C3'])
    bondsArray.append([6, 1, 1, '0.445',  '1500', 'P3 C1'])
    bondsArray.append([7, 2, 1, '0.345',  '4750', 'P5 C2'])
    bondsArray.append([-1, 'Phosphodiester and glycerol backbone'])
    bondsArray.append([1, 5, 1, '0.33 ',  '5000', 'C1 PO4'])
    bondsArray.append([5, 8, 1, '0.358',  '2000', 'PO4 GL1'])  # This links the head to the linker
    bondsArray.append([5, 9, 1, '0.518',  ' 900', 'PO4 GL2'])  # PIs have adjusted G-G links
    bondsArray.append([-1, 'Lipid tail (uses standard Martini tail rules)'])

    anglesArray.append([-1, 'Orient the headgroup'])
    anglesArray.append([2, 1, 5, 10,  '102.0',  '45.0', 'C2 C1 PO4'])
    anglesArray.append([1, 5, 8, 10,  ' 90.0',  ' 3.0', 'C1 PO4 GL1'])
    anglesArray.append([6, 3, 1, 10,  '100.0',  '20.0', 'P3 C3 C1'])
    anglesArray.append([7, 2, 3, 10,  ' 78.0',  '65.0', 'P5 C2 C3'])
    anglesArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    dihedralsArray.append([-1, 'Orient the headgroup'])
    dihedralsArray.append([3, 2, 1, 5,  2,  -142.0,  110.0, '',  'C3 C2 C1 PO4'])
    dihedralsArray.append([2, 1, 5, 8,  1,   180.0,    3.0,  2,  'C3 C2 C1 PO4'])
    dihedralsArray.append([6, 3, 1, 2,  2,  -162.0,   80.0, '',  'P3 C3 C1 C2'])
    dihedralsArray.append([7, 2, 3, 1,  2,   163.0,   80.0, '',  'P5 C2 C3 C1'])
    #dihedralsArray.append([-1, 'Tail part (uses standard Martini tail rules)'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.2910', 'C1 C2'])
    constraintsArray.append([1, 3, '0.3050', 'C1 C3'])
    constraintsArray.append([2, 3, '0.3150', 'C2 C3'])
    constraintsArray.append([-2, '#endif'])

    #Virtual site is just hardcoded through the constraints at the moment.
    #Not very elegant.
    constraintsArray.append([-2, ''])
    constraintsArray.append([-2, '[ virtual_sitesn ]'])
    constraintsArray.append([-2, '; site funct  constructing atom indices'])
    constraintsArray.append([-2, '   4     2     1 2 3'])

    exclusionsArray.append([-2, '  4 3 2 1 5 6 7'])
    exclusionsArray.append([-2, '  3 2 1'])
    exclusionsArray.append([-2, '  2 1'])

elif len(headsArray)>0 and headsArray[0]=='X' and len(linkersArray)==3 and linkersArray[0]=='G' and linkersArray[1]=='G' and linkersArray[2]=='G':
    # @TODO update to Martini3

    # Now add TGA glycerol (both head and linker - later no linker will be added). Has head as 'X' and linker as 'G G G'

    # This is from the TOG topology published in Vuorela et al.2010
    # WARNING this is always TAG so no head beads before or after - also no linker will be added
    beadArray.append([1, 'C1', 1, lipidName, 'GLY', 1, 0,    ''])
    beadArray.append([2, 'Na', 1, lipidName, 'GL1', 2, 0,    '; was called ES1'])
    beadArray.append([3, 'Na', 1, lipidName, 'GL2', 3, 0,    '; was called ES2'])
    beadArray.append([4, 'Na', 1, lipidName, 'GL3', 4, 0,    '; was called ES3'])
    index += 4
    linkersIndex.append(2) # Connects tail 1 - A
    linkersIndex.append(3) # Connects tail 2 - B
    linkersIndex.append(4) # Connects tail 3 - C
    beadArray.append([-1, 'Glycerol based on TOG topology from Vuorela et al. 2010'])
    beadArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    bondsArray.append([1, 2, defBlength, defBforce, ''])
    bondsArray.append([1, 3, defBlength, defBforce, ''])
    bondsArray.append([1, 4, defBlength, defBforce, ''])

    anglesArray.append([-1, 'Orient Glycerol'])
    anglesArray.append([-3, 2, 1, 4, ' 60.0', defAforce2, ''])
    anglesArray.append([-3, 2, 1, 3, '130.0', defAforce2, ''])
    anglesArray.append([-3, 3, 1, 4, '130.0', defAforce2, ''])
    anglesArray.append([-1, 'Orient Tails'])

elif len(headsArray)>0: # Else build head (works for simple PC, PE, PA, PG, PS, etc head groups)

    # Start head beads
    for i in range(len(headsArray)):
        beadArray.append([index, headMapp[headsArray[i]][0], 1, lipidName, headMapp[headsArray[i]][1], index, headMapp[headsArray[i]][2], ''])
        drawingArray[0].append(index)
        drawingHeadsize += 1
        lCharge += float(headMapp[headsArray[i]][2]) # Keep track of overall lipid charge
        if index > 1: # link head beads
            if headsArray[i-1] == 'PS1' and headsArray[i] == 'PS2':
                # Special case PS head group
                constraintsArray.append([index - 1, index, '0.248', 'PS1-PS2 PS head group constraint'])
                # Note used - All PS shoudl have [PS1 PS2 P] head group followed by linker add in PS dihedral
                # dihedralsArray.append([-3, index - 1, index, index + 1, index + 2, f"d_{headMapp[headsArray[i-1]][1]}_{headMapp[headsArray[i]][1]}_{headMapp[headsArray[i+1]][1]}_{linkerMapp[linkersArray[0]][3]}_cbt", ''])
                # and PS angle
                anglesArray.append([-3, index - 1, index, index + 1, f"a_{headMapp[headsArray[i-1]][1]}_{headMapp[headsArray[i]][1]}_{headMapp[headsArray[i+1]][1]}_{lipidType}", ''])
            else:
                # Normal head group bond
                bondsArray.append([-3, index - 1, index, f"b_{headMapp[headsArray[i-1]][1]}_{headMapp[headsArray[i]][1]}_def", ''])
        index += 1

    if len(linkersArray) > 0:
        # This links the head to the linker
        bondsArray.append([-3, index - 1, index, f"b_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[0]][3]}_{lipidType}", ''])
        if len(linkersArray) == 2 and len(headsArray) == 1 and headsArray[0] == "COH": # Add second bond to other linker bead case of CER lipids
            bondsArray.append([-3, index - 1, index + 1, f"b_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[1]][3]}_{lipidType}", ''])
        elif len(linkersArray) == 2 and lipidType == "sm": # Add second bond to other linker bead case of SM lipids
            bondsArray.append([-3, index - 1, index + 1, f"b_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[1]][3]}_{lipidType}", ''])
        elif len(linkersArray) == 2 and lipidType == "def": # Add second bond to other linker bead case of regular glycerolipids lipids
            bondsArray.append([-3, index - 1, index + 1, f"b_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[1]][3]}_{lipidType}_long", ''])

        if len(headsArray) > 1:
            # Add angle between head group (-2) to head group (-1) to linker to orient head group - only if at least x2 head group beads
            anglesArray.append([-3, index - 2, index - 1, index, f"a_{headMapp[headsArray[-2]][1]}_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[0]][3]}_{lipidType}", ''])
    else:
        print("Some type of head group to tail linker is needed", file=sys.stderr)
        sys.exit()

    if len(linkersArray) > 1:
        # Orient lipid head, add angles between linkers and head + head-first linker-first tail
        # Add angle between the last head bead and first x2 linkers
        if len(tailsArray[0]) == 1 and tailsArray[0][0]=='-':  # First tail is missing (like in LPC) keep <head - linker_1 - linker_2> straight
            #anglesArray.append([index - 1, index, index + 1, defAngle3, defAforce2, ''])
            print("@TODO lyso linker not supported yet", file=sys.stderr)
            sys.exit()
        elif len(headsArray) == 1 and headsArray[0] == "COH":
            # This is CER lipid then no head to linker linker angle
            pass
        else:  # All normal cases <head - linker_1 - linker_2> have a defAngle2
            # Angle removed in v17 as the x2 tails to first head now have triangle bonds so this angle not needed
            #anglesArray.append([-3, index - 1, index, index + 1, f"a_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[0]][3]}_{linkerMapp[linkersArray[1]][3]}_{linkerType}", ''])

            # no dihedral used - code here in chase something similar will be used later  
            #if len(headsArray) > 1:
            #    # Add dihedral between head group (-2) to head group (-1) to linker and other linker to orient head group - only if at least x2 head group beads
            #    dihedralsArray.append([-3, index - 2, index - 1, index, index + 1, f"d_{headMapp[headsArray[-2]][1]}_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[0]][3]}_{linkerMapp[linkersArray[1]][3]}_cbt", ''])
            pass

    # Add angle and dihedrals down to tail
    if len(linkersArray) > 1 and len(tailsArray[0]) == 1 and tailsArray[0][0]=='-':
        # Special case first tail is missing (like in LPC) in stead of aligning <last head, linker_1, first tail> allign <linker_1, linker_2, first tail>
        #anglesArray.append([index, index+1, index+len(linkersArray), defAngle3, defAforce2, ''])
        print("@TODO lyso linker not supported yet", file=sys.stderr)
        sys.exit()
    else:
        # Normal case, add angle between the last head bead, the first linker and the first tail
        # The first CER head-linker angle is named specially
        if lipidType == "sm" and headMapp[headsArray[-1]][1] == "COH":
            specialType = "cera"
        else:
            specialType = lipidType
        anglesArray.append([-3, index - 1, index, index+len(linkersArray), f"a_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[0]][3]}_C_{specialType}", ''])
        # no dihedral used - code here in chase something similar will be used later  
        #if len(headsArray) > 1:
        #    # Add dihedral between head group (-2) to head group (-1) to linker and first tail bead to orient head group - only if at least x2 head group beads
        #    dihedralsArray.append([-3, index - 2, index - 1, index, index + 2, f"d_{headMapp[headsArray[-2]][1]}_{headMapp[headsArray[-1]][1]}_{linkerMapp[linkersArray[0]][3]}_C_cbt", ''])

else: #If -alhead was empty len(headsArray)==0 then no head to add (DAG, CER etc)
    pass

# End Add heads
headIndex = index - 1 # To know what was the last headbead

# Add linker beads
for linkerBeadIndex in range(0, len(linkersArray)):
    cLinker = linkersArray[linkerBeadIndex]

    if (cLinker not in linkerMapp.keys()): # So far G, @TODO add A
        print(f"This script currently only supports {linkerMapp.keys()} linkers", file=sys.stderr)
        sys.exit()

    if len(headsArray)>0 and headsArray[0]=='X' and len(linkersArray)==3 and linkersArray[0]=='G' and linkersArray[1]=='G' and linkersArray[2]=='G':
        # This is TGA glycerol that was added above - don't add any linker
        break

    if (cLinker == "G" or cLinker == "E" or cLinker == "L"): # GLY (E and L are gly but for plasmalogen or ether lipids)
        if len(tailsArray[linkerBeadIndex]) == 1 and tailsArray[linkerBeadIndex][0]=='-':  # No tail given for current linker bead so free OH group Na changed to P1
            #beadArray.append([index, 'P1', 1, lipidName, "GL"+str(linkerBeadIndex+1), index, 0, ''])  # Example of this is in LPC (it as not tail A so GL1 changes to P1 type
            print("Not yet supported @TODO add", file=sys.stderr)
            sys.exit()
        elif len(headsArray)==0 and index == 1: # No head and this is the first Gly linker then change bead type to reflect the extra OH group.
            #beadArray.append([index, 'P1', 1, lipidName, "GL"+str(linkerBeadIndex+1), index, 0, '; bead type changed to P1 to account for the added OH group (mod. 2016.09.27)'])
            print("Not yet supported @TODO add", file=sys.stderr)
            sys.exit()
        else:  # For all other lipids
            beadArray.append([index, linkerMapp[cLinker][0], 1, lipidName, linkerMapp[cLinker][1].replace("?",str(linkerBeadIndex+1)), index, 0, ''])
    else: # A linker
        # Currently only expecting A1 or A2
        beadArray.append([index, linkerMapp[cLinker][0], 1, lipidName, linkerMapp[cLinker][1], index, 0, ''])

    if index > 1: # There have been beads before (don't do anything if this was the first linker bead and no head)
        if (index -1) > headIndex:  # Then this a linker - linker bond
            bondsArray.append([-3, index - 1, index, f"b_{linkerMapp[linkersArray[linkerBeadIndex-1]][3]}_{linkerMapp[cLinker][3]}_{linkerType}", ''])

    # Add to drawing
    if len(drawingArray) <= linkerBeadIndex:
        # new linker level add new line to drawing
        drawingArray.append([])
    drawingArray[linkerBeadIndex].append(index)

    #lCharge += 0 # Keep track of overall lipid charge, not needed as current linkers are all uncharged
    linkersIndex.append(index)
    index += 1
# End linkersArray loop

# Add linker angles and dihedrals @TODO extend this currently work only for x2 linker glycerol linked lipids with x2 tails
if len(linkersArray) == 2 and len(tailsArray) >= 2:
    # WARNING index is now +1 from last linker
    anglesArray.append([-3, linkersIndex[1], linkersIndex[0], index, f"a_{linkerMapp[linkersArray[1]][3]}_{linkerMapp[linkersArray[0]][3]}_C_{linkerType}", ''])
    anglesArray.append([-3, linkersIndex[0], linkersIndex[1], index+len(tailsArray[0]), f"a_{linkerMapp[linkersArray[0]][3]}_{linkerMapp[linkersArray[1]][3]}_C_{linkerType}", ''])
    # CBT dihedral not used - code here in chase something similar will be used later  
    #dihedralsArray.append([-3, linkersIndex[1], linkersIndex[0], index, index + 1, f"d_{linkerMapp[linkersArray[1]][3]}_{linkerMapp[linkersArray[0]][3]}_C_C_{linkerType}", ''])
    #dihedralsArray.append([-3, linkersIndex[0], linkersIndex[1], index+len(tailsArray[0]), index+len(tailsArray[0]) + 1, f"d_{linkerMapp[linkersArray[0]][3]}_{linkerMapp[linkersArray[1]][3]}_C_C_{linkerType}", ''])
    #dihedralsArray.append([-3, index, linkersIndex[0], linkersIndex[1], index+len(tailsArray[0]), f"d_C_{linkerMapp[linkersArray[0]][3]}_{linkerMapp[linkersArray[1]][3]}_C_{linkerType}", ''])
else:
    # @TODO add lyso and tryclyceroids
    print("This script currently only supports x2 linker lipids with x2 tails", file=sys.stderr)
    sys.exit()


# Add tail beads
indexToLetter = "A B C D E F G H I J K L M N".split()
for tailIndex in range(0, len(tailsArray)):
    cTail = tailsArray[tailIndex]
    cLinker = linkersArray[tailIndex] # assuming the same number of linkers and tails

    # If tail is empty (-) skip this section (like in LPC with no GL1 tail A)
    if len(cTail) == 1 and cTail[0]=='-':
        continue

    # Add bond from current tail to assosiated linker (to first tail bead)
    if tailMapp[cTail[0]][3][0] == 'S': # add 5long bond unless first tail bead is small
        bondsArray.append([-3, linkersIndex[tailIndex], index, f"b_{linkerMapp[cLinker][3]}_{tailMapp[cTail[0]][3]}_{linkerType}", ''])
    else:
        bondsArray.append([-3, linkersIndex[tailIndex], index, f"b_{linkerMapp[cLinker][3]}_{tailMapp[cTail[0]][3]}_{linkerType}_5long", ''])

    for cTailBeadIndex in range(0, len(cTail)):
        cTailBead = cTail[cTailBeadIndex]
        if (cTailBead not in tailMapp.keys()):
            print("Tail definition \"%s\" not recognized" % (cTailBead), file=sys.stderr)
            sys.exit()
        bType = tailMapp[cTailBead][0]
        atomName = tailMapp[cTailBead][1].replace("??", str(cTailBeadIndex+1) + indexToLetter[tailIndex])
        beadArray.append([index, bType, 1, lipidName, atomName, index, 0, ''])
        #lCharge += 0 # Keep track of overall lipid charge, not needed as current tails are all uncharged

        # Add bond between tail beads
        if cTailBeadIndex > 0: # linker to first tail alreaddy added
            if cTailBeadIndex == 1 and \
               tailMapp[cTail[cTailBeadIndex-1]][0][0] != 'S' and \
               tailMapp[cTail[cTailBeadIndex]][0][0] != 'S':
                # If first two tail beads are normal sized (not small) add long tail bond between them
                bondsArray.append([-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}_mid_5long", ''])
            else:
                #if len(tailsArray) >= 2 and \
                #   (cTailBeadIndex + 1) == len(cTail) and \
                #  cTail[cTailBeadIndex-1][0] == 'C' and \
                #   cTail[cTailBeadIndex][0] == 'C':
                #    # If last x2 beads are normal CC then add a longer end bond
                #    bondsArray.append([-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}_end", ''])
                # Simplified in v17 now all ends end with _end and defined in bond list b_C1_C1_end, b_C1_C4_end, b_C4_C1_end, b_C4_C4_end
                if len(tailsArray) >= 2 and (cTailBeadIndex + 1) == len(cTail):
                    bondsArray.append([-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}_end", ''])
                else:
                    # Add normal mid bond
                    bondsArray.append([-3, index-1, index, f"b_{tailMapp[cTail[cTailBeadIndex-1]][3]}_{tailMapp[cTail[cTailBeadIndex]][3]}_mid", ''])

        # Add angles to support the tails except not for the last tail (which can't have an angle)
        if (cTailBeadIndex + 1) < len(cTail):
            if cTailBeadIndex == 0: # top angle connecting to the tail linker [linker / head tail / head+1 tail]
                #if len(linkersArray)==3 and linkersArray[0]=='G' and linkersArray[1]=='G' and linkersArray[2]=='G': # WARNING ugly hack
                #    # this is only for TAG lipids, in the TOG lipid they are based on there is an extra glycerol bond 1 - linker bead - tail
                #    anglesArray.append([1, linkersIndex[tailIndex], index, cdefAngle, cdefAforce, ''])
                anglesArray.append([-3, linkersIndex[tailIndex], index, index + 1, f"a_{linkerMapp[cLinker][3]}_{tailMapp[cTail[cTailBeadIndex]][4]}_{tailMapp[cTail[cTailBeadIndex+1]][4]}_{linkerType}", ''])
            else: # regular angle connecting to the tail linker [current-1 tail / current tail / current+1 tail bead]
                # note, these angles are the came in glyco, sm and ceramides 
                anglesArray.append([-3, index - 1, index, index + 1, f"a_{tailMapp[cTail[cTailBeadIndex-1]][4]}_{tailMapp[cTail[cTailBeadIndex]][4]}_{tailMapp[cTail[cTailBeadIndex+1]][4]}_def", ''])
        # end angle stuff

        # Add dihedrals to support the tails except not for the last two tail (which can't have a dihedral)
        if (cTailBeadIndex + 2) < len(cTail):
            if cTailBeadIndex == 0: # top dihedrals connecting to the tail linker [linker / head tail / head+1 tail / head+2 tail]
                #dihedralsArray.append([-3, linkersIndex[tailIndex], index, index + 1, index + 2, f"d_{linkerMapp[cLinker][3]}_{tailMapp[cTail[cTailBeadIndex]][4]}_{tailMapp[cTail[cTailBeadIndex+1]][4]}_{tailMapp[cTail[cTailBeadIndex+2]][4]}_cbt", ''])
                # no dihedral used 
                pass
            else:
                # no dihedral used 
                # regular dihedrals connecting to the tail linker [current-1 tail / current tail / current+1 tail bead / current+2 tail bead]
                #dihedralsArray.append([linkersIndex[tailIndex], index, index + 1, f"d_{tailMapp[cTail[cTailBeadIndex-1]][4]}_{tailMapp[cTail[cTailBeadIndex]][4]}_{tailMapp[cTail[cTailBeadIndex+1]][4]}_{tailMapp[cTail[cTailBeadIndex+2]][4]}_cbt", ''])
                #dihedralsArray.append([-3, index - 1, index, index + 1, index + 2, f"d_{tailMapp[cTail[cTailBeadIndex-1]][4]}_{tailMapp[cTail[cTailBeadIndex]][4]}_{tailMapp[cTail[cTailBeadIndex+1]][4]}_{tailMapp[cTail[cTailBeadIndex+2]][4]}_cbt", ''])
                pass
        # end angle stuff

        # Add current bead to drawing @WARNINIG assumes linkers and tails add up (will fail for LPC lipids)
        drawingArray[tailIndex].append(index)

        index += 1
# End tailsArray loop


# Make .itp headder
itpFile = open(itpFileName,"w")
print(';;;;;; Martini lipid topology for ' + lipidCommonName + ', generated using:', file=itpFile)
print('; ' + progString, file=itpFile)
print('; WARNING: Autogenerated and not all lipid headgroup and tail combinations have been tested;', file=itpFile)
print(';          used with care and see K.B. Pedersen et al. 2024 for guidance.\n;', file=itpFile)
print('; Description:', file=itpFile)
print(';   ' + lipidDesc.replace('\\n','\n;  '), file=itpFile)
current = options["-keyw"].value
if current!=None and len(current) > 0:
    print(';@Keywords: '+current, file=itpFile)
print('; Parameterization:', file=itpFile)
print(';   ' + lipidParm.replace('\\n','\n;  '), file=itpFile)
current = options["-refs"].value
if current!=None and len(current) > 0:
    print('; Reference(s): ', file=itpFile)
    print(';   ' + current.replace('\\n','\n;  '), file=itpFile)
current = options["-crea"].value
if current!=None and len(current) > 0:
    print('; Created: ' + current, file=itpFile)
current = options["-auth"].value
if current!=None and len(current) > 0:
    print('; Author(s): ' + current, file=itpFile)
current = options["-modi"].value
if current!=None and len(current) > 0:
    print('; Modified:', file=itpFile)
    print(';   ' + current.replace('\\n','\n;  '), file=itpFile)
lArea = options["-area"].value
if lArea!=None and len(lArea) > 0:
    print('; Reference area per lipid: ' + lArea + ' nm^2', file=itpFile)
current = options["-warn"].value
if current!=None and len(current) > 0:
    print('; Warning(s)/Note(s):', file=itpFile)
    print(';   ' + current.replace('\\n','\n;  '), file=itpFile)
print(';', file=itpFile)
print('; Molecular topology and mapping of indices:', file=itpFile)
#print(';', file=itpFile)
for i in range(len(drawingArray)):
    cSringNames = ""
    if i > 0: # no head in this line
        cSringNames = cSringNames.rjust(drawingHeadsize * 4)
    for j in range(len(drawingArray[i])-1):
        cSringNames += f"{beadArray[drawingArray[i][j]-1][4]}-"
    cSringNames += f"{beadArray[drawingArray[i][-1]-1][4]}"
    print(f'; {cSringNames}', file=itpFile)
    if i+1 < len(drawingArray): # if not last row add a seperator line
        print(f'; {"|".rjust((drawingHeadsize * 4)+2)}', file=itpFile)
print(';', file=itpFile)
for i in range(len(drawingArray)):
    cSringIndex = ""
    if i > 0: # no head in this line
        cSringIndex = cSringIndex.rjust(drawingHeadsize * 4)
    for j in range(len(drawingArray[i])-1):
        cSringIndex += str(drawingArray[i][j]).center(3) + "-"
    cSringIndex += str(drawingArray[i][-1]).center(3)
    print(f'; {cSringIndex}', file=itpFile)
    if i+1 < len(drawingArray): # if not last row add a seperator line
        print(f'; {"|".rjust((drawingHeadsize * 4)+2)}', file=itpFile)
print(';', file=itpFile)

# Add @INSANE input string
cTail = lipidTail.replace('c','C')  # remove generator specific names
cTail = cTail.replace('F','D') # remove generator specific names
current = '@INSANE alhead='+lipidHead+', allink='+lipidLinker+', altail='+cTail+', alname='+lipidName
if lCharge!=None:
    current += ', charge='+str(lCharge)
if lArea!=None and len(lArea) > 0:
    current += ', area='+lArea
print(';' + current, file=itpFile)

# Add @RESNTEST, test if using x3 bead resnames how to fine the last letter (e.g. POP is it POPC, POPE, POPS etc)
resntest = ""
cutoflen = 3
if len(lipidName) < cutoflen: # fix test for short lipid names
    cutoflen = len(lipidName)
if len(headsArray)>0 and headsArray[0]=='PI': # Add PI head\
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[4]==GL1'
elif len(headsArray)>0 and headsArray[0]=='P1': # Add PIP_1 head
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[4]==P1'
elif len(headsArray)>0 and headsArray[0]=='P2': # Add PIP_2 head
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[5]==P2'
elif len(headsArray)>0 and headsArray[0]=='P3': # Add PIP_3 head
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[6]==P3'
elif len(headsArray)>0: # Else build head (works for simple PC, PE, PA, PG, etc heads)
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4]
else: #If -alhead was empty len(headsArray)==0 then no head to add (DAG, CER, TAG etc)
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4]
# Add test if using x3 bead resnames how to fine the last letter (e.g. POP is it POPC, POPE, POPS etc)
print(';@RESNTEST '+resntest, file=itpFile)

# Add all beads
sBeads = ""
beadNameDict = {}
for cBead in beadArray:
    if cBead[0] > 0:
        sBeads += cBead[4].strip()+" "
        beadNameDict[cBead[0]] = cBead[4]
print(';@BEADS '+sBeads, file=itpFile)

# List all bonds
sBonds = ""
for cBond in bondsArray:
    if cBond[0] > 0:
        # Regular bond
        sBonds += beadNameDict[cBond[0]].strip()+"-"+beadNameDict[cBond[1]].strip()+" "
    elif cBond[0] == -3:
        # Named bond then
        sBonds += beadNameDict[cBond[1]].strip()+"-"+beadNameDict[cBond[2]].strip()+" "
print(';@BONDS '+sBonds, file=itpFile)

print(';', file=itpFile)
print('', file=itpFile)
print('[moleculetype]', file=itpFile)
print('; molname      nrexcl', file=itpFile)
print('  ' + lipidName + '          1', file=itpFile)

# Write beads
print('\n[atoms]', file=itpFile)
print('; id 	type 	resnr 	residu 	atom 	cgnr 	charge    (mass)', file=itpFile)
for cBead in beadArray:
    if cBead[0] > 0:
        if len(cBead) > 8: # bead with mass and comment (no option of only mass without comment)
            print('  %2i  %5s  %2i  %4s  %3s  %2i \t%s \t%s \t; %s' % (cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6], cBead[7], cBead[8]), file=itpFile)
        elif len(cBead[7]) > 0: # also print comment
            print('  %2i  %5s  %2i  %4s  %3s  %2i \t%s \t; %s' % (cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6], cBead[7]), file=itpFile)
        else: # no comment no mass
            print('  %2i  %5s  %2i  %4s  %3s  %2i \t%s' % (cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6]), file=itpFile)
    elif cBead[0] == -1:  # Regular comment
        print('; ' + cBead[1], file=itpFile)
    elif cBead[0] == -2:  # gromacs system line
        print(cBead[1], file=itpFile)

# Write lipid bonds
print('\n[bonds]', file=itpFile)
print(';  i  j 	name 	(using named bondtypes from martini_v3.0.0_ffbonded_v2.itp)', file=itpFile)
print(';  i  j 	funct 	force.c.', file=itpFile)
for cBond in bondsArray:
    if cBond[0] > 0: ## Regular bond
        if len(cBond[5]) > 0: # also print comment
            print('  %2i %2i  %2i \t%s \t%s \t; %s' % (cBond[0], cBond[1], cBond[2], cBond[3], cBond[4], cBond[5]), file=itpFile)
        else: # no comment
            print('  %2i %2i  %2i \t%s \t%s' % (cBond[0], cBond[1], cBond[2], cBond[3], cBond[4]), file=itpFile)
    elif cBond[0] == -1:  # Regular comment
        print('; ' + cBond[1], file=itpFile)
    elif cBond[0] == -2:  # gromacs system line
        print(cBond[1], file=itpFile)
    elif cBond[0] == -3:  ## named value in ref file
        if len(cBond[4]) > 0: # also print comment
            print('  %2i %2i \t%s \t; %s' % (cBond[1], cBond[2], cBond[3], cBond[4]), file=itpFile)
        else: # no comment
            print('  %2i %2i \t%s' % (cBond[1], cBond[2], cBond[3]), file=itpFile)

# Write lipid angles
print('\n[angles]', file=itpFile)
print(';  i  j  k 	name 	(using named angletypes from martini_v3.0.0_ffbonded_v2.itp)', file=itpFile)
print(';  i  j  k 	funct 	angle 	force.c.', file=itpFile)
for cAngle in anglesArray:
    if cAngle[0] > 0:  ### Regular angle
        if len(cAngle[6]) > 0: # also print comment
            print('  %2i %2i %2i  %2i \t%s \t%s \t; %s' % (cAngle[0], cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5], cAngle[6]), file=itpFile)
        else: # no comment
            print('  %2i %2i %2i  %2i \t%s \t%s' % (cAngle[0], cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5]), file=itpFile)
    elif cAngle[0] == -1:  # Regular comment
        print('; ' + cAngle[1], file=itpFile)
    elif cAngle[0] == -2:  # gromacs system line
        print(cAngle[1], file=itpFile)
    elif cAngle[0] == -3:  ### named value in ref file.
        if len(cAngle[5]) > 0: # also print comment
            print('  %2i %2i %2i \t%s \t; %s' % (cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5]), file=itpFile)
        else: # no comment
            print('  %2i %2i %2i \t%s' % (cAngle[1], cAngle[2], cAngle[3], cAngle[4]), file=itpFile)

# Write lipid dihedrals
if len(dihedralsArray) > 0:
    print('\n[dihedrals]', file=itpFile)
    print(';  i  j  k  l 	name 	(using named dihedraltypes from martini_v3.0.0_ffbonded_v2.itp)', file=itpFile)
    print(';  i  j  k  l 	funct 	angle 	force.c.', file=itpFile)
    for cDihedral in dihedralsArray:
        if cDihedral[0] > 0:
            if len(cDihedral[8]) > 0: # also print comment
                print('  %2i %2i %2i %2i  %2i \t%i \t%s \t%s \t; %s' % (cDihedral[0], cDihedral[1], cDihedral[2], cDihedral[3], cDihedral[4], cDihedral[5], cDihedral[6], cDihedral[7], cDihedral[8]), file=itpFile)
            else: # no comment
                print('  %2i %2i %2i %2i  %2i \t%i \t%s \t%s' % (cDihedral[0], cDihedral[1], cDihedral[2], cDihedral[3], cDihedral[4], cDihedral[5], cDihedral[6], cDihedral[7]), file=itpFile)
        elif cDihedral[0] == -1:  # Regular comment
            print('; ' + cDihedral[1], file=itpFile)
        elif cDihedral[0] == -2:  # gromacs system line or manual entry
            print(cDihedral[1], file=itpFile)
        elif cDihedral[0] == -3:  # named value in ref file
            if len(cDihedral[6]) > 0: # also print comment
                print('  %2i %2i %2i %2i \t%s \t; %s' % (cDihedral[1], cDihedral[2], cDihedral[3], cDihedral[4],  cDihedral[5],  cDihedral[6]), file=itpFile)
            else: # no comment
                print('  %2i %2i %2i %2i \t%s' % (cDihedral[1], cDihedral[2], cDihedral[3], cDihedral[4],  cDihedral[5]), file=itpFile)

# Write lipid constraints
if len(constraintsArray) > 0:
    print('\n[constraints]', file=itpFile)
    print(';  i  j  k 	funct 	length', file=itpFile)
    for cConstraint in constraintsArray:
        if cConstraint[0] > 0:
            if len(cConstraint[3]) > 0: # also print comment
                print('  %2i %2i \t1 \t%s \t; %s' % (cConstraint[0], cConstraint[1], cConstraint[2], cConstraint[3]), file=itpFile)
            else: # no comment
                print('  %2i %2i \t1 \t%s' % (cConstraint[0], cConstraint[1], cConstraint[2]), file=itpFile)
        elif cConstraint[0] == -1:  # Regular comment
            print('; ' + cConstraint[1], file=itpFile)
        elif cConstraint[0] == -2:  # gromacs system line
            print(cConstraint[1], file=itpFile)

# Write lipid exclusions
if len(exclusionsArray) > 0:
    print('\n[exclusions]', file=itpFile)
    print('; i  j  k  ...', file=itpFile)
    for cExclusion in exclusionsArray:
        if cExclusion[0] > 0:
            print('  %2i %2i \t%s' % (cExclusion[0], cExclusion[1], cExclusion[2]), file=itpFile)
        elif cExclusion[0] == -1:  # Regular comment
            print('; ' + cExclusion[1], file=itpFile)
        elif cExclusion[0] == -2:  # gromacs system line
            print(cExclusion[1], file=itpFile)

print('\n\n', file=itpFile)
itpFile.close()
# End lipid-martini-itp
