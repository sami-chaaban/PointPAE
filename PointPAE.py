
#Dependencies

from pathlib import Path
import json
import os

####################################################################
### INPUT ########################


pdbfile = "mycomplex_rank_1_model_1.pdb"

paejson = "mycomplex_rank_1_model_1_scores.json"

chainletter = "D" #Chain of interest. Open the PDB in chimera to verify the letter.

residuenum = 29 #Residue of interest on that chain

sample = "y" #Choose "x" or "y" to sample along that axis. Choose "average" to average both axes.

chimerax = True #True for a ChimeraX-compatible attribute file, False for Chimera

####################################################################
##################################

directory=os.path.dirname(pdbfile)
scores = json.loads(Path(paejson).read_text())
pae = scores["pae"]

print("\n>There are " + str(len(pae)) + " residues in the stitched sequence")

##################################

#Identify the chain letters and their sequence length

lengths=[]
chain_names = []
previouschain = "*"
i = 0

for line in open(pdbfile,"r"):
    if line[:4] == "ATOM":
        chain = line[21:22]
        resnum = int(line[22:22+5])
        if resnum == 1 and chain != previouschain:
            if i != 0:
                lengths.append(previousresnum)
            chain_names.append(chain)
            previouschain = chain
        i +=1
        previousresnum=resnum
lengths.append(previousresnum)

print(">The identified chains in the pdb are " + str(chain_names))
print(">The lengths of those chains are " + str(lengths))

##################################

#Extract the PAE values for a given residue

residuenum_stitch = residuenum -1
chainnum = chain_names.index(chainletter)

for c in range(chainnum):
    residuenum_stitch = residuenum_stitch + lengths[c]

#print(pae)
#print(len(pae))



if sample == "y":
    pae_atresidue=[]
    for p in pae:
        pae_atresidue.append(p[residuenum_stitch])
elif sample == "x":
    pae_atresidue = pae[residuenum_stitch]
elif sample == "average":
    pae_atresidue_x=[]
    for p in pae:
        pae_atresidue_x.append(p[residuenum_stitch])
    pae_atresidue_y = pae[residuenum_stitch]
    pae_atresidue = [(x + y) / 2 for x, y in zip(pae_atresidue_x, pae_atresidue_y)]

##################################

#Generate individual attribute files for each chain

outputnames = []
attributenames = []

if chimerax:
    extension="_cx.defattr"
else:
    extension=".txt"

for chain,length in enumerate(lengths):
    
    outputname = os.path.join(directory, "pae_"+chain_names[chain]+"_relative_"+ chainletter + "_" +str(residuenum)+extension)
    outputnames.append(outputname)
    
    output = open(outputname, 'w')
    
    attributename = "pae_chain" + str(chain+1) + "_relative_" + chainletter + "_" + str(residuenum)
    attributenames.append(attributename)
    
    output.write("attribute: " + attributename + "\n")
    output.write("match mode: 1-to-1\n")
    output.write("recipient: residues\n")
    
    startingposition = 0
    for c in range(chain):
        startingposition = startingposition + lengths[c]
    
    for s in range(startingposition, startingposition+length):
        p = pae_atresidue[s]
        if chimerax:
            output.write("\t/" + chain_names[chain] + ":" + str(s-startingposition+1) + "\t" + str(p) + "\n")
        else:
            output.write("\t:" + str(s-startingposition+1) + "."+chain_names[chain]+"\t" + str(p) + "\n")
    
    output.close()
    
##################################

#Generate a script to open the attribute files and color them properly

if chimerax:
    extension=".cxc"
else:
    extension=".com"
    
scriptname = os.path.join(directory, "script"+"_relative_" + chainletter + "_" + str(residuenum) + extension)
    
script = open(scriptname, 'w')

if chimerax:
    script.write("open \"" + pdbfile + "\"\n")
else:
    script.write("open " + pdbfile + "\n") #Adding quotations to filename didn't work in Chimera for some reason

for outputname,attributename in zip(outputnames,attributenames):
    
    if chimerax:
        script.write("open " + outputname + "\n")
        script.write("color byattribute " + attributename + " range 0.75,35 palette cyanmaroon\n")
    else:
        script.write("defattr " + outputname + "\n")
        script.write("rangecol " + attributename + " 0.75 blue 17.875 white 35 red\n")
    
script.close()

print(">Wrote: " + scriptname)

if chimerax:
    print("-----------------------------------------------")
    print(">Open " + scriptname + " in ChimeraX.")
    print("-----------------------------------------------\n")
else:
    print("-----------------------------------------------")
    print(">Open " + scriptname + " in Chimera.")
    print("-----------------------------------------------\n")