import rdkit
from rdkit.Chem.rdMolAlign import GetBestRMS
geom_file = "../butane_NMS/butane_geom.xyz"

f = open(geom_file).readlines()
files = []

for line in f:
    if line.startswith("14"):
        files.append([])
    files[-1].append(line)

s = ""
for line in files[0]:
    s += line
print(s)
