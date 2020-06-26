#! /usr/local/bin/python3
import numpy as np
import os
from sys import argv


data           = np.genfromtxt("{}/output.txt".format(argv[1]),delimiter=", ")

conditions = {}
with open("{}/conditions.txt".format(argv[1])) as f:
    for line in f:
       (key, val) = line.split()
       conditions[key] = float(val)

for i in range(int(conditions["tmax"]/conditions["outputInterval"])):
    os.system("clear")
    print("Rendering {}/{}".format(i,int(conditions["tmax"]/conditions["outputInterval"])))
    outfile = open("{}/povrayTmp{:03d}.pov".format(argv[1],i),"w")
    outfile.write("#include \"colors.inc\"\n")
    outfile.write("camera {\n" )
    outfile.write("  sky <0,0,1>           \n")
    outfile.write("  direction <-1,0,0>      \n")
    outfile.write("  right <-4/3,0,0>      \n")
    outfile.write("  location <0,0,20> \n" )
    outfile.write("  look_at <0,0,0>     \n" )
    outfile.write("  angle 15      \n")
    outfile.write("}\n")
    outfile.write("global_settings { ambient_light White }\n")
    outfile.write("light_source {\n" )
    outfile.write("  <10,-10,20>   \n")
    outfile.write("  color White*2 \n")
    outfile.write("}\n")
    outfile.write("background { color White }\n" )
    subdata = data[i*int(conditions["Ncells"]):(i+1)*int(conditions["Ncells"]),:]
    for j in range(int(conditions["Ncells"])):
        outfile.write("sphere{{<{},{},0>,{} texture{{pigment{{color Green}}}}}}\n".format(subdata[j,0],subdata[j,1],conditions["σ"]/2.0))
    outfile.close()
    os.system("povray {}/povrayTmp{:03d}.pov > /dev/null 2>&1".format(argv[1],i))
    os.system("rm {}/povrayTmp{:03d}.pov".format(argv[1],i))

os.system("convert -delay 10 -loop 0 {}/povrayTmp*.png {}/animated.gif".format(argv[1],argv[1]))
os.system("rm {}/*.png".format(argv[1]))
