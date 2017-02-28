import os

def convert(file1, file2, c):
    fin = open(file1, 'r')
    fout = open(file2, 'a')
    
    lines = list()
    
    for line in fin:
         lines.append(line.replace("\n", ""))

    fout.write("V_"+ str(c) + " = [");

    first = True
    
    for i in range(1, len(lines)-1):
        if "--" not in lines[i]:
            fout.write(lines[i].replace(";", "").replace("\n", ""))
            if "--" not in lines[i+1]:
                fout.write(",")
            else:
                fout.write("]\n")
        else:
            if first:
                fout.write("E_"+ str(c) + " = [");
                first = False
            edge = lines[i].replace(";", "").split(" -- ")
            fout.write("[" + str(edge[0]) + "," + str(edge[1]) + "]")
            if "}" not in lines[i+1]:
                fout.write(",")
            else:
                fout.write("]\n")
                
    fout.write("G_name = \"" + file1.replace(".dot", "") + "\"\n")
    fout.write("\n")


    fin.close()
    fout.close()

c = 0

for file in os.listdir("."):
    if ".dot" in file:
        convert(file, "pytdlib_zoo.py", c)
        c += 1

