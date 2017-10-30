import os

def convert(dir, file1, file2, c):
    fin = open(dir+file1, 'r')
    fout = open(file2, 'a')

    lines = list()

    for line in fin:
         lines.append(line.replace("\n", ""))

    first = True

    num_v = 0

    name = dir.replace(".", "").replace("/", "")

    E = list()

    for i in range(1, len(lines)-1):
        if ";" not in lines[i]:
            continue

        if "label" in lines[i]:
            num_v +=1
        else:
            if first:
                fout.write("V_"+ str(c) + " = range(" + str(num_v) + ")\n")
                fout.write("E_"+ str(c) + " = [");
                first = False
            edge = lines[i].replace(";", "").split("->")

            st = "[" + str(edge[0]) + "," + str(edge[1]) + "]"
            if st not in E:
                E.append(st)

    if first:
        fout.write("V_"+ str(c) + " = range(" + str(num_v) + ")\n")
        fout.write("E_"+ str(c) + " = []\n")

    for i in range(0, len(E)-1):
        fout.write(E[i] + ",")

    if len(E) > 0:
        fout.write(E[len(E)-1] + "]\n")

    fout.write("name_"+str(c) + " = \"" + name+"_"+file1.replace(".dot", "") + "\"\n")
    fout.write("\n")


    fin.close()
    fout.close()

def to_it(dir, cnt):
    for file in os.listdir(dir):
        if ".dot" in file:
            convert(dir, file, "NewCFGs.py", cnt)
            cnt += 1

    return cnt

cnt = 0

cnt = to_it("./atomthreads/", cnt)
cnt = to_it("./contiki/", cnt)
cnt = to_it("./coremark/", cnt)
cnt = to_it("./dhrystone/", cnt)
cnt = to_it("./fuzix/", cnt)
cnt = to_it("./nuttx/", cnt)
cnt = to_it("./stdlib/", cnt)
cnt = to_it("./whetstone/", cnt)


