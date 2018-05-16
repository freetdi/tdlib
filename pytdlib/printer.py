def dump_treedecomposition_as_dot(T, fname):
    print(str(T.vertices()))
    print(str(T.edges()))

    fout = open(fname, 'w')

    fout.write("graph G {\n")

    for i in range(0, len(T.vertices())):
        fout.write(str(i) + "[label=\"" + str(T.vertices()[i]) + "\"];\n")

    for i in range(0, len(T.edges())-1, 2):
        fout.write(str(T.edges()[i]) + " --" + str(T.edges()[i+1]) + ";\n")

    fout.write("}\n")
    fout.close()

def dump_graph_as_dot(G, fname):

    print(str(G.vertices()))
    print(str(G.edges()))

    fout = open(fname, 'w')

    fout.write("graph G {\n")

    for i in range(0, len(G.vertices())):
        fout.write(str(i) + ";\n")

    for i in range(0, len(G.edges())-1, 2):
        fout.write(str(G.edges()[i]) + " --" + str(G.edges()[i+1]) + ";\n")

    fout.write("}\n")
    fout.close()












