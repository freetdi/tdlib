import sys

sys.path.append('./..')

LONG_VERBOSE = True

def print_graph_name(PREFIX, c):
    if not LONG_VERBOSE:
        return
    import Zoo
    print(str(eval(PREFIX+".G_"+str(c)+"_name")) + "[" + str(len(eval(PREFIX+".V_"+str(c)))) + "," + str(len(eval(PREFIX+".E_"+str(c))))+"]")

def skip(PREFIX, c, f):
    import Zoo
    if f(len(eval(PREFIX+".V_"+str(c))), len(eval(PREFIX+".E_"+str(c)))):
        return True
    return False
