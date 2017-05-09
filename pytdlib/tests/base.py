import sys

sys.path.append('./..')

LONG_VERBOSE = True

def print_graph_name(PREFIX, c):
    if not LONG_VERBOSE:
        return
    import Zoo
    import Dimacs
    import Networks
    print(str(eval(PREFIX+".name_"+str(c))) + "[" + str(len(eval(PREFIX+".V_"+str(c)))) + "," + str(len(eval(PREFIX+".E_"+str(c))))+"]")

def skip(PREFIX, c, f):
    import Zoo
    import Dimacs
    import Networks
    if f(len(eval(PREFIX+".V_"+str(c))), len(eval(PREFIX+".E_"+str(c)))):
        return True
    return False
