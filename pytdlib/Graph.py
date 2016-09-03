
class Graph:
    def __init__(self, V, E, backend="boost_graph_undirected"):
        self._V = V
        self._E = E

        if backend not in ["boost_graph_undirected", "boost_graph_directed"]:
            raise Exception

        self._backend = backend

    def vertices(self):
        return self._V

    def edges(self):
        return self._E

    def backend(self):
        return self._backend
