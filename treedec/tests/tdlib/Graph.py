# compatiblility hack, nothing to see here.

from treedec._graph import _balsvu, _balvvu

class Graph:
	def vertices(self):
		return self._v

	def __init__(self, V, E, type=1):
		self._v = V
		self._e = E

		m = {}
		map_e = []
		i=0
		for v in V:
			m[v] = i
			i+=1

		for e in E:
			map_e.append((m[e[0]], m[e[1]]))
		
		#  from pytdlib...
		#  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t; //type 0
		#  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> TD_graph_vec_t; //type 1
		#  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS> TD_graph_directed_t; //type 2
		#  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> TD_graph_directed_vec_t; //type 3

		if(type==0):
			self._graph = _balsvu(map_e, len(m))
		elif(type==1):
			self._graph = _balvvu(map_e, len(m))
		elif(type==4):
			from treedec._graph import _gsgvvu16
			self._graph = _gsgvvu16(map_e, len(m))
		elif(type==5):
			from treedec._graph import _gsgvvu32
			self._graph = _gsgvvu32(map_e, len(m))
		elif(type==6):
			from treedec._graph import _gsgvvu64
			self._graph = _gsgvvu64(map_e, len(m))
		else:
			raise ValueError("wrong type")

