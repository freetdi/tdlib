#include <tdlib/bucket_sorter.hpp>
#include <boost/property_map/property_map.hpp>


int main(){

	typedef boost::iterator_property_map<unsigned*,
			  boost::identity_property_map, unsigned, unsigned&> map_type;
	typedef boost::bucket_sorter<unsigned, unsigned,
			  map_type, boost::identity_property_map > container_type;

	std::vector<unsigned> V(10);
	map_type M(&V[0], boost::identity_property_map());

	container_type B(10, 10, M);

	V[0]=0;
	B.push(0);

	V[2]=1;
	B.push(2);

	V[0]=1;
	B.update(0); // kills 2

	unsigned c=0;
	for(auto i : B[1] ){
		++c;
		std::cout<< i<< "\n";
	}
	assert(c==2);

}
