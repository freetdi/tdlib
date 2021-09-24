#include <treedec/bucket_sorter.hpp>
#include <boost/property_map/property_map.hpp>


int main(){

	typedef boost::iterator_property_map<unsigned*,
			  boost::identity_property_map, unsigned, unsigned&> map_type;
	typedef boost::bucket_sorter<unsigned, unsigned,
			  map_type, boost::identity_property_map> container_type;

	std::vector<unsigned> V(10, 1);
	map_type M(&V[0], boost::identity_property_map());

	container_type B(10, 10, M);
	container_type Bb(10, 10, M);

 	B.push_back(4);
 	B.push_back(5);

	unsigned c = 0;
	for(auto i : B[1] ){
		++c;
		std::cout<< i<< "\n";
	}
	assert(c==2);

	B.push_front(0);
	B.push_front(1);
	B.push_front(2);
 	B.push_back(3);
//  	B.push_back(4);

	c = 0;
	for(auto i : B[1] ){
		++c;
		std::cout<< i<< "\n";
	}
	assert(c==6);

	Bb.push_back(3);
	Bb.push_back(4);
	Bb.push_front(0);
	Bb.push_front(1);
	V[1]=2;
	Bb.update(1);
	Bb.push_front(2);

	c = 0;
	for(auto i : Bb[1] ){
		++c;
		std::cout<< i<< "\n";
	}
	assert(c==4);

}
