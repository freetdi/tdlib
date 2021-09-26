#include <treedec/bucket_sorter.hpp>
#include <boost/property_map/property_map.hpp>


int main(){

    typedef boost::iterator_property_map<unsigned*,
        boost::identity_property_map, unsigned, unsigned&> map_type;
    typedef boost::bucket_sorter<unsigned, unsigned,
            map_type, boost::identity_property_map > container_type;

	 std::vector<unsigned> V(10);
	 map_type M(&V[0], boost::identity_property_map());


	 container_type B(10, 10, M);

	 for(auto i : B[0] ){ untested();
		 (void)i;
		 assert(false);
	 }

	 V[0]=1;
	 V[1]=1;
	 V[2]=1;
	 V[3]=1;
	 V[4]=1;

	 B.push(0);
	 B.push(1);

	 assert(B[0].begin() == B[0].end());
	 for(auto i : B[0] ){ untested();
		 (void)i;
		 assert(false);
	 }

	 assert(B[1].begin() != B[1].end());
	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }

	 B.remove(0);
	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }
	 assert(B[1].begin() != B[1].end());
	 B.remove(1);
	 assert(B[1].begin() == B[1].end());

	 B.push(1);
	 B.push(4);
	 B.remove(4);
	 B.push(2);
	 B.push(3);
	 V[2]=3;
	 B.update(2);

	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }

	 B.remove(0);
	 B.update(4);
	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }

	 std::cout << "up0\n";
	 B.update(0);
	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }

#if 0
	 B.push(0); // not allowed. crash
	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }
#endif
	 B.remove(0);
	 B.push(0);
	 for(auto i : B[1] ){
		 std::cout<< i<< "\n";
	 }

	 // try to bypass map update
	 B.remove(4);
	 B.remove(4);
	 B.remove(4);
	 B.remove(4);
	 B.update(4);
	 B.remove(1);
	 B.update(1);
	 for(auto i : B[1] ){
		 std::cout << ". " << i << "\n";
	 }

	 B.remove(4);
	 B.remove(1);
	 B.remove(0);
	 for(auto i : B[1] ){ untested();
		 std::cout << ".. " << i << "\n";
	 }

	 std::cout<< "done\n";

}
