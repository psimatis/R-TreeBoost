#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/index/detail/rtree/utilities/statistics.hpp>
/*#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
*/
#include <boost/foreach.hpp>

#include <bits/stdc++.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;
using namespace std::chrono;

#define CAPACITY 10

void createQuerySet(string fileName, vector<tuple<float, float, float, float>> &queryArray) {
	string line;
	int i = 0;
	ifstream file(fileName);
	if (file.is_open()) {
		getline(file, line);
		while (getline(file, line)){
			float xl, xh, yl, yh;
			istringstream buf(line);
			buf >> xl >> xh >> yl >> yh;
			queryArray.emplace_back(make_tuple(xl, xh, yl, yh));
			i++;
			if (i >= 100000) break;
		}
	}
	file.close();
} 

/*template <typename T, size_t I = 0, size_t S = tuple_size<T>::value>
struct print_tuple{
	template <typename Os>
	static inline Os & apply(Os & os, T const& t){
		os << get<I>(t) << ", ";
	    return print_tuple<T, I+1>::apply(os, t);
	}
};

template <typename T, size_t S>
struct print_tuple<T, S, S>{
	template <typename Os>
	static inline Os & apply(Os & os, T const&){
		return os;
	}
};*/

int main(){
	typedef bg::model::point<float, 2, bg::cs::cartesian> point;
	typedef bg::model::box<point> box;
	typedef pair<box, unsigned> value;

	// create the rtree using default constructor
	vector<tuple<float, float, float, float>> queryArray;        
	createQuerySet("aisSample", queryArray);	 

	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	//bgi::rtree<value, bgi::rstar<40,20>> rtree;
	bgi::rtree<value, bgi::rstar<CAPACITY>> rtree;
	for (int q = 0; q < queryArray.size(); q++){
		box b(point(get<0>(queryArray[q]), get<1>(queryArray[q])), point(get<0>(queryArray[q]), get<1>(queryArray[q]))); // create a box
		rtree.insert(make_pair(b, q)); // insert new value
	}
	double time = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
	cout << "Creation time: " << time << endl;

	queryArray.clear();
	createQuerySet("aisSampleQR1000", queryArray);	 

	// find values intersecting some area defined by a box
	startTime = high_resolution_clock::now();
	for (auto q: queryArray){
		box query_box(point(get<0>(q), get<1>(q)), point(get<2>(q), get<3>(q)));
		//box query_box(point(0,0), point(1,1));			
		vector<value> result_s;
		rtree.query(bgi::intersects(query_box), back_inserter(result_s));
	}
	time = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
	cout << "Range query time: " << time << endl;

	// find 5 nearest values to a point
	startTime = high_resolution_clock::now();
	for (auto q: queryArray){
		vector<value> result_n;
		rtree.query(bgi::nearest(point(get<0>(q), get<1>(q)), 5), back_inserter(result_n));
	}
	time = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
	cout << "KNN query time: " << time << endl;

	//typedef tuple<size_t, size_t, size_t, size_t, size_t, size_t> S;
	tuple<size_t, size_t, size_t, size_t, size_t, size_t> S;
	//S s;
	//print_tuple<S>::apply(cout, bgi::detail::rtree::utilities::statistics(rtree)) << endl;
	S = bgi::detail::rtree::utilities::statistics(rtree);

	cout << get<0>(S) << endl;



	// note: in Boost.Geometry WKT representation of a box is polygon
	// display results
	/*
	cout << "spatial query box:" << endl;
	cout << bg::wkt<box>(query_box) << endl;
	cout << "spatial query result:" << endl;
	BOOST_FOREACH(value const& v, result_s)
		cout << bg::wkt<box>(v.first) << " - " << v.second << endl;

	cout << "knn query point:" << endl;
	cout << bg::wkt<point>(point(0, 0)) << endl;
	cout << "knn query result:" << endl;
	BOOST_FOREACH(value const& v, result_n)
		cout << bg::wkt<box>(v.first) << " - " << v.second << endl;
	return 0;*/
}
