#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/index/detail/rtree/utilities/statistics.hpp>
#include <boost/foreach.hpp>

#include <bits/stdc++.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;
using namespace std::chrono;

#define CAPACITY 100

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
			if (i >= 500000) break;
		}
	}
	file.close();
} 

int main(){
	typedef bg::model::point<float, 2, bg::cs::cartesian> point;
	typedef bg::model::box<point> box;
	typedef pair<box, unsigned> value;

	vector<float> boundary = {180.0, -90.0, 180.0, 90.0};
	
	vector<tuple<float, float, float, float>> queryArray;        
	createQuerySet("aisSample", queryArray);	 

	// index construction
	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	//bgi::rtree<value, bgi::rstar<40,20>> rtree;
	bgi::rtree<value, bgi::rstar<CAPACITY>> rtree;
	for (int q = 0; q < queryArray.size(); q++){
		box b(point(get<1>(queryArray[q]), get<0>(queryArray[q])), point(get<1>(queryArray[q]), get<0>(queryArray[q]))); // create a box
		rtree.insert(make_pair(b, q)); // insert new value
	}
	double time = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
	cout << "Creation time: " << time << endl;

	queryArray.clear();
	createQuerySet("aisSampleQR1000", queryArray);	 

	// range queries
	vector<float> ranges = {0.0025, 0.005, 0.01, 0.02, 0.04, 0.08};
	for (auto rs: ranges){
		time = 0;
		for (auto q: queryArray){
			float pl[2], ph[2];
			pl[0] = get<2>(q) - 0.01;         
			pl[1] = get<1>(q) - 0.01;

  			ph[0] = min(boundary[2], pl[0] + rs*(boundary[2] + abs(boundary[0])));
  			ph[1] = min(boundary[3], pl[1] + rs*(boundary[3] + abs(boundary[1])));

			startTime = high_resolution_clock::now();
			box query_box(point(pl[0], pl[1]), point(ph[0], ph[1]));
			vector<value> result_s;
			rtree.query(bgi::intersects(query_box), back_inserter(result_s));
			time += duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();

		}
		cout << rs << " range query time: " << time << endl;
	}
	
	// kNN queries
	vector<int> kappas = {1, 4, 16, 64};
	for (auto k: kappas){
		startTime = high_resolution_clock::now();
		for (auto q: queryArray){
			vector<value> result_n;
			rtree.query(bgi::nearest(point(get<1>(q), get<0>(q)), k), back_inserter(result_n));
		}
		time = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
		cout << k << " NN query time: " << time << endl;
	}

	// get statistics and size
	tuple<size_t, size_t, size_t, size_t, size_t, size_t> S;
	S = bgi::detail::rtree::utilities::statistics(rtree);

	cout << "Number of levels: " << get<0>(S) << endl;  
	cout << "Number of internal nodes: " << get<1>(S) << endl;  
	cout << "Number of leaf nodes: " << get<2>(S) << endl;  
	cout << "Number of values: " << get<3>(S) << endl;  
	cout << "Min values per node: " << get<4>(S) << endl;  
	cout << "Max values per node: " << get<5>(S) << endl;  
	cout << "RTree size in MB: " << ((get<1>(S) + get<2>(S)) * 4 * sizeof(float) // rectangle size
									+ get<1>(S) * get<4>(S) * 8) / float(1e6)<< endl; // pointer size


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
