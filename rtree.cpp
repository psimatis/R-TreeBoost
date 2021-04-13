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

void parseDataFile(string fileName, vector<tuple<int, float, float>> &dataArray) {
	string line;
	ifstream file(fileName);
	if (file.is_open()) {
		getline(file, line);
		while (getline(file, line)){
			int id;
			float x, y;
			istringstream buf(line);
			buf >> id >> x >> y;
			dataArray.emplace_back(make_tuple(id, x, y));
		}
		file.close();
	} else cout << "Cant open file: " << fileName << endl;
} 

void parseQueryFile(string fileName, vector<tuple<char, float, float, float>> &queryArray) {
	string line;
	int i = 0;
	ifstream file(fileName);
	if (file.is_open()) {
		getline(file, line);
		while (getline(file, line)){
			char type;
			float x, y, size;
			istringstream buf(line);
			buf >> type >> x >> y >> size;
			queryArray.emplace_back(make_tuple(type, x, y, size));
			i++;
			if (i >= 500000) break;
		}
		file.close();
	} else cout << "Cant open file: " << fileName << endl;
} 

int main(){
	typedef bg::model::point<float, 2, bg::cs::cartesian> point;
	typedef bg::model::box<point> box;
	typedef pair<box, unsigned> value;

	vector<float> boundary = {180.0, -90.0, 180.0, 90.0};
	
	vector<tuple<int, float, float>> dataArray;        
	parseDataFile("data/aisCleanSample1e6.txt", dataArray);	 

	cout << "#Data: " << dataArray.size() << endl;

	vector<tuple<char, float, float, float>> queryArray;        
	parseQueryFile("test.txt", queryArray);	 

	// index construction
	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	//bgi::rtree<value, bgi::rstar<40,20>> rtree;
	bgi::rtree<value, bgi::rstar<CAPACITY>> rtree;
	for (auto q: dataArray){
		box b(point(get<2>(q), get<1>(q)), point(get<2>(q), get<1>(q))); // create a box
		rtree.insert(make_pair(b, get<0>(q))); // insert new value
	}
	double time = duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
	cout << "Creation time: " << time << endl;

	map<string, double> rangeLog, knnLog;
	
	for (auto q: queryArray){
		if (get<0>(q) == 'r'){
			float pl[2], ph[2], rs = get<3>(q);
			pl[0] = get<2>(q) - 0.01;         
			pl[1] = get<1>(q) - 0.01;

  			ph[0] = min(boundary[2], pl[0] + rs*(boundary[2] + abs(boundary[0])));
  			ph[1] = min(boundary[3], pl[1] + rs*(boundary[3] + abs(boundary[1])));

			startTime = high_resolution_clock::now();
			box query_box(point(pl[0], pl[1]), point(ph[0], ph[1]));
			vector<value> result_s;
			rtree.query(bgi::intersects(query_box), back_inserter(result_s));
			rangeLog["time " + to_string(rs)] += duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
			rangeLog["count " + to_string(rs)]++;
			cout << "spatial query box:" << endl;
			cout << bg::wkt<box>(query_box) << endl;
			//cout << "spatial query result:" << endl;
			int rCount = 0;
			BOOST_FOREACH(value const& v, result_s)
				rCount++;
				//cout << bg::wkt<box>(v.first) << " - " << v.second << endl;
			cout << rCount << endl;
		}
		else if (get<0>(q) == 'k'){
			startTime = high_resolution_clock::now();
			vector<value> result_n;
			rtree.query(bgi::nearest(point(get<2>(q), get<1>(q)), get<3>(q)), back_inserter(result_n));
			knnLog["time " + to_string( get<3>(q))] += duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
			knnLog["count " + to_string(get<3>(q))]++;
			
			//cout << "knn query point:" << endl;
			//cout << bg::wkt<point>(point(0, 0)) << endl;
			//cout << "knn query result:" << endl;
			//BOOST_FOREACH(value const& v, result_n)
			//	cout << bg::wkt<box>(v.first) << " - " << v.second << endl;
		}
	}

	for (auto it = rangeLog.cbegin(); it != rangeLog.cend(); ++it)  
		cout<< it->first << ": " << it->second << endl;   

	for (auto it = knnLog.cbegin(); it != knnLog.cend(); ++it)  
			cout<< it->first << ": " << it->second << endl;   
	

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
									
	return 0;
}
