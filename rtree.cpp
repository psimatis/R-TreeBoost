// Counters for node visits
// Warning! 
// 1. They need to be declared before the includes
// 2. Boost files spatial_query.hpp and distance_query.hpp need to be modified
int rangeInternalCount = 0;
int rangeLeafCount = 0;
int knnInternalCount = 0;
int knnLeafCount = 0;

#include <boost/geometry.hpp>
#include <boost/geometry/index/detail/rtree/utilities/statistics.hpp>
#include <boost/foreach.hpp>
#include <bits/stdc++.h>

#define CAPACITY 1000

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;
using namespace std::chrono;
namespace bgid = bgi::detail;

// Visitor that counts the number of pointers in the RTree
// it assists with size computation
int pointerCount = 0;
template <typename Value, typename Options, typename Translator, typename Box, typename Allocators>
class pointer_visitor
    : public bgid::rtree::visitor<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag, true>::type
{
    typedef typename bgid::rtree::internal_node<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type internal_node;
    typedef typename bgid::rtree::leaf<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type leaf;

public:
    void operator()(internal_node const& n){
        typedef typename bgid::rtree::elements_type<internal_node>::type elements_type;
        elements_type const& elements = bgid::rtree::elements(n);

		pointerCount += elements.size();

        for ( typename elements_type::const_iterator it = elements.begin(); it != elements.end() ; ++it){
            bgid::rtree::apply_visitor(*this, *(it->second));
        }
    }
};

template <typename Rtree>
inline void pointerVisitor(Rtree const& tree){
    typedef bgid::rtree::utilities::view<Rtree> RTV;
    RTV rtv(tree);

    pointer_visitor<
        typename RTV::value_type,
        typename RTV::options_type,
        typename RTV::translator_type,
        typename RTV::box_type,
        typename RTV::allocators_type
    > v;   
    rtv.apply_visitor(v);
}

void parseDataFile(string fileName, vector<tuple<int, float, float>> &dataArray, int limit) {
	string line;
	int i = 0;
	ifstream file(fileName);
	if (file.is_open()) {
		getline(file, line);
		while (getline(file, line)){
			int id;
			float x, y;
			istringstream buf(line);
			buf >> id >> x >> y;
			dataArray.emplace_back(make_tuple(id, x, y));
			i++;
			if (i >= limit) break;
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
			float x, y, info;
			istringstream buf(line);
			buf >> type >> x >> y >> info;
			queryArray.emplace_back(make_tuple(type, x, y, info));
		}
		file.close();
	} else cout << "Cant open file: " << fileName << endl;
} 

int main(int argc, char** argv){

	if (argc != 4){
		cout << "Usage: ./rtree dataFile limit queryFile" << endl;
		exit(1);
	}

	typedef bg::model::point<float, 2, bg::cs::cartesian> point;
	typedef bg::model::box<point> box;
	typedef pair<point, unsigned> value;
	
	vector<float> boundary = {180.0, -90.0, 180.0, 90.0};
	
	vector<tuple<int, float, float>> dataArray;        
	parseDataFile(argv[1], dataArray, atoi(argv[2]));	 

	vector<tuple<char, float, float, float>> queryArray;        	 
	parseQueryFile(argv[3], queryArray);	 

	vector<point> contourCenters;
	vector<value> cloud;
	for (auto q: dataArray){
		point p(get<2>(q), get<1>(q));
		contourCenters.push_back(p);
	}	
	size_t id_gen = 0;
	transform(contourCenters.begin(), contourCenters.end(), back_inserter(cloud), 
	          [&](point const& p) { return make_pair(p, id_gen++); });
	cout << "Done parsing input data" << endl;

	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	//bgi::rtree<value, bgi::rstar<CAPACITY>> rtree(cloud.begin(), cloud.end());
	//bgi::rtree<value, bgi::linear<CAPACITY>> rtree(cloud.begin(), cloud.end());
	bgi::rtree<value, bgi::quadratic<CAPACITY>> rtree(cloud.begin(), cloud.end());
	double time = duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
	cout << "Creation time: " << time << endl;

	map<string, double> rangeLog, knnLog, inLog;
	for (auto q: queryArray){
		if (get<0>(q) == 'r'){
			float pl[2], ph[2], rs = get<3>(q);
			pl[0] = get<2>(q) - 0.01;         
			pl[1] = get<1>(q) - 0.01;

  			ph[0] = min(boundary[2], pl[0] + rs*(boundary[2] + abs(boundary[0])));
  			ph[1] = min(boundary[3], pl[1] + rs*(boundary[3] + abs(boundary[1])));

			vector<value> result_s;
			box queryBox(point(pl[0], pl[1]), point(ph[0], ph[1]));
			startTime = high_resolution_clock::now();
			rtree.query(bgi::intersects(queryBox), back_inserter(result_s));
			rangeLog["time " + to_string(rs)] += duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
			rangeLog["count " + to_string(rs)]++;
			rangeLog["pages " + to_string(rs)] += rangeInternalCount + rangeLeafCount;
			rangeLog["leaf " + to_string(rs)] += rangeLeafCount;
			rangeLog["internal " + to_string(rs)] += rangeInternalCount;
			rangeLeafCount = 0;
			rangeInternalCount = 0;
			//cout << "Range query:" << bg::wkt<box>(queryBox) << endl;
			//BOOST_FOREACH(value const& v, result_s)
				//cout << bg::wkt<point>(v.first) << " - " << v.second << endl;
		}
		else if (get<0>(q) == 'k'){
			vector<value> result_n;
			point queryPoint(point(get<2>(q), get<1>(q)));
			startTime = high_resolution_clock::now();
			rtree.query(bgi::nearest(queryPoint, get<3>(q)), back_inserter(result_n));
			knnLog["time " + to_string(get<3>(q))] += duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
			knnLog["count " + to_string(get<3>(q))]++;
			knnLog["pages " + to_string(get<3>(q))] += knnInternalCount + knnLeafCount; 
			knnLog["leaves " + to_string(get<3>(q))] += knnLeafCount; 
			knnLog["internal " + to_string(get<3>(q))] += knnInternalCount; 
			knnLeafCount = 0;
			knnInternalCount = 0;
			//cout << "kNN query: " << bg::wkt<point>(point(queryPoint)) << endl;
			//BOOST_FOREACH(value const& v, result_n)
				//cout << bg::wkt<point>(v.first) << " - " << v.second << endl;
		}
		else if (get<0>(q) == 'i'){
			point p(point(get<2>(q), get<1>(q)));
			//auto pair = make_pair(p, get<3>(q));
			startTime = high_resolution_clock::now();
			rtree.insert(make_pair(p, get<3>(q)));
			inLog["time"] += duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
			inLog["count"]++;			
		}
	}

	cout << "---Insertions---" << endl;
	for (auto it = inLog.cbegin(); it != inLog.cend(); ++it)
		cout << it->first << ": " << it->second << endl;

	cout << "---Range---" << endl;
	for (auto it = rangeLog.cbegin(); it != rangeLog.cend(); ++it)  
		cout<< it->first << ": " << it->second << endl;   

	cout << "---KNN---" << endl;
	for (auto it = knnLog.cbegin(); it != knnLog.cend(); ++it)  
			cout<< it->first << ": " << it->second << endl; 

	cout << "---RTree Statistics---" << endl;
	tuple<size_t, size_t, size_t, size_t, size_t, size_t> S;
	S = bgi::detail::rtree::utilities::statistics(rtree);

	cout << "Number of levels: " << get<0>(S) << endl;  
	cout << "Number of internal nodes: " << get<1>(S) << endl;  
	cout << "Number of leaf nodes: " << get<2>(S) << endl;  
	cout << "Number of values: " << get<3>(S) << endl;  
	cout << "Min values per node: " << get<4>(S) << endl;  
	cout << "Max values per node: " << get<5>(S) << endl;  
	cout << "Capacity: " << CAPACITY << endl;
	pointerVisitor(rtree);
	cout << "Internal pointer count: " << pointerCount << endl;
	cout << "RTree size in MB (correct): " << ((get<1>(S) + get<2>(S)) * 4 * sizeof(float) // rectangle size
												+ ((get<1>(S) + get<2>(S)) * sizeof(int)) // int for element count
	 											+ pointerCount * 8) / float(1e6) << endl; // internal pointer size
}
