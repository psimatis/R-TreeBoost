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
#include <boost/geometry/index/detail/rtree/utilities/print.hpp>
#include <boost/foreach.hpp>
#include <bits/stdc++.h>

#define CAPACITY 512

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

void parseQueryFile(string fileName, vector<tuple<char, vector<float>, float>> &queryArray) {
    cout << "Begin query parsing for R-Tree" << endl;
    string line;

    ifstream file(fileName);
    if (file.is_open()) {
        // getline(file, line); // Skip the header line
        while (getline(file, line)) {
            char type = line[line.find_first_not_of(" ")];
            vector<float> q;
            line = line.substr(line.find_first_of(type) + 1);
            const char *cs = line.c_str();
            char *end;
            int params = (type == 'r') ? 4 : 2;
            for (int d = 0; d < params; d++) {
                q.emplace_back(strtof(cs, &end));
                cs = end;
            }
            float info = strtof(cs, &end);
            queryArray.emplace_back(make_tuple(type, q, info));
        }
        file.close();
    }
    cout << "Finish query parsing for R-Tree" << endl;
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

	vector<tuple<char, vector<float>, float>> queryArray;        	 
	parseQueryFile(argv[3], queryArray);	 

	vector<point> contourCenters;
	vector<value> cloud;
	for (auto q: dataArray){
		point p(get<1>(q), get<2>(q));
		contourCenters.push_back(p);
	}	
	size_t id_gen = 0;
	transform(contourCenters.begin(), contourCenters.end(), back_inserter(cloud), 
	          [&](point const& p) { return make_pair(p, id_gen++); });
	cout << "Done parsing input data" << endl;

	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	bgi::rtree<value, bgi::rstar<CAPACITY>> rtree(cloud.begin(), cloud.end());
	//bgi::rtree<value, bgi::linear<CAPACITY>> rtree(cloud.begin(), cloud.end());
	//bgi::rtree<value, bgi::quadratic<CAPACITY>> rtree(cloud.begin(), cloud.end());
	double time = duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
	cout << "Creation time: " << time << endl;

	map<string, double> rangeLog, knnLog, inLog;
	for (auto q: queryArray){
		if (get<0>(q) == 'r'){
			array<float, 4> query;
    		for (int i = 0; i < query.size(); i++)
        		query[i] = get<1>(q)[i];
			float rs = get<2>(q);
			
			vector<value> result_s;
			box queryBox(point(query[0], query[1]), point(query[2], query[3]));
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
			//cout << "results size: " << result_s.size() << endl;
		}
		else if (get<0>(q) == 'k'){
			array<float, 2> p;
    		for (uint i = 0; i < p.size(); i++)
        		p[i] = get<1>(q)[i];
    		int k = get<2>(q);
			
			vector<value> result_n;
			point queryPoint(point(p[0], p[1]));
			startTime = high_resolution_clock::now();
			rtree.query(bgi::nearest(queryPoint, k), back_inserter(result_n));
			knnLog["time " + to_string(k)] += duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
			knnLog["count " + to_string(k)]++;
			knnLog["pages " + to_string(k)] += knnInternalCount + knnLeafCount; 
			knnLog["leaves " + to_string(k)] += knnLeafCount; 
			knnLog["internal " + to_string(k)] += knnInternalCount; 
			knnLeafCount = 0;
			knnInternalCount = 0;
			//cout << "kNN query: " << bg::wkt<point>(point(queryPoint)) << endl;
			//BOOST_FOREACH(value const& v, result_n)
				//cout << bg::wkt<point>(v.first) << " - " << v.second << endl;			
		}
		else if (get<0>(q) == 'i'){
			array<float, 2> p;
    		for (uint i = 0; i < p.size(); i++)
        		p[i] = get<1>(q)[i];
    		int id = get<2>(q);
			
			point pNew(point(p[0], p[1]));
			auto pair = make_pair(pNew, id);
			startTime = high_resolution_clock::now();
			rtree.insert(pair);
			inLog["time"] += duration_cast<microseconds>(high_resolution_clock::now() - startTime).count();
			inLog["count"]++;
			//if (long(inLog["count"]) % long(1e4) == 0)
				//cerr << "Count: " << inLog["count"] << endl;			
		}
	}

	filebuf fb;
 	fb.open ("STR", ios::out);
 	ostream os(&fb);
	bgi::detail::rtree::utilities::print(os, rtree);
 	fb.close();

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
