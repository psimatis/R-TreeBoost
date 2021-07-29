#include <boost/geometry.hpp>
#include <boost/geometry/index/detail/rtree/utilities/statistics.hpp>
#include <boost/foreach.hpp>
#include <bits/stdc++.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;
using namespace std::chrono;
namespace bgid = bgi::detail;

#define CAPACITY 100

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

void parseDataFile(string fileName, vector<tuple<int, float, float>> &dataArray) {
	string line;
	//int i = 0;
	ifstream file(fileName);
	if (file.is_open()) {
		getline(file, line);
		while (getline(file, line)){
			int id;
			float x, y;
			istringstream buf(line);
			buf >> id >> x >> y;
			dataArray.emplace_back(make_tuple(id, x, y));
			//i++;
			//if (i >= 5000) break;
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
		}
		file.close();
	} else cout << "Cant open file: " << fileName << endl;
} 

int main(int argc, char** argv){

	if (argc != 3){
		cout << "Usage: ./rtree dataFile queryFile" << endl;
		exit(1);
	}
	
	typedef bg::model::point<float, 2, bg::cs::cartesian> point;
	typedef bg::model::box<point> box;
	typedef pair<point, unsigned> value;
	
	vector<float> boundary = {180.0, -90.0, 180.0, 90.0};
	
	vector<tuple<int, float, float>> dataArray;        
	parseDataFile(argv[1], dataArray);	 

	vector<tuple<char, float, float, float>> queryArray;        	 
	parseQueryFile(argv[2], queryArray);	 

	high_resolution_clock::time_point startTime = high_resolution_clock::now();
	bgi::rtree<value, bgi::rstar<CAPACITY>> rtree;
	for (auto q: dataArray){
		point p(get<2>(q), get<1>(q));
		rtree.insert(make_pair(p, get<0>(q)));
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

			vector<value> result_s;
			box queryBox(point(pl[0], pl[1]), point(ph[0], ph[1]));
			startTime = high_resolution_clock::now();
			rtree.query(bgi::intersects(queryBox), back_inserter(result_s));
			rangeLog["time " + to_string(rs)] += duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
			rangeLog["count " + to_string(rs)]++;
			
			cout << "Range query:" << bg::wkt<box>(queryBox) << endl;
			BOOST_FOREACH(value const& v, result_s)
				cout << bg::wkt<point>(v.first) << " - " << v.second << endl;
		}
		else if (get<0>(q) == 'k'){
			vector<value> result_n;
			point queryPoint(point(get<2>(q), get<1>(q)));
			startTime = high_resolution_clock::now();
			rtree.query(bgi::nearest(queryPoint, get<3>(q)), back_inserter(result_n));
			knnLog["time " + to_string(get<3>(q))] += duration_cast<duration<double>>(high_resolution_clock::now() - startTime).count();
			knnLog["count " + to_string(get<3>(q))]++;

			cout << "kNN query: " << bg::wkt<point>(point(queryPoint)) << endl;
			BOOST_FOREACH(value const& v, result_n)
				cout << bg::wkt<point>(v.first) << " - " << v.second << endl;
		}
	}

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

	pointerVisitor(rtree);
	cout << "Internal pointer count: " << pointerCount << endl;
	cout << "RTree size in MB (correct): " << ((get<1>(S) + get<2>(S)) * 4 * sizeof(float) // rectangle size
	 											+ pointerCount * 8) / float(1e6)<< endl; // internal pointer size
}
