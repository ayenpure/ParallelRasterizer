#include <string>
#include <iostream>
#include <unordered_map>

using namespace std;

struct vertex {
  double x,y,z;
  bool operator==(const vertex &other) const {
    return ( x==other.x &&
             y==other.y &&
             z==other.z );
  }
};

namespace std {
  template<> struct hash<vertex>
  {
    /*typedef vertex argument_type;
    typedef std::size_t result_type;*/
    size_t operator()(vertex const& v) const {
      size_t const h1 ( std::hash<double>{}(v.x) );
      size_t const h2 ( std::hash<double>{}(v.y) );
      size_t const h3 ( std::hash<double>{}(v.z) );
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };
}

int main(int argc, char *argv[]) {
  vertex v1 = {1,2,3};
  vertex v2 = {2,3,4};
  vertex v3 = {4,2,3};
  vertex v4 = {5,6,7};
  vertex v5 = {2,3,4};
  vertex v6 = {4,2,3};
  /*unordered_map<vertex, int> vertexMap;
  vertexMap.insert(std::make_pair(v1,1));
  vertexMap.insert(std::make_pair(v2,1));
  vertexMap.insert(std::make_pair(v3,1));
  vertexMap.insert(std::make_pair(v4,1));
  vertexMap.insert(std::make_pair(v5,1));
  vertexMap.insert(std::make_pair(v6,1));
  cout << "Final count : " << vertexMap.size() << endl;*/
  /*std::unordered_map<int, std::string> dict = {{1, "one"}, {2, "two"}};
    dict.insert({3, "three"});
    dict.insert(std::make_pair(4, "four"));
    dict.insert({{4, "another four"}, {5, "five"}});

    bool ok = dict.insert({1, "another one"}).second;
    std::cout << "inserting 1 -> \"another one\" "
              << (ok ? "succeeded" : "failed") << '\n';

    std::cout << "contents:\n";
    for(auto& p: dict)
        std::cout << " " << p.first << " => " << p.second << '\n';
    cout << "Final count : " << dict.size() << endl;*/

  unordered_map<vertex, int> vertexMap;
  vertexMap.insert(std::make_pair(v1,1));
  vertexMap.insert(std::make_pair(v2,2));
  vertexMap.insert(std::make_pair(v3,3));
  vertexMap.insert(std::make_pair(v4,4));
  vertexMap.insert(std::make_pair(v5,5));
  vertexMap.insert(std::make_pair(v6,6));
  cout << "Final count : " << vertexMap.size() << endl;
  cout << "Count : " << vertexMap.count(v2) << endl;
  auto search = vertexMap.find(v3);
  if( search != vertexMap.end() ) {
    cout << "Found value : " << search->second << endl;
  }
  return 0;
}
