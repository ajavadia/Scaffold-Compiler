#include <iostream>   //std::cout, std::cin
#include <string>
#include <stdlib.h>
#include <algorithm>  //std::random_shuffle
#include <vector>
#include <random>
#include <assert.h>
#include <cstring>    //strcmp

int n; // problem size
int P; // parallelism factor
int M; // 

using namespace std;

// draw 2M random numbers from the range [0...width-1]. 
vector<int> draw_2M (int C, int width) {
  assert(2*C <= width);
  vector<int> result (2*C);
  if (2*C < width/2) {  // this method is more efficient for small 2P
    for (int i=0; i<2*C; i++) {
      int r;
      do {
        r = rand()%width;
      } while (find(result.begin(), result.begin()+i+1, r)!=result.begin()+i+1);
      result[i] = r;
    }
  }
  else {                // this method is more efficient for large 2P
    vector<int> shuffle (width);
    for (int i=0; i<width; i++)
      shuffle[i] = i;
    random_shuffle(shuffle.begin(), shuffle.end());
    shuffle.resize(2*C);
    result = shuffle;
  }
  return result;
}

int main(int argc, char *argv[]) {

  n = -1; P = -1;
  
  for (int i = 0; i<argc; i++) {
    if (strcmp(argv[i],"--n")==0)
      n = atoi(argv[i+1]);
    if (strcmp(argv[i],"--P")==0)
      P = atoi(argv[i+1]);
  }

  int width = n*n;
  int depth = n*n*n*n;
  int repeat = n*n;
  int unique_depth = depth/repeat;

  random_device rd;
  mt19937 e(rd());    
  normal_distribution<> dist(P, P/4);

  cout << "int const sd_pairs[unique_depth][2*P][2] = " << endl;
  cout << "{";
  for (int i=0; i<unique_depth; i++) {  
    int C;  // how many CNOTs for this individual stage
    do {
      C = round(dist(e));    
    } while (2*C > width);
    vector<int> rand_pairs = draw_2M (C, width);      
    cout << "{";    
    for (int j=0; j<C; j++) {
      cout << "{" << rand_pairs[2*j] << "," << rand_pairs[2*j+1] << "}";
      if (j<2*P-1) cout << ",";
    }
    for (int k=C; k<2*P; k++) {
      cout << "{-1,-1}";
      if (k<2*P-1) cout << ",";      
    }
    cout << "}";                
    if (i<unique_depth-1) cout << ",\n";
  }
  cout << "};\n";
  
  return 0;
}
