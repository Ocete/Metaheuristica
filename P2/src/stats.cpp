#include <iostream>
#include <vector>
#include <string>

using namespace std;

/////////////////// INPUT //////////////////////

void readInput (vector<vector<double> > &mat) {
  unsigned i, j;
  double f;
  for (i=1; i<mat.size(); i++) {
    for (j=i+1; j<mat.size(); j++) {
      cin >> f >> f >> f;
      mat[i][j] = mat[j][i] = f;
    }
  }
}

void inputNames(vector<string> &input_names) {
  string str;
  int k;
  vector<string> input_names_bases = {"GKD-c_","SOM-b_","MDG-a_"};

  input_names.clear();
  for (int i=1; i<=30; i++) {
    k = (i-1)/10;
    str = input_names_bases[k] + to_string(i);
    input_names.push_back(str);
  }
}

//////////////////// MAIN /////////////////////

void computeStats(int iterations, string input) {
  double mean_fitness = 0, mean_time = 0, time_read, fitness_read, max_fitness = 0;
  int mean_steps = 0, steps_read;

  for (int i=0; i<iterations; i++) {
    cin >> fitness_read >> time_read >> steps_read;
    mean_fitness += fitness_read;
    max_fitness = max(max_fitness, fitness_read);
    mean_time += time_read;
    mean_steps += steps_read;
  }
  cout << mean_fitness/iterations << "\t";
  // cout << max_fitness << "\t";
  cout << mean_time/iterations << "\t";
  cout << (double) mean_steps/iterations << endl;
}

int main( int argc, char *argv[] ) {
  int iterations;
  string str;

  cin >> iterations;
  vector<string> input_names;
  inputNames(input_names);

  for (unsigned i=0; i<input_names.size(); i++) {
    computeStats(iterations, input_names[i]);
  }
}
