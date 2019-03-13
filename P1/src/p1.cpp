#include <iostream>
#include <vector>
#include <unordered_set>
#include <stdlib.h>     // srand y rand
#include <algorithm>    // sort

using namespace std;

struct solution {
  vector< pair<int, double> > v;
  double fitness;
};

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

/////////////////// GREEDY //////////////////////

void greedy(vector<int> &sol, vector<vector<double> > &mat) {

}

////////////////// EVALUATION ///////////////////////

// Computes the contribution of the element "elem" for the solution "sol" given
// AKA, the sum of the distances from elemet "elem" to each other element in sol
double singleContribution(solution &sol, vector<vector<double> > &mat, int elem) {
  double result = 0;
  for (unsigned i=0; i<sol.v.size(); i++) {
    result += mat[ elem ][ sol.v[i].first ];
  }
  return result;
}

// Updates sol.fitness and the contribution of each element (sol.v[i].second)
double evaluateSolution(solution &sol, vector<vector<double> > &mat) {
  sol.fitness = 0;

  for (unsigned i=0; i<sol.v.size(); i++) {
    sol.v[i].second = singleContribution(sol, mat, sol.v[i].first);
    sol.fitness += sol.v[i].second;
  }

  // Counting twice all the possible distances
  sol.fitness /= 2;
  return sol.fitness;
}

////////////////// LOCAL SEARCH ///////////////////////

// Creates a random solution
void randomSolution (solution &sol, int size, int choosen) {
  int currently_choosen = 0, random;
  unordered_set<int> s;
  // set seed
  srand (time(NULL));

  // Select 'choosen' elements
  while (currently_choosen  < choosen) {
    random = rand() % size;
    if ( s.find(random) == s.end() ) {
      s.insert(random);
      currently_choosen++;
    }
  }

  // Dump the set into the solution
  int i = 0;
  sol.v.resize(choosen);
  for (auto it = s.begin(); it != s.end(); it++) {
    sol.v[i].first = *it;
    i++;
  }
}

// Comparison operator for ordering the solution vector
// Keeps at the front the element with less contribution to the fitness
bool operator < (const pair<int,double> &p1, const pair<int,double> &p2) {
    return p1.second < p2.second;
}

// Computes a single step in the exploration, changing "sol"
bool stepInNeighbourhood (solution &sol, vector<vector<double> > &mat) {
  bool noImprov = true;
  unsigned i = 0, j, element_out;
  float newContribution, oldContribution;

  sort(sol.v.begin(), sol.v.end());

  // Fill hash with our used values
  unordered_set<int> s;
  for (unsigned i=0; i<sol.v.size(); i++) {
    s.insert( sol.v[i].first );
  }

  // Explore the neighbourhood and return the first better option found
  while (noImprov && i < sol.v.size()) {
    // Save data of the element we are trying to swap
    oldContribution = sol.v[i].second;
    element_out = sol.v[i].first;

    j = 0;
    while (noImprov && j < mat.size()) {
      // Try the swap if the element 'j' is not in the current solution
      if ( s.find(j) == s.end() ) {
        newContribution = singleContribution(sol, mat, j) - mat[j][element_out];
        if ( newContribution > oldContribution ) {
          sol.v[i].first = j;
          sol.v[i].second = newContribution;
          sol.fitness = sol.fitness + newContribution - oldContribution;
          noImprov = false;
        }
      }
      j++;
    }

    i++;
  }
  return noImprov;
}

// Computes the local search algorithm for a random starting solution
double localSearch(vector<vector<double> > &mat, int choosen) {
  solution sol;
  bool stop = false;
  int iterations = 0;

  randomSolution(sol, mat.size(), choosen);
  evaluateSolution(sol, mat);

  while (!stop && iterations < 10000) {
    stop = stepInNeighbourhood(sol, mat);
    iterations++;
  }

  cout << "Iterations: " << iterations << endl;
  return sol.fitness;
}

/////////////// TEST ////////////////////

void printMat (vector<vector<double> > &mat) {
  for (unsigned i=0; i<mat.size(); i++) {
    for (unsigned j=0; j<mat.size(); j++) {
      cout << mat[i][j] << " ";
    }
    cout << endl;
  }
}

template <class T>
void printVector (vector<T> &v) {
  for (unsigned i=0; i<v.size(); i++) {
    cout << v[i] << " ";
  }
  cout << endl;
}

void testEvaluation(vector<vector<double> > &mat, int size, int choosen) {
  solution sol;
  float fitness;
  for (unsigned i=0; i<10000; i++) {
    randomSolution(sol, size, choosen);
    fitness = evaluateSolution(sol, mat);
    cout << "Iteration " << i << " - fitness " << fitness << endl;
  }
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;
  double fitness;

  cin >> size >> choosen;
  vector<double> v (size, 0);
  vector<vector<double > > mat (size, v);
  readInput(mat);

  //testEvaluation(mat, size, choosen);

  fitness = localSearch(mat, choosen);
  cout << "Local search fitness: " << fitness << endl;
}
