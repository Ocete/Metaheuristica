#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <stdlib.h>     // srand y rand
#include <algorithm>    // sort
#include <time.h>

using namespace std;

vector<vector<double> > MAT;

/////////////////// INPUT //////////////////////

void readInput ( int size ) {
  vector<double> v (size, 0);
  vector<vector<double > > mat (size, v);
  MAT = mat;
  unsigned i, j;
  double f;
  for (i=1; i<MAT.size(); i++) {
    for (j=i+1; j<MAT.size(); j++) {
      cin >> f >> f >> f;
      MAT[i][j] = MAT[j][i] = f;
    }
  }
}

/////////////////// OUTPUT //////////////////////

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

////////////////// EVALUATION ///////////////////////

// Computes the contribution of the element "elem" for the solution "sol" given
// AKA, the sum of the distances from elemet "elem" to each other element in sol
/*double singleContribution(solution &sol, vector<vector<double> > &mat, int elem) {
  double result = 0;
  for (unsigned i=0; i<sol.v.size(); i++) {
    result += mat[ elem ][ sol.v[i] ];
  }
  return result;
}
*/
// Computes the contribution of the element "elem" for the solution "sol" given
// AKA, the sum of the distances from elemet "elem" to each other element in sol

double singleContribution(vector<bool> &v, int elem) {
  double result = 0;
  for (unsigned i=0; i<v.size(); i++) {
    if ( v[i] ) {
      result += v[i] * MAT[ elem ][ i ];
    }
  }
  return result;
}

// Returns the fitness of the whole solution
double evaluateContainer(vector<bool> &v) {
  double fitness = 0;
  for (unsigned i=0; i<v.size(); i++) {
    if ( v[i] ) {
      fitness += singleContribution(v, i);
    }
  }
  // Counting twice all the possible distances
  return fitness /= 2;
}

////////////////// GENETIC ///////////////////////

struct solution {
  vector<bool> v;
  double fitness;
  bool evaluated;
};

// Returns a random int in [a,b)
int random(int a, int b) {
  return a + rand() % b;
}

double evaluateSolution(solution &sol) {
  sol.fitness = evaluateContainer(sol.v);
  sol.evaluated = true;
  return sol.fitness;
}

// Creates a random solution
// Prec.: Random seed already set
void randomSolution (solution &sol, int choosen) {
  int size = MAT.size();
  int currently_choosen = 0, rand_n;

  // Set the flag and clear the solution
  sol.evaluated = false;
  sol.v = vector<bool> (size, false);

  // Select 'choosen' elements
  while (currently_choosen  < choosen) {
    rand_n = random(0, size);
    if ( !sol.v[rand_n] ) {
      sol.v[rand_n] = true;
      currently_choosen++;
    }
  }
}

// Comparison operator for ordering the solution vector
// Keeps at the front the element with less contribution to the fitness
bool operator < (const pair<int,double> &p1, const pair<int,double> &p2) {
    return p1.second < p2.second;
}

// Order sol.v by conribution to the solution in ascending order
void orderSolutionByContribution ( solution &sol, vector<int> &ordered ) {
  vector< pair<int, double> > pairs_v;
  double cont;

  // Initialize the auxiliar vector
  for (unsigned i=0; i<sol.v.size(); i++) {
    if ( sol.v[i] ) {
      cont = singleContribution(sol.v, i);
      pairs_v.push_back( pair<int,double> (i, cont) );
    }
  }

  // Order the vector by contribution
  sort(pairs_v.begin(), pairs_v.end());

  // Save the ordering
  ordered.resize( pairs_v.size() );
  for (unsigned i=0; i< pairs_v.size(); i++) {
    ordered[i] = pairs_v[i].first;
  }
}

// Computes a single step in the exploration, changing "sol"
bool stepInNeighbourhood ( solution &sol, int &tries, int max_tries ) {
  double percentage_studied;
  unsigned i, j, element_out, max_i, max_randoms, k;
  double newContribution, oldContribution;

  vector<int> ordered;
  orderSolutionByContribution(sol, ordered);

  percentage_studied = 0.2;

  max_i = max(percentage_studied * sol.v.size(), 1.0);
  max_randoms = max_tries / max_i;

  // Explore the neighbourhood and return the firstly found better option
  i = 0;
  while (i < max_i && tries < max_tries) {
    // Save data of the element we are trying to swap
    element_out = ordered[i];
    oldContribution = singleContribution(sol.v, element_out);

    k = 0;
    j = random(0, MAT.size());
    while ( k < max_randoms && tries < max_tries ) {
      // Try the swap if the element 'j' is not in the current solution
      if ( !sol.v[j] ) {
        newContribution = singleContribution(sol.v, j) - MAT[j][element_out];
        tries++;
        if ( newContribution > oldContribution ) {
          sol.v[element_out] = false;
          sol.v[j] = true;
          sol.fitness = sol.fitness + newContribution - oldContribution;
          return false;
        }
        k++;
      }
      j = random(0, MAT.size());
    }
    i++;
  }
  return true;
}

// Computes the local search algorithm for a random starting solution
double localSearch(int choosen, int &tries, int max_tries) {
  solution sol;
  bool stop = false;
  clock_t t_start, t_total;

  randomSolution(sol, choosen);
  evaluateSolution(sol);

  t_start = clock();
  while (!stop && tries < max_tries) {
    stop = stepInNeighbourhood(sol, tries, max_tries);
    // cout << sol.fitness << "\t" << iterations << endl;
  }
  t_total = clock() - t_start;

  // output: Fitness - Time - Iterations
  cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << tries << endl;
  return sol.fitness;
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;
  int tries = 0, max_tries = 50000;

  // set seed
  srand (time(NULL));

  cin >> size >> choosen;
  readInput(size);

  // testEvaluation(mat, size, choosen);
  // testFactorization(mat, size, choosen);

  localSearch(choosen, tries, max_tries);
}
