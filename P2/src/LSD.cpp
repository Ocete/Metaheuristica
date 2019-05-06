#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>
#include <stdlib.h>     // srand y rand
#include <algorithm>    // sort
#include <time.h>

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

template <class T>
double singleContribution(T &container, vector<vector<double> > &mat, int elem) {
  double result = 0;
  typename T::iterator it;
  for (it = container.begin(); it != container.end(); it++) {
    result += mat[ elem ][ *it ];
  }
  return result;
}

// Returns the fitness of the whole solution
template <class T>
double evaluateSolution(T &container, vector<vector<double> > &mat) {
  double fitness = 0;
  typename T::iterator it;
  for (it = container.begin(); it != container.end(); it++) {
    fitness += singleContribution(container, mat, *it);
  }
  // Counting twice all the possible distances
  return fitness /= 2;
}

////////////////// LOCAL SEARCH ///////////////////////

struct solution {
  vector<int> v;
  double fitness;
};

double updateSolution(solution &sol, vector<vector<double> > &mat) {
  sol.fitness = evaluateSolution(sol.v, mat);
  return sol.fitness;
}

// Creates a random solution
// Prec.: Random seed already set
void randomSolution (solution &sol, int size, int choosen) {
  int currently_choosen = 0, random;
  unordered_set<int> s;

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
    sol.v[i] = *it;
    i++;
  }
}

// Comparison operator for ordering the solution vector
// Keeps at the front the element with less contribution to the fitness
bool operator < (const pair<int,double> &p1, const pair<int,double> &p2) {
    return p1.second < p2.second;
}

// Order sol.v by conribution to the solution in ascending order
void orderSolutionByContribution(solution &sol, vector<vector<double> > &mat ) {
  pair<int, double> p (0, 0.0);
  vector< pair<int, double> > pairs_v ( sol.v.size(), p);

  // Initialize the auxiliar vector
  for (unsigned i=0; i< pairs_v.size(); i++) {
    pairs_v[i].first = sol.v[i];
    pairs_v[i].second = singleContribution(sol.v, mat, pairs_v[i].first);
  }

  // Order the vector by contribution
  sort(pairs_v.begin(), pairs_v.end());

  // Save the ordering
  for (unsigned i=0; i< pairs_v.size(); i++) {
    sol.v[i] = pairs_v[i].first;
  }
}

// Compute a good order to try the swaps
// Orders the element from (0 .. mat.size()) that are not in sol.v
// in ascending order of how they contribute to the solution
void obtainBestOrdering (vector<int> &best_ordering, solution &sol,
      vector<vector<double> > &mat) {
  set<int> sol_elements;

  // Put the elements the already appear on the solution into a set
  for (unsigned i=0; i<sol.v.size(); i++) {
    sol_elements.insert( sol.v[i] );
  }

  // Initialize the auxiliar vector
  pair<int, double> p (0, 0.0);
  int tam = mat.size() - sol.v.size();
  unsigned j = 0;
  vector< pair<int, double> > pairs_v (tam, p);

  for (unsigned i=0; i<mat.size(); i++) {
    if ( sol_elements.find(i) == sol_elements.end() ) {
      pairs_v[j].first = i;
      pairs_v[j].second = singleContribution(sol.v, mat, pairs_v[j].first);
      j++;
    }
  }

  // Order the vector by contribution
  sort(pairs_v.begin(), pairs_v.end());

  // Save the ordering
  best_ordering.resize(pairs_v.size());
  for (int i=pairs_v.size()-1; i>=0; i--) {
    best_ordering[i] = pairs_v[i].first;
  }
}

// Computes a single step in the exploration, changing "sol"
// The element chosen depending on how it contributes to the current solution
bool stepInNeighbourhoodDet (solution &sol, vector<vector<double> > &mat,
      int &evaluations, int MAX_EVALUATIONS) {
  unsigned i = 0, j, k, element_out, total_tries, max_i, max_k;
  double newContribution, oldContribution, percent_i;
  vector<int> best_ordering;
  int real_evaluations = 0;

  orderSolutionByContribution(sol, mat);
  obtainBestOrdering(best_ordering, sol, mat);

  percent_i = 0.1;
  total_tries = MAX_EVALUATIONS;

  // percent_i = 1;
  // total_tries = numeric_limits<int>::max();

  max_i = sol.v.size() * percent_i;
  max_k = min( total_tries / max_i, (unsigned) best_ordering.size());

  if ( max_k == best_ordering.size() ) {
    max_i = min(total_tries / max_k, (unsigned) sol.v.size());
  }

  // Fill hash with our used values
  unordered_set<int> s;
  for (unsigned i=0; i<sol.v.size(); i++) {
    s.insert( sol.v[i] );
  }

  // Explore the neighbourhood wisely and return the firstly found better option
  while (i < max_i) {
    // Save data of the element we are trying to swap
    oldContribution = singleContribution(sol.v, mat, sol.v[i]);
    element_out = sol.v[i];
    k = 0;
    while (k < max_k) {
      j = best_ordering[k];
      // Try the swap if the element 'j' is not in the current solution
      if ( s.find(j) == s.end() ) {
        newContribution = singleContribution(sol.v, mat, j) - mat[j][element_out];
        if ( newContribution > oldContribution ) {
          sol.v[i] = j;
          sol.fitness = sol.fitness + newContribution - oldContribution;
          s.erase(element_out);
          s.insert(j);
          return false;
        }
        real_evaluations++;
        if (real_evaluations % sol.v.size() == 0) {
          evaluations++;
        }
      }
      k++;
    }
    i++;
  }
  return true;
}

// Computes the local search algorithm for a random starting solution
double localSearchDet(vector<vector<double> > &mat, int choosen,
      int MAX_EVALUATIONS) {
  solution sol;
  bool stop = false;
  int evaluations = 0;
  clock_t t_start, t_total;

  // set seed
  srand (time(NULL));

  randomSolution(sol, mat.size(), choosen);
  updateSolution(sol, mat);

  t_start = clock();
  while (!stop && evaluations < MAX_EVALUATIONS) {
    stop = stepInNeighbourhoodDet(sol, mat, evaluations, MAX_EVALUATIONS-evaluations);
    // cout << sol.fitness << "\t" << iterations << endl;
  }
  t_total = clock() - t_start;

  // output: Fitness - Time - Iterations
  cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << " " << evaluations << endl;
  return sol.fitness;
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;
  int MAX_EVALUATIONS = 50000;

  cin >> size >> choosen;
  vector<double> v (size, 0);
  vector<vector<double > > mat (size, v);
  readInput(mat);

  // testEvaluation(mat, size, choosen);
  // testFactorization(mat, size, choosen);

  localSearchDet(mat, choosen, MAX_EVALUATIONS);
}
