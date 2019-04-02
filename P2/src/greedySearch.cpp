#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>
#include <stdlib.h>     // srand y rand
#include <algorithm>    // sort
#include <time.h>

using namespace std;

struct solution {
  vector<int> v;
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

/////////////////// GREEDY //////////////////////

// Returns the element from "non_selected" that is the farthest from "selected"
int farthestToSel(set<int> &non_selected, set<int> &selected, vector<vector<double> > &mat) {
  int farthest;
  double max_sum_dist, current_sum_dist;
  set<int>::iterator it;

  it = non_selected.begin();
  farthest = *it;
  max_sum_dist = singleContribution(selected, mat, farthest);

  for ( ; it!=non_selected.end(); it++) {
    current_sum_dist = singleContribution(selected, mat, *it);
    if (current_sum_dist > max_sum_dist) {
      max_sum_dist = current_sum_dist;
      farthest = *it;
    }
  }

  return farthest;
}

// Returns the element that is the farthest from the rest of elements
int farthestToAll(vector<vector<double> > &mat) {
  set<int> all_elements;
  for (unsigned i=0; i<mat.size(); i++) {
    all_elements.insert(i);
  }
  return farthestToSel( all_elements, all_elements, mat);
}

void greedy(solution &sol, vector<vector<double> > &mat, unsigned choosen) {
  set<int> non_selected, selected;
  int farthest;

  // Initialize selected with the farthestToAll and non_selected with the rest
  for (unsigned i=0; i<mat.size(); i++) {
    non_selected.insert(i);
  }
  farthest = farthestToAll(mat);
  non_selected.erase( farthest );
  selected.insert( farthest );

  while( selected.size() < choosen ) {
    farthest = farthestToSel(non_selected, selected, mat);
    non_selected.erase( farthest );
    selected.insert( farthest );
  }

  sol.fitness = evaluateSolution(selected, mat);

  sol.v.resize( selected.size() );
  set<int>::iterator it = selected.begin();
  for (unsigned i=0; i<sol.v.size(); i++) {
    sol.v[i] = *it;
    it++;
  }
}

////////////////// GREEDY SEARCH ///////////////////////

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

// Computes a single step in the exploration, changing "sol"
bool stepInNeighbourhood (solution &sol, vector<vector<double> > &mat) {
  double percentage_studied;
  unsigned i = 0, j, element_out, total_tries, max_i, max_randoms, k;
  double newContribution, oldContribution;

  orderSolutionByContribution(sol, mat);

  percentage_studied = 0.1;
  total_tries = 50000;

  // percentage_studied = 1;
  // total_tries = 1000000;

  max_i = max(percentage_studied * sol.v.size(), 1.0);
  max_randoms = total_tries / max_i;

  // Fill hash with our used values
  unordered_set<int> s;
  for (unsigned i=0; i<sol.v.size(); i++) {
    s.insert( sol.v[i] );
  }

  // Explore the neighbourhood and return the firstly found better option
  while (i < max_i) {
    // Save data of the element we are trying to swap
    element_out = sol.v[i];
    oldContribution = singleContribution(sol.v, mat, element_out);

    k = 0;
    j = rand() % mat.size();
    while (j < mat.size() && k < max_randoms) {
      // Try the swap if the element 'j' is not in the current solution
      if ( s.find(j) == s.end() ) {
        newContribution = singleContribution(sol.v, mat, j) - mat[j][element_out];
        if ( newContribution > oldContribution ) {
          sol.v[i] = j;
          sol.fitness = sol.fitness + newContribution - oldContribution;
          return false;
        }
        k++;
      }
      j = rand() % mat.size();
    }
    i++;
  }
  return true;
}

// Computes the local search algorithm for a random starting solution
double greedySearch(vector<vector<double> > &mat, int choosen) {
  solution sol;
  bool stop = false;
  int iterations = -1;
  clock_t t_start, t_total;

  // set seed
  srand (time(NULL));

  t_start = clock();
  // Initialize the solution with the greedy
  greedy(sol, mat, choosen);

  // Apply local search
  while (!stop && iterations < 10000) {
    stop = stepInNeighbourhood(sol, mat);
    iterations++;
    // cout << sol.fitness << "\t" << iterations << endl;
  }
  t_total = clock() - t_start;

  // output: Fitness - Time - Iterations
  cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << iterations<< endl;
  return sol.fitness;
}

/////////////// TEST ////////////////////

void testEvaluation(vector<vector<double> > &mat, int size, int choosen) {
  solution sol;
  for (unsigned i=0; i<10000; i++) {
    randomSolution(sol, size, choosen);
    updateSolution(sol, mat);
    cout << "Iteration " << i << " - fitness " << sol.fitness << endl;
  }
}

void testFactorization(vector<vector<double> > &mat, int size, int choosen) {
  solution sol;
  double fit_before, fit_after;
  bool OK = true;
  for (unsigned i=0; i<100; i++) {
    randomSolution(sol, size, choosen);
    updateSolution(sol, mat);

    for (unsigned i=0; i<1000; i++) {
      stepInNeighbourhood(sol, mat);
      fit_before = sol.fitness;
      fit_after = updateSolution(sol,mat);

      if ( abs(fit_after - fit_before) > 0.0001 ) {
        cerr << "ERROR in stepInNeighbourhood. Fit before: " << fit_before;
        cerr << " Fit after: " << fit_after << endl;
        OK = false;
      }
    }
  }
  if ( OK ) {
    cout << "No error found in stepInNeighbourhood" << endl;
  }
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;

  cin >> size >> choosen;
  vector<double> v (size, 0);
  vector<vector<double > > mat (size, v);
  readInput(mat);

  // testEvaluation(mat, size, choosen);
  // testFactorization(mat, size, choosen);

  greedySearch(mat, choosen);
}
