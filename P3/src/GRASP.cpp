#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <stdlib.h>     // srand y rand
#include <algorithm>    // sort y random_shuffle
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

////////////////// LOCAL SEARCH ///////////////////////

// Returns a random int in [a,b)
int random(int a, int b) {
  return a + rand() % b;
}

template <class T>
double singleContribution(T &container, int elem) {
  double result = 0;
  typename T::iterator it;
  for (it = container.begin(); it != container.end(); it++) {
    result += MAT[ elem ][ *it ];
  }
  return result;
}

// Returns the fitness of the whole solution
template <class T>
double evaluateSolution(T &container) {
  double fitness = 0;
  typename T::iterator it;
  for (it = container.begin(); it != container.end(); it++) {
    fitness += singleContribution(container, *it);
  }
  // Counting twice all the possible distances
  return fitness /= 2;
}

// Comparison operator for ordering the solution vector
// Keeps at the front the element with less contribution to the fitness
bool operator < (const pair<int,double> &p1, const pair<int,double> &p2) {
    return p1.second < p2.second;
}

struct solution_int {
  vector<int> v;
  double fitness;
};

double updateSolution(solution_int &sol) {
  sol.fitness = evaluateSolution(sol.v);
  return sol.fitness;
}

// Order sol.v by contribution to the solution in ascending order
void orderSolutionByContribution(solution_int &sol) {
  pair<int, double> p (0, 0.0);
  vector< pair<int, double> > pairs_v ( sol.v.size(), p);

  // Initialize the auxiliar vector
  for (unsigned i=0; i< pairs_v.size(); i++) {
    pairs_v[i].first = sol.v[i];
    pairs_v[i].second = singleContribution(sol.v, pairs_v[i].first);
  }

  // Order the vector by contribution
  sort(pairs_v.begin(), pairs_v.end());

  // Save the ordering
  for (unsigned i=0; i< pairs_v.size(); i++) {
    sol.v[i] = pairs_v[i].first;
  }
}

// Computes a single step in the exploration, changing "sol"
bool stepInNeighbourhood (solution_int &sol, int &evaluations, int MAX) {
  double percentage_studied;
  unsigned i = 0, j, element_out, total_tries, max_i, max_randoms, k;
  double newContribution, oldContribution;
  int real_evaluations = 0;

  orderSolutionByContribution(sol);

  // percentage_studied = 0.1;
  // total_tries = 50000;

  percentage_studied = 0.1;
  total_tries = MAX;

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
    oldContribution = singleContribution(sol.v, element_out);

    k = 0;
    j = rand() % MAT.size();
    while (j < MAT.size() && k < max_randoms) {
      // Try the swap if the element 'j' is not in the current solution
      if ( s.find(j) == s.end() ) {
        newContribution = singleContribution(sol.v, j) - MAT[j][element_out];
        if ( newContribution > oldContribution ) {
          sol.v[i] = j;
          sol.fitness = sol.fitness + newContribution - oldContribution;
          return false;
        }
        k++;

        real_evaluations++;
        if (real_evaluations % sol.v.size() == 0) {
          evaluations++;
          real_evaluations = 0;
        }
      }
      j = rand() % MAT.size();
    }
    i++;
  }
  return true;
}

struct solution {
  vector<bool> v;
  double fitness;
  bool evaluated;
};

double evaluateSolution(solution &sol) {
  sol.fitness = evaluateContainer(sol.v);
  sol.evaluated = true;
  return sol.fitness;
}

void BitsToInt(solution &sol_bits, solution_int &sol) {
  sol.v.clear();
  for(unsigned i=0; i<sol_bits.v.size(); i++) {
    if ( sol_bits.v[i] ) {
      sol.v.push_back(i);
    }
  }
}

void IntToBits(solution_int &sol, solution &sol_bits, int tam) {
  sol_bits.v = vector<bool> (tam, false);
  for(unsigned i=0; i<sol.v.size(); i++) {
    sol_bits.v[ sol.v[i] ] = true;
  }
}

// Computes the local search algorithm for a random starting solution
int localSearch(solution &sol_bits, int MAX_EVALUATIONS) {
  int tam_sol_bits = sol_bits.v.size();
  solution_int sol;
  bool stop = false;
  int evaluations = 0;

  BitsToInt(sol_bits, sol);
  updateSolution(sol);

  while (!stop && evaluations < MAX_EVALUATIONS) {
    stop = stepInNeighbourhood(sol, evaluations, MAX_EVALUATIONS-evaluations);
    // cout << sol.fitness << "\t" << iterations << endl;
  }

  // output: Fitness - Time - Iterations
  // cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << iterations<< endl;
  IntToBits(sol, sol_bits, tam_sol_bits);
  return evaluations;
}

/////////////////////////// END LOCAL SEARCH //////////////////////////////////
///////////////////////////   RAND  GREEDY     ///////////////////////////

// Returns the element from "non_selected" that is the farthest from "selected"
int farthestToSel(set<int> &non_selected, set<int> &selected) {
  int farthest;
  double max_sum_dist, current_sum_dist;
  set<int>::iterator it;

  it = non_selected.begin();
  farthest = *it;
  max_sum_dist = singleContribution(selected, farthest);

  for ( ; it!=non_selected.end(); it++) {
    current_sum_dist = singleContribution(selected, *it);
    if (current_sum_dist > max_sum_dist) {
      max_sum_dist = current_sum_dist;
      farthest = *it;
    }
  }

  return farthest;
}

// Returns an element from selected. It is quite randomly selected
// between the best possible elements.
int farthestToSelRandom(set<int> &non_selected, set<int> &selected, double alfa) {
  double current_sum_dist;
  set<int>::iterator it;
  vector<pair<int,double> > v;

  for ( it = non_selected.begin(); it!=non_selected.end(); it++) {
    current_sum_dist = singleContribution(selected, *it);
    v.push_back(pair<int,double> (*it, current_sum_dist));
  }

  // sorts in ascending order
  sort(v.begin(), v.end());
  // reverse the sorting
  reverse(v.begin(), v.end());
  double dist_min = v[ v.size()-1 ].second, dist_max = v[0].second;
  double umbral = dist_min + alfa*(dist_max - dist_min);

  unsigned i_umbral = 0;
  while ( i_umbral < v.size() && v[i_umbral].second >= umbral ) {
    i_umbral++;
  }

  return v[ random(0,i_umbral+1) ].first;
}

// Returns the element that is the farthest from the rest of elements
int farthestToAll() {
  set<int> all_elements;
  for (unsigned i=0; i<MAT.size(); i++) {
    all_elements.insert(i);
  }
  return farthestToSel( all_elements, all_elements);
}

void randGreedy( solution &sol, unsigned choosen, double alfa ) {
  set<int> non_selected, selected;
  int farthest;

  // Initialize selected with the farthestToAll and non_selected with the rest
  for (unsigned i=0; i<MAT.size(); i++) {
    non_selected.insert(i);
  }
  farthest = farthestToAll();
  non_selected.erase( farthest );
  selected.insert( farthest );

  while( selected.size() < choosen ) {
    farthest = farthestToSelRandom(non_selected, selected, alfa);
    non_selected.erase( farthest );
    selected.insert( farthest );
  }

  solution_int sol_i;
  for ( int s : selected ) {
    sol_i.v.push_back(s);
  }
  IntToBits(sol_i, sol, MAT.size());
}

///////////////////////////////////////////////////////////////////////////////

// Creates a random solution
// Prec.: Random seed already set
void randomSolution (solution &sol, int choosen) {
  int size = MAT.size();
  int currently_choosen = 0, rand_n;

  // Set the flag and clear the solution
  sol.evaluated = false;
  sol.fitness = 0;
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

void GRASP( int choosen, int MAX_EVALUATIONS ) {
  int evaluations = 0, max_evaluations = 50000, total_tries = 25;
  double best_fitness = 0, alfa = 0.3;
  clock_t t_start, t_total;

  t_start = clock();
  best_fitness = 0;
  for (int i=0; i<total_tries; i++) {
    solution sol;
    randGreedy(sol, choosen, alfa);
    evaluateSolution(sol);
    evaluations += localSearch( sol, max_evaluations );
    if ( sol.fitness > best_fitness ) {
      best_fitness = sol.fitness;
    }
  }
  t_total = clock() - t_start;

  // output: Fitness - Time - Iterations
  cout << best_fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << evaluations << endl;
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;
  int MAX_EVALUATIONS = 50000;

  // set seed
  srand (time(NULL));

  cin >> size >> choosen;
  readInput(size);

  GRASP(choosen, MAX_EVALUATIONS);
}
