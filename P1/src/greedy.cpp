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

/////////////////// GREEDY //////////////////////

// Returns the element that is the farthest from the rest of elements
int farthestToAll(vector<vector<double> > &mat) {
  int farthest;
  double max_sum_dist, current_sum_dist;

  vector<int> all_elements (0, mat.size());
  for (unsigned i=0; i<all_elements.size(); i++) {
    all_elements[i] = i;
  }
  farthest = 0;
  max_sum_dist = singleContribution(all_elements, mat, 0);

  for (unsigned i=1; i<mat.size(); i++) {
    current_sum_dist = singleContribution(all_elements, mat, i);
    if (current_sum_dist > max_sum_dist) {
      max_sum_dist = current_sum_dist;
      farthest = i;
    }
  }

  return farthest;
}

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

void greedy(vector<vector<double> > &mat, unsigned choosen) {
  set<int> non_selected, selected;
  int farthest;
  clock_t t_start, t_total;

  t_start = clock();
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
  t_total = clock() - t_start;

  double fitness = evaluateSolution(selected, mat);
  // output: Fitness - Time
  cout << fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << endl;
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

  greedy(mat, choosen);
}
