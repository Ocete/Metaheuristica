#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <stdlib.h>     // srand y rand
#include <algorithm>    // sort y random_shuffle
#include <time.h>
#include <math.h>

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

// Order sol.v by conribution to the solution in ascending order
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

// Mutate the solution by swapping two random elements
void mutateSolution(solution &sol) {
  int r_on, r_off;

  // Find elements to swap
  do {
    r_on = random(0, sol.v.size());
  } while ( !sol.v[r_on] );

  do {
    r_off = random(0, sol.v.size());
  } while ( sol.v[r_off] || r_off == r_on);

  // Swap elements
  sol.v[r_on] = false;
  sol.v[r_off] = true;

  // Update the fitness if possible
  if ( sol.evaluated ) {
    double oldContribution = singleContribution(sol.v, r_on) - MAT[r_on][r_off];
    double newContribution = singleContribution(sol.v, r_off);
    sol.fitness = sol.fitness + newContribution - oldContribution;
    // iterations++; AQUI (mirar la nota importante en el main)
  }
}

void abruptMutation(solution &sol ) {
  int n_mutations = sol.v.size()*0.1;
  for (int i=0; i<n_mutations; i++) {
    mutateSolution(sol);
  }
}

int ES( solution &sol_init, int max_evaluations ) {
  double best_fitness = 0;
  int evaluations = 0;
  solution sol, saved_sol, best_sol;;
  double mu = 0.3, phi = 0.3, beta;
  double temp, start_temp, final_temp = 0.001, diff;
  int max_vecinos = 10*MAT.size();
  int M = max_evaluations / max_vecinos;

  // Inicialización
  sol = sol_init;
  saved_sol = sol;
  best_sol = sol;
  best_fitness = sol.fitness;
  start_temp = (mu*sol.fitness)/(-log(phi));
  temp = start_temp;

  if (start_temp <= final_temp) {
    cerr << "Error en temperatura inicial" << endl;
    return 0;
  }

  while ( temp > final_temp && evaluations < max_evaluations ) {
    for (int i=0; i<max_vecinos; i++) {
      mutateSolution(sol);
      evaluations++;

      diff = sol.fitness - saved_sol.fitness;
      if ( diff > 0 || rand() < exp(diff/temp) ) {
        saved_sol = sol;
        if ( sol.fitness > best_fitness ) {
          best_fitness = sol.fitness;
          best_sol = sol;
        }
      } else {
        sol = saved_sol;
      }
    }
    beta = (start_temp - temp) / (M * temp * start_temp);
    temp = temp / (1 + beta*temp);
  }

  sol_init = best_sol;
  return evaluations;
}


void ILS_ES( int choosen ) {
  int evaluations = 0, max_evaluations = 50000, total_tries = 25;
  int max_eval_ES = max_evaluations / total_tries;
  double best_fitness = 0;
  clock_t t_start, t_total;
  solution sol, saved_sol;

  t_start = clock();
  randomSolution(sol, choosen);
  evaluateSolution(sol);
  saved_sol = sol;
  best_fitness = sol.fitness;
  for (int i=0; i<total_tries; i++) {
    abruptMutation(sol);
    evaluations++;
    evaluations += ES( sol, max_eval_ES );

    if ( saved_sol.fitness > sol.fitness ) {
      sol = saved_sol;
    } else {
      saved_sol = sol;
      if ( sol.fitness > best_fitness ) {
        best_fitness = sol.fitness;
      }
    }
  }
  t_total = clock() - t_start;

  // output: Fitness - Time - Iterations
  cout << best_fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << evaluations << endl;
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;

  // set seed
  srand (time(NULL));

  cin >> size >> choosen;
  readInput(size);

  ILS_ES(choosen);
  // NOTA IMPORTANTE de cara a la memoria: no estamos realizando
  // el aumento de evaluaciones en las mutaciones como haciamos en
  // los genéticos, comentar esto en la memoria.

}