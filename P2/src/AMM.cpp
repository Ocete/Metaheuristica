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

// Compute a good order to try the swaps
// Orders the element from (0 .. mat.size()) that are not in sol.v
// in ascending order of how they contribute to the solution
void obtainBestOrdering (vector<int> &best_ordering, solution_int &sol) {
  set<int> sol_elements;

  // Put the elements the already appear on the solution into a set
  for (unsigned i=0; i<sol.v.size(); i++) {
    sol_elements.insert( sol.v[i] );
  }

  // Initialize the auxiliar vector
  pair<int, double> p (0, 0.0);
  int tam = MAT.size() - sol.v.size();
  unsigned j = 0;
  vector< pair<int, double> > pairs_v (tam, p);

  for (unsigned i=0; i<MAT.size(); i++) {
    if ( sol_elements.find(i) == sol_elements.end() ) {
      pairs_v[j].first = i;
      pairs_v[j].second = singleContribution(sol.v, pairs_v[j].first);
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
bool stepInNeighbourhoodDet (solution_int &sol, int &evaluations, int MAX) {
  unsigned i = 0, j, k, element_out, total_tries, max_i, max_k;
  double newContribution, oldContribution, percent_i;
  vector<int> best_ordering;
  int real_evaluations = 0;

  orderSolutionByContribution(sol);
  obtainBestOrdering(best_ordering, sol);

  percent_i = 0.1;
  total_tries = MAX*sol.v.size();

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
    oldContribution = singleContribution(sol.v, sol.v[i]);
    element_out = sol.v[i];
    k = 0;
    while (k < max_k) {
      j = best_ordering[k];
      // Try the swap if the element 'j' is not in the current solution
      if ( s.find(j) == s.end() ) {
        newContribution = singleContribution(sol.v, j) - MAT[j][element_out];
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
          real_evaluations = 0;
        }
      }
      k++;
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
int localSearchDet(solution &sol_bits, int MAX_EVALUATIONS ) {
  solution_int sol;
  bool stop = false;
  int iterations = -1, evaluations = 0;
  int tam_sol_bits = sol_bits.v.size();

  BitsToInt(sol_bits, sol);
  updateSolution(sol);

  while (!stop && evaluations < MAX_EVALUATIONS) {
    stop = stepInNeighbourhoodDet(sol, evaluations, MAX_EVALUATIONS-evaluations);
    iterations++;
    // cout << sol.fitness << "\t" << iterations << endl;
  }

  // output: Fitness - Time - Iterations
  // cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << iterations<< endl;

  IntToBits(sol, sol_bits, tam_sol_bits);
  return evaluations;
}

///////////////////////////// GENETIC ///////////////////////////////////////

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

class population {
public:
  vector<solution> v;
  double max_fitness;
  int best_sol;
  int tam;

  population() {
    best_sol = 0;
    tam = 0;
    max_fitness = 0;
  }
};

void initializePop(population &pop, int tam_pob, int choosen) {
  pop.v.resize(tam_pob);
  for (int i=0; i<tam_pob; i++) {
    randomSolution(pop.v[i], choosen);
  }
  pop.tam = pop.v.size();
}

void evaluatePop(population &pop, int &iterations) {
  for (int i=0; i<pop.tam; i++) {
    if ( !pop.v[i].evaluated ) {
      evaluateSolution( pop.v[i] );
      iterations++;
    }
    if ( pop.v[i].fitness > pop.max_fitness ) {
      pop.best_sol = i;
      pop.max_fitness = pop.v[i].fitness;
    }
  }
}

// Prec: Assumes the population is evaluated
int binaryTournament(population &pop) {
  int r1 = random(0, pop.tam);
  int r2 = random(0, pop.tam);
  r1 = pop.v[r1].fitness > pop.v[r2].fitness ? r1 : r2;
  return  r1;
}

void selection(population &pop, population &new_pop) {
  new_pop.tam = pop.tam;
  new_pop.max_fitness = 0;
  int select;
  new_pop.v.resize(new_pop.tam);
  for (int i=0; i<pop.tam; i++) {
    select = binaryTournament ( pop );
    new_pop.v[i] = pop.v[select];
    if ( new_pop.v[i].fitness > new_pop.max_fitness ) {
      new_pop.best_sol = i;
      new_pop.max_fitness = pop.v[i].fitness;
    }
  }
}

// Mutate the solution by swapping two random elements
void mutateSolution(solution &sol, int &iterations) {
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
    iterations++;
  }
}

void mutatePop(population &pop, double &mut_prob,
      int choosen, int &iterations) {
  int r_sol=0, n_mut = mut_prob*pop.tam;
  if (n_mut > 0) {
    for (int i=0; i<n_mut; i++) {
      r_sol = random(0, pop.tam);
      mutateSolution( pop.v[r_sol], iterations );
    }
  } else {
    if ( rand() < mut_prob*pop.tam ) {
      r_sol = random(0, pop.tam);
      mutateSolution( pop.v[r_sol], iterations );
    }
  }
  if (pop.max_fitness < pop.v[r_sol].fitness) {
    pop.max_fitness = pop.v[r_sol].fitness;
    pop.best_sol = r_sol;
  }
}

// TODO implementar el operador que realmente hay que implementar
void repairSolution(solution &sol, int n_chosen) {
  int actualy_chosen = 0, r;

  for (unsigned i=0; i<sol.v.size(); i++) {
    if ( sol.v[i] ) {
      actualy_chosen++;
    }
  }

  while (actualy_chosen > n_chosen) {
    r = random(0, sol.v.size());
    if ( sol.v[r] ) {
      sol.v[r] = false;
      actualy_chosen--;
    }
  }

  while (actualy_chosen < n_chosen) {
    r = random(0, sol.v.size());
    if ( !sol.v[r] ) {
      sol.v[r] = true;
      actualy_chosen++;
    }
  }
}

// Operador de cruce uniforme
solution uniformCross(solution &p1, solution &p2) {
  solution child = p1;
  child.evaluated = false;
  int n_chosen = 0;

  for (unsigned i=0; i<p1.v.size(); i++) {
    if ( p1.v[i] ) {
      n_chosen++;
    }
    if ( p1.v[i] && p2.v[i] ) {
      child.v[i] = true;
    } else if ( !p1.v[i] && !p2.v[i] ) {
      child.v[i] = false;
    } else {
      child.v[i] = random(0,2) == 0;
    }
  }

  repairSolution (child, n_chosen);
  return child;
}

void cross(population &pop, double cross_prob) {
  int n_crosses = cross_prob * pop.tam / 2 ;
  int p1, p2;
  solution ch1, ch2;

  for (int i=0; i<n_crosses; i++) {
    p1 = random(0, pop.tam);
    do {
      p2 = random(0, pop.tam);
    } while (p1 == p2);
    ch1 = uniformCross( pop.v[p1], pop.v[p2]);
    ch2 = uniformCross( pop.v[p1], pop.v[p2]);
    pop.v[p1] = ch1;
    pop.v[p2] = ch2;
  }
}

void replace(population &pop, population &new_pop) {
  solution saved_best;
  bool elitist = false;
  int i_change;

  if ( new_pop.max_fitness < pop.max_fitness ) {
    elitist = true;
    i_change = new_pop.best_sol == 0 ? 1 : 0;
    saved_best = pop.v[ pop.best_sol ];
  }

  pop.v.swap(new_pop.v);

  if ( elitist ) {
    new_pop.v[ i_change ] = saved_best;
  }
}

// Tipo 0 = Poblacion entera cada 10 iteraciones
// Tipo 1 = Cada elemento con probabilidad p_mem
// Tipo 2 = Unicamente el mejor de la pob
void memetize(population &pop, int mem_type, int max_iterations,
    int &evaluations, double p_mem) {
  int pop_tam = pop.v.size();
  if ( mem_type == 0 ) {
    for ( int i=0; i<pop_tam; i++) {
      evaluations += localSearchDet( pop.v[i], max_iterations );
      evaluateSolution( pop.v[i] );
      if  ( pop.v[i].fitness > pop.max_fitness ) {
        pop.max_fitness = pop.v[i].fitness;
        pop.best_sol = i;
      }
    }
  } else if ( mem_type == 1) {
    int n_mems = pop_tam * p_mem;
    if (n_mems >= 1) {
      for (int i=0; i<n_mems; i++) {
        int r_rand = random(0, pop_tam);
        evaluations += localSearchDet( pop.v[r_rand], max_iterations );
        evaluateSolution( pop.v[ r_rand ] );
        if  ( pop.v[ r_rand ].fitness > pop.max_fitness ) {
          pop.max_fitness = pop.v[ r_rand ].fitness;
          pop.best_sol = r_rand;
        }
      }
    } else {
      if ( rand() < p_mem ) {
        int r_rand = random(0, pop_tam);
        evaluations += localSearchDet( pop.v[r_rand], max_iterations );
        if  ( pop.v[ r_rand ].fitness > pop.max_fitness ) {
          pop.max_fitness = pop.v[ r_rand ].fitness;
          pop.best_sol = r_rand;
        }
      }
    }
  } else if ( mem_type == 2 ) {
    evaluations += localSearchDet( pop.v[ pop.best_sol ], max_iterations);
    evaluateSolution( pop.v[ pop.best_sol ] );
    pop.max_fitness = pop.v[ pop.best_sol ].fitness;
  }
}

// Computes the local search algorithm for a random starting solution
// This implementation doesn't assume the pop is ordered
double AMM( int choosen, int MAX_EVALUATIONS, int mem_type) {
  int evaluations = 0, tam_pob = 50, generations = 0, max_iterations = 400;
  double mut_prob = 0.001, cross_prob = 0.7, p_mem = 0.1;
  clock_t t_start, t_total;
  population new_pop, pop;

  // set seed
  srand (time(NULL));
  t_start = clock();

  initializePop(pop, tam_pob, choosen);
  evaluatePop(pop, evaluations);

  while (evaluations < MAX_EVALUATIONS) {
    generations++;
    selection(pop, new_pop);
    cross(new_pop, cross_prob);
    mutatePop(new_pop, mut_prob, choosen, evaluations);
    evaluatePop(new_pop, evaluations);
    replace(pop, new_pop);
    if (generations % 10 == 0) {
      memetize(pop, mem_type, max_iterations, evaluations, p_mem);
    }
  }
  t_total = clock() - t_start;

  solution sol = pop.v[ pop.best_sol ];

  // output: Fitness - Time - Iterations
  cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << generations << endl;
  return sol.fitness;
}

//////////////////// MAIN /////////////////////

int main( int argc, char *argv[] ) {
  int size, choosen;
  int MAX_EVALUATIONS = 50000;

  // set seed
  srand (time(NULL));

  cin >> size >> choosen;
  readInput(size);

  AMM(choosen, MAX_EVALUATIONS, 2);
}
