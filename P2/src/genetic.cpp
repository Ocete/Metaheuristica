#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
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
double evaluateContainer(T &container, vector<vector<double> > &mat) {
  double fitness = 0;
  typename T::iterator it;
  for (it = container.begin(); it != container.end(); it++) {
    fitness += singleContribution(container, mat, *it);
  }
  // Counting twice all the possible distances
  return fitness /= 2;
}

////////////////// GENETIC ///////////////////////

struct solution {
  vector<int> v;
  double fitness;
  bool evaluated;
};

// Returns a random int in [a,b)
int random(int a, int b) {
  return a + rand() % b;
}

double evaluateSolution(solution &sol, vector<vector<double> > &mat) {
  sol.fitness = evaluateContainer(sol.v, mat);
  sol.evaluated = true;
  return sol.fitness;
}

// Creates a random solution
// Prec.: Random seed already set
void randomSolution (solution &sol, int size, int choosen) {
  int currently_choosen = 0, random;
  unordered_set<int> s;

  // Set the flag
  sol.evaluated = false;

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
void orderSolutionByContribution ( solution &sol, vector<vector<double> > &mat ) {
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

  // percentage_studied = 0.1;
  // total_tries = 50000;

  percentage_studied = 1;
  total_tries = 1000000;

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

void initializePop(vector<solution> &pop, int tam_pob, int size, int choosen) {
  pop.resize(tam_pob);
  for (unsigned i=0; i<tam_pob; i++) {
    randomSolution(pop[i], size, choosen);
  }
}

void evaluatePop(vector<solution> &pop, vector<vector<double> > &mat,
      int &iterations) {
  for (unsigned i=0; i<pop.size(); i++) {
    if ( !pop[i].evaluated ) {
      evaluateSolution(pop[i], mat);
      iterations++;
    }
  }
}

void selection(vector<solution> &pop, vector<solution> &new_pop) {
  unsigned pop_tam = pop.size();
  int r1, r2;
  new_pop.resize(pop_tam);
  for (unsigned i=0; i<n_crosses; i++) {
    r1 = random(0, pop_tam);
    r2 = random(0, pop_tam);
    r1 = min (r1,r2);
    new_pop[i] = pop[r1];
  }
}

bool mutateSolution(solution &sol, int choosen) {
  unordered_set<int> s;
  int r1, r2;
  bool evaluated;

  for ( i in sol.v ) {
    s.insert(i);
  }

  do {
    r_in = random(0, choosen);
    r_out = random(0, choosen);
  } while ( /* TODO COMPLETAR CON EL FIND DE LA CLASE VECTOR Y QUITAR EL SET */  );

  // Necesitamos la posicion para despues reemplazarlo. Se puede usar un map!!!!!!!

}

void mutatePop(vector<solution> &pop, double &mut_prob,
      int choosen, int &iterations) {
  int r_sol, n_mut = mut_prob*pop.size();
  bool evaluated;
  for (unsigned i=0; i<n_mut; i++) {
    r_sol = random(0,pop.size());
    evaluated = mutateSolution(pop[r_sol]);
    if ( evaluated ) {
      iterations++;
    }
  }
}

// Computes the local search algorithm for a random starting solution
double genetic(vector<vector<double> > &mat, int choosen) {
  int evaluations = 0, MAX_EVALUATIONS = 50000, tam_pob = 50, generations = 0;
  double mut_prob = 0.7, cross_prob = 0.001;
  clock_t t_start, t_total;
  vector<solution> new_pop, pop;

  // set seed
  srand (time(NULL));
  t_start = clock();

  initializePop(pop, tam_pob, mat.size(), choosen);
  evaluatePop(pop, mat, iterations);

  while (evaluations < MAX_EVALUATIONS) {
    // NOTE AHORA MISMO EL OP. DE SELECCION ASUME QUE HAY ORDEN, PERO NO HAY ORDEN
    // ESTUDIAR SI MERECE LA PENA
    generations++;
    selection(pop, new_pop);
    cross(new_pop, cross_prob);
    mutatePop(new_pop, mut_prob, choosen, iterations);
    replace(new_pop, pop);
    evaluatePop(new_pop, iterations);
  }
  t_total = clock() - t_start;

  solution sol;
  getBestSolution(pop, sol);

  // output: Fitness - Time - Iterations
  cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << generations<< endl;
  return sol.fitness;
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

  localSearch(mat, choosen);
}
