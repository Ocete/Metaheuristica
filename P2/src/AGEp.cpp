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

///////////////////////////// GENETIC ///////////////////////////////////////

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

// Comparison operator for ordering the population
// Keeps at the front the element with the highest fitness
bool operator < (const solution &s1, const solution &s2) {
    return s1.fitness > s2.fitness;
}

void initializePop(population &pop, int tam_pob, int choosen) {
  pop.v.resize(tam_pob);
  for (int i=0; i<tam_pob; i++) {
    randomSolution(pop.v[i], choosen);
  }
  pop.tam = pop.v.size();
}

void evaluatePop(population &pop, int &evaluations) {
  for (int i=0; i<pop.tam; i++) {
    if ( !pop.v[i].evaluated ) {
      evaluateSolution( pop.v[i] );
      if ( pop.v[i].fitness > pop.max_fitness ) {
        pop.best_sol = i;
        pop.max_fitness = pop.v[i].fitness;
      }
      evaluations++;
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

void selectPair(population &pop, solution &s1, solution &s2) {
  int i1, i2;
  i1 = binaryTournament ( pop );
  do {
    i2 = binaryTournament ( pop );
  } while (i1 == i2);
  s1 = pop.v[i1];
  s2 = pop.v[i2];
}

// Mutate the solution by swapping two random elements
void mutateSolution(solution &sol, int &evaluations) {
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
    evaluations++;
  }
}

void mutatePair(solution &c1, solution &c2, double mut_prob) {
  int aux = 0;
  if ( rand() < mut_prob*2 ) {
    mutateSolution( c1, aux );
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

void positionalCross(solution &p1, solution &p2,
                          solution &c1, solution &c2) {
  c1 = p1;
  c2 = p2;
  c1.evaluated = false;
  c2.evaluated = false;

  vector<bool> shuffled (p1.v.size(), false);
  vector<bool> to_shuffle, to_shuffle2;

  for (unsigned i=0; i<p1.v.size(); i++) {
    if ( p1.v[i] == p2.v[i] ) {
      c1.v[i] = p1.v[i];
      c2.v[i] = p1.v[i];
    } else {
      shuffled[i] = true;
      to_shuffle.push_back( p1.v[i] );
    }
  }

  random_shuffle( to_shuffle.begin(), to_shuffle.end() );
  to_shuffle2 = to_shuffle;
  random_shuffle( to_shuffle2.begin(), to_shuffle2.end() );

  int j=0;
  for (unsigned i=0; i<p1.v.size(); i++) {
    if ( shuffled[i] ) {
      c1.v[i] = to_shuffle[j];
      c2.v[i] = to_shuffle2[j];
      j++;
    }
  }
}

void crossPair(solution &s1, solution &s2, solution &c1, solution &c2) {
  positionalCross( s1, s2, c1, c2 );
}

void replace(population &pop, solution &s1, solution &s2) {
  pop.v.push_back(s1);
  pop.v.push_back(s2);
  sort( pop.v.begin(), pop.v.end() );
  pop.v.resize( pop.v.size() - 2);
}

// Computes the local search algorithm for a random starting solution
// This implementation doesn't assume the pop is ordered
double AGEp( int choosen, int MAX_EVALUATIONS ) {
  int evaluations = 0, tam_pob = 50, generations = 0;
  double mut_prob = 0.001;
  clock_t t_start, t_total;
  population new_pop, pop;
  solution p1, p2, c1, c2;

  // set seed
  srand (time(NULL));
  t_start = clock();

  initializePop(pop, tam_pob, choosen);
  evaluatePop(pop, evaluations);
  sort(pop.v.begin(), pop.v.end());
  pop.best_sol = 0;

  while (evaluations < MAX_EVALUATIONS) {
    generations++;
    selectPair(pop, p1, p2);
    crossPair(p1, p2, c1, c2);
    mutatePair(c1, c2, mut_prob*tam_pob);
    evaluateSolution(c1);
    evaluateSolution(c2);
    evaluations += 2;
    replace(pop, c1, c2);
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

  AGEp(choosen, MAX_EVALUATIONS);
}
