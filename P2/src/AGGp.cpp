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

// Comparison operator for ordering the solution vector
// Keeps at the front the element with less contribution to the fitness
bool operator < (const pair<int,double> &p1, const pair<int,double> &p2) {
    return p1.second < p2.second;
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
  int r_sol, n_mut = mut_prob*pop.tam;
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

// Operador de cruce uniforme
void cross(population &pop, double cross_prob) {
  int n_crosses = cross_prob * pop.tam / 2 ;
  int p1, p2;
  solution ch1, ch2;

  for (int i=0; i<n_crosses; i++) {
    p1 = random(0, pop.tam);
    do {
      p2 = random(0, pop.tam);
    } while (p1 == p2);
    positionalCross( pop.v[p1], pop.v[p2], ch1, ch2);
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

// Computes the local search algorithm for a random starting solution
// This implementation doesn't assume the pop is ordered
double AGGu( int choosen, int MAX_EVALUATIONS ) {
  int evaluations = 0, tam_pob = 50, generations = 0;
  double mut_prob = 0.001, cross_prob = 0.7;
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
  }
  t_total = clock() - t_start;

  solution sol = pop.v[ pop.best_sol ];

  // output: Fitness - Time - Iterations
  cout << sol.fitness << "\t" << (double) t_total / CLOCKS_PER_SEC << "\t" << generations<< endl;
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

  AGGu(choosen, MAX_EVALUATIONS);
}
