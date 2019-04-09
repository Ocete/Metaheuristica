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

///////////////////////////// GENETIC ///////////////////////////////////////

void initializePop(vector<solution> &pop, int tam_pob, int choosen) {
  int size = MAT.size();
  pop.resize(tam_pob);
  for (unsigned i=0; i<tam_pob; i++) {
    randomSolution(pop[i], choosen);
  }
}

void evaluatePop(vector<solution> &pop, int &iterations) {
  for (unsigned i=0; i<pop.size(); i++) {
    if ( !pop[i].evaluated ) {
      evaluateSolution( pop[i] );
      iterations++;
    }
  }
}

// Assumes the population is evaluated
int binaryTournament(vector<solution> &pop) {
  int r1 = random(0, pop.size());
  int r2 = random(0, pop.size());
  r1 = pop[r1].fitness > pop[r2].fitness ? r1 : r2;
  return  r1;
}

void selection(vector<solution> &pop, vector<solution> &new_pop) {
  unsigned pop_tam = pop.size();
  int select;
  new_pop.resize(pop_tam);
  for (unsigned i=0; i<n_crosses; i++) {
    select = binaryTournament ( pop );
    new_pop[i] = pop[select];
  }
}

// Mutate the solution by swapping two random elements
void mutateSolution(solution &sol, int &iterations) {
  int r_on, r_off;

  // Find elements to swap
  do {
    r_on = random(0, choosen);
  } while ( !sol.v[r_on] );

  do {
    r_out = random(0, choosen);
  } while ( sol.v[r_off] );

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

void mutatePop(vector<solution> &pop, double &mut_prob,
      int choosen, int &iterations) {
  int r_sol, n_mut = mut_prob*pop.size();
  for (unsigned i=0; i<n_mut; i++) {
    r_sol = random(0,pop.size());
    mutateSolution( pop[r_sol], iterations );
  }
}

void cross(vector <solution> &pop, double cross_prob) {

}

void replace(vector<solution> &new_pop, vector<solution> &pop) {
  pop.swap(new_pop);

  // TODO falta por hacer el elitismo. Para eso lo mas facil seria hacer una
  // struct de "pop" que mantenga guardado el fitness maximo de sus elementos
  // y quizas tambien la posicion del mejor elemento
}

// Computes the local search algorithm for a random starting solution
// This implementation doesn't assume the pop is ordered
double genetic( int choosen, int MAX_EVALUATIONS ) {
  int evaluations = 0, MAX_EVALUATIONS = 50000, tam_pob = 50, generations = 0;
  double mut_prob = 0.7, cross_prob = 0.001;
  clock_t t_start, t_total;
  vector<solution> new_pop, pop;

  // set seed
  srand (time(NULL));
  t_start = clock();

  initializePop(pop, tam_pob, choosen);
  evaluatePop(pop, iterations);

  while (evaluations < MAX_EVALUATIONS) {
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
  int MAX_EVALUATIONS = 50000;

  // set seed
  srand (time(NULL));

  cin >> size >> choosen;
  readInput(size);

  genetic(choosen, MAX_EVALUATIONS);
}
