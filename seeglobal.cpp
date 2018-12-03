// This is used for no noise. 
// All information and variables names are based on those in Doi 10.1038/ncomms7664
//Emergent long-range synchronization of oscillating
//ecological populations without external forcing
//described by Ising universality
//Independent study in fall 2018 semester under advisement of Prof Jon Machta 

#include <iostream>
//#include < omp.h > //http://www.bowdoin.edu/~ltoma/teaching/cs3225-GIS/fall17/Lectures/openmp.html
#include <time.h>
#include <fstream>
#include <sstream>
//#include <ppl.h>
#include <algorithm>                    
#include <cmath>                        
#include <cstdlib> 
#include <stdio.h>
//#include <armadillo>
#include <fstream>
#include <time.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <deque>
#include <time.h>
#include <random>
using namespace std;
//using namespace arma;


#define L 64  /* system size */
#define N L*L  /* number of spins */
#define D 4  /* number of neighbors */
int nn[N][D]; //neighbor table





double lattice[N]; //each spot holds a population density between 0 and 1. Not including 0 and 1. 
double last_lattice[N];
double neighbors[N]; //Is not directly related to the function neighor on line 82. 
double mjt[N]; //As named in cited paper in line 1. In the paper it is called the two-cycle amplitude
double mt = 0; //instantaneous synchronization. It is the sum of all mjt over the lattice for one discrete step in time
double order = 0; //Order is the time average of mt.





//double  pop = 0;
double sus_mt = 0;


int th_id, nthreads; //used in conjuction with omp.h... ideally

const double r = 3.3; //the r for the logistc map: rx(1-x)
double df = 0.1; //dispersal fraction for coupled logistica map. 0.9 stays and 0.1 goes to neighbors




//double order_squared;


//double ms=0;
//double mp=0;

//double mean_pop[range];
//double mean_square[range];
//double squares = 0;
//double logistic(int s) {
	//lattice[s] = r * last_lattice[s] * (1 - last_lattice[s]);


	//return 0;
//}



void neighbor(void) {   //two dimensional neighbor ... interesting note to remember. In c++ division of ints result in ints. Not floats! That's why this modular arith works. We need this because we have periodic boundary conditions.
	int disp;
	for (int k = 0; k < N; k++) {
		disp = L;
		nn[k][0] = (k / disp)*disp + (k + 1) % disp;  //+x
		nn[k][1] = (k / disp)*disp + (k - 1 + disp) % disp;  //-x
		disp = L * L;
		nn[k][2] = (k / disp)*disp + (k + L) % disp; //+y
		nn[k][3] = (k / disp)*disp + (k - L + disp) % disp; //-y
	}
}

//double remember{
	//for (int i = 0; i < (N - 1); i++) {
		//last_lattice[i] = lattice[i];

	//}
//}

long int const test = 1000; //Total number of interations for each local noise strengh
long int const burn = 0; //No data collected over the first burn# of steps
//int scope;
//long int to avoid overflow 

double mt_store[test - burn - 1];
double lattice_store[test - burn - 1][N];

int main() {
	clock_t start = clock();

	neighbor(); //create neighbor table




		//{
			//printf_s("%d\n", omp_get_max_threads());
		//}


		//{
			//printf_s("%d\n", omp_get_num_threads());
		//}
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < (N); i++) { //first lattice generation
		distribution.reset(); //Is this a bottleneck?
		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		lattice[i] = distribution(generator);
		//if (i>2048 & i< 2112) {
			//lattice[i] = 0.5;
		//}
	}
	



	for (long int year = 0; year < test; year++) {

		//if (year % burn == 0) {

			//std::cout << year << endl;
		//}











		for (int i = 0; i < (N); i++) {//remembering the last lattice
			last_lattice[i] = lattice[i];
		}





		for (int i = 0; i < (N); i++) {

			lattice[i] = r * last_lattice[i] * (1 - last_lattice[i]); //applying logistic map and local noise and global noise (if gn is included)
			if (lattice[i] >= 1)
			{
				lattice[i] = 0.999;
			}
			if (lattice[i] <= 0)
			{
				lattice[i] = 0.001;
				//	add lower bound check 
			}
		}

		//}

		for (int i = 0; i < (N); i++) {
			neighbors[i] = lattice[i];


		}

		for (int i = 0; i < (N); i++) {//adding neighbors 

			lattice[i] = (1 - df) * neighbors[i] + 0.25*df*(neighbors[nn[i][0]] + neighbors[nn[i][1]] + neighbors[nn[i][2]] + neighbors[nn[i][3]]);//no leaking of df (dispersal fraction). 1-.25*4=0
			if (lattice[i] >= 1)
			{
				lattice[i] = 0.999;
				//	add lower bound check 
			}
			if (lattice[i] <= 0)
			{
				lattice[i] = 0.001;
				//	add lower bound check 
			}

		}

		if (year > burn) {
			for (int i = 0; i < N; i++) {
				mjt[i] = .5 *  pow(-1, year) * (lattice[i] - last_lattice[i]);



			}

			for (int i = 0; i < N; i++) {
				mt = mt + mjt[i];

			}
			mt_store[year - burn-  1] = mt/ N ;

			mt = 0;
		}

		for (int i = 0; i < N; i++) {
			lattice_store[year - burn - 1][i] = lattice[i];
		}
	}
	

	
	



	ofstream myfile("mt_store.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int a = 0; a < test - burn - 1; a++) {

			myfile << mt_store[a] << " ";

		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("lattice_store.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int a = 0; a < test - burn - 1; a++) {
			for (int i = 0; i < N; i++) {
				myfile << lattice_store[a][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";





	clock_t end = clock();
	float time = (float)(end - start) / CLOCKS_PER_SEC;
	std::cout << "Time is ";
	std::cout << time << endl;
	//myfile.open("order.csv");
	//myfile << order_store;
	//myfile.close();


	//for (int i = 0; i < N; i++){
		//fout << lattice(i) << ",";

		//fout.close();
	myfile.open("notes.csv");//Recording testing conditions
	if (myfile.is_open())
	{


		//myfile << order_store[k] << " ";
		//myfile << ls << " ";
		//myfile << gs << " ";
		myfile << burn << " ";
		myfile << test << " ";
		myfile << time << " ";
		myfile << L << endl;
		myfile.close();
	}
	else cout << "Unable to open file";



	return 0;
}
