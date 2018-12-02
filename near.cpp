
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


#define L 10  /* system size */
#define N L*L  /* number of spins */
#define D 4  /* number of neighbors */
int nn[N][D]; //neighbor table
const int start_range = 0;
const int end_range = 100; //sets interval size of noise strength. It marks the last noise
const int range = end_range - start_range; //Counts how man different local noise levels are being checked


const int global_start = 10;
const int global_end = 30;
const int global_range = (global_end + 1) - global_start;

double noise_strength[range];
double global_strength[21] = { 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03 };

double lattice[N]; //each spot holds a population density between 0 and 1. Not including 0 and 1. 
double last_lattice[N];
double neighbors[N]; //Is not directly related to the function neighor on line 82. 
double mjt[N]; //As named in cited paper in line 1. In the paper it is called the two-cycle amplitude
double mt = 0; //instantaneous synchronization. It is the sum of all mjt over the lattice for one discrete step in time
double order = 0; //Order is the time average of mt.
double order_store[global_range][range]; //There is a local noise level for each order


double mt_squared = 0; //This is used to calculate the suspectilbity, ie the variance the variable order, aka the order parameter
double total_mt_squared = 0; //This is used to calculate the suspectilbity
double sus_mt_total = 0; //This is used to calculate the suspectilbity
double sus1[global_range][range];// Using the short cut for calculating variance <x^2> - <x>^2
double sus2[global_range][range];//

//double  pop = 0;
double sus_mt = 0;


int th_id, nthreads; //used in conjuction with omp.h... ideally

const double r = 3.3; //the r for the logistc map: rx(1-x)
double df = 0.1; //dispersal fraction for coupled logistica map. 0.9 stays and 0.1 goes to neighbors
double ln[N]; //local noise multiplied unto each spot on the lattice
double ls = 0;//local noise strength
double gs = 0;// global noise strength
double gn = 0;//
int k = 0; //Indexes order


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

long int test = 100; //Total number of interations for each local noise strengh
long int burn = 10; //No data collected over the first burn# of steps
//int scope;
//long int to avoid overflow 



int main() {
	clock_t start = clock();

	neighbor(); //create neighbor table

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0, 1.0);


	//#pragma omp parallel omp_set_num_threads(3);    
	//#pragma omp master  













	for (int i = start_range; i < end_range; i++) {

		noise_strength[i - start_range] = 0.001*i;// = 0.1 * noise_strength; "std of noise"... Setting up the array of local noise strength. This basically sets the mean. 
	}


	//for (int i = global_start; i < (global_end + 2); i++) {

		//global_strength[i - global_start] = global_start + i;// = 0.1 * noise_strength; "std of noise"... Setting up the array of local noise strength. This basically sets the mean. 
	//}





	//double measure_all[89999];//test-burn-1

	//gn = (1 + gs * randu()); //global noise




	//ls = 0.00; //local strength 
	//gs = 0.00; //global strength 
	//#pragma omp parallel //shared(a, b, c) private(i)

	//#pragma omp  parallel  for shared() private()


	for (int a = 0; a < 21; a++) {
		neighbor(); //create neighbor table
		for (int k = 0; k < range; k++)
		{
			//{
				//printf_s("%d\n", omp_get_max_threads());
			//}


			//{
				//printf_s("%d\n", omp_get_num_threads());
			//}
			for (int i = 0; i < (N); i++) { //first lattice generation

				distribution.reset(); //Is this a bottleneck?
				std::uniform_real_distribution<double> distribution(0.0, 1.0);
				lattice[i] = distribution(generator);
			}




			for (long int year = 0; year < test; year++) {

				//if (year % burn == 0) {

					//std::cout << year << endl;
				//}





				if (year > burn) {

					for (int i = 0; i < N; i++) {
						mjt[i] = .5 *  pow(-1, year) * (lattice[i] - last_lattice[i]);


						//pop = pop + lattice[i];

						//squares = squares + lattice[i] * lattice[i];
					}
					sus_mt = 0;
					mt = 0;
					mt_squared = 0;
					for (int i = 0; i < N; i++) {
						mt = mt + mjt[i];
						//mt_squared = mt_squared + mjt[i] * mjt[i];
					}

					sus_mt = abs(mt); //for variance, we do not need the absolute value of mt

					mt = abs(mt); //There is a typo in the paper, we need to do that abs here and not before. 


					total_mt_squared = total_mt_squared + sus_mt * sus_mt;//<x^2> - <x>^2 ...  Here we are collecting x^2
					sus_mt_total = sus_mt_total + sus_mt;//<x^2> - <x>^2 ...  Here we are collecting x

					order = order + mt;
					//order_squared = order_squared + order * order;

					mt = 0;
					//mp= mp + pop / N;

					//ms = ms + squares / N;

					//pop = 0;
					//squares = 0;

					//measure_all[year - burn - 1] = mt;

				}




				for (int i = 0; i < (N); i++) {//remembering the last lattice
					last_lattice[i] = lattice[i];
				}





				for (int i = 0; i < (N); i++) {//making multiplicattive random local noise
					distribution.reset();

					ln[i] = (1 + noise_strength[k] * distribution(generator));

					//ln[i] = (1 + noise_strength[k] * randn());
				}

				distribution.reset();

				gn = (1 + global_strength[a] * distribution(generator)); //global noise

	//			omp_set_num_threads(4);
	//#pragma omp	parallel shared(lattice, last_lattice, ln, r) private(i) 
				//{
					//{
						//printf_s("%d\n", omp_get_max_threads());
					//}
					//{
						//printf_s("%d\n", omp_get_num_threads());
					//}
	//#pragma omp for 

				for (int i = 0; i < (N); i++) {


					//{
						//printf_s("%d\n", omp_get_max_threads());
					//}
					//{
						//printf_s("%d\n", omp_get_num_threads());
					//}
					lattice[i] = r * last_lattice[i] * (1 - last_lattice[i]) * ln[i] * gn; //applying logistic map and local noise and global noise (if gn is included)
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



				//for (int m = 0; m < (N - 1);m++); {
					//lattice[i] = randu(1);
				//}
				//last_lattice = lattice.copy;
						//for (int i = 0; i < (N - 1); i++) {
						//	lattice[i] = r * last_lattice[i](1 - last_lattice[i]);

						//}
				//int n = tbb::task_scheduler_init::default_num_threads();
				//std::cout << lattice;

				//total_mt_squared = total_mt_squared + mt_squared;


			}
			//order(ls_step) = 7.228 * 1.1892 * accu(measure_all) / ((exper - burn) * 16);
			//for (int i = 0; i < (test - burn - 1); i++) {
				//measure_all[i];


			//}

			order_store[a][k] = 7.228 * pow(L, 0.125) * order / ((test - burn) * N); //The  7.228 * pow(L, 0.125) are there because they are finite scaling parameters, these were listed in the paper. 
			//The order is not order until we do the divsion we do here. Before the variable 'order' is just a preamble for what the order will be.
			//sus1[k] = order_squared / ((test - burn)* N);
			//sus2[k] = order / ((test - burn)* N); 
			//sus2[k] = order * order * N;

			order = 0; //Set it to zero for the next local noise str
			//order_squared = 0;

			sus1[a][k] = total_mt_squared / ((test - burn) * N * N);//N^2 happens here because mt should be divided by N to be truly mt. Now N^2 since it is x^2 for <x^2> - <x>^2
			//test - burn because it is really a time average over the time data is collected
			//sus1[k] = total_mt_squared / ((test - burn));
			total_mt_squared = 0;

			sus2[a][k] = sus_mt_total / ((test - burn) * N); //

			sus_mt_total = 0;

			//mean_pop[k] = mp;// / (test - burn);
			//mean_square[k] = ms;// / (test - burn);

			//mp = 0;
			//ms = 0;

			//std::cout << "Local noise is ";
			//std::cout << noise_strength[k] << endl;


			//std::cout << "Order is ";
			//std::cout << order_store[a][k] << endl;
			//343 through 348 tells the user the noise,order pair. 






			//ofstream myfile("lattice.csv");
			//if (myfile.is_open())
			//{

				//for (int count = 0; count < N; count++) {
					//myfile << lattice[count] << " ";
				//}
				//myfile.close();
			//}
			//else cout << "Unable to open file";






			//k = k + 1;
			//ls = ls + 0.001;

		}
		std::cout << "Global noise is ";
		std::cout << global_strength[a] << endl;
	}

	//ofstream myfile("mean_pop.csv");
	//if (myfile.is_open())
	//{

		//for (int i = 0; i < range; i++) {
			//myfile << mean_pop[i] << " ";
		//}
		//myfile.close();
	//}
	//else cout << "Unable to open file";

	//myfile.open("mean_square.csv");
	//if (myfile.is_open())
	//{

		//for (int i = 0; i < range; i++) {
			//myfile << mean_square[i] << " ";
		//}
		//myfile.close();
	//}
	//else cout << "Unable to open file";

	ofstream myfile("order.csv");//All the files are sent to matlab for visualization. In the matlab file I do variance = sus1 - sus2.^2 
	if (myfile.is_open())
	{
		for (int a = 0; a < 21; a++) {
			for (int i = 0; i < range; i++) {
				myfile << order_store[a][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";


	myfile.open("global_range.csv");
	if (myfile.is_open())
	{

		for (int i = 0; i < 21; i++) {
			myfile << global_strength[i] << " ";
		}
		myfile.close();
	}
	else cout << "Unable to open file";




	myfile.open("noise_range.csv");
	if (myfile.is_open())
	{

		for (int i = 0; i < range; i++) {
			myfile << noise_strength[i] << " ";
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("sus1.csv");
	if (myfile.is_open())
	{
		for (int a = 0; a < 21; a++) {
			for (int i = 0; i < range; i++) {
				myfile << sus1[a][i] << " ";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	myfile.open("sus2.csv");
	if (myfile.is_open())
	{
		for (int a = 0; a < 21; a++) {
			for (int i = 0; i < range; i++) {
				myfile << sus2[a][i] << " ";
			}
		}
		myfile << endl << endl;
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
		myfile << global_range - 1 << " ";
		myfile << L << endl;
		myfile.close();
	}
	else cout << "Unable to open file";



	return 0;
}