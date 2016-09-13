// start.cpp
//
// This program performs a series of mathematical operations on a set of integers read in from a file.
// 
//*****************
#include <iostream> // this is a standard C++ header, and should be included in all programs
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include "mprof.hh" // simple profile functions
using namespace std;
//*****************

// declaration of functions/subroutines

int polynomial(int i);	// a polynomial to be evaluated

template <typename T>
T part_i_slowest(const T x) {
  return 2 * pow(x, 6) - 3 * pow(x, 4) + 4 * pow(x, 2) - 3;
}

template <typename T>
T part_i_fast(const T x) {
  return 2 * x*x*x*x*x*x - 3 * x*x*x*x + 4 * x*x - 3;
}

template <typename T>
T part_i_fastest(const T x) {
  T x2 = x*x;
  T x4 = x2*x2;
  return 2 * x4*x2 - 3 * x4 + 4 * x2 - 3;
}

bool is_prime(int x);

// declaration of input and output stream

ifstream data_in("5input.txt"); // input data
ofstream data_out("5output.txt"); // output data

// the main part of the program - should always be type "int"

int main() {

  #define outs data_out // decide on output stream here

	// variable declarations
    int i,j,k; // loop counters
    const int upper = 8; // upper limit on calculation loop, corresponds to number of entries in the input file
    					  // declared as a "const int" so that arrays can be defined with this size
    int numbers[upper]; // array to store numbers, in C++, the array index starts at 0!
    
    bool sq; // boolean that indicates if the number is a perfect square
    int pol; // value of polynomial function
    int quotient, divisor, remainder;
    char square[200]; // strings for output
    
    // read in data from file
    for(i=0;i<upper;i++) data_in>>numbers[i];
    
    // calculation loop
    
    for(i=0;i<upper;i++){ // the loop will run over the elements in numbers[]
	   
	    // first, determine if the number is a perfect square
	   
	    sprintf(square," is not a perfect square"); //default
	    if(numbers[i]<0) sprintf(square," is negative so that the concept of a perfect square is undefined for real numbers");
		else {
			quotient = numbers[i];
			divisor = 2;
			sq = false;
			while(sq == false && quotient > divisor) {
				quotient  = numbers[i]/divisor;
				remainder = numbers[i]%divisor;
				if (quotient == divisor && remainder == 0) {
					sq = true;
					sprintf(square," is a perfect square"); 
				}	
				else divisor = divisor + 1;
			}
		}
	    
		
	    // second, evaluate the polynomial function
	    
	    pol = polynomial(numbers[i]);

      int resi;
      double resd;
      chrono::duration<double> elapsed_time;

      outs << "==========================================================================================================\n";
      outs << "calculating 2x^6 - 3x^4 + 4x^2 - 3\n";
      outs << "\nresult, time\n\n";
      outs << setw(35) << "`int part_i_slowest(int)`" 
           << setw(35) << "`int part_i_fast(int)`"   
           << setw(35) << "`int part_i_fastest(int)`" << '\n';
      outs << "----------------------------------------------------------------------------------------------------------\n";

	    tie(resi, elapsed_time) = profile(part_i_slowest<int>, numbers[i]);
      outs << setw(17) << resi << ", " << setw(15) << elapsed_time.count();

	    tie(resi, elapsed_time) = profile(part_i_fast<int>, numbers[i]);
      outs << setw(17) << resi << ", " << setw(15) << elapsed_time.count();

	    tie(resi, elapsed_time) = profile(part_i_fastest<int>, numbers[i]);
      outs << setw(17) << resi << ", " << setw(15) << elapsed_time.count();
      outs << "\n\n\n\n";

      outs << setw(35) << "`double part_i_slowest(double)`" 
           << setw(35) << "`double part_i_fast(double)`"   
           << setw(35) << "`double part_i_fastest(double)`" << '\n';
      outs << "-----------------------------------------------------------------------------------------------------------\n";

	    tie(resd, elapsed_time) = profile(part_i_slowest<double>, static_cast<double>(numbers[i]));
      outs << setw(17) << resd << ", " << setw(15) << elapsed_time.count();

	    tie(resd, elapsed_time) = profile(part_i_fast<double>, static_cast<double>(numbers[i]));
      outs << setw(17) << resd << ", " << setw(15) << elapsed_time.count();

	    tie(resd, elapsed_time) = profile(part_i_fastest<double>, static_cast<double>(numbers[i]));
      outs << setw(17) << resd << ", " << setw(15) << elapsed_time.count();
      outs << "\n\n";

      if (is_prime(numbers[i]))
        outs << numbers[i] << " is prime.\n";
      else
        outs << numbers[i] << " is not prime.\n";

	    //output the data
	    
	    outs<<"The number "<<numbers[i]<<square<<" and returns a value of "<<pol<<" when inserted into the function f"<<endl; 
     
    } 
}	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int polynomial(int i) {
    int f;
    
    int x2=i*i;
    f = 3*x2*x2+4*x2-3;
    
    return f;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool is_prime(int x) {
  if (x < 1) return false;
  if (x == 1 || x == 2) return true;
  if (x % 2 == 0) return false;
  
  for(int d = 3; d < x; d+=2)
    if (x % d == 0) return false;

  return true;
}
