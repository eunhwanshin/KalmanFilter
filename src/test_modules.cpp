#include <iostream>
#include <iomanip>

#include "sqrt_kf.h"

void PrintMatrix(const double A[], int rows, int cols, const char* name)
{
	std::cout << name << std::endl;

	for (int i = 0, i1 = 0; i < rows; ++i, i1 += cols)
	{
		for (int j = 0; j < cols; ++j)
			std::cout << " " << std::fixed << std::setprecision(6) << A[i1 + j];

		std::cout << std::endl;
	}
}

int main()
{
	double A[3 * 5] = {
	 0.927864,   0.519861,   0.916024,   0.981743,   0.939206,
     0.965873,   0.017677,   0.689787,   0.560441,   0.794960,
     0.545124,   0.677176,   0.933571,   0.256789,   0.798057
	};

	sqrt_kf::HouseholderTriangular(A, 3, 5);

	PrintMatrix(A, 3, 5, "A");

	double B[3 * 3] = {
	 0.927864,   0.519861,   0.916024,   
	 0.965873,   0.017677,   0.689787,  
	 0.545124,   0.677176,   0.933571,  
	};

	double C[3 * 2] = {
	  0.981743,   0.939206,
	  0.560441,   0.794960,
	  0.256789,   0.798057
	};

	sqrt_kf::SchmidtHouseholderPredict(B, C, 3, 2);

	PrintMatrix(B, 3, 3, "C(-)");
	PrintMatrix(C, 3, 2, "Zeroed");

	double x[3] = { 0.020503,  0.494263 , 0.819909 };
	double H[3] = { 0.36211,   0.23911,   0.16862 };
	double R = 2;
	double z = 1;

	sqrt_kf::CarlsonMeasUpdate(3, x, B, z, H, R);

	PrintMatrix(x, 1, 3, "x(+)");
	PrintMatrix(B, 3, 3, "C(+)");


	return 0;
}