#include <iostream>
#include <iomanip>

#include "sqrt_kf.h"

int UD_Index(int row, int col)
{
	return row + col * (col + 1) / 2;
}

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

void PrintUD(const double UD[], int rows, int cols, const char* name)
{
	std::cout << name << std::endl;

	for (int i = 0, i1 = 0; i < rows; ++i, i1 += cols)
	{
		for (int j = 0; j < i; ++j)
			std::cout << " " << std::fixed << std::setprecision(6) << 0.0;
		for (int j = i; j < cols; ++j)
			std::cout << " " << std::fixed << std::setprecision(6) << UD[UD_Index(i,j)];

		std::cout << std::endl;
	}
}

void TestUD(void)
{
	const int n = 5;
	const int p = 2;

	double P[n * n] = {
		0.88428,   0.51029,   0.98356,   1.10230,   0.69445,
		0.51029,   0.59402,   0.70002,   0.63662,   0.81465,
        0.98356,   0.70002,   1.98858,   1.93902,   1.61581,
        1.10230,   0.63662,   1.93902,   2.29528,   1.89578,
        0.69445,   0.81465,   1.61581,   1.89578,   2.34828
	};

	double U[n * n] = {0};
	double D[n];

	sqrt_kf::UD_Factor(n, P, U, D);

	PrintMatrix(U, n, n, "U");
	PrintMatrix(D, 1, n, "D");

	double PhiU[n * n];

	memcpy(PhiU, U, sizeof(U));

	PrintMatrix(PhiU, n, n, "PhiU");

	double GUq[n * p] = {
		0.2253,   0.6838,
        0.6340,   1.0639,
        0.9176,   0.6798,
        0.7892,   0.8583,
        0.6533,   0.3939	
	};

	double Dq[p] = { 0.2134 , 0.7688 };

	double Dm[5] = { 0 };

	sqrt_kf::UD_Predict(n, p, PhiU, D, U, Dm, GUq, Dq);

	PrintMatrix(U, n, n, "U(-)");
	PrintMatrix(Dm, 1, n, "D(-)");

	memcpy(D, Dm, sizeof(D));

	// Copy U and D into vectorized UD
	double UD[n * (n + 1) / 2];

	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < j; ++i)
		{
			UD[UD_Index(i, j)] = U[i * n + j];
		}

		UD[UD_Index(j, j)] = Dm[j];
	}

	PrintUD(UD, n, n, "UD(-)");

	double x[n] = { 0 };
	double z = 1.0;
	double H[5] = { 0.878975,   0.207268,   0.098002,   0.943587,   0.244249 };
	double R = 1.0;

	sqrt_kf::BiermanUpdate(n, x, U, D, z, H, R);

	PrintMatrix(x, 1, n, "x(+)");
	PrintMatrix(U, n, n, "U(+)");
	PrintMatrix(D, 1, n, "D(+)");

	memset(x, 0, sizeof(x));
	double a[5] = { 0.878975,   0.207268,   0.098002,   0.943587,   0.244249 };

	sqrt_kf::UD_Update(n, x, UD, UD_Index, z, a, R);

	PrintMatrix(x, 1, n, "x(+)");
	PrintUD(UD, n, n, "UD(+)");

}

void TestChol(void)
{
	const int n = 5;
	double P[n * n] = {
		0.88428,   0.51029,   0.98356,   1.10230,   0.69445,
		0.51029,   0.59402,   0.70002,   0.63662,   0.81465,
		0.98356,   0.70002,   1.98858,   1.93902,   1.61581,
		1.10230,   0.63662,   1.93902,   2.29528,   1.89578,
		0.69445,   0.81465,   1.61581,   1.89578,   2.34828
	};

	double C[n * n] = { 0 };

	if (sqrt_kf::CholUT(n, P, C))
		PrintMatrix(C, n, n, "Cholesky UT");
	else
		std::cout << "Not a positive difinte matrix." << std::endl;
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

	TestChol();

	TestUD();

	return 0;
}