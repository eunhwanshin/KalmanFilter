// Eun Hwan Shin, 2022
// Reference: Mohinder S. Grewal and Angus P. Andrews. (2001).
//   Kalman Filgering: Theorey and Practice, 2nd Ed., John-Wiley & Sons, Inc.

#include <vector>
#include <cmath>  // sqrt()

#include <iostream>

#include "sqrt_kf.h"

namespace sqrt_kf
{
	//-------------------------------------------------------------------------
	//
	void ToMatrix(
		double A[], 
		size_t rows, size_t cols, 
		std::vector<double*>& mA)
	{
		if (mA.size() != rows)
		{
			mA.resize(rows);
		}

		for (size_t i = 0, ii = 0; i < rows; ++i, ii += cols)
			mA[i] = &A[ii];
	}

	//-------------------------------------------------------------------------
	//
	void HouseholderTriangular(double A[], int rows, int cols)
	{
		std::vector<double*> mA(rows);

		const int r = cols - rows;

		ToMatrix(A, rows, cols, mA);

		std::vector<double> v(cols);

		for (int k = rows - 1; k >= 0; --k)
		{
			double sigma = 0.0;

			for (int j = 0; j <= r + k; ++j)
				sigma += mA[k][j] * mA[k][j];

			double alpha = sqrt(sigma);

			sigma = 0.0;

			for (int j = 0; j <= r + k; ++j)
			{
				if (j == r + k)
					v[j] = mA[k][j] - alpha;
				else
					v[j] = mA[k][j];

				sigma += v[j] * v[j];
			}

			alpha = 2.0 / sigma;

			for (int i = 0; i <= k; ++i)
			{
				sigma = 0.0;
				for (int j = 0; j <= r + k; ++j)
					sigma += mA[i][j] * v[j];
				double beta = alpha * sigma;

				for (int j = 0; j <= r + k; ++j)
					mA[i][j] -= beta * v[j];
			}
		}
	}

	//-------------------------------------------------------------------------
	//
	void SchmidtHouseholderPredict(
		double A[], double B[],
		int n, int r)
	{
		std::vector<double*> mA(n);
		std::vector<double*> mB(n);

		ToMatrix(A, n, n, mA);
		ToMatrix(B, n, r, mB);

		std::vector<double> v(n);
		std::vector<double> w(r);

		for (int k = n - 1; k >= 0; --k)
		{
			double sigma = 0.0;

			for (int j = 0; j < r; ++j)
				sigma += mB[k][j] * mB[k][j];

			for (int j = 0; j <= k; ++j)
				sigma += mA[k][j] * mA[k][j];

			double alpha = sqrt(sigma);
			sigma = 0.0;

			for (int j = 0; j < r; ++j)
			{
				w[j] = mB[k][j];
				sigma += w[j] * w[j];
			}

			for (int j = 0; j <= k; ++j)
			{
				if (j == k)
					v[j] = mA[k][j] - alpha;
				else
					v[j] = mA[k][j];

				sigma += v[j] * v[j];
			}

			alpha = 2.0 / sigma;

			for (int i = 0; i <= k; ++i)
			{
				sigma = 0.0;
				for (int j = 0; j < r; ++j)
					sigma += mB[i][j] * w[j];
				for (int j = 0; j <= k; ++j)
					sigma += mA[i][j] * v[j];

				double beta = alpha * sigma;

				for (int j = 0; j < r; ++j)
					mB[i][j] -= beta * w[j];
				for (int j = 0; j <= k; ++j)
					mA[i][j] -= beta * v[j];
			}
		}

	}

	//-------------------------------------------------------------------------
	//
	void CarlsonMeasUpdate(
		int n,
		double x[], double C[],
		double z, double H[], double R)
	{
		std::vector<double*> mC(n);
		ToMatrix(C, n, n, mC);

		std::vector<double> w(n);
		double alpha = R;
		double inno = z;

		for (int j = 0; j < n; ++j)
		{
			inno -= H[j] * x[j];

			double sigma = 0.0;

			for (int i = 0; i <= j; ++i)
				sigma += mC[i][j] * H[i];

			double beta = alpha;

			alpha += sigma * sigma;

			double gamma = sqrt(alpha * beta);
			double eta = beta / gamma;
			double zeta = sigma / gamma;

			w[j] = 0.0;

			for (int i = 0; i <= j; ++i)
			{
				double tau = mC[i][j];
				mC[i][j] = eta * mC[i][j] - zeta * w[i];
				w[i] += tau * sigma;
			}

		}

		double epsilon = inno / alpha;

		for (int i = 0; i < n; ++i)
			x[i] += w[i] * epsilon;
	}
}