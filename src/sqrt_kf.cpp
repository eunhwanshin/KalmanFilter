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
	template <class T>
	void ToMatrix(
		T A[], 
		size_t rows, size_t cols, 
		std::vector<T*>& mA)
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
	bool CholUT(int n, const double P[], double S[])
	{
		std::vector<const double*> mP(n);
		ToMatrix(P, n, n, mP);

		std::vector<double*> mS(n);
		ToMatrix(S, n, n, mS);

		for (int j = n - 1; j >= 0; --j)
		{
			for (int i = j; i >= 0; --i)
			{
				double sigma = mP[i][j];

				for (int k = j + 1; k < n; ++k)
				{
					sigma -= mS[i][k] * mS[j][k];
				}

				if (i == j)
				{
					if (sigma <= 0.0) return false;

					mS[i][j] = sqrt(sigma);
				}
				else
					mS[i][j] = sigma / mS[j][j];
			}
			
		}

		return true;
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

	//-------------------------------------------------------------------------
    //
	void UD_Factor(int n, const double P[], double U[], double D[])
	{
		std::vector<const double*> mP(n);
		ToMatrix(P, n, n, mP);

		std::vector<double*> mU(n);
		ToMatrix(U, n, n, mU);

		for (int j = n - 1; j >= 0; --j)
		{
			for (int i = j; i >= 0; --i)
			{
				double sigma = mP[i][j];

				for (int k = j + 1; k < n; ++k)
					sigma -= mU[i][k] * D[k] * mU[j][k];

				if (i == j)
				{
					D[j] = sigma;
					mU[j][j] = 1.0;
				}
				else
				{
					mU[i][j] = sigma / D[j];
				}
			}
		}
	}

	//-------------------------------------------------------------------------
	//
	void UD_Predict(
	    int n, int p,
		double PhiU[], const double D[],
		double Um[], double Dm[],
		double GUq[], double Dq[])
	{
		std::vector<double*> mPhiU(n);
		ToMatrix(PhiU, n, n, mPhiU);

		std::vector<double*> mU(n);
		ToMatrix(Um, n, n, mU);

		std::vector<double*> mGUq(n);
		ToMatrix(GUq, n, p, mGUq);

		for (int i = n - 1; i >= 0; --i)
		{
			double sigma = 0.0;

			for (int j = 0; j < n; ++j)
				sigma += mPhiU[i][j] * mPhiU[i][j] * D[j];

			for (int j = 0; j < p; ++j)
			{
				sigma += mGUq[i][j] * mGUq[i][j] * Dq[j];
			}
				

			Dm[i] = sigma;
			double inv_D = 1.0 / Dm[i];
			mU[i][i] = 1.0;

			for (int j = 0; j < i; ++j)
			{
				sigma = 0;

				for (int k = 0; k < n; ++k)
					sigma += mPhiU[i][k] * D[k] * mPhiU[j][k];

				for (int k = 0; k < p; ++k)
					sigma += mGUq[i][k] * Dq[k] * mGUq[j][k];

				mU[j][i] = sigma * inv_D;

				for (int k = 0; k < n; ++k)
					mPhiU[j][k] -= mU[j][i] * mPhiU[i][k];

				for (int k = 0; k < p; ++k)
					mGUq[j][k] -= mU[j][i] * mGUq[i][k];
			}
		}
	}

	//-------------------------------------------------------------------------
	//
	void BiermanUpdate(
		int n, double x[],
		double U[], double D[],
		double z, double H[], double R)
	{
		std::vector<double*> mU(n);
		ToMatrix(U, n, n, mU);

		std::vector<double> v(n, 0);
		std::vector<double> w(n, 0);

		double delta = z;

		for (int j = 0; j < n; ++j)
		{
			delta -= H[j] * x[j];
			v[j] = H[j];

			for (int i = 0; i < j; ++i)
				v[j] += mU[i][j] * H[i];
		}

		double sigma = R;

		for (int j = 0; j < n; ++j)
		{
			double nu = v[j];
			v[j] *= D[j];
			w[j] = v[j];
			for (int i = 0; i < j; ++i)
			{
				double tau = mU[i][j] * v[j];
				mU[i][j] -= nu * w[i] / sigma;
				w[i] += tau;
			}

			D[j] *= sigma;
			sigma += nu * v[j];
			D[j] /= sigma;
		}

		delta /= sigma;
		for (int i = 0; i < n; ++i)
			x[i] += w[i] * delta;
			
	}

	//-------------------------------------------------------------------------
	//
	void UD_Update(
		int n, double x[], double UD[], 
		std::function<int(int, int)> f_idx,
		double z, double a[], double R)
	{
		// Compute innovation: z= z - a*x
		for (int j = 0; j < n; ++j)
		{
			z -= a[j] * x[j];
		}

		// Unweighted Kalman gain: K = b/alpha
		std::vector<double> b(n);

		// b = D*U'*a; a = U'*a
		for (int j = n - 1; j > 0; --j)
		{
			for (int k = 0; k < j; ++k)
			{
				a[j] += UD[f_idx(k, j)] * a[k];
			}
			b[j] = UD[f_idx(j, j)] * a[j];
		}

		b[0] = UD[f_idx(0, 0)] * a[0];

		// variance of the innovation
		double alpha = R + b[0] * a[0];
		double gamma = 1. / alpha;

		UD[f_idx(0, 0)] *= R * gamma;

		for (int j = 1; j < n; ++j)
		{
			double beta = alpha;
			alpha += b[j] * a[j];
			double lambda = -a[j] * gamma;

			gamma = 1. / alpha;

			UD[f_idx(j, j)] *= beta * gamma;

			for (int i = 0; i < j; ++i)
			{
				double& uij = UD[f_idx(i, j)];
				beta = uij;

				uij = beta + b[i] * lambda;
				b[i] += b[j] * beta;

			}
		}

		z *= gamma;

		for (int j = 0; j < n; ++j)
		{
			x[j] += b[j] * z;
		}
	}
}