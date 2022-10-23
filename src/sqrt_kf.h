#pragma once

#include <functional>
//
// Eun Hwan Shin, 2022
// 
// References: 
// [1] Mohinder S. Grewal and Angus P. Andrews. (2001).
//     Kalman Filgering: Theorey and Practice, 2nd Ed., John-Wiley & Sons, Inc.
// [2] Gerald J. Bierman. (1977).
//     Factorization Methods for Discrete Sequential Estimation, Academic Press, Inc.
//
namespace sqrt_kf
{
	//-------------------------------------------------------------------------
	// @brief Computes an upper-triangular matrix S such that P = S * S'
	//
	// @param[in]  n   Dimension of the matrix
	// @param[in]  P   (n x n) symmetric positive definite matrix
	// @param[out] S   Upper-triangular Cholesky factor
	//
	// @return true if positive definite
	//
	bool CholUT(int n, const double P[], double S[]);

	//-------------------------------------------------------------------------
	// @brief This function can be used to implement squre-root Kalman prediction with
	//     A = [ Phi*C(+),   G*Cq] 
	// where C(+) C(+)' = P(+) and Cq Cq'= Q.
	// The resulting C(-) is stored in the upper-right triangular part.
	//
	void HouseholderTriangular(double A[], int rows, int cols);

	//-------------------------------------------------------------------------
	// Inputs:  A(n x n) = Phi*C(+), B(n x r) = G*Cq
	// Outputs: A = C(-), B = zeroed
	void SchmidtHouseholderPredict(
		double A[], double B[],
		int n, int r);

	//-------------------------------------------------------------------------
	// @brief Carlson's Squar-Root Measurement Update
	//
	// @param[in]     n  Number of states
	// @param[in/out] x  State vector (n x 1)
	// @param[in/out] C  Cholesky factor of P (n x n)
	// @param[in]     z  Measurement
	// @param[in]     H  Sensitivity matrix (1 x n)
	// @param[in]     R  Measurement uncertainty
	//
	void CarlsonMeasUpdate(
		int n,
		double x[], double C[], 
		double z, double H[], double R);

	//-------------------------------------------------------------------------
	// @brief Computes U and D matrix such that P = U * D * U'
	// 
	// @param[in]  n  Dimension of matrix
	// @param[in]  P  (n x n) symmetric positive definite matrix
	// @param[out] U  Upper-triangular factor
	// @param[out] D  Diagonal factor; store only diagonal elements
	//
	void UD_Factor(int n, const double P[], double U[], double D[]);

	//-------------------------------------------------------------------------
	// @brief Kalman covariance prediction in UD-factored form:
	//
	// @param[in]     n    Number of the states
	// @param[in]     p    Number of system noise parameters
	// @param[in]     PhiU Phi * U
	// @param[in]     D    Diagonal factor of the covariance
	// @param[out]    Um   Predicted U: U(-)
	// @param[out]    Dm   Predicted D: D(-)
	// @param[in/out] GUq  G * Uq, where Q = Uq * Dq * Uq'
	// @param[in/out] Dq   Diagonal factor of Q
	// 
	void UD_Predict(
		int n, int p,
		double PhiU[], const double D[],
		double Um[], double Dm[],
		double GUq[], double Dq[]);

	//-------------------------------------------------------------------------
	// @brief Kalman measurement update in UD-factored form:
	//
	// @param[in]     n Number of states
	// @param[in/out] x state vector
	// @param[in/out] U factor of the covariance
	// @param[in/out] D factor of the covariance
	// @param[in]     z Measurement
	// @param[in]     H Sensitivity matrix
	// @param[in]     R Measurement uncertainty
	//
	void BiermanUpdate(
		int n, double x[],
		double U[], double D[],
		double z, double H[], double R);

	//-------------------------------------------------------------------------
   // @brief Kalman measurement update in UD-factored form with UD storage in
	//        generic form.
	// Reference: [2], pp. 100-101.
	// 
	// @param[in]     n     Number of states
	// @param[in/out] x     State vecor
	// @param[in/out] UD    UD-factored covariance with diagonals storing D
	// @param[in]     f_idx UD element indexing function
	// @param[in]     z     Measurement
	// @param[in]     H     Measurement sensitivity matrix
	// @param[in]     R     Measurement uncertainty variance
	//
	void UD_Update(
		int n, double x[], double UD[],
		std::function<int(int, int)> f_idx,
		double z, double H[], double R);
}