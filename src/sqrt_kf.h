#pragma once

// Eun Hwan Shin, 2022
// Reference: Mohinder S. Grewal and Angus P. Andrews. (2001).
//   Kalman Filgering: Theorey and Practice, 2nd Ed., John-Wiley & Sons, Inc.

namespace sqrt_kf
{
	//-------------------------------------------------------------------------
	// This function can be used to implement squre-root Kalman prediction with
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
	// Carlson's Squar-Root Measurement Update
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
}