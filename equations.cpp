#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <complex>
#include <chrono>
#include "Matrix.h"

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using namespace std::literals::chrono_literals;

void Jacobi (Matrix&, Matrix&, Matrix&);
void Gauss_Seidel (Matrix&, Matrix&, Matrix&);
void LU (Matrix&, Matrix&);

int main ()
{
	std::cout << std::showpos << std::fixed << std::setprecision (2);      //   Setting the way of printing numbers

	constexpr int N = 918;
	constexpr int a1 = 5 + 3;
	constexpr int a2 = -1;
	constexpr int a3 = -1;

	Matrix A (N, N);
	Matrix b (1, N);
	Matrix x (1, N);

	for (int i = 0; i < A.y; i++)
		for (int j = 0; j < A.x; j++)
		{
			if (i == j)
				A.matrix[i][j] = a1;
			else if (i == j + 1)
				A.matrix[i][j] = a2;
			else if (j == i + 1)
				A.matrix[i][j] = a2;
			else if (i == j + 2)
				A.matrix[i][j] = a3;
			else if (j == i + 2)
				A.matrix[i][j] = a3;
		}

	for (int i = 0; i < b.y; i++)
		for (int j = 0; j < b.x; j++)
			b.matrix[i][j] = sin (i * 6);

	Jacobi (A, x, b);

	for (int i = 0; i < x.y; i++)
		for (int j = 0; j < x.x; j++)
			x.matrix[i][j] = 0;

	Gauss_Seidel (A, x, b);

// ========================================================================================   Next  ================================== //

	constexpr int a1c = 3;

	for (int i = 0; i < A.y; i++)
		for (int j = 0; j < A.x; j++)
			if (i == j)
				A.matrix[i][j] = a1c;

	Matrix L (4, 4);
	Matrix U (4, 4);
	
	U.matrix[0][0] = 2;
	U.matrix[0][1] = 1;
	U.matrix[0][2] = 1;
	U.matrix[0][3] = 0;
	U.matrix[1][0] = 4;
	U.matrix[1][1] = 3;
	U.matrix[1][2] = 3;
	U.matrix[1][3] = 1;
	U.matrix[2][0] = 8;
	U.matrix[2][1] = 7;
	U.matrix[2][2] = 9;
	U.matrix[2][3] = 5;
	U.matrix[3][0] = 6;
	U.matrix[3][1] = 7;
	U.matrix[3][2] = 9;
	U.matrix[3][3] = 8;

	U.print ();

	for (int i = 0; i < L.y; i++)
		for (int j = 0; j < L.x; j++)
			if (i == j)
				L.matrix[i][j] = 1;

	LU (L, U);
	L.print ();
	U.print ();


	
	return 0;
}


void Jacobi (Matrix& A, Matrix& x, Matrix& b)
{
	time_point<Clock> start = Clock::now ();
	std::cout << "\nJacobi method below!\n\n";
	double stop = 0.000000001;
	double residuum_norm = 1;
	double sum_to_i;
	double sum_from_i;
	Matrix prev (x.x, x.y);
	Matrix residuum (x.x, x.y);
	int count = 0;
	while (residuum_norm > stop)
	{
		count++;
		for (int i = 0; i < prev.y; i++)				// Remember last iteration
			prev.matrix[i][0] = x.matrix[i][0];			// of Jacobi

		for (int i = 0; i < x.y; i++)
		{
			sum_from_i = 0;
			sum_to_i = 0;

			for (int j = 0; j <= i - 1; j++)
				sum_to_i += A.matrix[i][j] * prev.matrix[j][0];

			for (int j = i + 1; j < x.y; j++)
				sum_from_i += A.matrix[i][j] * prev.matrix[j][0];

			x.matrix[i][0] = b.matrix[i][0] - sum_to_i - sum_from_i;
			x.matrix[i][0] /= A.matrix[i][i];
		}

 		residuum = A*x;
		residuum = residuum - b;
		residuum_norm = 0;
		for (int j = 0; j < residuum.y; j++)
			residuum_norm += residuum.matrix[j][0] * residuum.matrix[j][0];
		residuum_norm = sqrt (residuum_norm);
	}
	std::cout << "in iteration number: " << count << "\n";
	time_point<Clock> end = Clock::now ();
	milliseconds time = duration_cast<milliseconds>(end - start);
	std::cout << time.count () << "ms\n\n";
}

void Gauss_Seidel (Matrix& A, Matrix& x, Matrix& b)
{
	time_point<Clock> start = Clock::now ();
	std::cout << "\n\nGauss-Seidel method below!\n\n";
	double stop = 0.000000001;
	double residuum_norm = 1;
	double sum_to_i;
	double sum_from_i;
	Matrix prev (x.x, x.y);
	Matrix residuum (x.x, x.y);
	int count = 0;
	while (residuum_norm > stop)
	{
		count++;
		for (int i = 0; i < prev.y; i++)				// Remember last iteration
			prev.matrix[i][0] = x.matrix[i][0];			// of Gauss-Seidel

		for (int i = 0; i < x.y; i++)
		{
			sum_from_i = 0;
			sum_to_i = 0;

			for (int j = 0; j <= i - 1; j++)
				sum_to_i += A.matrix[i][j] * x.matrix[j][0];

			for (int j = i + 1; j < x.y; j++)
				sum_from_i += A.matrix[i][j] * prev.matrix[j][0];

			x.matrix[i][0] = b.matrix[i][0] - sum_to_i - sum_from_i;
			x.matrix[i][0] /= A.matrix[i][i];
		}

		residuum = A*x;
		residuum = residuum - b;
		residuum_norm = 0;
		for (int j = 0; j < residuum.y; j++)
			residuum_norm += residuum.matrix[j][0] * residuum.matrix[j][0];
		residuum_norm = sqrt (residuum_norm);
	}
	std::cout << "in iteration number: " << count << "\n";
	time_point<Clock> end = Clock::now ();
	milliseconds time = duration_cast<milliseconds>(end - start);
	std::cout << time.count () << "ms\n\n";
}

void LU (Matrix& L, Matrix& U)
{
	for(int k = 0; k < L.y - 1; k++)
		for (int j = k + 1; j < L.x; j++)
		{
			L.matrix[j][k] = U.matrix[j][k] / U.matrix[k][k];
			for (int i = k; i < L.x; i++)
				U.matrix[j][i] = U.matrix[j][i] - (L.matrix[j][k] * U.matrix[k][i]);
		}
}