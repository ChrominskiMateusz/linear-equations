#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <complex>
#include "Matrix.h"

void Jacobi (Matrix&, Matrix&, Matrix&);

int main ()
{
	constexpr int N = 10;
	constexpr int a1 = 5 + 3;
	constexpr int a2 = -1;
	constexpr int a3 = -1;

	Matrix A (N, N);
	Matrix b (1, N);
	Matrix x (1, N);

	for (int i = 0; i < A.y; i++)
	{
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
	}

	for (int i = 0; i < b.y; i++)
	{
		for (int j = 0; j < b.x; j++)
		{
			b.matrix[i][j] = sin (i * 6);
		}
	}


	Matrix tmp_a (3, 3);
	tmp_a.matrix[0][0] = 1.5;
	tmp_a.matrix[0][1] = 1.1;
	tmp_a.matrix[0][2] = -0.3;
	tmp_a.matrix[1][0] = 1.1;
	tmp_a.matrix[1][1] = 1.5;
	tmp_a.matrix[1][2] = 1.1;
	tmp_a.matrix[2][0] = 0.3;
	tmp_a.matrix[2][1] = 1.1;
	tmp_a.matrix[2][2] = 2.0;
	
	Matrix tmp_b (1, 3);
	tmp_b.matrix[0][0] = 1;
	tmp_b.matrix[1][0] = 2;
	tmp_b.matrix[2][0] = 1;

	Matrix tmp_x (1, 3);
	tmp_x.matrix[0][0] = 1;
	tmp_x.matrix[1][0] = 1;
	tmp_x.matrix[2][0] = 1;

	

	Jacobi (tmp_a, tmp_x, tmp_b);
	
	
	return 0;
}


void Jacobi (Matrix& A, Matrix& x, Matrix& b)
{
	std::cout << "Jacobi method below!\n\n";
	double stop = 0.00000001;
	double residuum_norm = 1;
	double sum_to_i;
	double sum_from_i;
	Matrix prev (x.x, x.y);
	Matrix residuum (x.x, x.y);
	std::cout << "A matrix:\n";
	A.print ();
	std::cout << "b matrix:\n";
	b.print ();
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
	std::cout << "x calculated by Jacobi method: \n";
	x.print ();
	std::cout << "in iteration number: " << count << "\n\n";
}