#include "Matrix.h"
#include <cmath>
#include <iomanip>
#include <iostream>

void Matrix::print ()
{
	for (int i = 0; i < x; i++)
		if (i == 0)
			std::cout << " _      ";
		else
			std::cout << "       ";

	std::cout << " _\n";
	for (int i = 0; i < y - 1; i++)
	{
		std::cout << "|  ";
		for (int j = 0; j < x; j++)
			std::cout << matrix[i][j] << "  ";
		std::cout << "|\n";
	}
	std::cout << "|_ ";
	for (int j = 0; j < x - 1; j++)
		std::cout << matrix[y - 1][j] << "  ";
	std::cout << matrix[y - 1][x - 1] << " _|\n\n\n";
}

Matrix Matrix::operator*(const Matrix& secound)
{
	Matrix to_ret (secound.x, y);
	if (x != secound.y)

	{
		std::cout << "Size not correct\n";
		return to_ret;
	}

	for (int i = 0; i < y; i++)
		for (int j = 0; j < secound.x; j++)
			for (int k = 0; k < secound.y; k++)
				to_ret.matrix[i][j] += matrix[i][k] * secound.matrix[k][j];

	return to_ret;
}

Matrix Matrix::operator+(const Matrix& secound)
{
	Matrix to_ret (x, y);
	if (x != secound.x || y != secound.y)
	{
		std::cout << "Size not correct\n";
		return to_ret;
	}
	
	for (int i = 0; i < y; i++)
		for (int j = 0; j < x; j++)
			to_ret.matrix[i][j] = matrix[i][j] + secound.matrix[i][j];
	
	return to_ret;
}

Matrix Matrix::operator-(const Matrix& secound)
{
	Matrix to_ret (x, y);
	if (x != secound.x && y != secound.y)
	{
		std::cout << "Size not correct\n";
		return to_ret;
	}

	for (int i = 0; i < y; i++)
		for (int j = 0; j < x; j++)
			to_ret.matrix[i][j] = matrix[i][j] - secound.matrix[i][j];
	
	return to_ret;
}

Matrix& Matrix::operator=(const Matrix& secound)
{
	double** k = new double*[secound.y];
	for (int i = 0; i < secound.y; i++)
		k[i] = new double[secound.x];

	for (int i = 0; i < secound.y; i++)
		for (int j = 0; j < secound.x; j++)
			k[i][j] = secound.matrix[i][j];

	this->~Matrix ();

	matrix = k;
	y = secound.y;
	x = secound.x;
	return *this;
}

Matrix::Matrix (Matrix&& m)
	:matrix{ m.matrix },
	x{ m.x },
	y{ m.y }
{
	m.matrix = nullptr;
}

Matrix::Matrix (const Matrix& m)
	:x{ m.x },
	 y{ m.y }
{
	matrix = new double *[y];
	for (int i = 0; i < y; i++)
	{
		matrix[i] = new double[x];
		for (int j = 0; j < x; j++)
			matrix[i][j] = m.matrix[i][j];
	}
}

Matrix::Matrix (int x, int y) :
	x (x), y (y)
{
	matrix = new double *[y];
	for (int i = 0; i < y; i++)
	{
		matrix[i] = new double[x];
		for (int j = 0; j < x; j++)
			matrix[i][j] = 0;
	}
}

Matrix::~Matrix ()
{
	if (matrix != nullptr)
		for (int i = 0; i < y; i++)
			if (matrix[i] != nullptr)
			{
				delete[] matrix[i];
				matrix[i] = nullptr;
			}
	delete[] matrix;
}