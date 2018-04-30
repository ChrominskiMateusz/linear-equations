#pragma once

class Matrix {
public:
	int x;
	int y;
	double **matrix;
	void print ();
	Matrix operator*(const Matrix&);
	Matrix operator+(const Matrix&);
	Matrix operator-(const Matrix&);
	Matrix& operator=(const Matrix&);
	Matrix (int, int);
	~Matrix ();
};

