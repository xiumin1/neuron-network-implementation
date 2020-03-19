#pragma once
#pragma once
#include<iostream>
#include<time.h>
#include<stdlib.h>
#include<vector>
using namespace std;
class Matrix
{
private:
	int rowNum;
	int colNum;
	int size;
	vector<double> mat;
public:
	Matrix() { rowNum = 0; colNum = 0; size = 0; }
	Matrix(int row, int col, bool rand_init);
	Matrix(const Matrix &a); //copy construnctor

	vector<double> getMat() const { return mat; }
	int getRow() const { return rowNum; }
	int getCol() const { return colNum; }
	int getSize() const { return size; }

	void setV(int row, int col, double v) { mat[row*colNum + col] = v; }
	double getV(int row, int col) const { return mat[row*colNum + col]; }
	void print();

	//operator overloading
	void Add(const Matrix& a, const Matrix&b);
	void Add(const Matrix& a, double b);
	//void Add(const Matrix& a, vector<double>b);
	void Times(const Matrix& a, const Matrix& b);
	void Times(const Matrix& a, double b);
	void Divide(const Matrix& a, const Matrix& b);
	bool operator==(const Matrix& a) const;
	void operator=(const Matrix& a);

	void Multiply(const Matrix& a, const Matrix& b);
	void Transpose(const Matrix& a);

	double Exp(double x) const;
	void Log(const Matrix& a);
	//relu, sigmoid, tanh, 3 activation functions for neuron network layer parameters calculation
	void Relu(const Matrix& a);
	void Sigmoid(const Matrix& a);
	void Tanh(const Matrix& a);
	void DRelu(const Matrix& a);
	void DSigmoid(const Matrix& a);
	void DTanh(const Matrix& a);

	// for bias b to average the cololumns and then averaged to the same value, then expand the same value back to all cololums
	//example, b=[2 3 5 6 8]==> b=[4.8 4.8 4.8 4.8 4.8]
	void Reassign(const Matrix& a);

	void SplitCol(const Matrix& a, int start, int end);  //[start, end]

	~Matrix();
};

