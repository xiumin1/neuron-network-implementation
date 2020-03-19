#include "Matrix.h"

Matrix::Matrix(int row, int col, bool rand_init)
{
	rowNum = row;
	colNum = col;
	size = row * col;
	//mat.reserve(size);
	mat.resize(size);
	srand((int)time(0));
	if (rand_init)
	{
		for (int i = 0; i < size; i++)
		{
			//mat.emplace_back((double)rand() / (double)RAND_MAX);
			mat[i] = (double)rand() / (double)RAND_MAX;
			mat[i] /= 1000.0;
		}
	}
	else
	{
		for (int i = 0; i<size; i++)
		{
			//mat.emplace_back(0.0);
			mat[i] = 0;
		}
	}
}
Matrix::Matrix(const Matrix &a)

{
	cout << "copy called !!!!!!!!" << endl;
	rowNum = a.rowNum;
	colNum = a.colNum;
	size = a.getSize();
	mat = a.mat;
}
void Matrix::print()
{
	cout << "the matrix is: " << endl;
	for (int i = 0; i<rowNum; i++)
	{
		for (int j = 0; j<colNum; j++)
		{
			cout << mat[i*colNum + j] << " ";
		}
		cout << endl;
	}
}
void Matrix::Add(const Matrix& a, const Matrix&b)
{
	if (a.getRow() == b.getRow())
	{
		this->size = a.getSize();
		this->rowNum = a.getRow();
		this->colNum = a.getCol();
		//this->mat.reserve(this->size);
		this->mat.resize(this->size);

		if (a.getCol() == b.getCol())
		{
			for (int i = 0; i < this->size; i++)
			{
				//this->mat.emplace_back(a.mat[i]+b.mat[i]);
				this->mat[i] = a.mat[i] + b.mat[i];
			}
		}
		else if (b.getCol() == 1)
		{
			for (int i = 0; i < this->size; i++)
			{

				this->mat[i] = a.mat[i] + b.mat[i / this->colNum];

			}
		}
		else
		{
			cout << "Matrix dimensions do not match!" << endl;
		}

	}
	else
	{
		cout << "Matrix dimensions do not match!" << endl;
	}
}
void Matrix::Add(const Matrix& a, double b)
{
	this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();

	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	for (int i = 0; i < this->size; i++)
	{
		//this->mat.emplace_back(a.mat[i] + b);
		this->mat[i] = a.mat[i] + b;
	}
}
void Matrix::Times(const Matrix& a, const Matrix& b)
{
	if (a.getRow() == b.getRow() && a.getCol() == b.getCol())
	{
		this->size = a.getSize();
		this->rowNum = a.getRow();
		this->colNum = a.getCol();
		//this->mat.reserve(this->size);
		this->mat.resize(this->size);
		for (int i = 0; i < size; i++)
		{
			//this->mat.emplace_back(a.mat[i] * b.mat[i]);
			this->mat[i] = a.mat[i] * b.mat[i];
		}
	}
	else
	{
		cout << "Matrix dimensions do not match!" << endl;
	}
}
void Matrix::Times(const Matrix& a, double b)
{
	//this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();
	this->size = this->rowNum*this->colNum;
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	for (int i = 0; i < size; i++)
	{
		//this->mat.emplace_back(a.mat[i] * b);
		this->mat[i] = a.mat[i] * b;
	}
}
void Matrix::Divide(const Matrix& a, const Matrix& b)
{
	if (a.getRow() == b.getRow() && a.getCol() == b.getCol())
	{
		this->size = a.getSize();
		this->rowNum = a.getRow();
		this->colNum = a.getCol();
		//this->mat.reserve(this->size);
		this->mat.resize(this->size);
		for (int i = 0; i < size; i++)
		{
			//this->mat.emplace_back(a.mat[i] / b.mat[i]);
			this->mat[i] = a.mat[i] / b.mat[i];
		}
	}
	else
	{
		cout << "Matrix dimensions do not match!" << endl;
	}
}
bool Matrix::operator==(const Matrix& a) const
{
	if (this->getRow() == a.getRow() && this->getCol() == a.getCol())
	{
		for (int i = 0; i < this->size; i++)
		{
			if (this->mat[i] != a.mat[i])
			{
				return false;
			}
		}
	}
	return true;
}
void Matrix::operator=(const Matrix& a)
{
	this->size = a.getSize();
	this->colNum = a.getCol();
	this->rowNum = a.getRow();
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	for (int i = 0; i < this->size; i++)
	{
		//this->mat.emplace_back(a.mat[i]);
		this->mat[i] = a.mat[i];
	}
}

void Matrix::Multiply(const Matrix& a, const Matrix& b)
{
	int a_row = a.getRow();
	int a_col = a.getCol();
	int b_row = b.getRow();
	int b_col = b.getCol();

	if (a_col == b_row)
	{
		this->rowNum = a_row;
		this->colNum = b_col;
		this->size = a_row * b_col;
		this->mat.resize(this->size);

		double temp = 0.0;
		for (int i = 0; i<a_row; i++)
		{
			for (int j = 0; j<b_col; j++)
			{
				temp = 0;
				for (int k = 0; k<b_row; k++)
				{
					//tempv += this->getV(i, k)*other.getV(k, j);
					temp += a.getV(i, k)* b.getV(k, j);
				}
				this->setV(i, j, temp);
			}
		}
	}
	else
	{
		cout << "Matrix dimensions do not match!" << endl;
	}
}
void Matrix::Transpose(const Matrix& a)
{
	int row = a.getRow();
	int col = a.getCol();
	this->rowNum = col;
	this->colNum = row;
	this->size = row * col;
	this->mat.resize(this->size);

	for (int i = 0; i < col; i++)
	{
		for (int j = 0; j < row; j++)
		{
			this->setV(i, j, a.getV(j, i));
		}
	}
}

double Matrix::Exp(double x) const
{
	double sum = 1.0; // initialize sum of series
	int item = 10; //to add up the first 10 items
				   //e^x=1 + (x/1) (1 + (x/2) (1 + (x/3) (........) ) )
	for (int i = 10 - 1; i > 0; --i)
		sum = 1.0 + x * sum / (double)i;
	return sum;
}
void Matrix::Log(const Matrix& a)
{
	this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	for (int i = 0; i < this->size; i++)
	{
		//this->mat.emplace_back(a.mat[i]);
		this->mat[i] = log(a.mat[i]);
	}
}
void Matrix::Relu(const Matrix& a)
{
	this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	double temp;
	for (int i = 0; i < this->size; i++)
	{
		if (a.mat[i] < 0.0) temp = 0.0;
		else temp = a.mat[i];
		//this->mat.emplace_back(temp);
		this->mat[i] = temp;
	}
}

void Matrix::Sigmoid(const Matrix& a)
{
	this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	double temp, temp1;
	for (int i = 0; i < this->size; i++)
	{
		temp1 = 1.0 + exp(a.mat[i] * (-1.0));
		temp = 1.0 / temp1;
		//this->mat.emplace_back(temp);
		this->mat[i] = temp;
	}
}
void Matrix::Tanh(const Matrix& a)
{
	this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	double temp, temp1;
	for (int i = 0; i < this->size; i++)
	{
		temp1 = 1.0 + exp(2.0*a.mat[i] * (-1.0));
		temp = 2.0 / temp1 - 1.0;
		//this->mat.emplace_back(temp);
		this->mat[i] = temp;
	}
}
void Matrix::DRelu(const Matrix& a)
{
	this->size = a.getSize();
	this->rowNum = a.getRow();
	this->colNum = a.getCol();
	//this->mat.reserve(this->size);
	this->mat.resize(this->size);
	double temp;
	for (int i = 0; i < this->size; i++)
	{
		if (a.mat[i] < 0.0) temp = 0.0;
		else temp = 1.0;
		//this->mat.emplace_back(temp);
		this->mat[i] = temp;
	}
}
void Matrix::DSigmoid(const Matrix& a)
{
	this->Sigmoid(a);
	double temp;
	for (int i = 0; i < this->size; i++)
	{
		temp = this->mat[i] * (1.0 - this->mat[i]);
		//this->mat.emplace_back(temp);
		this->mat[i] = temp;
	}
}
void Matrix::DTanh(const Matrix& a)
{
	this->Tanh(a);
	double temp;
	for (int i = 0; i < this->size; i++)
	{
		temp = 1.0 - this->mat[i] * this->mat[i];
		//this->mat.emplace_back(temp);
		this->mat[i] = temp;
	}
}
void Matrix::Reassign(const Matrix& a)
{
	/*
	int row = a.getRow();
	int col = a.getCol();

	this->rowNum = row;
	this->colNum = col;
	this->size = row * col;
	this->mat.resize(this->size);

	double temp;
	for (int i = 0; i < row; i++)
	{
	temp = 0;
	for (int j = 0; j < col; j++)
	{
	temp += a.mat[i*col + j];
	}
	temp /= (double)col;

	for (int j = 0; j < col; j++)
	{
	this->mat[i*col + j] = temp;
	}
	}
	*/
	int row = a.getRow();
	int col = a.getCol();

	this->rowNum = row;
	this->colNum = 1;
	this->size = row * 1;
	this->mat.resize(this->size);

	double temp;
	for (int i = 0; i < row; i++)
	{
		temp = 0;
		for (int j = 0; j < col; j++)
		{
			temp += a.mat[i*col + j];
		}
		temp /= (double)col;

		this->mat[i] = temp;
	}
}
void Matrix::SplitCol(const Matrix& a, int start, int end)
{
	int row = a.getRow();
	int col = a.getCol();

	if (start >= 0 && end >= 0 && start <= end && end < col)
	{
		this->rowNum = row;
		this->colNum = end - start + 1;
		this->size = this->rowNum*this->colNum;
		this->mat.resize(this->size);

		for (int i = 0; i < this->rowNum; i++)
		{
			for (int j = 0; j < this->colNum; j++)
			{
				this->mat[i*this->colNum + j] = a.mat[i*col + j + start];
			}
		}
	}
	else
	{
		cout << "Split range is wrong!!" << endl;
	}
}
Matrix::~Matrix()
{
	//delete mat;
	cout << "called" << endl;
	//delete temp;
}
