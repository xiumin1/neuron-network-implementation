//#include "MyRL_QLearning.h"
#include "Matrix.h"
//#include "NeuronNetwork.h"
#include "NeuronNetwork.h"
#include<iostream>
#include<fstream>
#include<time.h>
#include<stdlib.h>

using namespace std;
int reverseInt(int i)
{
	unsigned char c1, c2, c3, c4;

	c1 = i & 255;
	c2 = (i >> 8) & 255;
	c3 = (i >> 16) & 255;
	c4 = (i >> 24) & 255;

	return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}
void read_mnist_image(string filename, vector<vector<double>> &vec)
{
	//cout << "in" << endl;
	//ifstream file("D:\\Box Sync\\Research Projects\\7.07.2019-projects\\RL_c++\\Test_c++\\data\\train-images-idx3-ubyte\\train-images.idx3-ubyte");
	ifstream file(filename, ios::binary);
	//ifstream file(/*full_path*/"t10k-images-idx3-ubyte.gz");
	if (file.is_open())
	{
		cout << "file open" << endl;
		int magic_number = 0;
		int number_of_images = 0;
		int n_rows = 0;
		int n_cols = 0;
		file.read((char*)&magic_number, sizeof(magic_number));
		magic_number = reverseInt(magic_number);
		file.read((char*)&number_of_images, sizeof(number_of_images));
		number_of_images = reverseInt(number_of_images);
		cout << number_of_images << endl;
		file.read((char*)&n_rows, sizeof(n_rows));
		n_rows = reverseInt(n_rows);
		cout << n_rows << endl;
		file.read((char*)&n_cols, sizeof(n_cols));
		n_cols = reverseInt(n_cols);
		cout << n_cols << endl;
		for (int i = 0; i<number_of_images; ++i)
		{
			vector<double> tp;
			for (int r = 0; r<n_rows; ++r)
			{
				for (int c = 0; c<n_cols; ++c)
				{
					unsigned char temp = 0;
					file.read((char*)&temp, sizeof(temp));
					tp.push_back((double)temp);


				}
			}
			vec.push_back(tp);
		}
	}
}
void read_mnist_label(string filename, vector<double> &vec)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		int n_rows = 0;
		int n_cols = 0;
		file.read((char*)&magic_number, sizeof(magic_number));
		magic_number = reverseInt(magic_number);
		file.read((char*)&number_of_images, sizeof(number_of_images));
		number_of_images = reverseInt(number_of_images);
		for (int i = 0; i < number_of_images; ++i)
		{
			unsigned char temp = 0;
			file.read((char*)&temp, sizeof(temp));
			vec.push_back((double)temp);
			//vec[i] = (double)temp;
		}
	}
}
vector<int> sampleNum(int num)
{
	vector<int> indexs;
	srand((int)time(0));

	for (int i = 0; i < num; i++)
	{
		indexs.push_back(rand() % 59999);
	}
	return indexs;
}
int main(int argc, const char * argv[])
{
	vector<int> indexs;
	int num = 1000;
	indexs = sampleNum(num);

	string filename_image = "..\\test_nn\\data\\train-images-idx3-ubyte\\train-images.idx3-ubyte";
	string filename_label = "..\\test_nn\\data\\train-labels-idx1-ubyte\\train-labels.idx1-ubyte";
	vector<vector<double>> image;
	read_mnist_image(filename_image, image);

	vector < double> label;
	read_mnist_label(filename_label, label);

	vector < vector < double >> image_sample;
	vector<double> label_sample;
	for (int i = 0; i < num; i++)
	{
		image_sample.push_back(image[indexs[i]]);
		label_sample.push_back(label[indexs[i]]);

		for (int j = 0; j < image_sample[0].size(); j++)
		{
			image_sample[i][j] /= (double)255; // here normalize the input data to (0-1)
		}
	}

	int sample = image_sample.size();
	// prediction the training samples
	Matrix* x = new Matrix(sample*0.7, image_sample[0].size(), false);
	Matrix* y = new Matrix(sample*0.7, 10, false);
	Matrix* x_input = new Matrix();
	Matrix* y_input = new Matrix();

	for (int i = 0; i < x->getRow(); i++)
	{
		for (int j = 0; j < x->getCol(); j++)
		{
			x->setV(i, j, image_sample[i][j]);
		}
	}
	// here the output has 10 neurons, only the position = the digit value equals 1, all others 0.
	// exp:  0 0 0 1 0 0 ==>3
	for (int i = 0; i < y->getRow(); i++)
	{
		//int digit = label[i];
		for (int j = 0; j < 10; j++)
		{
			y->setV(i, j, 0);
		}
		y->setV(i, label_sample[i], 1);
	}
	x_input->Transpose(*x);
	y_input->Transpose(*y);

	// prepare the evaluation samples
	Matrix* x1 = new Matrix(sample*0.2, image_sample[0].size(), false);
	Matrix* y1 = new Matrix(sample*0.2, 1, false);
	Matrix* x_test = new Matrix();
	Matrix* y_test = new Matrix();

	for (int i = 0; i < x1->getRow(); i++)
	{
		for (int j = 0; j < x1->getCol(); j++)
		{
			x1->setV(i, j, image_sample[i + x->getRow()][j]);
		}
	}

	for (int i = 0; i < y1->getRow(); i++)
	{
		y1->setV(i, 0, label_sample[i + x->getRow()]);
	}
	x_test->Transpose(*x1);
	y_test->Transpose(*y1);

	// prepare the prediction samples
	Matrix* x2 = new Matrix(sample*0.1, image_sample[0].size(), false);
	Matrix* x_predict = new Matrix();

	for (int i = 0; i < x2->getRow(); i++)
	{
		for (int j = 0; j < x2->getCol(); j++)
		{
			x2->setV(i, j, image_sample[i + x->getRow() + x1->getRow()][j]);
		}
	}
	x_predict->Transpose(*x2);

	// start training
	int x_dim = 784;// image[0].size();
	int batchsize = 32;
	double learningrate = 0.3f;
	int epoches = 25;
	int layernum = 3;

	int* neuronnums = new int[layernum]; // to change the neuron layers and nueron numbers
	neuronnums[0] = 128;
	neuronnums[1] = 128;
	neuronnums[2] = 10;

	int* actfunnums = new int[layernum]; // to change the activation function
	actfunnums[0] = 1;
	actfunnums[1] = 1;
	actfunnums[2] = 2;

	NeuronNetwork* nn = new NeuronNetwork(x_dim, layernum, neuronnums, actfunnums, learningrate, batchsize);
	cout << "training start" << endl;
	nn->Training(epoches, *x_input, *y_input);
	cout << "evaluation start" << endl;
	nn->Evaluation(*x_test, *y_test);
	cout << "prediction start" << endl;
	nn->Prediction(*x_predict);

	nn->LoadModel();

	cout << "done" << endl;

	getchar();

	return 0;
}