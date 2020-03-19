#pragma once
#include "Matrix.h"
#include<vector>
#define relu 1
#define sigmoid 2
#define tanh 3

class Layer
{
	//MatMath matMath;
	// parameters need for forward propagation
	Matrix* weight;
	Matrix* bias;
	Matrix* linearValue;
	Matrix* activatedValue;
	// parameters need for backward propagation
	Matrix* d_weight;
	Matrix* d_bias;
	Matrix* d_linearValue;
	Matrix* d_activatedValue;
	int activateFun;
	int batchNum; // sample numbers passed into the neuronnetwork


public:
	Layer(const Matrix& w, const Matrix& b, int actfun);

	void setWeight(const Matrix& w) { *weight = w; }
	void setBias(const Matrix& b) { *bias = b; }
	void setLinearValue(const Matrix& lv) { *linearValue = lv; }
	void setActivatedValue(const Matrix& av) { *activatedValue = av; }
	Matrix* getWeight() { return weight; }
	Matrix* getBias() { return bias; }
	Matrix* getLinearValue() { return linearValue; }
	Matrix* getActivatedValue() { return activatedValue; }

	int getActFun() { return activateFun; }
	int getBatchNum() { return batchNum; }
	void calcValues(const Matrix& batchinput);
	// for backward propagation
	void setDWeight(const Matrix& dw) { *d_weight = dw; }
	void setDBias(const Matrix& db) { *d_bias = db; }
	void setDLinearValue(const Matrix& dlv) { *d_linearValue = dlv; }
	void setDActivatedValue(const Matrix& dav) { *d_activatedValue = dav; }
	Matrix* getDWeight() { return d_weight; }
	Matrix* getDBias() { return d_bias; }
	Matrix* getDLinearValue() { return d_linearValue; }
	Matrix* getDActivatedValue() { return d_activatedValue; }
};

class NeuronNetwork
{
	int* neuronNums; // use to store neuron numbers for each layer {2,3, 4,5, ...}, the input layer will not counted as a layer, so neuronNums[0] means the first hidden layer
	int* actfunNums;
	vector<Layer*> Layers; // use to store all the layer information
	int layerNum; //number of layers, except the first input layer
	double learningRate;
	int batchSize;
	int x_dim;
	//Matrix* cost;
	vector<double> cost;
	vector<Matrix*> temp;
	//vector<Matrix> savemodel;
	//Matrix  test;
	// 0 layer is the first hidden layer
public:
	NeuronNetwork();
	NeuronNetwork(int x_dim, int layernum, int neuronnums[], int actfunnums[], double learningrate, int batchsize);
	Layer* AddLayer(int neuron, int neuron_next, int actFun);
	void FeedForward(const Matrix& x_train); //here call the calactivated value function for each layer to calculated activatedValue for each layer
	double CrossEntropyCost(const Matrix& y_train, const Matrix& y_activate);// |real value-calculated value|
	void BackPropagation(const Matrix& x_train, const Matrix& y_train);
	void ParameterUpdate();

	void Training(int n_iteration, const Matrix& input, const Matrix& ouput);
	void Evaluation(const Matrix& x_test, const Matrix& y_test);
	void Prediction(const Matrix& x_data);
	void PrintCost();
	void SaveModel();
	void LoadModel();
	//~NeuronNetwork();
};




