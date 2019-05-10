#pragma once
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

class Data
{
public:
	long double H;
	long double L;
	long double nH;
	long double nL;
	long double K;
	long double T;
	long double C;
	long double RO;
	long double A;
	long double dT;
	long double T0;
	long double universal_ksi[4][4];
	long double universal_eta[4][4];
	long double  **Global_H;
	long double  **Global_C;
	long double  **Global_HB;
	long double  **H_final;
	long double *temp;
	long double *GlobalP;
	Data();
	~Data();
};

