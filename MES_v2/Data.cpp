#include "Data.h"




Data::Data()
{
	fstream plik;
	plik.open("dane.txt", ios::in);
	if (!plik.good()) {
		cout << "plik nie dzia³a" << endl;
		return;
	}

	string reader;
	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	this->H = stod(reader);


	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	this->L = stod(reader);


	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	reader[2] = '0';
	this->nH = stod(reader);

	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	reader[2] = '0';
	this->nL = stod(reader);

	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	this->K = stod(reader);

	
	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	T = stof(reader);

	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	C = stof(reader);

	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	reader[2] = '0';
	RO = stof(reader);
	
	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	A = stof(reader);
	
	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	reader[2] = '0';
	dT = stof(reader);

	plik >> reader;
	reader[0] = '0';
	reader[1] = '0';
	reader[2] = '0';
	T0 = stof(reader);


	plik.close();

	long double wymiar = nH * nL;

	Global_H = new long double *[wymiar];    //deklaracja GlobalH
	for (int i = 0; i < wymiar; i++)
		Global_H[i] = new long double[wymiar];

	for (int i = 0; i < wymiar; i++) {
		for (int j = 0; j < wymiar; j++) {
			Global_H[i][j] = 0;
		}
	}

	Global_C = new long double *[wymiar];    //deklaracja GlobalC
	for (int i = 0; i < wymiar; i++)
		Global_C[i] = new long double[wymiar];

	for (int i = 0; i < wymiar; i++) {
		for (int j = 0; j < wymiar; j++) {
			Global_C[i][j] = 0;
		}
	}

	Global_HB = new long double *[wymiar];    //deklaracja GlobalHB
	for (int i = 0; i < wymiar; i++)
		Global_HB[i] = new long double[wymiar];

	for (int i = 0; i < wymiar; i++) {
		for (int j = 0; j < wymiar; j++) {
			Global_HB[i][j] = 0;
		}
	}

	H_final = new long double *[wymiar];    //deklaracja H_final
	for (int i = 0; i < wymiar; i++)
		H_final[i] = new long double[wymiar];

	for (int i = 0; i < wymiar; i++) {
		for (int j = 0; j < wymiar; j++) {
			H_final[i][j] = 0;
		}
	}


	//deklaracja tablicy z temperaturami
	temp = new long double[wymiar];
	for (int i = 0; i <wymiar; i++) {  
		temp[i] = T0;
	}

	//zerowanie wektora P
	GlobalP = new long double[wymiar];
	for (int i = 0; i < wymiar; i++) {
		GlobalP[i] = 0;
	}

}


Data::~Data()
{
}
