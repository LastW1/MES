#include "Grid.h"
#include <iostream>
#include "Data.h"
#include "cstdlib"

using namespace std;


Grid::Grid()
{
	Data data;
	long double nodeWidth = data.L / (data.nL - 1);
	long double nodeHight = data.H / (data.nH - 1);
	long int id = 1;

	//tworzenie i nadawanie wartoœci nodom
	for (int i = 0; i < data.nL; i++) {
		for (int j = 0; j < data.nH; j++) {
			if (i == 0 || j == 0 || i == data.nL - 1 || j == data.nH - 1) {      //tu definiowanie brzegu dla warunku
				gridNodes.push_back(Node(i*nodeWidth, j*nodeHight, id, 1));
				id++;
			}
			else {
				gridNodes.push_back(Node(i*nodeWidth, j*nodeHight, id, 0));
				id++;
			}
		}
		
	}


	//tworzenie elementów i przypisywanie do nich nodów
	id = 1;
	for ( int i = 0; i < data.nL - 1; i++)
	{
		for ( int j = 0; j < data.nH - 1; j++)
		{

			gridElements.push_back(Element((j + data.nH * i)+1,(j + data.nH * (i + 1))+1,(j + 1 + data.nH * (i + 1))+1, (j + 1 + data.nH * i)+1, id));
			id++;
		}
	}
	for (int i = 0; i < (data.nH - 1)*(data.nL - 1); i++)
	//cout << gridNodes[gridElements[8].node[2]-1].id << endl;

	//Brzegi dla elementu
	for (int i = 0; i < (data.nL - 1)*(data.nH - 1); i++) {
		if (gridNodes[gridElements[i].node[0]-1].boundary * gridNodes[gridElements[i].node[3]-1].boundary) {         //brzeg 4
			gridElements[i].boundary[3] = 1;
		}
		else
			gridElements[i].boundary[3] = 0;

		if (gridNodes[gridElements[i].node[2]-1].boundary * gridNodes[gridElements[i].node[3]-1].boundary) {    //brzeg 3
			gridElements[i].boundary[2] = 1;
		}
		else
			gridElements[i].boundary[2] = 0;

		if (gridNodes[gridElements[i].node[1]-1].boundary * gridNodes[gridElements[i].node[2]-1].boundary) {   //brzeg 2
			gridElements[i].boundary[1] = 1;
		}
		else
			gridElements[i].boundary[1] = 0;

		if (gridNodes[gridElements[i].node[0]-1].boundary * gridNodes[gridElements[i].node[1]-1].boundary) {   //brzeg 1
			gridElements[i].boundary[0] = 1;
		}
		else
			gridElements[i].boundary[0] = 0;

	}


	//wyœwietlanie siatki
	for (int i = 0; i < (data.nH - 1)*(data.nL - 1); i++) {

		cout << "id elementu:" << gridElements[i].id 
			<<" brzeg1: "<<gridElements[i].boundary[0] 
			<< " brzeg2: " << gridElements[i].boundary[1]
			<< " brzeg3: " << gridElements[i].boundary[2]
			<< " brzeg4: " << gridElements[i].boundary[3]
			<<"        wezly: " << gridElements[i].node[0] <<  ", " << gridElements[i].node[1] << ", " << gridElements[i].node[2] << ", " << gridElements[i].node[3] << endl;
		
	}

    //element uniwersalny
		universal(&data);

		//cout << "\nksi";
		//for (int i = 0; i < 4; i++) {
		//	cout << endl;
		//	for (int j = 0; j < 4; j++) {
		//		cout << data.universal_ksi[i][j] << " ";
		//	}
		//}

		//cout << "\n\neta";
		//for (int i = 0; i < 4; i++) {
		//	cout << endl;
		//	for (int j = 0; j < 4; j++) {
		//		cout << data.universal_eta[i][j] << " ";
		//	}
		//}
		//cout << endl;
		//cout << endl;
	
		//macierz H
		matrixH(&data);
}

void Grid::show_coordinate(long int id) {
	cout << "polozenie noda o id:" << id  << " x:" << gridNodes[id-1].x << " y:" << gridNodes[id-1].y << endl; //id-1 poniewa¿ id liczone jest od 1, ale w rzeczywistoœci, w kodzie od 0



}

void Grid::universal(Data *data) {
	long double pc[4][2] = { { -1 / sqrt(3), -1 / sqrt(3)},            //wartoœci wspó³rzêdnych w lokalnym uk³adzie wspó³rzêdnych
						{ 1 / sqrt(3), -1 / sqrt(3)},
						{ 1 / sqrt(3), 1 / sqrt(3)},
						{ -1 / sqrt(3), 1 / sqrt(3)}};
	//cout << pc[0][0]<<" "<<pc[2][1] << endl;
	
	long double universal_eta[4][4];       //wartoœci N/eta i N/ksi
	long double universal_ksi[4][4];



	//cout << "\ntest punktow calkowania";
	//for (int i = 0; i < 4; i++) {
	//	cout << endl;
	//	for (int j = 0; j < 2; j++) {
	//		cout << pc[i][j] << " ";
	//	}
	//}
	cout << endl;

	for (int i = 0; i < 4; i++) {
		switch (i) {
		case 0:

			for (int j = 0; j < 4; j++) {
				universal_ksi[j][i] = -0.25*(1 - pc[j][1]);
				universal_eta[j][i] = -0.25*(1 - pc[j][0]);
			}

			break;
		case 1:

			for (int j = 0; j < 4; j++) {
				universal_ksi[j][i] = 0.25*(1 - pc[j][1]);
				universal_eta[j][i] = -0.25*(1 + pc[j][0]);
			}

			break;
		case 2:

			for (int j = 0; j < 4; j++) {
				universal_ksi[j][i] = 0.25*(1 + pc[j][1]);
				universal_eta[j][i] = 0.25*(1 + pc[j][0]);
			}

			break;
		case 3:

			for (int j = 0; j < 4; j++) {
				universal_ksi[j][i] = -0.25*(1 + pc[j][1]);
				universal_eta[j][i] = 0.25*(1 - pc[j][0]);
			}

			break;
		}

	}


	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			data->universal_eta[i][j] = universal_eta[i][j];
			data->universal_ksi[i][j] = universal_ksi[i][j];
		}
	}

}

void Grid::matrixH(Data* data) {

	
	long double jacobian[4];
	long double detJ;
	long double final_J[4][4];

	long double pc[4][2] = { { -1 / sqrt(3), -1 / sqrt(3)},            //wartoœci wspó³rzêdnych w lokalnym uk³adzie wspó³rzêdnych
					  	{ 1 / sqrt(3), -1 / sqrt(3)},
						{ 1 / sqrt(3), 1 / sqrt(3)},
						{ -1 / sqrt(3), 1 / sqrt(3)} };
	

	
	//pêtla po wszystkich elementach
	for (int i = 0; i < (data->nH - 1) * (data->nL - 1); i++) {

 //   long double x[4] = { gridElements[i].node[0].x,
	//					 gridElements[i].node[1].x,
	//					 gridElements[i].node[2].x, 
	//					 gridElements[i].node[3].x };

	//long double y[4] = { gridElements[i].node[0].y,
	//					 gridElements[i].node[1].y,
	//					 gridElements[i].node[2].y,
	//					 gridElements[i].node[3].y };

	long double x[4] = { gridNodes[gridElements[i].node[0]-1].x,
						 gridNodes[gridElements[i].node[1]-1].x,
						 gridNodes[gridElements[i].node[2]-1].x,
						 gridNodes[gridElements[i].node[3]-1].x };

	long double y[4] = { gridNodes[gridElements[i].node[0]-1].y,
						 gridNodes[gridElements[i].node[1]-1].y,
						 gridNodes[gridElements[i].node[2]-1].y,
						 gridNodes[gridElements[i].node[3]-1].y };



	long double pc_N_x[4][4];
	long double pc_N_y[4][4];
	long double tmp_H_x[4][4];
	long double tmp_H_y[4][4];

	//wyœwietlanie wspó³rzêdnych
	//cout << "\n\n\n\n\nwspolrzedne";
	//cout << "\nx: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
	//cout << "y: " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << "\n" << endl; 

	//pêtla dla nodów
	for (int j = 0; j < 4; j++) {
		long double x_ksi = x[0] * data->universal_ksi[j][0] + x[1] * data->universal_ksi[j][1] + x[2] * data->universal_ksi[j][2] + x[3] * data->universal_ksi[j][3];
		long double x_eta = x[0] * data->universal_eta[j][0] + x[1] * data->universal_eta[j][1] + x[2] * data->universal_eta[j][2] + x[3] * data->universal_eta[j][3];
		long double y_ksi = y[0] * data->universal_ksi[j][0] + y[1] * data->universal_ksi[j][1] + y[2] * data->universal_ksi[j][2] + y[3] * data->universal_ksi[j][3];
		long double y_eta = y[0] * data->universal_eta[j][0] + y[1] * data->universal_eta[j][1] + y[2] * data->universal_eta[j][2] + y[3] * data->universal_eta[j][3];

		if (x_ksi < 0.000001)
			x_ksi = 0;
		if (x_eta < 0.000001)
			x_eta = 0;
		if (y_ksi < 0.000001)
			y_ksi = 0;
		if (y_eta < 0.000001)
			y_eta = 0;


		detJ = (x_ksi * y_eta) - (y_ksi * x_eta);

		//cout << x_ksi <<" "<< x_eta <<" " << y_ksi <<" "<< y_eta << endl;
		//cout << detJ << endl;
		final_J[3][j] = x_ksi/detJ;
		final_J[2][j] = x_eta/detJ;
		final_J[1][j] = y_ksi/detJ;
		final_J[0][j] = y_eta/detJ;

		//N/x i N/y w 4 pkt ca³kowania j odpowiada za pkt ca³kowania, k za funkcjê kszta³tu
		for (int k = 0; k<4 ; k++){ 
			pc_N_x[j][k] = 1/detJ *(y_eta * data->universal_ksi[j][k]) + (-y_ksi * data->universal_ksi[j][k]);  
			pc_N_y[j][k] = 1 / detJ * (x_ksi * data->universal_eta[j][k]) + (-x_eta * data->universal_eta[j][k]);
		}

	}
	/*
	cout << "pc N/x i N/y" << endl;
	for (int x = 0; x < 4; x++) {
		cout << endl;
		for (int y = 0; y < 4; y++)
			cout << pc_N_x[x][y] << " ";
	}
	cout << endl;

	for (int x = 0; x < 4; x++) {
		cout << endl;
		for (int y = 0; y < 4; y++)
			cout << pc_N_y[x][y] << " ";
	}
	*/

	/*cout << "Jacobian" << endl;
	for (int x = 0; x < 4; x++) {
		cout << endl;
		for (int y = 0; y < 4; y++)
			cout << final_J[x][y] << " ";
	}
		cout << endl;*/

	//namna¿anie N/x i N/y przez siebie 
	for(int j = 0 ; j<4 ; j++)
		for (int k = 0; k < 4; k++) {

			tmp_H_x[j][k] = pc_N_x[0][k] * pc_N_x[0][j] + pc_N_x[1][k] * pc_N_x[1][j] + pc_N_x[2][k] * pc_N_x[2][j] + pc_N_x[3][k] * pc_N_x[3][j];
			tmp_H_y[j][k] = pc_N_y[0][k] * pc_N_y[0][j] + pc_N_y[1][k] * pc_N_y[1][j] + pc_N_y[2][k] * pc_N_y[2][j] + pc_N_y[3][k] * pc_N_y[3][j];
		}

	//przemno¿enie macierzy przez wsp K oraz detJ
	for (int j = 0; j < 4; j++)
		for (int k = 0; k < 4; k++) {
			tmp_H_x[j][k] *= data->K* detJ;
			tmp_H_y[j][k] *= data->K* detJ;

		}

	//sumowani macierzy x oraz y do wspólnej macierzy H
	for (int j = 0; j < 4; j++)
		for (int k = 0; k < 4; k++)
			gridElements[i].matrixH[j][k] = tmp_H_x[j][k] + tmp_H_y[j][k];


	//wyœwitlanie macierzy H lokalnej
	//cout << "macierz H" << endl;
	//for (int x = 0; x < 4; x++) {
	//	cout << endl;
	//	for (int y = 0; y < 4; y++)
	//		cout << gridElements[i].matrixH[x][y] << " ";
	//}
	//agregacja do global H
	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++) {
			data->Global_H[gridElements[i].node[x] - 1][gridElements[i].node[y] - 1] += gridElements[i].matrixH[x][y];
		}

	//CZAS NA TABLICÊ C!!!!!!!!!!

	//funkcje kszta³tu
	long double tabN[4][4];
	long double tabC[4][4];

	for (int j = 0; j < 4; j++) {
		switch (j) {
		case 0:

			for (int k = 0; k < 4; k++) {
				tabN[k][j] = 0.25*((1 - pc[k][0])*(1 - pc[k][1]));

			}

			break;
		case 1:

			for (int k = 0; k < 4; k++) {
				tabN[k][j] = 0.25*((1 + pc[k][0])*(1 - pc[k][1]));

			}

			break;
		case 2:

			for (int k = 0; k < 4; k++) {
				tabN[k][j] = 0.25*((1 + pc[k][0])*(1 + pc[k][1]));

			}

			break;
		case 3:

			for (int k = 0; k < 4; k++) {
				tabN[k][j] = 0.25*((1 - pc[k][0])*(1 + pc[k][1]));

			}

			break;
		}

	}
	/*cout << "Tablica funkcji kszta³tow" << endl;
	for (int x = 0; x < 4; x++) {
		cout << endl;
		for (int y = 0; y < 4; y++)
			cout << tabN[x][y] << " ";
	}*/


	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++) {
			tabC[x][y] = ((tabN[0][x] * tabN[0][y]) + (tabN[1][x] * tabN[1][y]) + (tabN[2][x] * tabN[2][y]) + (tabN[3][x] * tabN[3][y])) * detJ * data->C * data->RO;
		}

	/*cout << endl;
	for (int x = 0; x < 4; x++) {
		cout << endl;
		for (int y = 0; y < 4; y++) {
			cout << tabC[x][y] << " ";
		}
	}*/
	//agregacja do global C
	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++) {
			data->Global_C[gridElements[i].node[x] - 1][gridElements[i].node[y] - 1] += tabC[x][y];
		}


	//obliczanie H z warunkami brzegowymi

	//zerowanie lokalnej tablicy HB
	long double matrixHB[4][4];
	long double d[4]; //d³ugoœæ boku
	d[0]= abs(gridNodes[gridElements[i].node[0] - 1].x - gridNodes[gridElements[i].node[1] - 1].x);
	d[1]= abs(gridNodes[gridElements[i].node[1] - 1].y - gridNodes[gridElements[i].node[2] - 1].y);
	d[2]= abs(gridNodes[gridElements[i].node[2] - 1].x - gridNodes[gridElements[i].node[3] - 1].x);
	d[3]= abs(gridNodes[gridElements[i].node[0] - 1].y - gridNodes[gridElements[i].node[3] - 1].y);
	
	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++)
			matrixHB[x][y] = 0;
	//wektor P
	long double tabP[4] = { 0,0,0,0 };

	//sprawdzenie czy ten brzeg posiada warunek
	if (gridElements[i].boundary[0]) {
		calculate_HB(matrixHB,-0.57735, 0.57735,-1,-1,data->A, d[0],tabP);
	}
	if (gridElements[i].boundary[1]) {
		calculate_HB(matrixHB, 1, 1, -0.57735, 0.57735, data->A, d[1],tabP);
	}
	if (gridElements[i].boundary[2]) {
		calculate_HB(matrixHB, 0.57735, -0.57735, 1, 1, data->A, d[2],tabP);
	}
	if (gridElements[i].boundary[3]) {
		calculate_HB(matrixHB, -1, -1, 0.57735, -0.57735, data->A, d[3],tabP);
	}

	//for (int x = 0; x < 4; x++) {
	//	cout << endl;
	//	for (int y = 0; y < 4; y++)
	//		cout << matrixHB[x][y] << " ";
	//}

   //agregacja do global HB
	for (int x = 0; x < 4; x++)
		for (int y = 0; y < 4; y++) {
			data->Global_HB[gridElements[i].node[x] - 1][gridElements[i].node[y] - 1] += matrixHB[x][y];
		}



	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			switch (j) {
			case 0: {
				tabP[i] += (tabC[i][j] / data->dT)*100;
				break;
			}
			case 1: {
				tabP[i] += (tabC[i][j] / data->dT)*100;
				break;
			}
			case 2: {
				tabP[i] += (tabC[i][j] / data->dT)*100;
				break;
			}
			case 3: {
				tabP[i] += (tabC[i][j] / data->dT)*100;
				break;
			}
			}
		}


	//agregacja do globalnego P
	for (int x = 0; x < 4; x++)
			data->GlobalP[gridElements[i].node[x] - 1] += tabP[x];

	




	}
	

	show(data);


}


void Grid::show(Data*data) {

	cout << "macierz globalna H" << endl;
	for (int i = 0; i < data->nH  * data->nL; i++) {
		cout << endl;
		for (int j = 0; j < data->nH * data->nL; j++)
			cout << data->Global_H[i][j] << " ";

	}

	cout << "\n\nmacierz globalna HB" << endl;
	for (int i = 0; i < data->nH * data->nL; i++) {
		cout << endl;
		for (int j = 0; j < data->nH * data->nL; j++)
			cout << data->Global_HB[i][j] << " ";

	}

	cout << "\n\nmacierz globalna C" << endl;
	for (int i = 0; i < data->nH * data->nL; i++) {
		cout << endl;
		for (int j = 0; j < data->nH * data->nL; j++)
			cout << data->Global_C[i][j] << " ";

	}

	for (int i = 0; i < data->nH * data->nL; i++)
		for (int j = 0; j < data->nH * data->nL; j++)
			data->H_final[i][j] = data->Global_H[i][j] + data->Global_HB[i][j] + (data->Global_C[i][j] / data->dT);


	cout << "\n\nmacierz H (H+C/dt)" << endl;
	for (int i = 0; i < data->nH * data->nL; i++) {
		cout << endl;
		for (int j = 0; j < data->nH * data->nL; j++)
			cout << data->H_final[i][j] << " ";

	}

	cout << "\n\nwektor P" << endl;
	for (int i = 0; i < data->nH * data->nL; i++)
			cout << data->GlobalP[i] << " ";



}


void Grid::calculateP(){

}


void Grid::calculate_HB(long double tabHB[][4],long double ksi1, long double ksi2, long double eta1, long double eta2,long double a,long double d,long double p[]) {
	

	long double tab[4][2];

				tab[0][0] = 0.25*((1 - ksi1)*(1 - eta1));
				tab[0][1] = 0.25*((1 - ksi2)*(1 - eta2));

				tab[1][0] = 0.25*((1 + ksi1)*(1 - eta1));
				tab[1][1] = 0.25*((1 + ksi2)*(1 - eta2));

				tab[2][0] = 0.25*((1 + ksi1)*(1 + eta1));
				tab[2][1] = 0.25*((1 + ksi2)*(1 + eta2));

				tab[3][0] = 0.25*((1 - ksi1)*(1 + eta1));
				tab[3][1] = 0.25*((1 - ksi2)*(1 + eta2));



	//cout << "Macierz pcHB" << endl;
	//for (int x = 0; x < 4; x++) {
	//	cout << endl;
	//	for (int y = 0; y < 2; y++)
	//		cout << tab[x][y] << " ";
	//}


				for(int x = 0 ; x<4 ; x++){
					p[x] += (tab[x][0] + tab[x][1])*a* (d / 2) * 100;
					cout << p[x] << endl;
				}


				for (int i = 0; i < 2; i++) {
					for (int x = 0; x < 4; x++)
						for (int y = 0; y < 4; y++) {
							tabHB[x][y] += tab[x][i] * tab[y][i] * a* (d / 2);
						}
				}


}

Grid::~Grid()
{
}
