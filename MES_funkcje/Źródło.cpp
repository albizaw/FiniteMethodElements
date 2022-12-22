#include<iostream>
#include <string>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "funkcjeKsztaltow.h"
#include "macierzH.h"

using namespace std;

struct node {
	double x;
	double y;
	double t;
	bool BC = 0;
};

struct element {
	int ID[4];
};

struct grid {
	int nN;
	int nE;
	vector<node> ND;
	vector<element>EL;
};

struct GlobalData {
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
};

struct kwadratury {
	double PC2p[2] = {-sqrt(1/3.), sqrt(1/3.)};
	double W2p[2] = { 1,1 };

	double PC3p[3] = { -sqrt(3/5.),0,sqrt(3/5.) };
	double W3p[3] = { 5 / 9.,8 / 9.,5 / 9. };

	double PC4p[4] = { -(sqrt((3/7.)+((2/7.)*(sqrt(6/5.))))),-(sqrt((3 / 7.) - ((2 / 7.) * (sqrt(6 / 5.))))), (sqrt((3 / 7.) - ((2 / 7.) * (sqrt(6 / 5.))))), (sqrt((3 / 7.) + ((2 / 7.) * (sqrt(6 / 5.))))) };
	double W4p[4] = { (18-sqrt(30))/36., (18 + sqrt(30)) / 36., (18 + sqrt(30)) / 36., (18 - sqrt(30)) / 36. };
};

//przekazujemy flagê, która tworzy odpowiedni wymiar tablicy
struct Elem4 {
	int numberOfPoints;
	int punktyNaPowierzchni;

	double** tabKsi;
	double** tabEta;

	double** tabNPowierzchnia;
	double** tabFunkcjiKsztaltow;
	
	Elem4(int numberOfPoints, int punktyNaPowierzchni) : numberOfPoints(numberOfPoints),punktyNaPowierzchni(punktyNaPowierzchni)
	{
		
		int k = numberOfPoints * numberOfPoints;
		tabKsi = new double* [k];
		tabEta = new double* [k];
		tabFunkcjiKsztaltow = new double* [k];

		tabNPowierzchnia = new double* [punktyNaPowierzchni * 4]; //bo 4 œciany
		for (int i = 0; i < k; i++)
		{
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
			tabFunkcjiKsztaltow[i] = new double[4];
		}

		for (int i = 0; i < punktyNaPowierzchni * 4; i++)
		{
			tabNPowierzchnia[i] = new double[4];
		}
	}
	
};

void elem4(kwadratury& kwadratury, Elem4& elem) {

	int counter = 0;
	double eta;
	double ksi;

	if (elem.numberOfPoints == 2) {
		for (int i = 0; i < elem.numberOfPoints * elem.numberOfPoints; i++)
		{
			if (i % 2 == 0)
			{
				counter = 0;
			}

			ksi = kwadratury.PC2p[counter];
			eta = kwadratury.PC2p[i / 2];

			//cout << ksi << " " << eta << endl;

			elem.tabKsi[i][0] = N1ksi(eta);
			elem.tabKsi[i][1] = N2ksi(eta);
			elem.tabKsi[i][2] = N3ksi(eta);
			elem.tabKsi[i][3] = N4ksi(eta);

			elem.tabEta[i][0] = N1eta(kwadratury.PC2p[counter]);
			elem.tabEta[i][1] = N2eta(kwadratury.PC2p[counter]);
			elem.tabEta[i][2] = N3eta(kwadratury.PC2p[counter]);
			elem.tabEta[i][3] = N4eta(kwadratury.PC2p[counter]);

			elem.tabFunkcjiKsztaltow[i][0] = N1ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][1] = N2ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][2] = N3ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][3] = N4ksieta(ksi, eta);

			cout << N1ksieta(ksi, eta) << " " << N2ksieta(ksi, eta) << " " << N3ksieta(ksi, eta) << " " << N4ksieta(ksi, eta) << " " << endl;
			counter++;
		}
		cout << endl;
	}

	if (elem.numberOfPoints == 3)
	{
		for (int i = 0; i < elem.numberOfPoints * elem.numberOfPoints; i++)
		{
			if (i % 3 == 0)
			{
				counter = 0;
			}

			ksi = kwadratury.PC3p[counter];
			eta = kwadratury.PC3p[i / 3];


			elem.tabKsi[i][0] = N1ksi(eta);
			elem.tabKsi[i][1] = N2ksi(eta);
			elem.tabKsi[i][2] = N3ksi(eta);
			elem.tabKsi[i][3] = N4ksi(eta);

			elem.tabEta[i][0] = N1eta(ksi);
			elem.tabEta[i][1] = N2eta(ksi);
			elem.tabEta[i][2] = N3eta(ksi);
			elem.tabEta[i][3] = N4eta(ksi);

			elem.tabFunkcjiKsztaltow[i][0] = N1ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][1] = N2ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][2] = N3ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][3] = N4ksieta(ksi, eta);

			counter++;

		}
	}

	if (elem.numberOfPoints == 4) {
		for (int i = 0; i < elem.numberOfPoints * elem.numberOfPoints; i++)
		{

			if (i % 4 == 0)
			{
				counter = 0;
			}

			ksi = kwadratury.PC4p[counter];
			eta = kwadratury.PC4p[i / 4];

			elem.tabKsi[i][0] = N1ksi(eta);
			elem.tabKsi[i][1] = N2ksi(eta);
			elem.tabKsi[i][2] = N3ksi(eta);
			elem.tabKsi[i][3] = N4ksi(eta);

			elem.tabEta[i][0] = N1eta(ksi);
			elem.tabEta[i][1] = N2eta(ksi);
			elem.tabEta[i][2] = N3eta(ksi);
			elem.tabEta[i][3] = N4eta(ksi);

			elem.tabFunkcjiKsztaltow[i][0] = N1ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][1] = N2ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][2] = N3ksieta(ksi, eta);
			elem.tabFunkcjiKsztaltow[i][3] = N4ksieta(ksi, eta);

			counter++;
		}
	}


	
	if (elem.punktyNaPowierzchni == 2)
	{
		//dolnasciana
		ksi = kwadratury.PC2p[0];
			elem.tabNPowierzchnia[0][0] = N1ksieta(ksi,-1);
			elem.tabNPowierzchnia[0][1] = N2ksieta(ksi, -1);
			elem.tabNPowierzchnia[0][2] = N3ksieta(ksi, -1);
			elem.tabNPowierzchnia[0][3] = N4ksieta(ksi, -1);

			ksi = kwadratury.PC2p[1];
			elem.tabNPowierzchnia[1][0] = N1ksieta(ksi, -1);
			elem.tabNPowierzchnia[1][1]= N2ksieta(ksi, -1);
			elem.tabNPowierzchnia[1][2] = N3ksieta(ksi, -1);
			elem.tabNPowierzchnia[1][3] = N4ksieta(ksi, -1);

			//prawa sciana
			ksi = kwadratury.PC2p[0];
			elem.tabNPowierzchnia[2][0] = N1ksieta(1, ksi);
			elem.tabNPowierzchnia[2][1] = N2ksieta(1, ksi);
			elem.tabNPowierzchnia[2][2] = N3ksieta(1, ksi);
			elem.tabNPowierzchnia[2][3] = N4ksieta(1, ksi);

			ksi = kwadratury.PC2p[1];
			elem.tabNPowierzchnia[3][0] = N1ksieta(1, ksi);
			elem.tabNPowierzchnia[3][1] = N2ksieta(1, ksi);
			elem.tabNPowierzchnia[3][2] = N3ksieta(1, ksi);
			elem.tabNPowierzchnia[3][3] = N4ksieta(1, ksi);

			//gorna sciana
			ksi = kwadratury.PC2p[0];
			elem.tabNPowierzchnia[4][0] = N1ksieta(ksi,1);
			elem.tabNPowierzchnia[4][1] = N2ksieta(ksi,1);
			elem.tabNPowierzchnia[4][2] = N3ksieta(ksi,1);
			elem.tabNPowierzchnia[4][3] = N4ksieta(ksi,1);

			ksi = kwadratury.PC2p[1];
			elem.tabNPowierzchnia[5][0] = N1ksieta(ksi, 1);
			elem.tabNPowierzchnia[5][1] = N2ksieta(ksi, 1);
			elem.tabNPowierzchnia[5][2] = N3ksieta(ksi, 1);
			elem.tabNPowierzchnia[5][3] = N4ksieta(ksi, 1);

			//lewa
			ksi = kwadratury.PC2p[0];
			elem.tabNPowierzchnia[6][0] = N1ksieta(-1, ksi);
			elem.tabNPowierzchnia[6][1] = N2ksieta(-1, ksi);
			elem.tabNPowierzchnia[6][2] = N3ksieta(-1, ksi);
			elem.tabNPowierzchnia[6][3] = N4ksieta(-1, ksi);

			ksi = kwadratury.PC2p[1];
			elem.tabNPowierzchnia[7][0] = N1ksieta(-1, ksi);
			elem.tabNPowierzchnia[7][1] = N2ksieta(-1, ksi);
			elem.tabNPowierzchnia[7][2] = N3ksieta(-1, ksi);
			elem.tabNPowierzchnia[7][3] = N4ksieta(-1, ksi);
	}

	if (elem.punktyNaPowierzchni == 3)
	{
		//dolna sciana elementu
		//1p
		ksi = kwadratury.PC3p[0];
		elem.tabNPowierzchnia[0][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[0][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[0][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[0][3] = N4ksieta(ksi, -1);
		//2p
		ksi = kwadratury.PC3p[1];
		elem.tabNPowierzchnia[1][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[1][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[1][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[1][3] = N4ksieta(ksi, -1);
		//3p
		ksi = kwadratury.PC3p[2];
		elem.tabNPowierzchnia[2][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[2][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[2][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[2][3] = N4ksieta(ksi, -1);

		//prawa sciana elementu
		//1p
		ksi = kwadratury.PC3p[0];
		elem.tabNPowierzchnia[3][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[3][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[3][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[3][3] = N4ksieta(1, ksi);
		//2p
		ksi = kwadratury.PC3p[1];
		elem.tabNPowierzchnia[4][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[4][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[4][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[4][3] = N4ksieta(1, ksi);
		//3p
		ksi = kwadratury.PC3p[2];
		elem.tabNPowierzchnia[5][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[5][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[5][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[5][3] = N4ksieta(1, ksi);

		//gorna sciana elementu
		//1p
		ksi = kwadratury.PC3p[0];
		elem.tabNPowierzchnia[6][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[6][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[6][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[6][3] = N4ksieta(ksi, 1);
		//2p
		ksi = kwadratury.PC3p[1];
		elem.tabNPowierzchnia[7][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[7][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[7][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[7][3] = N4ksieta(ksi, 1);
		//3p
		ksi = kwadratury.PC3p[2];
		elem.tabNPowierzchnia[8][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[8][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[8][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[8][3] = N4ksieta(ksi, 1);

		//lewa sciana elementu
		//1p
		ksi = kwadratury.PC3p[0];
		elem.tabNPowierzchnia[9][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[9][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[9][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[9][3] = N4ksieta(-1, ksi);
		//2p
		ksi = kwadratury.PC3p[1];
		elem.tabNPowierzchnia[10][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[10][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[10][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[10][3] = N4ksieta(-1, ksi);
		//3p
		ksi = kwadratury.PC3p[2];
		elem.tabNPowierzchnia[11][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[11][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[11][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[11][3] = N4ksieta(-1, ksi);

	}

	if (elem.punktyNaPowierzchni == 4)
	{
		//dolna sciana elementu
		//1p
		ksi = kwadratury.PC4p[0];
		elem.tabNPowierzchnia[0][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[0][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[0][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[0][3] = N4ksieta(ksi, -1);
		//2p
		ksi = kwadratury.PC4p[1];
		elem.tabNPowierzchnia[1][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[1][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[1][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[1][3] = N4ksieta(ksi, -1);
		//3p
		ksi = kwadratury.PC4p[2];
		elem.tabNPowierzchnia[2][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[2][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[2][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[2][3] = N4ksieta(ksi, -1);

		//4p
		ksi = kwadratury.PC4p[3];
		elem.tabNPowierzchnia[3][0] = N1ksieta(ksi, -1);
		elem.tabNPowierzchnia[3][1] = N2ksieta(ksi, -1);
		elem.tabNPowierzchnia[3][2] = N3ksieta(ksi, -1);
		elem.tabNPowierzchnia[3][3] = N4ksieta(ksi, -1);

		//prawa sciana elementu
		//1p
		ksi = kwadratury.PC4p[0];
		elem.tabNPowierzchnia[4][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[4][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[4][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[4][3] = N4ksieta(1, ksi);
		//2p
		ksi = kwadratury.PC4p[1];
		elem.tabNPowierzchnia[5][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[5][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[5][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[5][3] = N4ksieta(1, ksi);
		//3p
		ksi = kwadratury.PC4p[2];
		elem.tabNPowierzchnia[6][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[6][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[6][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[6][3] = N4ksieta(1, ksi);
		//4p
		ksi = kwadratury.PC4p[3];
		elem.tabNPowierzchnia[7][0] = N1ksieta(1, ksi);
		elem.tabNPowierzchnia[7][1] = N2ksieta(1, ksi);
		elem.tabNPowierzchnia[7][2] = N3ksieta(1, ksi);
		elem.tabNPowierzchnia[7][3] = N4ksieta(1, ksi);

		//gorna sciana elementu
		//1p
		ksi = kwadratury.PC4p[0];
		elem.tabNPowierzchnia[8][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[8][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[8][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[8][3] = N4ksieta(ksi, 1);
		//2p
		ksi = kwadratury.PC4p[1];
		elem.tabNPowierzchnia[9][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[9][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[9][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[9][3] = N4ksieta(ksi, 1);
		//3p
		ksi = kwadratury.PC4p[2];
		elem.tabNPowierzchnia[10][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[10][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[10][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[10][3] = N4ksieta(ksi, 1);
		//4p
		ksi = kwadratury.PC4p[3];
		elem.tabNPowierzchnia[11][0] = N1ksieta(ksi, 1);
		elem.tabNPowierzchnia[11][1] = N2ksieta(ksi, 1);
		elem.tabNPowierzchnia[11][2] = N3ksieta(ksi, 1);
		elem.tabNPowierzchnia[11][3] = N4ksieta(ksi, 1);

		//lewa sciana elementu
		//1p
		ksi = kwadratury.PC4p[0];
		elem.tabNPowierzchnia[12][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[12][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[12][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[12][3] = N4ksieta(-1, ksi);
		//2p
		ksi = kwadratury.PC4p[1];
		elem.tabNPowierzchnia[13][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[13][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[13][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[13][3] = N4ksieta(-1, ksi);
		//3p
		ksi = kwadratury.PC4p[2];
		elem.tabNPowierzchnia[14][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[14][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[14][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[14][3] = N4ksieta(-1, ksi);

		//4p
		ksi = kwadratury.PC4p[3];
		elem.tabNPowierzchnia[15][0] = N1ksieta(-1, ksi);
		elem.tabNPowierzchnia[15][1] = N2ksieta(-1, ksi);
		elem.tabNPowierzchnia[15][2] = N3ksieta(-1, ksi);
		elem.tabNPowierzchnia[15][3] = N4ksieta(-1, ksi);

	}

}

double **JacobiMatrix(Elem4& elem, grid& newGrid, kwadratury& kwadratura, int _nrEl, GlobalData& globe)
{
	vector <double> x;
	vector <double> y;

	vector <int> BC;

	int nrEl = _nrEl;
	
	for (int i = 0; i < 4; i++)
	{
		int k = newGrid.EL[nrEl].ID[i];
		
		x.push_back(newGrid.ND[k-1].x);
		y.push_back(newGrid.ND[k - 1].y);
		
	}

	int liczbaScianBC = 0;
	for (int i = 0; i < 4; i++)
	{
		int k = newGrid.EL[nrEl].ID[i];
		BC.push_back(newGrid.ND[k - 1].BC);

	}

	vector <int> ktoraSciana;
	int cnt = 0;
	for (int i = 0; i < 3; i++)
	{
		if (BC[i] == 1 && BC[i + 1] == 1)
		{
			liczbaScianBC++;
			ktoraSciana.push_back(cnt);
		}
		cnt++;
	}
	
	if (BC[3]==1 && BC[0]==1)
	{
		liczbaScianBC++;
		ktoraSciana.push_back(cnt);
	}

	//double Hb[8][4][4]{};
	int liczba = elem.punktyNaPowierzchni;
	int pkt = liczba * 4;


	double*** Hbnew;
	Hbnew = new double**[pkt];
	for (size_t i = 0; i < pkt; i++)
	{
		Hbnew[i] = new double* [4];
		for (size_t j = 0; j < 4; j++)
		{
			Hbnew[i][j]=new double[4];
		}
	}

	//cout << "Liczba scian z warunkiem: " << liczbaScianBC << endl;
	//tutaj bede przekazywal do funkcji liczbe punktow na powierzchni
	int counterWagi;
	double waga;
	double cond = globe.Alfa;
	//double cond = 25;
	for (int i = 0; i < pkt; i++) //liczba bokow
	{
		
		if (i % elem.punktyNaPowierzchni == 0)
		{
			counterWagi = 0;
		}

		if (elem.punktyNaPowierzchni == 2)
		{
			waga = kwadratura.W2p[counterWagi];
		}
		else if (elem.punktyNaPowierzchni == 3)
		{
			waga = kwadratura.W3p[counterWagi];
		}
		else if (elem.punktyNaPowierzchni == 4)
		{
			waga = kwadratura.W4p[counterWagi];
		}


		for (int j = 0; j <4;j++)
		{
			for (int k = 0; k < 4; k++)
			{
				double wyn = cond * elem.tabNPowierzchnia[i][j] * elem.tabNPowierzchnia[i][k];
				wyn = wyn * waga;
				Hbnew[i][j][k] = wyn;

			}

		}

		counterWagi++;
		
	}

	//musimy utworzyæ tablice Hbc tylko na tych œcianach gdzie jest warunek brzegowy
	//double Hbc[2][4][4];
	double*** Hbc;
	
		Hbc = new double** [liczbaScianBC];
		for (size_t i = 0; i < liczbaScianBC; i++)
		{
			Hbc[i] = new double* [4];
			for (size_t j = 0; j < 4; j++)
			{
				Hbc[i][j] = new double[4];
			}
		}

		//dodawanie do siebie na danym punktow na bokach
		double detJBok;
		for (int i = 0; i < liczbaScianBC; i++)
		{

			

			if (ktoraSciana[i] + 1 == 4)
			{
				detJBok = sqrt((x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0]));
			}
			else
			{
				detJBok = sqrt((x[ktoraSciana[i] + 1] - x[ktoraSciana[i]]) * (x[ktoraSciana[i] + 1] - x[ktoraSciana[i]]) + (y[ktoraSciana[i] + 1] - y[ktoraSciana[i]]) * (y[ktoraSciana[i] + 1] - y[ktoraSciana[i]]));
			};

			detJBok = detJBok / 2;
		
			for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 4; k++)
						{
							//Hbc[i][j][k] = (Hbnew[ktoraSciana[i]*2][j][k] + Hbnew[ktoraSciana[i]*2 + 1][j][k]) * 0.0125;
							if (elem.punktyNaPowierzchni == 2)
							{
								Hbc[i][j][k] = (Hbnew[ktoraSciana[i] * 2][j][k] + Hbnew[ktoraSciana[i] * 2 + 1][j][k]) * detJBok;
							}
							else if (elem.punktyNaPowierzchni == 3)
							{
								Hbc[i][j][k] = (Hbnew[ktoraSciana[i] * 3][j][k] + Hbnew[ktoraSciana[i] * 3 + 1][j][k] + Hbnew[ktoraSciana[i] * 3 + 2][j][k]) * detJBok;
							}
							else if (elem.punktyNaPowierzchni == 4)
							{
								Hbc[i][j][k] = (Hbnew[ktoraSciana[i] * 4][j][k] + Hbnew[ktoraSciana[i] * 4 + 1][j][k] + Hbnew[ktoraSciana[i] * 4 + 2][j][k] + +Hbnew[ktoraSciana[i] * 4 + 3][j][k]) * detJBok;
							}
							
						}

					}
		}

	int k = elem.numberOfPoints;
	int k2 = k * k;
	
	double** allJacobis;
	double** Nidx;
	double** Nidy;
	double** H;

		allJacobis = new double* [k2];
		Nidx = new double* [k2];
		Nidy = new double* [k2];
		H = new double* [k2];
		for (int i = 0; i < k2; i++)
		{
			allJacobis[i] = new double[4];
			Nidx[i] = new double[4];
			Nidy[i] = new double[4];
			H[i] = new double[4];
			
		}

		double*** Hpci;
		double*** Hnx;
		double*** Hny;
		

		Hnx = new double** [k2];
		Hny = new double** [k2];
		Hpci = new double** [k2];
		

		for (int j = 0; j < k2; j++)
		{
			Hnx[j] = new double* [k2];
			Hny[j] = new double* [k2];
			Hpci[j] = new double* [k2];

			for (int k = 0; k < k2; k++)
			{
				Hnx[j][k] = new double[4];
				Hny[j][k] = new double[4];
				Hpci[j][k] = new double[4];
			}
		}

	
		int counter;
		double eta;
		double ksi;

	for (int i = 0; i < k2; i++)
	{

		if (i % k == 0)
		{
			counter = 0;
		}

		if (k == 2)
		{
			ksi = kwadratura.W2p[counter];

			eta = kwadratura.W2p[i / k];
		}
		else if (k == 3)
		{
			ksi = kwadratura.W3p[counter];

			eta = kwadratura.W3p[i / k];
		}
		else if (k == 4)
		{
			ksi = kwadratura.W4p[counter];
			eta = kwadratura.W4p[i / k];
		}

		

		//double JacobiMat[2][2];
		double sumxksi = 0.0;
		double sumyksi = 0.0;
		double sumxeta = 0.0;
		double sumyeta = 0.0;

		for (int j = 0; j < 4; j++)
		{
			sumxksi += elem.tabKsi[i][j] * x[j] / 1.0;
			sumyksi += elem.tabKsi[i][j] * y[j] / 1.0;
			sumxeta += elem.tabEta[i][j] * x[j] / 1.0;
			sumyeta += elem.tabEta[i][j] * y[j] / 1.0;

		}

		allJacobis[i][0] = sumxksi / 1.0;
		allJacobis[i][1] = sumyksi / 1.0;
		allJacobis[i][2] = sumxeta / 1.0;
		allJacobis[i][3] = sumyeta / 1.0;


		//1) dNi / dx
		double detJ = (allJacobis[i][0] * allJacobis[i][3] - allJacobis[i][1] * allJacobis[i][2]) / 1.0;
		//cout << "DetJ : pkt nr:" << i << " === " << detJ << " | " << endl;

		//odwrócenie macierzy
		double pom;
		pom = allJacobis[i][0];
		allJacobis[i][0] = allJacobis[i][3];
		allJacobis[i][3] = pom;
		allJacobis[i][1] = -1 *allJacobis[i][1];
		allJacobis[i][2] = -1 * allJacobis[i][2];

		//wyznacznik 1/det
		double reverseDetJ = 1.0 / (detJ);
		for (int k = 0; k < 4; k++)
		{
			allJacobis[i][k] *= reverseDetJ;
		}

		//slajd13 macierz tabeladwuwymiarowazmacierza
		//Ni/dx macierz

		for (int j = 0; j < 4; j++)
		{
			Nidx[i][j] = allJacobis[i][0] * elem.tabKsi[i][j] + allJacobis[i][1] * elem.tabEta[i][j];
			Nidy[i][j] = allJacobis[i][2] * elem.tabKsi[i][j] + allJacobis[i][3] * elem.tabEta[i][j];
			//cout << Nidy[i][j] << " ";
		}
		//cout << endl;
		
		//macierzH dla punktów ca³kowania
		double kt = globe.Conductivity;
		double den = globe.Density;
		double spec = globe.SpecificHeat;
		
		//cout << "WagaKsi: " << ksi << "   WagaEta: " << eta << endl;

		for (int j = 0; j < 4; j++)
		{

			for (int k = 0; k < 4; k++)
			{
				Hnx[i][j][k] = Nidx[i][j] * Nidx[i][k];
				Hny[i][j][k] = Nidy[i][j] * Nidy[i][k];
				Hpci[i][j][k] = kt * (Hnx[i][j][k] + Hny[i][j][k]) * detJ;
				Hpci[i][j][k] *= ksi * eta;
			}

		}
		counter++;

	}

	//dodawanie macierzy H; zmiany wprowadzic w schemacie punktow
	double suma;
	//i - kolumny
	//j - wiersz
	for (int i = 0; i < 4; i++)
	{
		suma = 0;

		for (int j = 0; j < 4; j++)
		{
			suma = Hpci[0][i][j];
		
			
			for (int k = 1;  k < k2;  k++)
			{
				suma += Hpci[k][i][j];
				
			}

			H[i][j] = suma;

		}
	}

	// bez Hbc
	if (liczbaScianBC > 0)
	{
		for (int i = 0; i < 4; i++)
		{
			suma = 0;
			for (size_t j = 0; j < 4; j++)
			{
				suma = Hbc[0][i][j];

				for (size_t k = 1; k < liczbaScianBC; k++)
				{
					suma += Hbc[k][i][j];
				}
				H[i][j] += suma;
			}
		}
	}
	
	delete[] allJacobis;
	delete[] Nidx;
	delete[] Nidy;
	delete[] Hpci;
	delete[] Hnx;
	delete[] Hny;

	return H;
}

double ***HforAllElements(Elem4& elem, grid& newGrid, kwadratury& kwadratura, GlobalData& globe)
{
	double*** HforAll;
	int nE = newGrid.nE;
	HforAll = new double** [nE];

	for (int j = 0; j < nE; j++)
	{
		
		HforAll[j] = new double* [4];

		for (int k = 0; k < 4; k++)
		{
			HforAll[k] = new double* [4];
		}
	}

	//int i = 0; i < newGrid.nE; i++
	for (int i = 0; i < newGrid.nE; i++)
	{

		HforAll[i] = JacobiMatrix(elem, newGrid, kwadratura, i, globe);
		
	}

	//int i = 0; i < newGrid.nE; i++
	HLocalForElements(HforAll, newGrid.nE);

	return HforAll;
	delete[] HforAll;
	
}

double** CforElement(Elem4& elem, grid& newGrid, kwadratury& kwadratura, int _nrEl, GlobalData& globe)
{
	vector <double> x;
	vector <double> y;

	int nrEl = _nrEl;

	for (int i = 0; i < 4; i++)
	{
		int k = newGrid.EL[nrEl].ID[i];

		x.push_back(newGrid.ND[k - 1].x);
		y.push_back(newGrid.ND[k - 1].y);

	}

	int numberOfPoints = elem.numberOfPoints;
	int allPointsElem = numberOfPoints * numberOfPoints;

	double*** CNew;
	CNew = new double** [allPointsElem];
	for (size_t i = 0; i < allPointsElem; i++)
	{
		CNew[i] = new double* [4];
		for (size_t j = 0; j < 4; j++)
		{
			CNew[i][j] = new double[4];
		}
	}

	int counter;
	double eta;
	double ksi;
	int k = numberOfPoints;

	double density = globe.Density;
	double specificHeat = globe.SpecificHeat;

	double** allJacobis;
	double** C;
	allJacobis = new double* [allPointsElem];
	C = new double* [allPointsElem];

	for (int i = 0; i < allPointsElem; i++)
	{
		allJacobis[i] = new double[4];
		C[i] = new double[4];

	}

	for (int i = 0; i < allPointsElem; i++)
	{
		if (i % k == 0)
		{
			counter = 0;
		}

		if (k == 2)
		{
			ksi = kwadratura.W2p[counter];

			eta = kwadratura.W2p[i / k];
		}
		else if (k == 3)
		{
			ksi = kwadratura.W3p[counter];

			eta = kwadratura.W3p[i / k];
		}
		else if (k == 4)
		{
			ksi = kwadratura.W4p[counter];
			eta = kwadratura.W4p[i / k];
		}

		//double JacobiMat[2][2];
		double sumxksi = 0.0;
		double sumyksi = 0.0;
		double sumxeta = 0.0;
		double sumyeta = 0.0;

		for (int j = 0; j < 4; j++)
		{
			sumxksi += elem.tabKsi[i][j] * x[j] / 1.0;
			sumyksi += elem.tabKsi[i][j] * y[j] / 1.0;
			sumxeta += elem.tabEta[i][j] * x[j] / 1.0;
			sumyeta += elem.tabEta[i][j] * y[j] / 1.0;

		}

		allJacobis[i][0] = sumxksi / 1.0;
		allJacobis[i][1] = sumyksi / 1.0;
		allJacobis[i][2] = sumxeta / 1.0;
		allJacobis[i][3] = sumyeta / 1.0;

		double detJ = (allJacobis[i][0] * allJacobis[i][3] - allJacobis[i][1] * allJacobis[i][2]) / 1.0;
		//cout << detJ << " ";

		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				double wyn = density * elem.tabFunkcjiKsztaltow[i][j] * elem.tabFunkcjiKsztaltow[i][k] * specificHeat;
				CNew[i][j][k] = wyn * ksi * eta * detJ;
			}
				
		}
		counter++;

	}

	double suma;

	for (int i = 0; i < 4; i++)
	{
		suma = 0;

		for (int j = 0; j < 4; j++)
		{
			suma = CNew[0][i][j];

			for (int k = 1; k < allPointsElem; k++)
			{
				suma += CNew[k][i][j];
			}

			C[i][j] = suma;
		}

	}

	return C;

}

double*** CforAllElements(Elem4& elem, grid& newGrid, kwadratury& kwadratura, GlobalData& globe)
{
	double*** CforAll;
	int nE = newGrid.nE;
	CforAll = new double** [nE];

	for (int j = 0; j < nE; j++)
	{

		CforAll[j] = new double* [4];

		for (int k = 0; k < 4; k++)
		{
			CforAll[k] = new double* [4];
		}
	}

	//int i = 0; i < newGrid.nE; i++
	for (int i = 0; i < newGrid.nE; i++)
	{

		CforAll[i] = CforElement(elem, newGrid, kwadratura, i, globe);

	}

	//int i = 0; i < newGrid.nE; i++
	for (int i = 0; i < newGrid.nE; i++)
	{
		cout << "Macierz lokalna C dla elementu nr [" << i + 1 << "]" << endl;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				cout << CforAll[i][j][k] << " | ";
			}
			cout << endl;
		}

		cout << endl << endl;
	}

	return CforAll;

	delete[] CforAll;
}

double** CGlobal(double*** HforAllElements, grid& newGrid)
{
	int countOfNodes = newGrid.nN;

	double** CGlobal;
	CGlobal = new double* [countOfNodes];

	for (int i = 0; i < countOfNodes; i++)
	{
		CGlobal[i] = new double[countOfNodes];
	};

	for (int i = 0; i < countOfNodes; i++)
	{
		for (int j = 0; j < countOfNodes; j++)
		{
			CGlobal[i][j] = 0;
		}
	}



	//int el = 0; el < newGrid.nE; el++
	for (int el = 0; el < newGrid.nE; el++)
	{

		vector <int> ID;

		//dla kazdego elementu wczytujemy ID
		for (int i = 0; i < 4; i++)
		{
			int k = newGrid.EL[el].ID[i];
			//cout << k << " | ";
			ID.push_back(k);
		}

		//cout << endl;

		for (int j = 0; j < 4; j++)
		{

			for (int k = 0; k < 4; k++)
			{

				CGlobal[ID[j] - 1][ID[k] - 1] += HforAllElements[el][j][k];
				//cout << HGlobal[ID[j]-1][ID[k]-1] << endl;
			}
			//cout << endl;

		}

	}


	cout << "\n\nMacierz globalna C: " << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;

	for (int i = 0; i < countOfNodes; i++)
	{
		//cout << i << "[] ";
		for (int j = 0; j < countOfNodes; j++)
		{
			cout << CGlobal[i][j] << "  |";

		}
		cout << endl;
	}
	cout << endl;
	return CGlobal;

	delete[] CGlobal;
	delete[] HforAllElements;

}

double *PforElement(Elem4& elem, grid& newGrid, kwadratury& kwadratura, int _nrEl, GlobalData& globe)
{
	vector <double> x;
	vector <double> y;

	//x = { 0,0.025,0.025,0 };
	//y = { 0,0,0.025,0.025 };

	vector <int> BC;

	int nrEl = _nrEl;
	for (int i = 0; i < 4; i++)
	{
		int k = newGrid.EL[nrEl].ID[i];

		x.push_back(newGrid.ND[k - 1].x);
		y.push_back(newGrid.ND[k - 1].y);

	}

	int liczbaScianBC = 0;
	for (int i = 0; i < 4; i++)
	{
		int k = newGrid.EL[nrEl].ID[i];
		BC.push_back(newGrid.ND[k - 1].BC);
		//cout << k << ":";
		//cout << BC[i] << " | ";


	}
	//cout << endl;

	vector <int> ktoraSciana;
	int cnt = 0;
	for (int i = 0; i < 3; i++)
	{
		if (BC[i] == 1 && BC[i + 1] == 1)
		{
			liczbaScianBC++;
			ktoraSciana.push_back(cnt);
		}
		cnt++;
	}

	if (BC[3] == 1 && BC[0] == 1)
	{
		liczbaScianBC++;
		ktoraSciana.push_back(cnt);
	}

	/*cout << "BC na scianie nr: ";
	for (int i = 0; i < liczbaScianBC; i++)
	{
		cout << ktoraSciana[i] + 1 << ", ";
	}*/

	
	int liczba = elem.punktyNaPowierzchni;
	int pkt = liczba * 4;

	double** Pnew;
	Pnew = new double* [pkt];
	for (size_t i = 0; i < pkt; i++)
	{
		Pnew[i] = new double[4];
	}

	//cout << "Liczba scian z warunkiem: " << liczbaScianBC << endl;
	//tutaj bede przekazywal do funkcji liczbe punktow na powierzchni
	int counterWagi;
	double waga;
	double cond = globe.Alfa;
	double tempOt = globe.Tot;
	for (int i = 0; i < pkt; i++) //liczba bokow
	{

		if (i % elem.punktyNaPowierzchni == 0)
		{
			counterWagi = 0;
		}

		if (elem.punktyNaPowierzchni == 2)
		{
			waga = kwadratura.W2p[counterWagi];
		}
		else if (elem.punktyNaPowierzchni == 3)
		{
			waga = kwadratura.W3p[counterWagi];
		}
		else if (elem.punktyNaPowierzchni == 4)
		{
			waga = kwadratura.W4p[counterWagi];
		}

		//cout << waga << " || ";
		//wektor ma jeden wiersz

		for (int k = 0; k < 4; k++)
		{
			//cout << elem.tabNPowierzchnia[i][j] << " " << elem.tabNPowierzchnia[i][k] << " | ";
			//cout << "Mno¿e: " << elem.tabNPowierzchnia[i][k] << " * " << tempOt << " || " << endl;
			double wyn = cond * elem.tabNPowierzchnia[i][k];

			wyn = wyn * waga * tempOt;
			Pnew[i][k] = wyn;
			//cout << Pnew[i][k] << " | ";

		}
		//cout << endl << endl;

		counterWagi++;

	}

	//musimy utworzyæ tablice Hbc tylko na tych œcianach gdzie jest warunek brzegowy
	////double Hbc[2][4][4];
	double** PileBC;
	PileBC = new double* [liczbaScianBC];
	for (size_t i = 0; i < liczbaScianBC; i++)
	{
		PileBC[i] = new double[4];
	};


	
	//dodawanie do siebie wektorów na danym  na boku
	double detJBok;
	for (int i = 0; i < liczbaScianBC; i++)
	{
		if (ktoraSciana[i] + 1 == 4)
		{
			detJBok = sqrt((x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0]));
		}
		else
		{
			detJBok = sqrt((x[ktoraSciana[i] + 1] - x[ktoraSciana[i]]) * (x[ktoraSciana[i] + 1] - x[ktoraSciana[i]]) + (y[ktoraSciana[i] + 1] - y[ktoraSciana[i]]) * (y[ktoraSciana[i] + 1] - y[ktoraSciana[i]]));
		};

		detJBok = detJBok / 2;
		//cout << "DetJBok sciana[" << i + 1 << "] element[" << nrEl << "] = " << detJBok << endl;
		//cout << "DetJbok = " << detJBok << endl;
		//cout << "Hbc el:[" << nrEl + 1 << "] sciana[" << ktoraSciana[i] + 1 << "]" << endl;

		for (int j = 0; j < 4; j++)
		{
			if (elem.punktyNaPowierzchni == 2)
			{
				PileBC[i][j] = (Pnew[ktoraSciana[i] * 2][j] + Pnew[ktoraSciana[i] * 2 + 1][j]) * detJBok;
				//cout << PileBC[i][j] << " | ";
			}
			else if (elem.punktyNaPowierzchni == 3)
			{
				PileBC[i][j] = (Pnew[ktoraSciana[i] * 3][j] + Pnew[ktoraSciana[i] * 3 + 1][j] + Pnew[ktoraSciana[i] * 3 + 2][j]) * detJBok;
				//cout << PileBC[i][j] << " | ";
			}
			else if (elem.punktyNaPowierzchni == 4)
			{
				PileBC[i][j] = (Pnew[ktoraSciana[i] * 4][j] + Pnew[ktoraSciana[i] * 4 + 1][j] + Pnew[ktoraSciana[i] * 4 + 2][j] + Pnew[ktoraSciana[i] * 4 + 3][j]) * detJBok;
				//cout << PileBC[i][j] << " | ";
			}
			
		}
		//cout << endl;

	}

	//tablica zwracajaca P dla kazdego elementu
	//double tabWholeElement[4];
	double* tabWholeElement;
	tabWholeElement = new double[4];
	for (int i = 0; i < 4; i++)
	{
		tabWholeElement[i] = 0;
	}

	if (liczbaScianBC > 0)
	{
		for (int i = 0; i < liczbaScianBC; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tabWholeElement[j] += PileBC[i][j];
			}
			
		}
	}
	

	for (int i = 0; i < 4; i++)
	{
		//cout << tabWholeElement[i] << " || ";
	}


	return tabWholeElement;
}

double** PforAllElements(Elem4& elem, grid& newGrid, kwadratury& kwadratura, GlobalData& globe)
{
	double** PforAll;
	//double *pom = PforElement(elem, newGrid, kwadratura, 0, globe);

	int nE = newGrid.nE;
	PforAll = new double* [nE];

	for (int j = 0; j < nE; j++)
	{
		PforAll[j] = new double[4];
	}

	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			PforAll[i][j] = 0.0;
		}
	}

	//int i = 0; i < newGrid.nE; i++
	for (int i = 0; i < newGrid.nE; i++)
	{
		PforAll[i] = PforElement(elem, newGrid, kwadratura, i, globe);
		//PforAll[i] = PforElement(elem, newGrid, kwadratura, i, globe);
	
		//cout << "\nMacierz lokalna wektora P dla elementu nr [" << i + 1 << "]" << endl;
		for (int j = 0; j < 4; j++)
		{
			
			//cout << PforAll[i][j] << " | ";
		}
		//cout << endl;

	}

	/*cout << "P for All Elements: \n";
	for (int i = 0; i < newGrid.nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << PforAll[i][j] << " | ";
		}
		cout << endl;
	}*/
	

	return PforAll;

	delete[] PforAll;

}

double** HGlobal(double ***HforAllElements, grid& newGrid) {
	int countOfNodes = newGrid.nN;

	double** HGlobal;
	HGlobal = new double* [countOfNodes];

	for (int i = 0; i < countOfNodes; i++)
	{
		HGlobal[i] = new double [countOfNodes];
	};

	for (int i = 0; i < countOfNodes; i++)
	{
		for (int j = 0; j < countOfNodes; j++)
		{
			HGlobal[i][j] = 0;	
		}
	}

	
	
	//int el = 0; el < newGrid.nE; el++
	for (int el = 0; el < newGrid.nE; el++)
	{

		vector <int> ID;

		//dla kazdego elementu wczytujemy ID
		for (int i = 0; i < 4; i++)
		{
			int k = newGrid.EL[el].ID[i];
			//cout << k << " | ";
			ID.push_back(k);
		}

		//cout << endl;

		for (int j = 0; j < 4; j++)
		{

			for (int k = 0; k < 4; k++)
			{
				
				HGlobal[ID[j]-1][ID[k]-1] += HforAllElements[el][j][k];
				//cout << HGlobal[ID[j]-1][ID[k]-1] << endl;
			}
			//cout << endl;

		}

	}


	cout << "Macierz globalna H: " << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;

	for (int i = 0; i < countOfNodes; i++)
	{
		//cout << i << "[] ";
		for (int j = 0; j < countOfNodes; j++)
		{
				
				cout << HGlobal[i][j] << "  |";
			
		}
		cout << endl;
	}

	return HGlobal;

	delete[] HGlobal;
	delete[] HforAllElements;

}

double** HCdT(double** HGlobal, double** CGlobal, grid& newGrid, GlobalData& globe)
{
	int k = newGrid.nN;
	double** HCdt;
	HCdt = new double* [k];

	double dT = globe.SimulationStepTime;

	for (int i = 0; i < k; i++)
	{
		HCdt[i] = new double[k];
	}

	cout << "\nMacierz [H] + [C]/dT\n===============================================\n\n";
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j <k; j++)
		{
			//cout << HGlobal[i][j] << " | " << CGlobal[i][j];
			HCdt[i][j] = HGlobal[i][j] + CGlobal[i][j] / dT;
			cout << HCdt[i][j] <<  " | ";
		}
		cout << endl;
	}
	cout << endl;

	return HCdt;
}

double* PGlobal(double **PforAllElements, grid& newGrid)
{

	int countOfNodes = newGrid.nN;

	double *PGlobal;
	PGlobal = new double [countOfNodes];

	for (int i = 0; i < countOfNodes; i++)
	{
		PGlobal[i] = 0;
	}



	//int el = 0; el < newGrid.nE; el++
	for (int el = 0; el < newGrid.nE; el++)
	{

		vector <int> ID;

		//dla kazdego elementu wczytujemy ID
		for (int i = 0; i < 4; i++)
		{
			int k = newGrid.EL[el].ID[i];
			//cout << k << " | ";
			ID.push_back(k);
			//cout << k << " | ";
		}

		//cout << endl;


		for (int j = 0; j < 4; j++)
		{
			PGlobal[ID[j]-1] += PforAllElements[el][j];
			
		}
		
		

	}


	cout << "\nMacierz globalna wektora P: " << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl;

	for (int i = 0; i < countOfNodes; i++)
	{
		
			cout << PGlobal[i] << "  |";

		
		
	}

	cout << endl;

	return PGlobal;
}

double* tVector(double** HGlobalArray, double* PGlobalArray, grid& newGrid)
{
	double** CopyGlobal = HGlobalArray;
	CopyGlobal = new double* [newGrid.nN];
	for (int i = 0; i < newGrid.nN; i++)
	{
		CopyGlobal[i] = new double[newGrid.nN];
	}

	for (int i = 0; i < newGrid.nN; i++)
	{
		for (int j = 0; j < newGrid.nN; j++)
		{
			CopyGlobal[i][j] = HGlobalArray[i][j];
		}
	}

	double* CopyPGlobal;
	CopyPGlobal = new double[newGrid.nN];

	for (int i = 0; i < newGrid.nN; i++)
	{
		CopyPGlobal[i] = PGlobalArray[i];
	}

	int numberOfNodes = newGrid.nN;
	//vector<double> t(numberOfNodes);
	double* t;
	t = new double[numberOfNodes];
	
	//cout << HGlobalArray[0][0] << " Global Array \n\n";
	for (int i = 0; i < numberOfNodes; i++)
	{
		double pivot = CopyGlobal[i][i];

		for (int j = 0; j < numberOfNodes; j++)
		{
			CopyGlobal[i][j] /= pivot;
		}
		CopyPGlobal[i] /= pivot;

		for (int k = 0; k < numberOfNodes; k++)
		{
			if (k != i)
			{
				double factor = CopyGlobal[k][i];
				for (int j = 0; j < numberOfNodes; j++)
				{
					CopyGlobal[k][j] -= factor * CopyGlobal[i][j];
				}
				CopyPGlobal[k] -= factor * CopyPGlobal[i];
			}
		}
	}

	for (int i = numberOfNodes -1; i >=0; i--)
	{
		t[i] = CopyPGlobal[i];
		for (int j = i+1; j < numberOfNodes; j++)
		{
			t[i] -= CopyGlobal[i][j] * t[j];
		}
	}

//cout << "Wektor {t}" << endl;
//	for (int i = 0; i < numberOfNodes; i++)
//	{
//		cout << t[i] << " ";
//	}

	return t;
}

void vectorP(double** C, double* P,double** HC, grid& newGrid, GlobalData& globe)
{
	int N = newGrid.nN;
	double* t;
	double* Pnew;
	t = new double[N];
	Pnew = new double[N];

	for (int i = 0; i < N; i++)
	{
		t[i] = globe.InitialTemp;
		Pnew[i] = 0;
	}

		double** Copy;
	Copy = new double* [N];
	for (int i = 0; i < N; i++)
	{
		Copy[i] = new double[N];
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Copy[i][j] = C[i][j] / globe.SimulationStepTime;
		}
	}

	
	for (int iteration = globe.SimulationStepTime; iteration <= globe.SimulationTime; iteration += globe.SimulationStepTime)
	{
		for (int el = 0; el < N; el++)
		{
			Pnew[el] = 0;
		}
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Pnew[i] += Copy[i][j] * t[j];
		}
	}

	for (int i = 0; i < N; i++)
	{
		Pnew[i] += P[i];
		//cout << Pnew[i] << " | ";
	}
	//cout << endl;

	t = tVector(HC, Pnew, newGrid);

	//zadac pytanie o poprawnosc test case'ow
	/*for (int k = 0; k < N; k++)
	{
		cout << t[k] << " | ";
	}
	cout << endl;*/

	double min = t[0];
	double max = t[0];
	for (int k = 1; k < N; k++)
	{
		min = std::min(min, t[k]);
		max = std::max(max, t[k]);
	}

	cout << iteration << "\t" << min << "\t\t" << max << endl;
	
	}

}

//lab1  
void createNodes(grid& nowa, GlobalData& globe, int nodesNumber, ifstream& odczyt) {
	for (int i = 0; i < nodesNumber; i++)
	{
		string a, b, c;
		odczyt >> a >> b >> c;

		a.pop_back();
		b.pop_back();

		double num1 = stod(a);
		double num2 = stod(b);
		double num3 = stod(c);

		node newNode;
		newNode.x = num2;
		newNode.y = num3;
		newNode.t = globe.InitialTemp;

		nowa.ND.push_back(newNode);

	}
}

void displayNodes(grid &newGrid,int nodesNumber) {
	cout << "\nWyswietlenie wezlow: " << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;
	for (int i = 0; i < nodesNumber; i++)
	{
		if (newGrid.ND[i].x > 0 && newGrid.ND[i].x < 0.1) {
			cout << i + 1 << ":";
			cout << newGrid.ND[i].x << "\t";
			cout << newGrid.ND[i].y << "\t\t";
			cout << newGrid.ND[i].BC << endl;
		}
		else {
			cout << i + 1 << ":";
			cout << newGrid.ND[i].x << "\t\t";
			cout << newGrid.ND[i].y << "\t\t";
			cout << newGrid.ND[i].BC << endl;
		}
	
		
	}

}

void createElements(grid &newGrid, int elementsNumber, ifstream &odczyt) {
	for (int i = 0; i < elementsNumber; i++)
	{
		string a, b, c, d, e;
		odczyt >> a >> b >> c >> d >> e;
		b.pop_back();
		c.pop_back();
		d.pop_back();

		int num1 = stoi(b);
		int num2 = stoi(c);
		int num3 = stoi(d);
		int num4 = stoi(e);


		element newElement;

		newElement.ID[0] = num1;
		newElement.ID[1] = num2;
		newElement.ID[2] = num3;
		newElement.ID[3] = num4;

		newGrid.EL.push_back(newElement);

	}
}

void displayElements(grid &newGrid, int elementsNumber) {
	cout << "\nWyswietlenie elementow: " << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;
	for (int i = 0; i < elementsNumber; i++)
	{
		cout << i + 1 << ".   ";
		cout << newGrid.EL[i].ID[0] << " " << newGrid.EL[i].ID[1] << " " << newGrid.EL[i].ID[2] << " " << newGrid.EL[i].ID[3] << endl;

	}
}

int main()
{
	
	cout << "\LABORATORIUM NR 1" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;

	int wybor;
	string plik;
	cout << "Wybierz plik:" << endl;
	cout << "1. Test1_4_4.txt" << endl;
	cout << "2. Test2_4_4_MixGrid.txt" << endl;
	cout << "3. Test3_31_31_kwadrat.txt" << endl;
	cin >> wybor;

	if (wybor==1)
	{
		plik = "Test1_4_4.txt";
	}
	else if(wybor==2)
	{
		plik = "Test2_4_4_MixGrid.txt";
	}
	else if (wybor ==3)
	{
		plik = "Test3_31_31_kwadrat.txt";
	}
	else {
		plik = "";
	}

	system("cls");

	cout << "\LABORATORIUM NR 1" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;

	GlobalData global;
	ifstream odczyt(plik);

	if (!odczyt) {
		cout << "Nie mozna otworzyc pliku!" << endl;
		system("pause");
	}

	{
		odczyt >> global.SimulationTime;
		odczyt >> global.SimulationStepTime;
		odczyt >> global.Conductivity;
		odczyt >> global.Alfa;
		odczyt >> global.Tot;
		odczyt >> global.InitialTemp;
		odczyt >> global.Density;
		odczyt >> global.SpecificHeat;
	}

	//odczyt iloœci wêz³ów
	int nodesNumber;
	odczyt >> nodesNumber;

	//odczyt iloœci elementów
	int elementsNumber;
	odczyt >> elementsNumber;

	//definicja struktury grid
	grid newGrid;
	newGrid.nN = nodesNumber;
	newGrid.nE = elementsNumber;



	//tworenie nodów
	createNodes(newGrid, global, nodesNumber, odczyt);

	//wyswietlanie nodów
	displayNodes(newGrid, nodesNumber);

	//tworzenie elementów
	createElements(newGrid, elementsNumber, odczyt);

	//wyœwietlanie elementów
	displayElements(newGrid, elementsNumber);
	
	//operacje potrzebne do mojego sposobu wpisywania do pliku
	//obliczenie ile bc==1
	vector <int> tab;
	cout << endl;
	int counter = 0;
	int wymiarSiatki = sqrt(nodesNumber);
	for (int i = 1; i <= wymiarSiatki; i++)
	{
		for (int j = 1; j <= wymiarSiatki; j++)
		{
			if (i == 1 || i == wymiarSiatki || j == 1 || j == wymiarSiatki)
			{
				counter++;
			}
		}
	}

	//cout << endl << counter << endl;
	//zapisanie numerów wezlow z BC==1
	for (int i = 0; i < counter; i++)
	{
		string a;
		odczyt >> a;
		//usuwanie przecinków
		if (i != counter - 1)
		{

			a.pop_back();
			int num = stoi(a);
			tab.push_back(num);
		}
		else
		{
			int num = stoi(a);
			//numery nodów, które maj¹ mieæ bc==1 zapisuje do wektora
			tab.push_back(num);
		}

	}

	for (int i = 0; i < counter; i++)
	{
		newGrid.ND[tab[i] - 1].BC = 1;
	}

	//wyswietlanie nodów z BC
	displayNodes(newGrid, nodesNumber);

	//lab2
	cout << endl;
	/*cout << "\LABORATORIUM NR 2" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;*/
	kwadratury kwadratura;
	//lab4
	cout << "\nLABORATORIUM NR 4" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;
	Elem4 elem(2,2);
	elem4(kwadratura, elem);
	//JacobiMatrix(elem, newGrid, kwadratura, 0);
	double*** HforAllElementss;
	HforAllElementss = HforAllElements(elem, newGrid, kwadratura, global);
	double** HGlobalArray = HGlobal(HforAllElementss, newGrid);
	//PforElement(elem, newGrid, kwadratura,0, global);
	double** PforAllElementss;
	PforAllElementss = PforAllElements(elem, newGrid, kwadratura, global);
	double* PGlobalArrray = PGlobal(PforAllElementss, newGrid);

	//wektor T
	double* t0;
	t0 = tVector(HGlobalArray, PGlobalArrray, newGrid);

	//CforElement
	double*** CforAllElementss;
	CforAllElementss = CforAllElements(elem, newGrid, kwadratura, global);
	double** CGlobalArray = CGlobal(CforAllElementss, newGrid);


	//HCdT(double** HGlobal, double** CGlobal, grid& newGrid, GlobalData& globe)
	double** HCArray;
	HCArray = HCdT(HGlobalArray, CGlobalArray, newGrid, global);

	vectorP(CGlobalArray, PGlobalArrray, HCArray, newGrid, global);

}