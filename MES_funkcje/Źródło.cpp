#include<iostream>
#include <string>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

//BC - bumper condition 1-jesli zachodzi 0-jesli nie

//lab1
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

//lab2
struct kwadratury {
	double PC2p[2] = { -1/sqrt(3),1/sqrt(3)};
	double W2p[2] = { 1,1 };

	double PC3p[3] = { -sqrt(0.6),0,sqrt(0.6)};
	double W3p[3] = { 5/9.,8/9.,5/9. };

};



double function1d(double x) {
	return (2*x*x+3*x-8);
}

double function2d(double ksi, double eta) {
	return (-5 * ksi * ksi * eta + 2*ksi*eta*eta+10);
}

//system - 1d/2d
//version - 2points, 3points
double integrationLab2(kwadratury &kwadratura,int version,int system) {
	
	double sum = 0;
	if (system==1)
	{
		if (version==2)
		{
			for (int i = 0; i < version; i++)
			{
				sum += function1d(kwadratura.PC2p[i]) * kwadratura.W2p[i];
				
			}
			
		}
		 
		if (version == 3)
		{
			for (int i = 0; i < version; i++)
			{
				sum += function1d(kwadratura.PC3p[i]) * kwadratura.W3p[i];
				
			}
			return sum;
		}
	
	}
	
	if (system==2)
	{
		if (version == 2)
		{
			for (int i = 0; i < version; i++)
			{
				for (int j = 0; j < version; j++) {
					sum += function2d(kwadratura.PC2p[i], kwadratura.PC2p[j]) * kwadratura.W2p[i] * kwadratura.W2p[j];
				}
			}
			
		}
		else if (version == 3)
		{
			for (int i = 0; i < version; i++)
			{
				for (int j = 0; j < version; j++) {
					sum += function2d(kwadratura.PC3p[i], kwadratura.PC3p[j]) * kwadratura.W3p[i] * kwadratura.W3p[j];
				}
			}
			
		}
	}
	return sum;
	
};

//lab3

double N1(double ksi) {
	double N1 = (1 - ksi) / (2);
	return N1;
}

double N2(double ksi) {
	double N2 = (1 + ksi) / (2);
	return N2;
}

double fun(double x) {
	return (3*x*x - 6*x + 1);
}


double integrationLab3(int version,double x1, double x2, kwadratury &kwadratura) {
	double Pcx1, Pcx2, Pcx3, detJ, integral;
	if (version == 2) {
		 Pcx1 = (N1(kwadratura.PC2p[0]) * x1) + (N2(kwadratura.PC2p[0]) * x2);

		 Pcx2 = (N1(kwadratura.PC2p[1]) * x1) + (N2(kwadratura.PC2p[1]) * x2);

		//cout << "Pcx1 = " << Pcx1 << endl;
		//cout << "Pcx1 = " << Pcx2 << endl;

		 detJ = (x2 - x1) / 2;

		 integral = (fun(Pcx1) * kwadratura.W2p[0] + fun(Pcx2) * kwadratura.W2p[1]) * detJ;

		cout << "Wynik calki 2p wynosi: " << integral << endl;
		return integral;

	}
	else if (version==3)
	{
		 Pcx1 = (N1(kwadratura.PC3p[0]) * x1) + (N2(kwadratura.PC3p[0]) * x2);

		 Pcx2 = (N1(kwadratura.PC3p[1]) * x1) + (N2(kwadratura.PC3p[1]) * x2);

		 Pcx3 = (N1(kwadratura.PC3p[2]) * x1) + (N2(kwadratura.PC3p[2]) * x2);

		double detJ = (x2 - x1) / 2;

		integral = (fun(Pcx1) * kwadratura.W3p[0] + fun(Pcx2) * kwadratura.W3p[1] + fun(Pcx3) * kwadratura.W3p[2]) * detJ;

		cout << "Wynik calki 3p wynosi: " << integral << endl;
		return integral;
	}
	
}

//lab4

//lab4
//tutaj wzór funkcji kszta³tu

double N1ksi(double eta)
{
	return (-0.25 * (1 - eta));
}

double N2ksi(double eta)
{
	return (0.25 * (1 - eta));
};

double N3ksi(double eta)
{
	return (0.25 * (1 + eta));
}

double N4ksi(double eta)
{
	return (-0.25 * (1 + eta));
};

double N1eta(double ksi)
{
	return (-0.25 * (1 - ksi));
}

double N2eta(double ksi)
{
	return (-0.25 * (1 + ksi));
}

double N3eta(double ksi)
{
	return (0.25 * (1 + ksi));
}

double N4eta(double ksi)
{
	return (0.25 * (1 - ksi));
}



struct Elem2D {
	double tabKsi[4][4];
	double tabEta[4][4];
};

struct Elem3D {
	double tabKsi[9][4];
	double tabEta[9][4];
};

//przekazujemy flagê, która tworzy odpowiedni wymiar tablicy
struct Elem4 {
	int numberOfPoints;

	double** tabKsi;
	double** tabEta;
	
	Elem4(int numberOfPoints) : numberOfPoints(numberOfPoints) {
		int k = numberOfPoints * numberOfPoints;
		tabKsi = new double* [k];
		tabEta = new double* [k];
		for (int i = 0; i < k; i++)
		{
			tabKsi[i] = new double[4];
			tabEta[i] = new double[4];
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
			//cout << ksi << endl;

			/*if (i < elem.numberOfPoints)
			{
				eta = kwadratury.PC2p[0];
			}
			else
			{
				eta = kwadratury.PC2p[1];
			}*/
			eta = kwadratury.PC2p[i / 2];


			elem.tabKsi[i][0] = N1ksi(eta);
			elem.tabKsi[i][1] = N2ksi(eta);
			elem.tabKsi[i][2] = N3ksi(eta);
			elem.tabKsi[i][3] = N4ksi(eta);

			elem.tabEta[i][0] = N1eta(kwadratury.PC2p[counter]);
			elem.tabEta[i][1] = N2eta(kwadratury.PC2p[counter]);
			elem.tabEta[i][2] = N3eta(kwadratury.PC2p[counter]);
			elem.tabEta[i][3] = N4eta(kwadratury.PC2p[counter]);


			//cout << eta << " " << N1ksi(eta) << " " << N2ksi(eta) << " " << N3ksi(eta) << " " << N4ksi(eta) << endl;
			//cout << ksi << " " << N1eta(kwadratury.PC2p[counter]) << " " << N2eta(kwadratury.PC2p[counter]) << " " << N3eta(kwadratury.PC2p[counter]) << " " << N4eta(kwadratury.PC2p[counter]) << endl;

			counter++;
		}
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


			//cout << eta << " " << N1ksi(eta) << " " << N2ksi(eta) << " " << N3ksi(eta) << " " << N4ksi(eta) << endl;
			//cout << ksi << " " << N1eta(ksi) << " " << N2eta(ksi) << " " << N3eta(ksi) << " " << N4eta(ksi) << endl;

			counter++;


		}
	}
}

double **JacobiMatrix(Elem4& elem, grid& newGrid, kwadratury& kwadratura)
{
	vector <double> x;
	vector <double> y;

	x = { 0,0.025,0.025,0 };
	y = { 0,0,0.025,0.025 };
	
	/*
	for (int i = 0; i < 4; i++)
	{
		int k = newGrid.EL[0].ID[i];
		cout << k << endl;
		x.push_back(newGrid.ND[k-1].x);
		y.push_back(newGrid.ND[k - 1].y);
		cout << x[i] << " " << y[i] << endl;
	}
	*/
	
	int k = elem.numberOfPoints;
	int k2 = k * k;
	
	double** allJacobis;
	double** Nidx;
	double** Nidy;
	double** H;
;
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

		

		//double JacobiMat[2][2];
		double sumxksi = 0;
		double sumyksi = 0;
		double sumxeta = 0;
		double sumyeta = 0;

		for (int j = 0; j < 4; j++)
		{
			sumxksi += elem.tabKsi[i][j] * x[j];
			sumyksi += elem.tabKsi[i][j] * y[j];
			sumxeta += elem.tabEta[i][j] * x[j];
			sumyeta += elem.tabEta[i][j] * y[j];

		}

		//cout << "[" << sumxksi << " " << sumyksi << " " << sumxeta << " " << sumyeta << "]" << endl;

		allJacobis[i][0] = sumxksi;
		allJacobis[i][1] = sumyksi;
		allJacobis[i][2] = sumxeta;
		allJacobis[i][3] = sumyeta;

		//1) dNi / dx
		double detJ = allJacobis[i][0] * allJacobis[i][3] - allJacobis[i][1] * allJacobis[i][4];

		//odwrócenie macierzy
		double pom;
		pom = allJacobis[i][0];
		allJacobis[i][0] = allJacobis[i][3];
		allJacobis[i][3] = pom;
		allJacobis[i][1] = -allJacobis[i][1];
		allJacobis[i][2] = -allJacobis[i][2];

		//cout << "DetJ:" << detJ << endl;

		//wyznacznik 1/det
		double reverseDetJ = 1. / detJ;
		//cout << reverseDetJ;
		for (int k = 0; k < 4; k++)
		{
			allJacobis[i][k] *= reverseDetJ;
			cout << allJacobis[i][k] << " | ";
		}

		cout << endl;

		//slajd13 macierz tabeladwuwymiarowazmacierza

		//Ni/dx macierz

		for (int j = 0; j < 4; j++)
		{
			Nidx[i][j] = allJacobis[i][0] * elem.tabKsi[i][j] + allJacobis[i][1] * elem.tabEta[i][j];
			Nidy[i][j] = allJacobis[i][2] * elem.tabKsi[i][j] + allJacobis[i][3] * elem.tabEta[i][j];
			
		}
		
		//macierzH dla punktów ca³kowania

		double kt = 30;

		cout << "Ksi: " << ksi << endl;
		cout << "Eta: " << eta << endl;

		for (int j = 0; j < k2; j++)
		{
			
			for (int k = 0; k < 4; k++)
			{
				Hnx[i][j][k] = Nidx[i][j] * Nidx[i][k];
				Hny[i][j][k] = Nidy[i][j] * Nidy[i][k];
				Hpci[i][j][k] = kt * (Hnx[i][j][k] + Hny[i][j][k]) * detJ;
				//obliczone z wagami
				Hpci[i][j][k] *= ksi * eta;
				//H[j][k] += Hpci[i][j][k];
				//cout << Hpci[i][j][k] << " ";
				
			}
			cout << endl;
			
		}
		cout << endl;

		counter++;

	}

	//dodawanie macierzy H; zmiany wprowadzic w schemacie punktow
	double suma;
	//i - kolumny
	//j - wiersz
	//
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

			//cout << suma << " | ";
			
			H[i][j] = suma;
			cout << H[i][j] << " ";
		}
		cout << endl;
	}
	

	return H;
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
	cout << "\LABORATORIUM NR 2" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;
	kwadratury kwadratura;
	cout << "2punkty 1D wynik: " << integrationLab2(kwadratura, 2, 1) << endl;
	cout << "3punkty 1D wynik: " << integrationLab2(kwadratura, 3, 1) << endl;
	cout << "2punkty 2D wynik: " << integrationLab2(kwadratura, 2, 2) << endl;
	cout << "3punkty 2D wynik: " << integrationLab2(kwadratura, 3, 2) << endl;

	//lab3
	cout << "\nLABORATORIUM NR 3" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;
	double integral2p = integrationLab3(2, -3, 6.5, kwadratura);
	double integral3p = integrationLab3(3, -3, 6.5, kwadratura);

	//lab4
	cout << "\nLABORATORIUM NR 4" << endl;
	cout << "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|" << endl << endl;
	Elem4 elem(2);
	elem4(kwadratura, elem);
	JacobiMatrix(elem, newGrid, kwadratura);





}