
#include "macierzH.h"

void HLocalForElements(double*** HforAll, double nE)
{
	for (int i = 0; i < nE; i++)
	{
		cout << "Macierz lokalna H dla elementu nr [" << i + 1 << "]" << endl;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				cout << HforAll[i][j][k] << " | ";
			}
			cout << endl;
		}

		cout << endl << endl;
	}
}