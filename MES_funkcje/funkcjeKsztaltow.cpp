#include "funkcjeKsztaltow.h"

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

double N1ksieta(double ksi, double eta)
{
	return (0.25 * (1 - ksi) * (1 - eta));
}

double N2ksieta(double ksi, double eta)
{
	return (0.25 * (1 + ksi) * (1 - eta));
}

double N3ksieta(double ksi, double eta)
{
	return (0.25 * (1 + ksi) * (1 + eta));
}

double N4ksieta(double ksi, double eta)
{
	return (0.25 * (1 - ksi) * (1 + eta));
}