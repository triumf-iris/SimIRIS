#ifndef ENERGYLOSSMANAGER_H
#define ENERGYLOSSMANAGER_H

#include <TMath.h>
#include <TRandom3.h>
#include "nucleus.h"

class EnergyLossManager()
{
public:
    void EnergyLossManager();
    double eval(double, double[100], double[100]);
    double eloss(nucleus, double, double, double, double[100], double[100]);
    double elossFi(double, double, double[100], double[100]);
};

#endif
