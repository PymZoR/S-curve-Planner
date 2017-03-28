#include "motionPlanning.h"
#include "bbpr.h"

void computePeriods(double Xpeak[], double T[]) {
    int n = 0;
    double coeffs[3] = {0, 0, 0};
    double Xmax = 0;

    n = 2;
    coeffs[0] = Xpeak[2];
    coeffs[1] = T[1]*Xpeak[2];
    coeffs[2] = -Xpeak[0];
    T[2] = findRoot(coeffs, n);

    Xmax = T[2]*Xpeak[2];

    if (Xmax > Xpeak[1]) {
        n = 1;
        coeffs[0] = Xpeak[2];
        coeffs[1] = -Xpeak[1];
        T[2] = findRoot(coeffs, n);
    } else {
        Xpeak[1] = Xmax;
    }

    n = 1;
    coeffs[0] = T[2]*Xpeak[2];
    coeffs[1] = pow(T[2], 2)*Xpeak[2] - Xpeak[0];
    T[1] = findRoot(coeffs, n);

}

double getSetpoint(double Xpeak[], double T[], double t) {
    if (T[0] <= t && t <= T[2]) {
        return pow((t-(T[0])), 2)*Xpeak[2]/2;
    }
    if (T[2] <= t && t <= T[1] + T[2]) {
        return (t-(T[2]))*Xpeak[1] + pow(T[2], 2)*Xpeak[2]/2;
    }
    if (T[1] + T[2] <= t && t <= T[1] + 2*T[2]) {
        return -pow((t-(T[1] + T[2])), 2)*Xpeak[2]/2 + (t-(T[1] + T[2]))*T[2]*Xpeak[2] + T[1]*Xpeak[1] + pow(T[2], 2)*Xpeak[2]/2;
    }
    else {
        return Xpeak[0];
    }
}
