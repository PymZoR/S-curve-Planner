#include "motionPlanning.h"
#include "bbpr.h"

void computePeriods(double Xpeak[], double T[]) {
    int n = 0;
    double coeffs[4] = {0, 0, 0, 0};
    double Xmax = 0;

    n = 3;
    coeffs[0] = 2*Xpeak[3];
    coeffs[1] = T[1]*Xpeak[3] + 3*T[2]*Xpeak[3];
    coeffs[2] = T[1]*T[2]*Xpeak[3] + pow(T[2], 2)*Xpeak[3];
    coeffs[3] = -Xpeak[0];
    T[3] = findRoot(coeffs, n);

    Xmax = (2*T[2] + 2*T[3])*T[3]*Xpeak[3]/2;

    if (Xmax > Xpeak[1]) {
        n = 2;
        coeffs[0] = Xpeak[3];
        coeffs[1] = T[2]*Xpeak[3];
        coeffs[2] = -Xpeak[1];
        T[3] = findRoot(coeffs, n);
    } else {
        Xpeak[1] = Xmax;
    }

    Xmax = T[3]*Xpeak[3];

    if (Xmax > Xpeak[2]) {
        n = 1;
        coeffs[0] = Xpeak[3];
        coeffs[1] = -Xpeak[2];
        T[3] = findRoot(coeffs, n);
    } else {
        Xpeak[2] = Xmax;
    }

    n = 2;
    coeffs[0] = T[3]*Xpeak[3];
    coeffs[1] = T[1]*T[3]*Xpeak[3] + 3*pow(T[3], 2)*Xpeak[3];
    coeffs[2] = T[1]*pow(T[3], 2)*Xpeak[3] + 2*pow(T[3], 3)*Xpeak[3] - Xpeak[0];
    T[2] = findRoot(coeffs, n);

    Xmax = (2*T[2] + 2*T[3])*T[3]*Xpeak[3]/2;

    if (Xmax > Xpeak[1]) {
        n = 1;
        coeffs[0] = T[3]*Xpeak[3];
        coeffs[1] = pow(T[3], 2)*Xpeak[3] - Xpeak[1];
        T[2] = findRoot(coeffs, n);
    } else {
        Xpeak[1] = Xmax;
    }

    n = 1;
    coeffs[0] = T[2]*T[3]*Xpeak[3] + pow(T[3], 2)*Xpeak[3];
    coeffs[1] = pow(T[2], 2)*T[3]*Xpeak[3] + 3*T[2]*pow(T[3], 2)*Xpeak[3] + 2*pow(T[3], 3)*Xpeak[3] - Xpeak[0];
    T[1] = findRoot(coeffs, n);

}

double getSetpoint(double Xpeak[], double T[], double t) {
    if (T[0] <= t && t <= T[3]) {
        return pow((t-(T[0])), 3)*Xpeak[3]/6;
    }
    if (T[3] <= t && t <= T[2] + T[3]) {
        return pow((t-(T[3])), 2)*Xpeak[2]/2 + (t-(T[3]))*pow(T[3], 2)*Xpeak[3]/2 + pow(T[3], 3)*Xpeak[3]/6;
    }
    if (T[2] + T[3] <= t && t <= T[2] + 2*T[3]) {
        return -pow((t-(T[2] + T[3])), 3)*Xpeak[3]/6 + pow((t-(T[2] + T[3])), 2)*T[3]*Xpeak[3]/2 + (t-(T[2] + T[3]))*(T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2) + pow(T[2], 2)*Xpeak[2]/2 + T[2]*pow(T[3], 2)*Xpeak[3]/2 + pow(T[3], 3)*Xpeak[3]/6;
    }
    if (T[2] + 2*T[3] <= t && t <= T[1] + T[2] + 2*T[3]) {
        return (t-(T[2] + 2*T[3]))*Xpeak[1] + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2)*T[3] + pow(T[2], 2)*Xpeak[2]/2 + T[2]*pow(T[3], 2)*Xpeak[3]/2 + pow(T[3], 3)*Xpeak[3]/2;
    }
    if (T[1] + T[2] + 2*T[3] <= t && t <= T[1] + T[2] + 3*T[3]) {
        return -pow((t-(T[1] + T[2] + 2*T[3])), 3)*Xpeak[3]/6 + (t-(T[1] + T[2] + 2*T[3]))*(T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]) + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2)*T[3] + T[1]*Xpeak[1] + pow(T[2], 2)*Xpeak[2]/2 + T[2]*pow(T[3], 2)*Xpeak[3]/2 + pow(T[3], 3)*Xpeak[3]/2;
    }
    if (T[1] + T[2] + 3*T[3] <= t && t <= T[1] + 2*T[2] + 3*T[3]) {
        return -pow((t-(T[1] + T[2] + 3*T[3])), 2)*Xpeak[2]/2 + (t-(T[1] + T[2] + 3*T[3]))*(T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2) + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2)*T[3] + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3])*T[3] + T[1]*Xpeak[1] + pow(T[2], 2)*Xpeak[2]/2 + T[2]*pow(T[3], 2)*Xpeak[3]/2 + pow(T[3], 3)*Xpeak[3]/3;
    }
    if (T[1] + 2*T[2] + 3*T[3] <= t && t <= T[1] + 2*T[2] + 4*T[3]) {
        return pow((t-(T[1] + 2*T[2] + 3*T[3])), 3)*Xpeak[3]/6 - pow((t-(T[1] + 2*T[2] + 3*T[3])), 2)*T[3]*Xpeak[3]/2 + (t-(T[1] + 2*T[2] + 3*T[3]))*pow(T[3], 2)*Xpeak[3]/2 + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2)*T[2] + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3]/2)*T[3] + (T[2]*Xpeak[2] + pow(T[3], 2)*Xpeak[3])*T[3] + T[1]*Xpeak[1] + T[2]*pow(T[3], 2)*Xpeak[3]/2 + pow(T[3], 3)*Xpeak[3]/3;
    }
    else {
        return Xpeak[0];
    }
}
