#!/usr/bin/python

import sys, os, getopt
from sympy import *
import re

def makePow(expr):
    '''Convert Python x**y syntax to C math.pow(x, y)'''
    return re.sub(r'([\w]+(?:\[\d\])?)\*\*([\d])', r'pow(\1, \2)', str(expr))

def genSigns(n):
    '''Generate the signs of XPeak[n] for each period'''
    if (n == 1):
        return [1]
    else:
        return genSigns(n-1) + [1] + [-i for i in genSigns(n-1)]

def genOrders(n, val):
    '''Generate which Xpeak is constant for each period'''
    if (val == 1):
        return [n]
    else:
        previous = genOrders(n, val-1)
        new =  previous + [ n-val+1 ] + previous

        return new

def genTemplates(n):
    '''Generate the intervals and the functions generating every X[p]'''
    Xpeak = IndexedBase('Xpeak')
    T     = IndexedBase('T')
    i, j, t  = symbols('i, j, t', integer=True)

    signs       = genSigns(n)

    orders      = genOrders(n, n)

    timePeriods = [ T[order] for order in orders]

    timeFrames  = [T[0]] + [ sum(timePeriods[:i+1]) for i in range(len(timePeriods)) ]

    intervals = [[timeFrames[i], timeFrames[i+1]] for i in range(len(timeFrames)-1)]
 
    templates = [ [ 0 for inter in intervals] for i in range(n+1) ]

    for i in range(n, -1, -1):
        for j, order in enumerate(orders):
            peakValue = signs[j] * Xpeak[order]

            if (i == order):
                templates[i][j] = peakValue

            elif (i < n):
                func = integrate(templates[i+1][j], t)

                for k in range(j+1, len(orders)):
                    templates[i][k] += func.subs({ t: T[order] })

                templates[i][j] += func

    return templates, intervals

def genPolyA(n, p):
    Xpeak = IndexedBase('Xpeak')
    T    = IndexedBase('T')
    i, j = symbols('i, j', integer=True)
    
    A = ((Xpeak[n]/2**(n)) * Product( Sum(2**j * T[n+1-i+j], (j, 0, i-1)).doit() + T[n+1-i], (i, 1, n)).doit() - Xpeak[0]).expand()
    
    return Poly(A, T[p])

def genExprB(n, q):
    Xpeak = IndexedBase('Xpeak')
    T     = IndexedBase('T')
    i, j  = symbols('i, j', integer=True)

    B = (Xpeak[n]/2**(n-q)) * Product( Sum(2**j * T[n+1-i+j], (j, 0, i-1)).doit() + T[n+1-i], (i, 1, n-q)).doit()
    
    return B

def genPolyC(n, p, q):
    Xpeak = IndexedBase('Xpeak')
    T    = IndexedBase('T')
    i, j = symbols('i, j', integer=True)
    
    C = ((Xpeak[n]/2**(n-q)) * Product( Sum(2**j * T[n+1-i+j], (j, 0, i-1)).doit() + T[n+1-i], (i, 1, n-q)).doit() - Xpeak[q])
    
    return Poly(C, T[p])


def generateFiles(n, outDir):
    Xpeak   = IndexedBase('Xpeak')
    T       = IndexedBase('T')
    i, j, t = symbols('i, j, t', integer=True)

    templates, intervals = genTemplates(n)

    try:
        with open(os.path.join(outDir, 'motionPlanning.h'), 'wr') as file:
            file.write('#include <math.h>\n\n')
            file.write('void computePeriods(double Xpeak[], double T[]);\n') 
            file.write('double getSetpoint(double Xpeak[], double T[], double t);\n') 
    except Exception as err: 
        error(err)

    try:
        with open(os.path.join(outDir, 'motionPlanning.c'), 'w')  as file: 
            file.write('#include "motionPlanning.h"\n')
            file.write('#include "bbpr.h"\n\n')
            file.write('void computePeriods(double Xpeak[], double T[]) {\n')
            file.write('    int n = 0;\n')
            file.write('    double coeffs[%d] = {%s};\n' % (n+1, ', '.join([ '0' for k in range(n+1)])))
            file.write('    double Xmax = 0;\n\n')

            for p in range(n, 1-1, -1):
                polyA   = genPolyA(n, p)
                degreeA = polyA.degree()
                coeffsA = polyA.all_coeffs()

                file.write('    n = %d;\n' % degreeA)
                for k, coeff in enumerate(coeffsA):
                    file.write('    coeffs[%d] = %s;\n' % (k, makePow(coeff)))
                file.write('    T[%d] = findRoot(coeffs, n);\n\n' % p)

                for q in range(1, p-1+1):
                    B   = genExprB(n, q)
                    polyC   = genPolyC(n, p, q)
                    degreeC = polyC.degree()
                    coeffsC = polyC.all_coeffs();

                    file.write('    Xmax = %s;\n\n' % makePow(B))
                    file.write('    if (Xmax > Xpeak[%d]) {\n' % (q))
                    file.write('        n = %d;\n' % (degreeC))
                    for k, coeff in enumerate(coeffsC):
                        file.write('        coeffs[%d] = %s;\n' % (k, makePow(coeff)))
                    file.write('        T[%d] = findRoot(coeffs, n);\n' % p)
                    file.write('    } else {\n')
                    file.write('        Xpeak[%d] = Xmax;\n' % q)
                    file.write('    }\n\n')
            file.write('}\n\n')
            
            file.write('double getSetpoint(double Xpeak[], double T[], double t) {\n')
            
            for i, template in enumerate(templates[0]):
                tA = intervals[i][0]
                tB = intervals[i][1]

                file.write('    if (%s <= t && t <= %s) {\n' % (tA, tB))
                file.write('        return %s;\n' % makePow(template).replace('t', '(t-(%s))' % tA))
                file.write('    }\n')
            file.write('    else {\n')
            file.write('        return Xpeak[0];\n')
            file.write('    }\n')
            file.write('}\n')
    
    except Exception as err: 
        error(err)

def error(err):
    print err
    sys.exit(1)

def exit():
    print 'usage: generator.py -n <order> -o <outputDirectory>'
    sys.exit(2)

def main(argv):
    opts = []

    try:
        opts = getopt.getopt(argv, "n:o:h", ["order", "out", "help"])[0]
    except getopt.GetoptError:
        exit()

    if len(opts) < 2:
        exit()

    order  = int(opts[0][1])
    outDir = opts[1][1]
    generateFiles(order, outDir)

if __name__ == "__main__":
   main(sys.argv[1:])
