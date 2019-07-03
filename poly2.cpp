/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019
///
/// poly2.cpp
/// Sample polynomial operations

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "poly.h"

int main ()
{
    int n;

    Polynomial a, b, c;

    a=Polynomial({1})+Term(1,1); /// $x+1$ // Works
    //a=Polynomial({1,1}); /// x+1$ // Does not work
    //a=Polynomial({1,1,1}); /// $x^2+x+1$ // Works
    std::cout << "a: " << a << std::endl;

    b=Polynomial({1,0,1}); /// $x^2+1$
    std::cout << "b: " << b << std::endl;

    c=a+b;
    std::cout << "(" << a << ")" <<  " + " << "(" << b << ")" << " = " << c << std::endl;

    c=a-b;
    std::cout << "(" << a << ")" <<  " - " << "(" << b << ")" << " = " << c << std::endl;

    c=a*b;
    std::cout << "(" << a << ")" <<  " * " << "(" << b << ")" << " = " << c << std::endl;

    n=c.degree();
    std::cout << "n: " << n << std::endl;

    return(0);
} /// of main()
