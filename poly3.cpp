/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019
///
/// poly3.cpp
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
    int m,n;

    Polynomial a, b, c;

    a=Polynomial({1,1,0,1,0,1}); /// $x^5 + x^3 + x + 1$
    std::cout << "a: " << a << std::endl;

    b=Polynomial({-2,0,0,0,0,1}); /// $x^5-2$
    std::cout << "b: " << b << std::endl;

    c=a-b;
    std::cout << "(" << a << ")" <<  " - " << "(" << b << ")" << " = " << c << std::endl;

    n=c.degree();
    std::cout << "n: " << n << std::endl;

    m=c.size();
    std::cout << "m: " << m << std::endl;

    return(0);
} /// of main()
