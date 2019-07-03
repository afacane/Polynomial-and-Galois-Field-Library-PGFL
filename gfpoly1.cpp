/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019
///
/// gfpoly1.cpp
/// Sample polynomial operations over GF(p)

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
    std::ofstream ToF;

    int p;
    Polynomial g,f,r,q;

    ToF.open("gfpoly1.txt");

    p=11; /// Galois field GF(11)

    g=Polynomial({2,5,0,0,3}); /// g = q*f + r /// g = 3 x^4 + 5x + 2
    f=Polynomial({5,0,1,2}); /// f = 2x^3+ x^2 + 5
    r=gfdeconv(g,f,p);
    q=gfdeconv_q(g,f,p);

    ToF << "g: " << g << std::endl;
    ToF << "q: " << q << std::endl;
    ToF << "f: " << f << std::endl;
    ToF << "r: " << r << std::endl;
    ToF << "g-(q*f+r): " << poly_sub( g, poly_add( poly_mul(q,f,p),r,p ), p ) << std::endl;

    ToF.close();
    return(0);
} /// of main()
