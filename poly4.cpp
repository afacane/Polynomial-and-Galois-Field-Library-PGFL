/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019
///
/// poly4.cpp
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

std::ofstream ToF;

int main ()
{
    Polynomial g,f,r,q;

    ToF.open("poly_sample.txt");

    g=Polynomial({4,0,5,-7,1}); /// g = q*f + r
    f=Polynomial({3,0,1});
    r=deconv(g,f);
    q=deconv_q(g,f);

    ToF << "g: " << g << std::endl;
    ToF << "q: " << q << std::endl;
    ToF << "f: " << f << std::endl;
    ToF << "r: " << r << std::endl;
    ToF << "g-(q*f+r): " << g-(q*f+r) << std::endl;

    ToF.close();
    return(0);
} /// of main()
