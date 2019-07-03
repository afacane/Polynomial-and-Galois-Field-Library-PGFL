/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019
///
/// poly1.cpp
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
    Polynomial a;

    a=Polynomial({1,2,3}); /// 3 x^2 + 2 x + 1
    std::cout << "a: " << a << std::endl;

    a=Polynomial({0,0,1}); /// x^2
    std::cout << "a: " << a << std::endl;

    a=(Term(1, 2)); /// x^2
    std::cout << "a: " << a << std::endl;

    a=Polynomial({0}); /// 0
    std::cout << "a: " << a << std::endl;

    a=Polynomial({1}); /// 1
    std::cout << "a: " << a << std::endl;

    /// Does not work
    /// a=Polynomial({0,1}); /// x
    /// std::cout << "a: " << a << std::endl;

    a=(Term(1, 1)); /// x
    std::cout << "a: " << a << std::endl;

    /// Does not work
    /// a=Polynomial({1,1}); /// x+1
    /// std::cout << "a: " << a << std::endl;

    a=Polynomial({1})+(Term(1, 1)); /// x+1
    std::cout << "a: " << a << std::endl;

    a=Polynomial({1,1,1}); /// x^2 + x + 1
    std::cout << "a: " << a << std::endl;

    a=Polynomial({1,0,0,2,0,3}); /// 3 x^5+ 2 x^3 + 1
    std::cout << "a: " << a << std::endl;

    a=(Term(5, 3)); /// 5 x^3
    std::cout << "a: " << a << std::endl;

    a=Polynomial({0,0,0,5}); /// 5 x^3
    std::cout << "a: " << a << std::endl;

    return(0);
} /// of main()
