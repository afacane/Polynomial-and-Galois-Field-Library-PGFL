/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019

#ifndef Poly_Included
#define Poly_Included

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

struct Term
{
    int poly;
    unsigned exp;
    Term(int c, unsigned e) : poly(c), exp(e) {};
};

/// Term({1,0}) --> 1
/// Term({1,1}) --> x
/// Term({1,2}) --> x^2
/// Term({1,3}) --> x^3
/// Term({2,5}) --> 2x^5

class Polynomial
{
    friend bool operator==(Polynomial lhs, Polynomial rhs);

    public:
        Polynomial();
        Polynomial(const Polynomial&);
        Polynomial(const std::vector<int>& c) : poly(c) {};
        Polynomial(const Term t);
        ~Polynomial();

        int size();
        int degree();
        int iszero();
        void zero();
        void resize(int n);
        std::string print() const;

        Polynomial& operator=(const Polynomial& other);
        Polynomial& operator*=(const Polynomial& other);
        Polynomial& operator+=(const Polynomial& other);
        Polynomial& operator-=(const Polynomial& other);

        Polynomial reduction(int p);
        int Get_Term(int n);

    private:
        std::vector<int> poly;
};

std::ostream& operator<<(std::ostream& lhs, const Polynomial& rhs);

Polynomial poly_mul(Polynomial f, Polynomial g, int p);
Polynomial poly_add(Polynomial f, Polynomial g, int p);
Polynomial poly_sub(Polynomial f, Polynomial g, int p);
const Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs);
const Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs);
const Polynomial operator-(const Polynomial& lhs, const Polynomial& rhs);
int gcd(int x, int y);
int phi(int i);
int isprime(int num);
void primefactor(int number,int *p,int *k);
int Zp(int a, int p);
int power(int x, int n, int p);

Polynomial deconv(Polynomial a, Polynomial b);
Polynomial deconv_q(Polynomial a, Polynomial b);
Polynomial gfdeconv(Polynomial f, Polynomial g, int p);
Polynomial gfdeconv_q(Polynomial a, Polynomial b, int p);
int gfprimck(Polynomial a, int p);

#endif
