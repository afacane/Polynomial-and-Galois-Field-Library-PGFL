/// Polynomial and Galois Field Library
///
/// Erkan Afacan, Gazi University
///
/// 03/07/2019

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "poly.h"

///---------------------------------------------------------------------------------------

Polynomial::Polynomial() {
    this->poly.push_back(0);
}

Polynomial::Polynomial(const Polynomial& polynomial) {
    this->poly = polynomial.poly;
}

Polynomial::Polynomial(const Term t) : poly(t.exp + 1)
{
    poly[t.exp] = t.poly;
}

Polynomial::~Polynomial() {
    this->poly.clear();
}

int Polynomial::size()
{
    return poly.size();
}

int Polynomial::degree()
{
    int n;

    for (n = poly.size()-1; n >= 0; n--) {
        if (poly[n] != 0) {
            break;
        }
    }
    return(n);
}

int Polynomial::iszero()
{
    int i, iszero=1;

    for (i = poly.size()-1; i>=0; --i) {
        if (poly[i] != 0) {
            iszero=0;
            break;
        }
    }
    return(iszero);
}

void Polynomial::zero()
{
    for (int i = poly.size()-1; i>=0; --i) {
        poly[i]=0;
    }
}

void Polynomial::resize(int n)
{
    poly.resize(n);
}

std::string Polynomial::print() const
{
    std::string temp;
    signed long long i;
    int iszero=1;

    for (i = poly.size()-1; i>=0; --i) {
        if (poly[i] != 0) {
            iszero=0;
            break;
        }
    }

    if (!iszero) {
        for (i = poly.size() - 1; i >= 0; --i) {
            if (poly[i] == 0)
                continue;
            if (temp != "")
                temp += " + ";
            if (poly[i] != 1 || i == 0)
                temp += patch::to_string(poly[i]);
            if (i != 0) {
                if (i==1)
                    temp += "x";
                else
                    temp += "x^" + patch::to_string(i);
            }
        }
    }
    else {
        temp="0";
    }
    return temp;
}

Polynomial& Polynomial::operator=(const Polynomial& other)
{
     if (this == &other) {
        return *this;
    }
    this->poly= other.poly;
    return *this;
}

Polynomial& Polynomial::operator*=(const Polynomial& other)
{
    int i,j,n,term;
    Polynomial h;

    n=poly.size()+other.poly.size();

    h.poly.resize(n-1); /// h[0,1,...,n-2]
    for (i=0;i<h.size();i++) {
        h.poly[i]=0;
    }

    for (i=0;i<poly.size();i++) {
        term=0;
        for (j=0;j<other.poly.size();j++) {
            term = poly[i]*other.poly[j];
            h.poly[i+j] += term;
        }
    }

    poly.resize(n-1);
    for (i=0;i<h.size();i++) {
        poly[i]=0;
    }

    for (i=0;i<h.size();i++) {
        poly[i]=h.poly[i];
    }

    return *this;
}

Polynomial& Polynomial::operator+=(const Polynomial& other)
{
    if (poly.size() < other.poly.size())
        poly.resize(other.poly.size());
    for (unsigned i = 0; i < std::min(poly.size(), other.poly.size()); ++i) {
        poly[i] += other.poly[i];
    }
    return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& other)
{
    if (poly.size() < other.poly.size())
        poly.resize(other.poly.size());
    for (unsigned i = 0; i < std::min(poly.size(), other.poly.size()); ++i) {
        poly[i] -= other.poly[i];
    }
    return *this;
}

bool operator==(Polynomial lhs, Polynomial rhs)
{
    if (lhs.degree() != rhs.degree()) {
        return false;
    }

    for (int i = 0; i <= std::min(lhs.degree(),rhs.degree()); i++) {
        if (lhs.poly[i] != rhs.poly[i]) {
            return false;
        }
    }
    return true;
}

Polynomial Polynomial::reduction(int p)
{
    for (int i = poly.size()-1; i>=0; --i) {
        poly[i]= Zp(poly[i],p);
    }
    return poly;
}

int Polynomial::Get_Term(int n)
/// Term x^n
{
    return poly[n];
}

const Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs)
{
    Polynomial tmp(lhs);
    tmp *= rhs;
    return tmp;
}

const Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs)
{
    Polynomial tmp(lhs);
    tmp += rhs;
    return tmp;
}

const Polynomial operator-(const Polynomial& lhs, const Polynomial& rhs)
{
    Polynomial tmp(lhs);
    tmp -= rhs;
    return tmp;
}

std::ostream& operator<<(std::ostream& lhs, const Polynomial& rhs)
{
    return lhs << rhs.print();
}

///---------------------------------------------------------------------------------------

Polynomial poly_mul(Polynomial f, Polynomial g, int p)
/// h(x)=f(x)*g(x) % p
{
    int n;
    Polynomial h;

    h.resize(f.size()+g.size()-1);

    h=f*g;
    h.reduction(p);

    n=h.degree();
    h.resize(n+1);

    return(h);
} /// of Polynomial poly_mul()

Polynomial poly_add(Polynomial a, Polynomial b, int p)
/// c(x) = a(x) + b(x) % p
{
    Polynomial c;

    c=a+b;
    c.reduction(p);

    return(c);
} /// of Polynomial poly_add()

Polynomial poly_sub(Polynomial a, Polynomial b, int p)
/// c(x) = a(x) - b(x) % p
{
    Polynomial c;

    c=a-b;
    c.reduction(p);

    return(c);
} /// of Polynomial poly_sub()

int phi(int i)
{
	int res; /// Result
	int j;

	if (i==1) return 1;

    res=i;

    /// Check for divisibility by every prime number below the square root.
    /// Start with 2.
    if (i%2==0) {
        res-=res/2;
		do i/=2; while (i%2==0) ;
    }

    /// Since this doesn't use a list of primes, check every odd number.
    /// Ideally, skip past composite numbers.
	for (j=3; j*j<=i; j+=2) {
		if (i%j==0) {
			res-=res/j;
			do i/=j; while (i%j==0) ;
		}
	}

    /// If i>1, then it's the last factor at this point.
	if (i>1) res-=res/i;

    /// Return the result.
	return(res);
} /// of int phi()

int isprime(int num)
/// If the number is prime returns 1,
/// If the number is nonprime returns 0
{
   int j,flag=1;
   for (j=2;j<=num/2;++j){
        if (num%j==0){
            flag=0;
            break;
        }
   }
   return flag;
} /// int isprime(int num)

void primefactor(int number,int *p,int *k)
/// Factorizes the number and returns it as p^k.
{
    int div = 2;
    int i = 2;

    while (number!=0) {
        if (number%div!=0)
            div = div + 1;
        else {
            number = number / div;
            i++;
            if (number==1)
                break;
        }
    }
    *p = div;
    *k = i-2;
} /// of void primefactor()

int Zp(int a, int p)
/// Returns a positive representative for a between [0,..,p-1]
{
    int n;

    n = a%p;
    while (n < 0) {
        n += p;
    }

    return(n);
} /// of inf Zp()

int power(int x, int n, int p)
/// Calculates x^n (mod p)
{
    int tmp;

    if (n == 0) return 1;
    if (n % 2 == 0) {
        tmp = power(x, n >> 1, p);
        return(tmp * tmp % p);
    }
    else {
        return(x * power(x, n - 1, p) % p);
    }
} /// int power()

int gcd(int m,int n)
/// m!=0 or n!=0
/// gcd(a,b)=gcd(b,a mod b)
{
    int a,b,c;
    int d; // d=gcd(a,b)
    int q,r;

    a=m;
    b=n;
    q=a/b;
    r=a-q*b;

    if (b==0) {
        d=abs(a);
    }
    else {
        while (b!=0) {
            c=b;
            b=a%b;
            a=c;
            if (b!=0) {
                q=c/b;
                r=a-q*b;
            }
        }
        d=a;
    }

    return(d);
} /// of int gcd()

int euclid(int m,int n,int *x,int *y)
/// m!=0 or n!=0.
/// Returns x,y,d where d=gcd(m,n)=mx+ny.
{
    int a,b;
    int d0,d1,d2,x0,x1,x2,y0,y1,y2,q;
    int d;

    a=m;
    b=n;

    d0=a;
    x0=1;
    y0=0;
    d1=b;
    x1=0;
    y1=1;
    while (d1!=0) {
        q=d0/d1;
        d2=d1;
        x2=x1;
        y2=y1;
        d1=d0-q*d1;
        x1=x0-q*x1;
        y1=y0-q*y1;
        d0=d2;
        x0=x2;
        y0=y2;
    }
    d=d0;
    *x=x0;
    *y=y0;
    return(d);
} /// of int euclid()

int mod_inverse(int a, int m)
/// m>0
/// Returns a^{-1} mod m or message that gcd(a,m)>1.
{
    int d0,d1,d2,x0,x1,x2,q;

    d0=a;
    d1=m;
    x0=1;
    x1=0;
    while (d1!=0) {
        q=d0/d1;
        d2=d0-q*d1;
        x2=x0-q*x1;
        x0=x1;
        x1=x2;
        d0=d1;
        d1=d2;
    }
    if (d0==1) {
        while (x0<0) {
            x0+=m;
        }
        return(x0); /// a^{-1} = x0 mod m
    }
    else {
        printf("Inverse does not exist\n");
    }
} /// of int mod_inverse()

Polynomial deconv(Polynomial a, Polynomial b)
/// a(x) = b(x)*q(x) + r(x). It returns the remaining part r(x).
{
    int i,j,n,t;
    double term1;
    Polynomial r,temp,q,c;
    int m;

    r=a;

    i=r.degree();
    j=b.degree();

    /// if the degree of a.poly[] is smaller than b.poly[]
    if (i<j) {
        return(r);
    }
    else {
        temp.resize(a.size());
        q.resize(a.size());
        c.resize(a.size());
        q.zero();

        for (n = a.degree();n > 0; n--) {
            if ((r.Get_Term(n) != 0) && (n>=j)) {
                temp.zero(); c.zero();
                m=r.Get_Term(n);
                term1 = (r.Get_Term(n))/(b.Get_Term(j));
                t = floor(term1);
                q += Term({t,n-j}); /// t x^{n-j}
                temp += Term({t,n-j});
                c = temp*b;
                r = r - c;
            }
        }
        /// quotient q(x)
        return(r); /// remainder r(x)
    }
} /// of Polynomial deconv()

Polynomial deconv_q(Polynomial a, Polynomial b)
/// a(x) = b(x)*q(x) + r(x). It returns the quotient q(x).
{
    int i,j,n,t;
    double term1;
    Polynomial r,temp,q,c;
    int m;

    r=a;

    i=r.degree();
    j=b.degree();

    /// if the degree of a.poly[] is smaller than b.poly[]
    if (i<j) {
        return(r);
    }
    else {
        temp.resize(a.size());
        q.resize(a.size());
        c.resize(a.size());
        q.zero();

        for (n = a.degree();n > 0; n--) {
            if ((r.Get_Term(n) != 0) && (n>=j)) {
                temp.zero(); c.zero();
                m=r.Get_Term(n);
                term1 = (r.Get_Term(n))/(b.Get_Term(j));
                t = floor(term1);
                q += Term({t,n-j}); /// t x^{n-j}
                temp += Term({t,n-j});
                c = temp*b;
                r = r - c;
            }
        }
        /// remainder r(x)
        return(q);  /// quotient q(x)
    }
} /// of Polynomial deconv_q()

Polynomial gfdeconv(Polynomial a, Polynomial b, int p)
/// a(x) = b(x)*q(x) + r(x) in GF(p). It returns the remaining part r(x).
{
    int i,j,n,t;
    int term1;
    Polynomial r,temp,q,c;

    r=a;

    i=r.degree();
    j=b.degree();

    /// if the degree of a.poly[] is smaller than b.poly[]
    if (i<j) {
        return(r);
    }
    else {
        temp.resize(a.size());
        q.resize(a.size());
        c.resize(a.size());
        q.zero();

        /// The inverse of the first (highest) coefficient of b.
        term1 = mod_inverse(b.Get_Term(j),p); /// term1*b.Get_Term(j) % p = 1

        for (n = a.degree();n > 0; n--) {
            if ((r.Get_Term(n) != 0) && (n>=j)) {
                t = r.Get_Term(n)*term1;
                t = Zp(t,p);
                temp = Term({t,n-j}); /// t x^{n-j}
                q += temp;
                c = poly_mul(temp, b, p);
                r = poly_sub(r, c, p);
            }
        }
        /// quotient q(x)
        return(r); /// remainder r(x)
    }
} /// of Polynomial gfdeconv()

Polynomial gfdeconv_q(Polynomial a, Polynomial b, int p)
/// a(x) = b(x)*q(x) + r(x) in GF(p). It returns the quotient q(x).
{
    int i,j,n,t;
    double term1;
    Polynomial r,temp,q,c;
    int m;

    r=a;

    i=r.degree();
    j=b.degree();

    /// if the degree of a.poly[] is smaller than b.poly[]
    if (i<j) {
        return(r);
    }
    else {
        temp.resize(a.size());
        q.resize(a.size());
        c.resize(a.size());
        q.zero();

        for (n = a.degree();n > 0; n--) {
            if ((r.Get_Term(n) != 0) && (n>=j)) {
                temp.zero(); c.zero();
                m=r.Get_Term(n);
                r = r - Term({m,n}); /// m x^n
                while (m < b.Get_Term(j)) {
                    m = m + p;
                }
                r = r + Term({m,n}); /// m x^n
                term1 = mod_inverse(b.Get_Term(j),p); /// term1*b.Get_Term(j) % p = 1
                t = r.Get_Term(n)*term1;
                t = Zp(t,p);
                q += Term({t,n-j}); /// t x^{n-j}
                temp += Term({t,n-j});
                c = poly_mul(temp, b, p);
                r = poly_sub(r, c, p);
            }
        }
        /// remainder r(x)
        return(q);  /// quotient q(x)
    }
} /// of Polynomial gfdeconv_q()

int gfprimck(Polynomial a, int p)
/// From Matlab: function ck = gfprimck(a, p)
///
/// GFPRIMCK Check whether a polynomial over a Galois field is primitive.
///
///     CK = GFPRIMCK(A, P) checks whether the degree-M GF(P) polynomial
///     A is a primitive polynomial for GF(P^M). P is a prime number.
///     The output CK is as follows:
///         CK = -1   A is not an irreducible polynomial;
///         CK =  0   A is irreducible but not a primitive polynomial;
///         CK =  1   A is a primitive polynomial.
{
    Polynomial test_poly,r;
    int i,j,m,n,tmp,y;
    int ck,test_ord;

    /// Allocate space for the result, assume primitive.
    ck = 1;
    m=a.degree();

    /// The polynomial is divisible by x, hence is reducible.
    /// The only exception is when the polynomial is x ...
    if (a.Get_Term(0) == 0) {
        if (m == 1)
            ck = 0;
        else
            ck = -1;
    }
    /// This polynomial is actually a constant.
    else if (m == 0)
        ck = 1;
    /// The typical case.
    else {
        /// First test if the current polynomial is irreducible.
        n = (int)floor(pow(p,m))-p-1;
        if (n<=0) n=1;

        /// 'test_dec' is a vector containing the decimal(scalar) representations of
        /// the polynomials that could be divisors of 'a'.

        int test_dec[n];
        for (i=0;i<n;i++) {
            test_dec[i] = i+p+1;
        }

        test_poly.resize(a.size());
        r.resize(a.size());

        i=0;
        /// Loop through all polynomials that could be divisors of 'a'.
        while (i<n) {
            /// Expand the scalar value to a polynomial in GF(p).
            test_poly.zero();
            r.zero();
            tmp = test_dec[i];
            for (j=0;j<m;j++) {
                test_poly += Term(tmp%p,j);
                tmp = floor(tmp/p);
            }
            r = gfdeconv(a,test_poly,p);
            if (r.iszero()==1) {
                ck = -1;
                break;
            }
            i++;
        }
        if (ck == 1) {
            /// If the current polynomial is irreducible then check if it is primitive.
            /// To be primitive, the polynomial must not be a factor of another
            /// polynomial of the form x^n + 1 for any value of n in the range
            /// m < n < p^m - 1
            /// To check for this we check to see if the polynomial divides x^n
            /// with a remainder of 1 for all values of n in this range.

            test_poly.zero();
            test_poly = Term({1,m}); /// x^m

            y = pow(p,m)-1;
            test_ord = m;
            while (test_ord<y) {
                r.zero();

                /// calculate the remainder
                r = gfdeconv(test_poly,a,p);

                if ((r.Get_Term(0)==1) && (r.degree()==0)) { /// if (r == 1)
                    /// If we find a value of n in this range for which the remainder is
                    /// 1, we can then conclude the test and declare that the polynomial
                    /// is not primitive.
                    ck = 0;
                    break;
                }
                else {
                    /// To reduce the computational load, on each successive test we
                    /// simply need to test against the remainder of the previous test
                    /// multiplied by x (i.e., a shifted version of the previous remainder).
                    test_poly.zero();
                    for (i=0;i<r.degree()+1;i++) {
                        test_poly += Term(r.Get_Term(i),i+1); /// r.Get_Term(i) x^{i+1}
                    }
                    test_ord++;
                }
            }
        }
    }
    return(ck);
} /// of int gfprimck()
