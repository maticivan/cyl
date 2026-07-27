//*************************************************************************************************
//*************************************************************************************************
//* The MIT License (MIT)                                                                         *
//* Copyright (C) 2026 Ivan Matic and Rados Radoicic                                              *
//*                                                                                               *
//* Permission is hereby granted, free of charge, to any person obtaining a copy of this          *
//* software and associated documentation files (the "Software"), to deal in the Software         *
//* without restriction, including without limitation the rights to use, copy, modify, merge,     *
//* publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons    *
//* to whom the Software is furnished to do so, subject to the following conditions:              *
//*                                                                                               *
//* The above copyright notice and this permission notice shall be included in all copies or      *
//* substantial portions of the Software.                                                         *
//*                                                                                               *
//* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,           *
//* INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR      *
//* PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE     *
//* FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR          *
//* OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        *
//* DEALINGS IN THE SOFTWARE.                                                                     *
//*************************************************************************************************

#ifndef _INCL_FA_CPP
#define _INCL_FA_CPP
namespace FA{
int GL_OVERFLOW=0;
int GL_DIVBY0=0;
__int128 multiplyBigNumbers(__int128 a, __int128 b){
    __int128 c;
    if(__builtin_mul_overflow(a, b, &c)){ GL_OVERFLOW = 1; }
    return c;
}
__int128 addBigNumbers(__int128 a, __int128 b) {
    __int128 c;
    if (__builtin_add_overflow(a, b, &c)) { GL_OVERFLOW = 1; }
    return c;
}
__int128 subBigNumbers(__int128 a, __int128 b) {
    __int128 c;
    if (__builtin_sub_overflow(a, b, &c)) { GL_OVERFLOW = 1; }
    return c;
}
__int128 absBig(__int128 a){
    if(a<0){return subBigNumbers(0,a);}
    return a;
}
__int128 gcdBig(__int128 a, __int128 b){
    if( (a==0)&&(b==0)){GL_DIVBY0=1; return 1;}
    a=absBig(a);
    b=absBig(b);
    __int128 t;
    while(b != 0){
        t = a % b;
        a = b;
        b = t;
    }
    return a;
}
__int128 leastCommonMultiple(__int128 a, __int128 b){
    a=absBig(a);
    b=absBig(b);
    if(a == 0 || b == 0){ return 0; }
    return multiplyBigNumbers(a / gcdBig(a, b), b);
}
int discreteLog2(__int128 a){
    //returns -1 if a is not a positive power of two, otherwise the exponent
    if(a<1){return -1;}
    if(a==1){return 0;}
    int counter=0;
    while(a>1){
        if(a%2){return -1;}
        ++counter;
        a/=2;
    }
    return counter;
}
// Materializes 2^e. USED ONLY on the non-dyadic paths, never in dyadic mode.
__int128 powerOfTwo(int e){
    if(e<0){return 1;}
    if(e>126){GL_OVERFLOW=1;return 1;}
    __int128 res=1;
    for(int i=0;i<e;++i){res=addBigNumbers(res,res);}
    return res;
}
std::string writeAsPowerOfTwoIfPossible(__int128 a){
    int dl2=discreteLog2(a);
    if(dl2<0){return "";}
    return "2^{"+std::to_string(dl2)+"}";
}
std::string toString(__int128 a){
    if(a==0){return "0";}
    std::string res;
    std::string sgn;
    if(a<0){sgn+="-";a = subBigNumbers(0, a);}
    while(a>0){
        res=std::to_string(static_cast<int>(a%10))+res;
        a/=10;
    }
    return sgn+res;
}
}

class Frac{
    // REPRESENTATION (two modes, selected by log2Den):
    //
    //   log2Den >= 0  ("dyadic mode"):  value = num / 2^log2Den .
    //                 The field `den` is NOT PART OF THE VALUE and is never read.
    //                 It is held at 1 purely so that the object is deterministic.
    //                 2^log2Den is NEVER materialized, so log2Den may be arbitrarily
    //                 large (thousands) with no overflow.
    //                 Canonical form: num odd, or log2Den==0.  Zero is (0,1,0).
    //
    //   log2Den == -1 ("general mode"): value = num / den, with den >= 3 not a power
    //                 of two, gcd(|num|,den)==1, den > 0, sign carried by num.
    //
    // Any new code path that writes num/den/log2Den MUST restore this by calling
    // reduce().  No path may read `den` without first checking log2Den == -1.
private:
    __int128 num;
    __int128 den;
    int log2Den;
    void normalizeDyadic();
    void normalizeGeneral();
    int cmp(const Frac&) const;   // -1, 0, or +1
public:
    Frac(const __int128& =0, const __int128& =1);
    __int128 getNum() const;
    int isDyadic() const;
    int getLog2Den() const;
    __int128 getDenGeneral() const;   // valid only when !isDyadic()
    void setNum(const __int128 &);
    void setDen(const __int128 &);
    int operator>(const Frac& ) const;
    int operator<(const Frac& ) const;
    int operator==(const Frac& ) const;
    void reduce();
    Frac& operator+=(const Frac& );
    Frac& operator-=(const Frac& );
    Frac& operator*=(const Frac& );
    Frac& operator=(const __int128 &);
    Frac& operator=(const Frac &);
    std::string toString() const;
};

// ---------------------------------------------------------------- normalization

void Frac::normalizeDyadic(){
    // assumes log2Den >= 0; strips common factors of two WITHOUT touching den
    if(num==0){den=1;log2Den=0;return;}
    if(log2Den<0){log2Den=0;}
    __int128 a=FA::absBig(num);
    int neg=(num<0);
    while( (log2Den>0) && (a%2==0) ){
        a/=2;--log2Den;
    }
    num = neg ? FA::subBigNumbers(0,a) : a;
    den = 1;
}
void Frac::normalizeGeneral(){
    // assumes the value is num/den with den != 0; picks the right mode
    if(den==0){den=1;FA::GL_DIVBY0=1;}
    if(den<0){ num=FA::subBigNumbers(0,num); den=FA::subBigNumbers(0,den); }
    if(num==0){den=1;log2Den=0;return;}
    __int128 g=FA::gcdBig(num,den);
    if(g>1){ num/=g; den/=g; }
    int e=FA::discreteLog2(den);
    if(e>=0){
        log2Den=e;
        normalizeDyadic();     // drops into dyadic mode; den released
        return;
    }
    log2Den=-1;                // genuinely non-dyadic; den is meaningful
}
void Frac::reduce(){
    if(log2Den>-1){normalizeDyadic();}
    else{normalizeGeneral();}
}

// ---------------------------------------------------------------- construction

Frac::Frac(const __int128 & n, const __int128 &d){
    num=n;den=d;log2Den=-1;
    if(den==0){den=1;FA::GL_DIVBY0=1;}
    normalizeGeneral();
}
__int128 Frac::getNum() const{return num;}
int Frac::isDyadic() const{return (log2Den>-1);}
int Frac::getLog2Den() const{return log2Den;}
__int128 Frac::getDenGeneral() const{return den;}
void Frac::setNum(const __int128 &n){
    num=n;
    reduce();
}
void Frac::setDen(const __int128 &d){
    den=d;log2Den=-1;
    if(den==0){den=1;FA::GL_DIVBY0=1;}
    normalizeGeneral();
}
Frac& Frac::operator=(const __int128 & x){
    num=x;den=1;log2Den=0;
    return *this;
}
Frac& Frac::operator=(const Frac & b){
    num=b.num;den=b.den;log2Den=b.log2Den;
    return *this;
}

// ---------------------------------------------------------------- comparison
// All comparisons are overflow-free for dyadics of arbitrary exponent: instead of
// scaling one side by 2^delta, we double it step by step and stop the moment it
// exceeds the other side, at which point the answer is already determined.

int compareDyadicMagnitudes(__int128 a1, int m1, __int128 a2, int m2){
    // a1,a2 > 0; compares a1/2^m1 with a2/2^m2; returns -1,0,1
    while( (m1>m2) && (a2<=a1) ){ a2=FA::addBigNumbers(a2,a2); --m1; }
    while( (m2>m1) && (a1<=a2) ){ a1=FA::addBigNumbers(a1,a1); --m2; }
    if(m1>m2){return -1;}   // a2 outgrew a1 and still has scaling left
    if(m2>m1){return 1;}
    if(a1<a2){return -1;}
    if(a2<a1){return 1;}
    return 0;
}
int compareMixedMagnitudes(__int128 a1, int e, __int128 a2, __int128 d){
    // a1,a2>0, d>0; compares a1/2^e with a2/d  <=>  a1*d with a2*2^e
    __int128 lhs=FA::multiplyBigNumbers(a1,d);
    while( (e>0) && (a2<=lhs) ){ a2=FA::addBigNumbers(a2,a2); --e; }
    if(e>0){return -1;}
    if(lhs<a2){return -1;}
    if(a2<lhs){return 1;}
    return 0;
}
int Frac::cmp(const Frac& b) const{
    if( (num==0) && (b.num==0) ){return 0;}
    if(num==0){return (b.num>0)?-1:1;}
    if(b.num==0){return (num>0)?1:-1;}
    int s1=(num>0)?1:-1;
    int s2=(b.num>0)?1:-1;
    if(s1!=s2){return (s1>s2)?1:-1;}
    __int128 a1=FA::absBig(num);
    __int128 a2=FA::absBig(b.num);
    int c;
    if( (log2Den>-1) && (b.log2Den>-1) ){
        c=compareDyadicMagnitudes(a1,log2Den,a2,b.log2Den);
    }
    else if( (log2Den>-1) && (b.log2Den==-1) ){
        c=compareMixedMagnitudes(a1,log2Den,a2,b.den);
    }
    else if( (log2Den==-1) && (b.log2Den>-1) ){
        c=-compareMixedMagnitudes(a2,b.log2Den,a1,den);
    }
    else{
        __int128 l=FA::multiplyBigNumbers(a1,b.den);
        __int128 r=FA::multiplyBigNumbers(a2,den);
        c = (l<r) ? -1 : ((r<l)?1:0);
    }
    return (s1>0)?c:(-c);
}
int Frac::operator<(const Frac & b) const{ return (cmp(b)<0)?1:0; }
int Frac::operator>(const Frac & b) const{ return (cmp(b)>0)?1:0; }
int Frac::operator==(const Frac & b) const{
    // canonical forms are unique within each mode, and a canonical dyadic is
    // never equal to a canonical general value, so mode mismatch means unequal
    if(log2Den!=b.log2Den){return 0;}
    if(!(num==b.num)){return 0;}
    if(log2Den==-1){return (den==b.den)?1:0;}
    return 1;                     // dyadic mode: den is not part of the value
}

// ---------------------------------------------------------------- addition

Frac addDyadicAndGeneral(__int128 dyNum, int dyExp, __int128 gNum, __int128 gDen){
    // dyNum/2^dyExp + gNum/gDen, with gDen>0 not a power of two.
    // The result is generally non-dyadic, so its denominator 2^dyExp*gDen must be
    // formed; FA::powerOfTwo flags an overflow if dyExp is too large to allow it.
    __int128 p=FA::powerOfTwo(dyExp);
    __int128 numNew=FA::addBigNumbers(FA::multiplyBigNumbers(dyNum,gDen),
                                      FA::multiplyBigNumbers(gNum,p));
    __int128 denNew=FA::multiplyBigNumbers(p,gDen);
    return Frac(numNew,denNew);
}
Frac& Frac::operator+=(const Frac& b) {
    if(b.num==0){return *this;}
    if(num==0){ *this=b; return *this; }
    if( (log2Den>-1) && (b.log2Den>-1) ){
        // pure dyadic path: align exponents by doubling numerators only
        __int128 n1=num, n2=b.num;
        int m1=log2Den, m2=b.log2Den;
        int m=(m1>m2)?m1:m2;
        while(m1<m){++m1;n1=FA::addBigNumbers(n1,n1);}
        while(m2<m){++m2;n2=FA::addBigNumbers(n2,n2);}
        num=FA::addBigNumbers(n1,n2);
        log2Den=m;
        normalizeDyadic();       // no 2^m anywhere
        return *this;
    }
    if( (log2Den>-1) && (b.log2Den==-1) ){
        *this=addDyadicAndGeneral(num,log2Den,b.num,b.den);
        return *this;
    }
    if( (log2Den==-1) && (b.log2Den>-1) ){
        *this=addDyadicAndGeneral(b.num,b.log2Den,num,den);
        return *this;
    }
    __int128 lcm=FA::leastCommonMultiple(den,b.den);
    __int128 numNew=FA::addBigNumbers(FA::multiplyBigNumbers(num,lcm/den),
                                      FA::multiplyBigNumbers(b.num,lcm/b.den));
    num=numNew;
    den=lcm;
    normalizeGeneral();
    return *this;
}
Frac& Frac::operator-=(const Frac& b) {
    Frac mb=b;
    mb.num=FA::subBigNumbers(0,mb.num);
    return operator+=(mb);
}

// ---------------------------------------------------------------- multiplication

Frac& Frac::operator*=(const Frac& b) {
    if(b.num==0){num=0;den=1;log2Den=0;return *this;}
    if(num==0){den=1;log2Den=0;return *this;}
    if( (log2Den>-1) && (b.log2Den>-1) ){
        // pure dyadic path: exponents add, nothing is materialized
        int e=log2Den+b.log2Den;
        num=FA::multiplyBigNumbers(num,b.num);
        log2Den=e;
        normalizeDyadic();
        return *this;
    }
    if( (log2Den>-1) || (b.log2Den>-1) ){
        // one dyadic, one general: cancel twos against the general numerator first,
        // then the result's denominator 2^e * d has to be formed (never happens for
        // the deposited certificate, in which every quantity is dyadic)
        __int128 dyNum, gNum, gDen; int e;
        if(log2Den>-1){ dyNum=num; e=log2Den; gNum=b.num; gDen=b.den; }
        else          { dyNum=b.num; e=b.log2Den; gNum=num; gDen=den; }
        while( (e>0) && (gNum%2==0) && (gNum!=0) ){ gNum/=2; --e; }
        __int128 g=FA::gcdBig(dyNum,gDen);
        if(g>1){ dyNum/=g; gDen/=g; }
        __int128 p=FA::powerOfTwo(e);
        num=FA::multiplyBigNumbers(dyNum,gNum);
        den=FA::multiplyBigNumbers(p,gDen);
        log2Den=-1;
        normalizeGeneral();
        return *this;
    }
    __int128 gcdDenBNum=FA::gcdBig(den,b.num);
    __int128 gcdNumBDen=FA::gcdBig(num,b.den);
    __int128 numNew=FA::multiplyBigNumbers(num/gcdNumBDen,b.num/gcdDenBNum);
    __int128 denNew=FA::multiplyBigNumbers(den/gcdDenBNum,b.den/gcdNumBDen);
    num=numNew;
    den=denNew;
    log2Den=-1;
    normalizeGeneral();
    return *this;
}
Frac operator+(Frac a, const Frac& b){ a+=b; return a; }
Frac operator-(Frac a, const Frac& b){ a-=b; return a; }
Frac operator*(Frac a, const Frac& b){ a*=b; return a; }

// ---------------------------------------------------------------- output

std::string Frac::toString() const{
    if(num==0){return "0";}
    std::string sgn="";
    __int128 absNum=FA::absBig(num);
    if(num<0){sgn="-";}
    std::string numStr=FA::toString(absNum);
    if(log2Den==0){return sgn+numStr;}
    if(log2Den>0){
        return sgn+"\\frac{"+numStr+"}{2^{"+std::to_string(log2Den)+"}}";
    }
    if(den==1){return sgn+numStr;}
    std::string denStr=FA::writeAsPowerOfTwoIfPossible(den);
    if(denStr==""){denStr=FA::toString(den);}
    return sgn+"\\frac{"+numStr+"}{"+denStr+"}";
}
#endif
