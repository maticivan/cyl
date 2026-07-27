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

#ifndef _INCL_PA_H
#define _INCL_PA_H
namespace PA{ 
  template<typename TTT>
  struct Monomial{
  public:
    TTT c;
    std::map<std::string, long> vExp;
    Monomial(const TTT& =0, const std::string& = "", long =0);
    std::string toString(int = 0) const;
    int operator<(const Monomial<TTT>& ) const;// lexicographic comparison so we can create sets and maps
    Monomial<TTT> evaluate(const std::map<std::string,TTT>& varVal) const;
    Monomial<TTT>& operator*=(const Monomial& );
    Monomial<TTT>& operator*=(const TTT& );
  };
  template<typename TTT>
  struct Polynomial{
  public:
    ssm::set<Monomial<TTT> > terms;
    Polynomial(TTT = 0); 
    Polynomial<TTT>& operator+=(const Polynomial<TTT>& );
    Polynomial<TTT>& operator*=(const TTT& );
    Polynomial<TTT>& operator*=(const Monomial<TTT>& );
    Polynomial<TTT>& operator*=(const Polynomial<TTT>& );
    Polynomial<TTT>& operator-=(const Polynomial<TTT>& );
    Polynomial<TTT>& operator+=(const Monomial<TTT>& );
    Polynomial<TTT>& operator-=(const Monomial<TTT>& );
    int operator<(const Polynomial<TTT>& )const; //total ordering: lexicographic tiebraking
    Polynomial<TTT> evaluate(const std::map<std::string,TTT>& ) const;
    std::string toString() const;
  };       
}

#endif
