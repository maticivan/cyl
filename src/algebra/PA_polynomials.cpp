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

#ifndef _INCL_PA_CPP
#define _INCL_PA_CPP
namespace PA{ 
    std::string printDiadicRational(const Frac& x){
        return x.toString();
    }
  template<typename TTT>
  Monomial<TTT>::Monomial(const TTT& coefficient, const std::string& variable, long exponent){
    c=coefficient;
    if(exponent<0){exponent=0;}
    if(exponent>0){ 
      vExp[variable]=exponent;
    }
  }
  template<typename TTT> 
  std::string Monomial<TTT>::toString(int skipCoefficient) const{
      
      std::string res;
      if(c==0){return res;}
      res+=printDiadicRational(c);
    if(skipCoefficient){res="";}
    std::map<std::string,long>::const_iterator it=vExp.begin(),itE=vExp.end();
    while(it!=itE){
      res+="\\cdot "+it->first+"^{"+std::to_string(it->second)+"}";
      ++it;
    }
    return res;
  }
  template<typename TTT>
  int Monomial<TTT>::operator<(const Monomial<TTT>& oth) const{
    if(vExp.size()<oth.vExp.size()){return 1;}
    if(vExp.size()>oth.vExp.size()){return 0;}
    std::map<std::string,long>::const_iterator itM=vExp.begin(),itO=oth.vExp.begin();
    while(itM!=vExp.end()){
      if(itM->first<itO->first){return 1;}
      if(itM->first>itO->first){return 0;}
      if(itM->second<itO->second){return 1;}
      if(itM->second>itO->second){return 0;}
      ++itM;++itO;
    }
    return 0;
  } 
  template<typename TTT>
  TTT power(TTT x, long e){
    TTT res=1;
    for(long i=0;i<e;++i){
      res*=x;
    }
    return res;
  }
  template<typename TTT>
  Monomial<TTT> Monomial<TTT>::evaluate(const std::map<std::string,TTT>& varVal) const{
    Monomial<TTT> zero,res;
    if(c==0){return zero;}
    res.c=c; 
    typename std::map<std::string,TTT>::const_iterator it=varVal.begin();
    std::map<std::string,long>::const_iterator itVEFinder;
    while(it!=varVal.end()){
      if(it->second==0){
        if(vExp.find(it->first)!=vExp.end()){
          return zero;
        }
      }
      else{
        itVEFinder=vExp.find(it->first);
        if(itVEFinder!=vExp.end()){
          res.c*=power(it->second,itVEFinder->second);
        }
      }
      ++it;
    }
    itVEFinder=vExp.begin();
    while(itVEFinder!=vExp.end()){
      if(varVal.find(itVEFinder->first)==varVal.end()){
        res.vExp[itVEFinder->first]=itVEFinder->second;
      }
      ++itVEFinder;
    }
    return res;
  }
  template<typename TTT>
  Monomial<TTT>& Monomial<TTT>::operator*=(const Monomial<TTT>& oth){
    c*=oth.c;
    if(c==0){
      vExp.clear();
    }
    else{
      std::map<std::string,long>::const_iterator it=oth.vExp.begin(),itThis;
      while(it!=oth.vExp.end()){
        itThis=vExp.find(it->first);
        if(itThis==vExp.end()){
          vExp[it->first]=it->second;
        }
        else{
          vExp[it->first]+=it->second;
        }
        ++it;
      }
    }
    return *this;
  }
  template<typename TTT>
  Monomial<TTT>& Monomial<TTT>::operator*=(const TTT& oth){
    Monomial<TTT> othM(oth);
    return operator*=(othM);
  }
  template<typename TTT>
  Monomial<TTT> operator*(const Monomial<TTT>& a, const Monomial<TTT>& b){
    Monomial<TTT> res=a;
    res*=b;
    return res;
  } 
  template<typename TTT> Polynomial<TTT>::Polynomial(TTT c){
    if(c.getNum()!=0){
      Monomial<TTT> m(c);
      terms.insert(m);
    }
  }
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator+=(const Polynomial<TTT>& oth){
    for(long i=0;i<(oth.terms).size();++i){
      long j=terms.find((oth.terms)[i]);
      if(j<0){
        terms.insert((oth.terms)[i]);
      }
      else{
        Monomial<TTT> tOld=terms[j];
        Monomial<TTT> tToAdd=(oth.terms)[i];
        terms.erase(tOld);
        tOld.c+=tToAdd.c;
        if(tOld.c.getNum()!=0){
          terms.insert(tOld);
        }
      }
    }
    return *this;
  }
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator*=(const TTT& m){
    if(m==0){
      terms.clear();
      return *this;
    }
    ssm::set<Monomial<TTT> > termsNew;
    for(long i=0;i<terms.size();++i){
      Monomial<TTT> t=terms[i];
      t.c *= m;
      termsNew.insert(t);
    }
    terms=std::move(termsNew);
    return *this;
  }
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator-=(const Polynomial<TTT>& oth){
    Polynomial<TTT> q=oth;
    q*=-1;
    (*this)+=q;
    return *this;
  }  
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator+=(const Monomial<TTT>& m){
    if(m.c==0){return *this;}
    Polynomial<TTT> q;
    q.terms.insert(m);
    return operator+=(q);
  }  
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator-=(const Monomial<TTT>& m){
    if(m.c==0){return *this;}
    Polynomial<TTT> q;
    q.terms.insert(m);
    return operator-=(q);
  }  
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator*=(const Monomial<TTT>& m){
    if(m.c==0){
      terms.clear();
      return *this;
    }
    ssm::set<Monomial<TTT> > termsNew;
    for(long i=0;i<terms.size();++i){
      Monomial<TTT> t=terms[i];
      t *= m;
      termsNew.insert(t); 
    }
    terms=std::move(termsNew);
    return *this; 
  }
  template<typename TTT>
  Polynomial<TTT>& Polynomial<TTT>::operator*=(const Polynomial<TTT>& oth){
    Polynomial<TTT> result;
    Polynomial<TTT> current;
    for(long i=0;i<oth.terms.size();++i){
      current.terms=terms;
      current*=oth.terms[i];
      result+=current;
    }
    terms=result.terms;
    return *this;
  }
  template<typename TTT>
  int Polynomial<TTT>::operator<(const Polynomial<TTT>& oth) const {
    if(terms.size()<oth.terms.size()){return 1;}
    if(terms.size()>oth.terms.size()){return 0;}
    for(long i=0;i<terms.size();++i){
      Monomial<TTT> tMI=terms[i];
      Monomial<TTT> othI=oth.terms[i];
      if(tMI<othI){return 1;}
      if(othI<tMI){return 0;}
    }
    return 0;
  } 
  template<typename TTT>
  Polynomial<TTT> Polynomial<TTT>::evaluate(const std::map<std::string,TTT>& varVal) const{
    Polynomial<TTT> res;
    for(long i=0;i<terms.size();++i){
      res+=terms[i].evaluate(varVal);
    }
    return res;
  }
  template<typename TTT>
  std::string Polynomial<TTT>::toString() const{
    std::string res;
    int printed=0;long charactersPrintedInLine=0;std::string newestAddition;
    for(long i=0;i<terms.size();++i){
      Monomial<TTT> t=terms[i];
      if(charactersPrintedInLine>100){
        res+="\\\\\n&\\quad ";
        charactersPrintedInLine=0;
      }
      if((printed)&&(t.c>0)){
        res+="+";
      }
      printed=1;
      newestAddition=t.toString();
      res+=newestAddition;
      charactersPrintedInLine+=newestAddition.size();
    }
    return res;
  } 
  template<typename TTT>
  long maxVarDegree(const PA::Monomial<TTT>& m){
    long mdeg=0;
    std::map<std::string,long>::const_iterator it=m.vExp.begin();
    while(it!=m.vExp.end()){
      if(it->second>mdeg){mdeg=it->second;}
      ++it;
    }
    return mdeg;
  }
  template<typename TTT>
  long maxVarDegree(const PA::Polynomial<TTT>& p){
    long mdeg=0;long cd;
    for(long i=0;i<p.terms.size();++i){
      cd=maxVarDegree(p.terms[i]);
      if(cd>mdeg){mdeg=cd;}
    }
    return mdeg;
  }
  template<typename TTT> ssm::set<std::string> allVariables(const PA::Monomial<TTT>& m){
    ssm::set<std::string> avr;
    std::map<std::string,long>::const_iterator it=m.vExp.begin();
    while(it!=m.vExp.end()){
      if(it->second>0){
        avr.insert(it->first);
      }
      ++it;
    }
    return avr;
  }
  template<typename TTT>
  ssm::set<std::string> allVariables(const PA::Polynomial<TTT>& p){
    ssm::set<std::string> avr;
    for(long i=0;i<p.terms.size();++i){
      avr+=allVariables(p.terms[i]);
    }
    return avr;
  }
  template<typename TTT>
  long degreeOf(const PA::Monomial<TTT>& m,const std::string& vN){
    std::map<std::string,long>::const_iterator it=m.vExp.find(vN);
    if(it==m.vExp.end()){return 0;}
    return it->second;
  }
  template<typename TTT>
  Polynomial<TTT> polynomialToPower(const Polynomial<TTT>& base, long exponent){
    Polynomial<TTT> res(1);
    for(long i=0;i<exponent;++i){
      res*=base;
    }
    return res;
  }
  template<typename TTT>
  Polynomial<TTT> substitute(const Monomial<TTT>& m, const std::map<std::string,Polynomial<TTT> > substituteMap){
    Polynomial<TTT> res(1);
    if(m.c==0){return res;}
    std::map<std::string,long>::const_iterator it=m.vExp.begin();
    typename std::map<std::string,Polynomial<TTT> >::const_iterator itS;
    while(it!=m.vExp.end()){ 
      itS=substituteMap.find(it->first);
      if(itS==substituteMap.end()){
        Monomial<TTT> tmp(1,it->first,it->second);
        res*=tmp;
      }
      else{
        res*=polynomialToPower(itS->second,it->second);
      }
      ++it;
    }
    res*=m.c;
    return res;
  }
  template<typename TTT>
  Polynomial<TTT> substitute(const Polynomial<TTT>& p, const std::map<std::string,Polynomial<TTT> > substituteMap){
    Polynomial<TTT> res;
    for(long i=0;i<p.terms.size();++i){ 
      res+=substitute(p.terms[i],substituteMap);
    }
    return res;
  }
  template<typename TTT>
  int equal(const Polynomial<TTT>& a, const Polynomial<TTT>& b){
    if(a.terms.size()!=b.terms.size()){return 0;}
    for(long i=0;i<a.terms.size();++i){
      Monomial<TTT> aM=a.terms[i];
      Monomial<TTT> bM=b.terms[i];
      if( !(aM.c==bM.c) ){return 0;}
      if(aM.vExp.size()!=bM.vExp.size()){return 0;}
      std::map<std::string,long>::const_iterator itA=aM.vExp.begin(),itB=bM.vExp.begin();
      while( itA!=aM.vExp.end() ){
        if(itA->first!=itB->first){return 0;}
        if(itA->second!=itB->second){return 0;}
        ++itA;++itB;
      }
    }
    return 1;
  }
  template<typename TTT>
  std::vector<std::vector<std::vector<TTT> > > polMatrix(const PA::Polynomial<TTT>& p){
    long mdeg=maxVarDegree(p);
    std::vector<TTT> p1;
    std::vector<std::vector<TTT> > p2;
    std::vector<std::vector<std::vector<TTT> > > p3;
    p1.resize(mdeg+1);
    p2.resize(mdeg+1);
    p3.resize(mdeg+1);
    for(long i=0;i<mdeg+1;++i){
      p2[i]=p1;
    }
    for(long i=0;i<mdeg+1;++i){
      p3[i]=p2;
    }
    ssm::set<std::string> avr=allVariables(p);
    long indSearch;
    for(long i=0;i<p.terms.size();++i){
      PA::Monomial<TTT> mi=p.terms[i];
      p3[degreeOf(mi,avr[0])][degreeOf(mi,avr[1])][degreeOf(mi,avr[2])]=mi.c;
    }
    return p3;
  }

}

#endif
