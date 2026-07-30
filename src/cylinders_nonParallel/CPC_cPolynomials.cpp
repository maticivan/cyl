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


#ifndef _INCL_CPC_CPP
#define _INCL_CPC_CPP
namespace CPC{
    long GL_skipMultVerification=0;
    long GL_skipOrthVerification=0;
    long GL_stringSizeLimit=1000;
    std::string GL_openTag="{";
    std::string GL_closeTag="}";
    template<typename TT>
    std::vector<TT> getColumn(const std::vector<std::vector<TT> > & A, long j){
        std::vector<TT> res;
        long numRows=A.size();
        if(numRows<1){return res;}
        long numCols=A[0].size();
        if(numCols<1){return res;}
        if((j<0)||(j>=numCols)){return res;}
        res.resize(numRows);
        for(long i=0;i<numRows;++i){
            res[i]=A[i][j];
        }
        return res;
    }
    template<typename TT>
    std::vector<std::vector<TT> > transpose(const std::vector<std::vector<TT> >& A){
        std::vector<std::vector<TT> > AT;
        long numRows=A.size();
        if(numRows<1){return AT;}
        long numCols=A[0].size();
        if(numCols<1){return AT;}
        AT.resize(numCols);
        for(long j=0;j<numCols;++j){
            AT[j]=getColumn(A,j);
        }
        return AT;
    }
std::vector<std::vector<PA::Polynomial<Frac> > > createR(){
    std::vector<std::vector<PA::Polynomial<Frac> > > res;
    res.resize(3);
    std::vector<PA::Polynomial<Frac> > row;
    row.resize(3);
    PA::Monomial<Frac> q0(1,"q_0",1);
    PA::Monomial<Frac> q1(1,"q_1",1);
    PA::Monomial<Frac> q2(1,"q_2",1);
    PA::Monomial<Frac> q3(1,"q_3",1);
    PA::Monomial<Frac> q0Sq=q0; q0Sq*=q0;
    PA::Monomial<Frac> q1Sq=q1; q1Sq*=q1;
    PA::Monomial<Frac> q2Sq=q2; q2Sq*=q2;
    PA::Monomial<Frac> q3Sq=q3; q3Sq*=q3;
    PA::Monomial<Frac> q0q1=q0; q0q1*=q1;
    PA::Monomial<Frac> q0q2=q0; q0q2*=q2;
    PA::Monomial<Frac> q0q3=q0; q0q3*=q3;
    PA::Monomial<Frac> q1q2=q1; q1q2*=q2;
    PA::Monomial<Frac> q1q3=q1; q1q3*=q3;
    PA::Monomial<Frac> q2q3=q2; q2q3*=q3;
    PA::Polynomial<Frac> tmp;
    tmp.terms.clear();tmp+=q0Sq;tmp+=q1Sq;tmp-=q2Sq;tmp-=q3Sq;row[0]=tmp;
    tmp.terms.clear();tmp+=q1q2;tmp-=q0q3;tmp*=2;row[1]=tmp;
    tmp.terms.clear();tmp+=q1q3;tmp+=q0q2;tmp*=2;row[2]=tmp;
    res[0]=row;
    tmp.terms.clear();tmp+=q1q2;tmp+=q0q3;tmp*=2;row[0]=tmp;
    tmp.terms.clear();tmp+=q0Sq;tmp-=q1Sq;tmp+=q2Sq;tmp-=q3Sq;row[1]=tmp;
    tmp.terms.clear();tmp+=q2q3;tmp-=q0q1;tmp*=2;row[2]=tmp;
    res[1]=row;
    tmp.terms.clear();tmp+=q1q3;tmp-=q0q2;tmp*=2;row[0]=tmp;
    tmp.terms.clear();tmp+=q2q3;tmp+=q0q1;tmp*=2;row[1]=tmp;
    tmp.terms.clear();tmp+=q0Sq;tmp-=q1Sq;tmp-=q2Sq;tmp+=q3Sq;row[2]=tmp;
    res[2]=row;
    return res;
}
template<typename TT>
TT dotProduct(const std::vector<TT>& a,const std::vector<TT>& b){
    TT defaultT;
    if((a.size()<1)||(a.size()!=b.size())){return defaultT;}
    TT result=a[0];
    result*=b[0];
    for(long i=1;i<a.size();++i){
        TT intermediate=a[i];
        intermediate*=b[i];
        result+=intermediate;
    }
    return result;
}
template<typename TT>
std::vector<TT> crossProduct(const std::vector<TT>& a, const std::vector<TT>& b){
    std::vector<TT> res;
    if( (a.size()!=3)||(b.size()!=3) ){return res;}
    res.resize(3);
    {
        res[0]=a[1]; res[0]*=b[2];
        TT oth=a[2]; oth*=b[1];
        res[0]-=oth;
    }
    {
        res[1]=a[2]; res[1]*=b[0];
        TT oth=a[0]; oth*=b[2];
        res[1]-=oth;
    }
    {
        res[2]=a[0]; res[2]*=b[1];
        TT oth=a[1]; oth*=b[0];
        res[2]-=oth;
    }
    return res;
}
}

template<typename TT>
std::vector<TT > operator*(const std::vector<std::vector<TT> >& M, const std::vector<TT>& col){
    std::vector<TT> result;
    long numRows=M.size();
    if(numRows<1){return result;}
    long numCols=M[0].size();
    if(numCols<1){return result;}
    if(numCols!=col.size()){return result;}
    result.resize(numRows);
    for(long i=0;i<numRows;++i){
        result[i]=CPC::dotProduct(M[i],col);
    }
    return result;
}
template<typename TT>
std::vector<TT> operator-(const std::vector<TT>& a, const std::vector<TT>& b){
    std::vector<TT> res;
    if((a.size()<1)||(a.size()!=b.size())){return res;}
    res=a;
    for(long i=0;i<b.size();++i){res[i]-=b[i];}
    return res;
}
template<typename TT>
std::vector<std::vector<TT > >operator*(const std::vector<std::vector<TT> >& A,
                                        const std::vector<std::vector<TT> >& B){
    std::vector<std::vector<TT > > ABT;
    long numRowsA=A.size();
    if(numRowsA<1){return ABT;}
    long numColsA=A[0].size();
    if(numColsA<1){return ABT;}
    long numRowsB=B.size();
    if(numColsA!=numRowsB){return ABT;}
    long numColsB=B[0].size();
    ABT.resize(numColsB);
    for(long i=0;i<numColsB;++i){
        ABT[i]=A*CPC::getColumn(B,i);
    }
    return CPC::transpose(ABT);
}

namespace CPC{
    PA::Polynomial<Frac> createQSquare(){
        PA::Polynomial<Frac> res;
        for(long i=0;i<4;++i){
            PA::Monomial<Frac> tmp(1,"q_"+std::to_string(i),1);
            tmp*=tmp;
            res+=tmp;
        }
        return res;
    }
    struct Components{
    public:
        std::vector<PA::Polynomial<Frac> > paper_u0;
        std::vector<PA::Polynomial<Frac> > paper_v0;
        std::vector<std::vector<PA::Polynomial<Frac> > > paper_R;
        std::vector<PA::Polynomial<Frac> > paper_u;
        std::vector<PA::Polynomial<Frac> > paper_gv0;
        std::vector<PA::Polynomial<Frac> > paper_w;
        std::vector<PA::Polynomial<Frac> > paper_n;
        PA::Polynomial<Frac> paper_N2;
        PA::Polynomial<Frac> paper_ww;
        PA::Polynomial<Frac> paper_wv0;
        PA::Polynomial<Frac> paper_wgv0;
        PA::Polynomial<Frac> paper_v0gv0;
        PA::Polynomial<Frac> paper_G6;
        PA::Polynomial<Frac> paper_U1Sq;
        PA::Polynomial<Frac> paper_Q;
        PA::Polynomial<Frac> paper_Uh;
        PA::Polynomial<Frac> paper_Qh;
        Components();
    } GL_C;
    Components::Components(){
        PA::Polynomial<Frac> pol_0(0);
        PA::Polynomial<Frac> pol_1(1);
        paper_u0.resize(3); paper_u0[0]=pol_0;paper_u0[1]=pol_0;paper_u0[2]=pol_1;
        paper_v0.resize(3); paper_v0[0]=pol_1;paper_v0[1]=pol_0;paper_v0[2]=pol_0;
        paper_R=CPC::createR();
        paper_u=paper_R*paper_u0;
        paper_gv0=paper_R*paper_v0;
        paper_w=paper_u-paper_u0;
        for(long i=0;i<paper_w.size();++i){paper_w[i]*=2;}
        paper_n=CPC::crossProduct(paper_v0,paper_gv0);
        paper_N2=CPC::dotProduct(paper_n,paper_n);
        paper_ww=CPC::dotProduct(paper_w,paper_w);
        paper_wv0=CPC::dotProduct(paper_w,paper_v0);
        paper_wgv0=CPC::dotProduct(paper_w,paper_gv0);
        paper_v0gv0=CPC::dotProduct(paper_v0,paper_gv0);
        {
            paper_G6=paper_v0gv0;
            paper_G6*=paper_wgv0;
            paper_G6*=(-2);
            paper_G6+=paper_wv0;
            paper_G6*=paper_wv0;
            PA::Polynomial<Frac> tmp=paper_wgv0;
            tmp*=paper_wgv0;
            paper_G6+=tmp;
        }
        {
            paper_U1Sq=paper_wv0;
            paper_U1Sq*=paper_wv0;
            paper_U1Sq*=-1;
            paper_U1Sq+=paper_ww;
        }
        {
            paper_Q=paper_ww;
            PA::Polynomial<Frac> four(4);
            paper_Q-=four;
            paper_Q*=paper_N2;
            paper_Q-=paper_G6;
        }
        PA::Polynomial<Frac> qSq=createQSquare();
        {
            std::vector<PA::Polynomial<Frac> > WHalf;
            WHalf.resize(3);
            WHalf[0]=paper_R[0][2];
            WHalf[1]=paper_R[1][2];
            WHalf[2]=paper_R[2][2];
            WHalf[2]-=qSq;
            PA::Polynomial<Frac> wHwH=dotProduct(WHalf,WHalf);
            paper_Uh=wHwH;
            PA::Polynomial<Frac> R13Sq;
            R13Sq=paper_R[0][2];
            R13Sq*=R13Sq;
            paper_Uh-=R13Sq;
            PA::Polynomial<Frac> q4=qSq;
            q4*=q4;
            paper_Uh-=q4;
            paper_Uh*=4;

            std::vector<PA::Polynomial<Frac> > D2;
            D2.resize(3);
            D2[0]=paper_R[0][0];D2[1]=paper_R[1][0];D2[2]=paper_R[2][0];
            PA::Polynomial<Frac> WHalfv0=dotProduct(WHalf,paper_v0);
            PA::Polynomial<Frac> WHalfD2=dotProduct(WHalf,D2);
            paper_Qh=wHwH;
            paper_Qh-=q4;
            PA::Polynomial<Frac> t12=paper_R[1][0];
            t12*=t12;
            PA::Polynomial<Frac> t122=paper_R[2][0];
            t122*=t122;
            t12+=t122;
            paper_Qh*=t12;
            PA::Polynomial<Frac> t2=q4;
            t2*=WHalfv0;t2*=WHalfv0;
            PA::Polynomial<Frac> t2oth=WHalfv0;
            t2oth*=WHalfD2;
            t2oth*=paper_R[0][0];
            t2oth*=2;
            t2-=t2oth;
            t2oth=WHalfD2;
            t2oth*=t2oth;
            t2+=t2oth;
            paper_Qh-=t2;
            paper_Qh*=4;
        }
    }
    PA::Polynomial<Frac> createP(const std::vector<std::vector<std::vector<Frac> > >& array_C,
                                   const std::vector<std::vector<std::vector<PA::Polynomial<Frac> > > >& array_TPolynomials,
                                   const std::vector<std::vector<Frac> >& array_invWDiagonals){
        PA::Polynomial<Frac> res;
        PA::Polynomial<Frac> qSq=createQSquare();
        PA::Polynomial<Frac> multiplier(1);
        long L4=array_C.size();//must be 4
        for(long ell=4;ell>0;--ell){
            PA::Polynomial<Frac> FracSummation;
            for(long i=0;i<array_C[ell-1].size();++i){
                for(long j=0;j<array_C[ell-1].size();++j){
                    PA::Polynomial<Frac> Th_lqij=(array_TPolynomials[ell-1])[i][j];
                    Frac C_lij=(array_C[ell-1])[i][j];
                    if( (Th_lqij.terms.size()>0) && ( !( C_lij==0 ) ) ){
                        Th_lqij *= C_lij*((array_invWDiagonals[ell-1])[j]);
                        FracSummation+=Th_lqij;
                    }
                }
            }
            FracSummation*=multiplier;
            res+=FracSummation;
            multiplier*=qSq;
        }
        return res;
    }
    std::vector<Frac> divideVectors(const std::vector<long>& num, const std::vector<long>& den){
        std::vector<Frac> res;
        if(num.size()<1){return res;}
        if(num.size()!=den.size()){return res;}
        res.resize(num.size());
        for(long i=0;i<num.size();++i){
            if(den[i]!=0){
                res[i].setNum(num[i]);
                res[i].setDen(den[i]);
            }
        }
        return res;
    }
    std::vector<PA::Polynomial<Frac> > hamiltonProductBasic(){
        std::vector<PA::Polynomial<Frac> > res;
        res.resize(4);
        std::vector<PA::Polynomial<Frac> > p,q;
        p.resize(4);q.resize(4);
        for(long i=0;i<4;++i){
            PA::Monomial<Frac> mP(1,"p_"+std::to_string(i),1);
            PA::Monomial<Frac> mQ(1,"q_"+std::to_string(i),1);
            p[i]+=mP;q[i]+=mQ;
        }
        PA::Polynomial<Frac> tmp;
        //                       (p*q)_0 = p0q0 - p1q1 - p2q2 - p3q3
        tmp.terms.clear(); tmp+=p[0];tmp*=q[0];
        res[0]+=tmp;
        tmp.terms.clear(); tmp+=p[1];tmp*=q[1];
        res[0]-=tmp;
        tmp.terms.clear(); tmp+=p[2];tmp*=q[2];
        res[0]-=tmp;
        tmp.terms.clear(); tmp+=p[3];tmp*=q[3];
        res[0]-=tmp;
        //                       (p*q)_1 = p0q1 + p1q0 + p2q3 - p3q2
        tmp.terms.clear(); tmp+=p[0];tmp*=q[1];
        res[1]+=tmp;
        tmp.terms.clear(); tmp+=p[1];tmp*=q[0];
        res[1]+=tmp;
        tmp.terms.clear(); tmp+=p[2];tmp*=q[3];
        res[1]+=tmp;
        tmp.terms.clear(); tmp+=p[3];tmp*=q[2];
        res[1]-=tmp;
        //                       (p*q)_2 = p0q2 - p1q3+ p2q0 + p3q1
        tmp.terms.clear(); tmp+=p[0];tmp*=q[2];
        res[2]+=tmp;
        tmp.terms.clear(); tmp+=p[1];tmp*=q[3];
        res[2]-=tmp;
        tmp.terms.clear(); tmp+=p[2];tmp*=q[0];
        res[2]+=tmp;
        tmp.terms.clear(); tmp+=p[3];tmp*=q[1];
        res[2]+=tmp;
        //                       (p*q)_3 = p0q3 + p1q2 - p2q1 + p3q0
        tmp.terms.clear(); tmp+=p[0];tmp*=q[3];
        res[3]+=tmp;
        tmp.terms.clear(); tmp+=p[1];tmp*=q[2];
        res[3]+=tmp;
        tmp.terms.clear(); tmp+=p[2];tmp*=q[1];
        res[3]-=tmp;
        tmp.terms.clear(); tmp+=p[3];tmp*=q[0];
        res[3]+=tmp;
        return res;
    }
    int checkMultiplicativity(const std::vector<std::vector<PA::Polynomial<Frac> > >& TL){
        if(GL_skipMultVerification){return 1;}
        std::map<std::string,PA::Polynomial<Frac> > hamiltonSubstitution,qpSubstitution;
        std::vector<PA::Polynomial<Frac> > hamiltonProduct=hamiltonProductBasic();
        for(long i=0;i<4;++i){
            hamiltonSubstitution["q_"+std::to_string(i)]=hamiltonProduct[i];
            PA::Polynomial<Frac> pI;
            PA::Monomial<Frac> mI(1,"p_"+std::to_string(i),1);
            pI+=mI;
            qpSubstitution["q_"+std::to_string(i)]=pI;
        }
        std::vector<std::vector<PA::Polynomial<Frac> > > TLH=TL, TP=TL;//just to get the formats right

        for(long i=0;i<TLH.size();++i){
            for(long j=0;j<TLH[i].size();++j){
                TLH[i][j]=PA::substitute(TL[i][j],hamiltonSubstitution);
                TP[i][j]=PA::substitute(TL[i][j],qpSubstitution);
            }
        }
        std::vector<std::vector<PA::Polynomial<Frac> > > TProd=TP*TL;
        for(long i=0;i<TLH.size();++i){
            for(long j=0;j<TLH[i].size();++j){
                if(equal(TLH[i][j],TProd[i][j])==0){
                    return 0;
                }
            }
        }
        return 1;
    }
    int checkOrthogonality(const std::vector<std::vector<PA::Polynomial<Frac> > >& TL,
                const std::vector<long>& array_dNum,
                const std::vector<long>& array_dDen){
        if(GL_skipOrthVerification){return 1;}
        std::vector<std::vector<PA::Polynomial<Frac> > > matrixD;
        std::vector<PA::Polynomial<Frac> > placeHolderRow;
        long dim=TL.size();
        placeHolderRow.resize(dim);
        matrixD.resize(dim);
        for(long i=0;i<dim;++i){matrixD[i]=placeHolderRow;}
        long lcmDen=1;
        for(long i=0;i<dim;++i){
            lcmDen=FA::leastCommonMultiple(lcmDen,array_dDen[i]);
        }
        for(long i=0;i<dim;++i){
            PA::Monomial<Frac> mI(static_cast<Frac>(array_dNum[i]*lcmDen/array_dDen[i]));
            matrixD[i][i]+=mI;
        }
        long ell=(dim-1)/2;
        PA::Polynomial<Frac> qSq=createQSquare();
        PA::Polynomial<Frac> qFourL(1);
        for(long i=0;i<ell;++i){
            qFourL*=qSq;
        }
        qFourL*=qFourL;
        std::vector<std::vector<PA::Polynomial<Frac> > > LHS;
        LHS=transpose(TL) * matrixD;
        LHS= LHS * TL;
        std::vector<std::vector<PA::Polynomial<Frac> > > RHS=matrixD;
        for(long i=0;i<dim;++i){
            RHS[i][i]*=qFourL;
        }
        for(long i=0;i<LHS.size();++i){
            for(long j=0;j<LHS[i].size();++j){
                if(equal(LHS[i][j],RHS[i][j])==0){return 0;}
            }
        }
        return 1;
    }
    
    int checkMultiplicativityAndOrthogonality(const std::vector<std::vector<std::vector<PA::Polynomial<Frac> > > >& array_TPolynomials,
                               const std::vector<std::vector<long> >& array_dNums,
                               const std::vector<std::vector<long> >& array_dDens){
        int success=1;
        long l=0;
        while(l<4){
            std::cout<<"Multiplicativity and orthogonality checks for ell="<<l+1<<".";
            if(l==2){std::cout<<" Please be patient.";}
            if(l==3){std::cout<<" Please be very patient.";}
            std::cout<<"\n";
            success*=checkMultiplicativity(array_TPolynomials[l])*checkOrthogonality(array_TPolynomials[l],array_dNums[l],array_dDens[l]);
            if(success==1){
                std::cout<<"Passed. \n";
            }
            if(success==0){
                std::cout<<"Failed. \n";
                return 0;
            }
            ++l;
        }
        return success;
    }
    std::vector<std::vector<PA::Polynomial<Frac> > > createTPolynomial(const std::vector<std::vector<long> >& tRaw,
                                                                         long tDenominator, long ell){
        std::vector<std::vector<PA::Polynomial<Frac> > > res;
        res.resize(2*ell+1);
        std::vector<PA::Polynomial<Frac> > placeHolder;
        placeHolder.resize(2*ell+1);
        for(long i=0;i<res.size();++i){res[i]=placeHolder;}
        for(long k=0;k<tRaw.size();++k){
            PA::Monomial<Frac> mn;
            mn.c.setNum(tRaw[k][6]);
            mn.c.setDen(tDenominator);
            for(long j=0;j<4;j++){
                if(tRaw[k][j+2]>0){
                    mn.vExp["q_"+std::to_string(j)]=tRaw[k][j+2];
                }
            }
            res[tRaw[k][0]][tRaw[k][1]]+=mn;
        }
        return res;
    }
    std::vector<std::vector<long> > matrixLongFromString(const std::string& _in){
        std::vector<std::string> rows=SF::stringToVector(_in,GL_openTag,GL_closeTag);
        std::vector<std::vector<long> > res;
        res.resize(rows.size());
        for(long i=0;i<rows.size();++i){
            res[i]=SF::stringToVectorLong(rows[i],GL_openTag,GL_closeTag);
        }
        return res;
    }
    int checkHomogeneity(const std::vector<std::vector<long> >& T, long twoEll){
        int res=1;
        long i=T.size();
        while((i>0)&&(res==1)){
            --i;
            if(T[i][2]+T[i][3]+T[i][4]+T[i][5]!=twoEll){res=0;}
        }
        return res;
    }
    int areWsPositive(const std::vector<long>& dens, const std::vector<long>& nums){
        int res=1;
        long i=dens.size();
        if(nums.size()!=i){return 0;}
        while((i>0)&&(res==1)){
            --i;
            if((dens[i]<1)||(nums[i]<1)){res=0;}
        }
        return res;
    }
    int getFromString(const std::string& _in,
                      std::vector<std::vector<std::vector<PA::Polynomial<Frac> > > >& array_TPolynomials,
                      std::vector<std::vector<Frac> >& array_invWDiagonals
                     ){
        long L4=4;
        array_TPolynomials.resize(L4);
        array_invWDiagonals.resize(L4);
        int indicatorHomogeneity=1;
        int indicatorWPositive=1;
        std::vector<std::vector<long> > array_dNums,array_dDens;
        array_dNums.resize(L4);array_dDens.resize(L4);
        std::map<std::string,std::string>::const_iterator it;
        std::map<std::string,std::string> mainMap=SF::stringToMap(_in,"[key]","[/key]","[value]","[/value]");
        it=mainMap.find("skipMultVerification");
        if(it!=mainMap.end()){
            if(it->second=="yes"){GL_skipMultVerification=1;}
        }
        it=mainMap.find("skipOrthVerification");
        if(it!=mainMap.end()){
            if(it->second=="yes"){GL_skipOrthVerification=1;}
        }
        it=mainMap.find("stringSizeLimit");
        if(it!=mainMap.end()){
            GL_stringSizeLimit=BF::stringToInteger(it->second);
        }
        for(long l=0;l<L4;++l){
            long ell=l+1;
            it=mainMap.find("Tden000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            long tDenominator=SF::stringToVectorLong(SF::stringToVector(it->second,GL_openTag,GL_closeTag)[0],GL_openTag,GL_closeTag)[0];
            it=mainMap.find("WNum000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            array_dNums[l]=SF::stringToVectorLong(SF::stringToVector(it->second,GL_openTag,GL_closeTag)[0],GL_openTag,GL_closeTag);
            it=mainMap.find("WDen000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            array_dDens[l]=SF::stringToVectorLong(SF::stringToVector(it->second,GL_openTag,GL_closeTag)[0],GL_openTag,GL_closeTag);
            it=mainMap.find("Tpoly000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            std::vector<std::vector<long> > tRaw=matrixLongFromString(it->second);
            indicatorHomogeneity*=checkHomogeneity(tRaw,2*ell);
            indicatorWPositive*=areWsPositive(array_dDens[l],array_dNums[l]);
            array_invWDiagonals[l]=divideVectors(array_dDens[l],array_dNums[l]);
            array_TPolynomials[l]=createTPolynomial(tRaw,tDenominator,ell);
        }
        return indicatorHomogeneity*indicatorWPositive*checkMultiplicativityAndOrthogonality(array_TPolynomials,array_dNums,array_dDens);
    }
    std::vector<std::vector<Frac> > divideByScalar(const std::vector<std::vector<long> >& M, long den){
        std::vector<std::vector<Frac> > res;
        std::vector<Frac> placeHolder;
        if(M.size()<1){return res;}
        res.resize(M.size());
        placeHolder.resize(M.size());
        for(long i=0;i<M.size();++i){res[i]=placeHolder;}
        for(long i=0;i<M.size();++i){
            for(long j=0;j<M.size();++j){
                res[i][j].setNum(M[i][j]);
                res[i][j].setDen(den);
            }
        }
        return res;
    }
    int checkPsdWitness(const std::vector<std::vector<long> >& C,
                            const std::vector<std::vector<long> >& W,
                            long sWitt){
        //number sWitt must be positive for this test to be meaningful
        if(sWitt<=0){return 0;}
        std::vector<std::vector<long> > E=C;
        std::vector<std::vector<long> > WTW=transpose(W)*W;
        for(long i=0;i<E.size();++i){
            for(long j=0;j<E.size();++j){
                E[i][j]=sWitt*C[i][j]-WTW[i][j];
            }
        }
        long absVal;
        for(long i=0;i<E.size();++i){
            for(long j=i+1;j<E.size();++j){
                if(E[i][j]!=E[j][i]){return 0;}
            }
            long aw=2*E[i][i];
            for(long j=0;j<E.size();++j){
                absVal=E[i][j];
                if(absVal<0){absVal*=-1;}
                aw-=absVal;
            }
            if(aw<0){return 0;}
        }
        return 1;
    }
    int getFromString(const std::string& _in,
                      std::vector<std::vector<std::vector<Frac> > >& array_C,
                      const std::vector<std::vector<Frac> >& array_invWDiagonals){
        long L4=4;
        array_C.resize(L4);
        long twoTo20=1048576;
        std::map<std::string,std::string>::const_iterator it;
        std::map<std::string,std::string> mainMap=SF::stringToMap(_in,"[key]","[/key]","[value]","[/value]");
        for(long l=0;l<L4;++l){
            long ell=l+1;
            it=mainMap.find("Swit000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            long sWitt=SF::stringToVectorLong(SF::stringToVector(it->second,GL_openTag,GL_closeTag)[0],GL_openTag,GL_closeTag)[0];
            it=mainMap.find("Cscaled000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            std::vector<std::vector<long> > cLong=matrixLongFromString(it->second);
            it=mainMap.find("Wwit000"+std::to_string(ell));
            if(it==mainMap.end()){return 0;}
            std::vector<std::vector<long> > wLong=matrixLongFromString(it->second);
            if(!checkPsdWitness(cLong,wLong,sWitt)){
                return 0;
            }
            array_C[l]=divideByScalar(cLong,twoTo20);
        }
        Frac wTracesTotal=0;
        for(long l=0;l<L4;++l){
            for(long i=0;i<array_C[l].size();++i){
                wTracesTotal += (array_C[l])[i][i] * array_invWDiagonals[l][i];
            }
        }
        return (wTracesTotal*static_cast<Frac>(4096) == static_cast<Frac>(23017));
    }
    PA::Polynomial<Frac> create20G(const std::vector<std::vector<std::vector<Frac> > >& array_C,
                                     const std::vector<std::vector<std::vector<PA::Polynomial<Frac> > > >& array_TPolynomials,
                                     const std::vector<std::vector<Frac> >& array_invWDiagonals){
        PA::Polynomial<Frac> twentyG=createP(array_C,array_TPolynomials,array_invWDiagonals);
        PA::Polynomial<Frac> tmp=createQSquare();
        tmp *= tmp;// now tmp=|q|^4
        tmp *= tmp;// now tmp=|q|^8
        tmp*=-19;
        twentyG*=-20;
        twentyG+=tmp;
        return twentyG;
    }


    struct TestSetup{
    public:
        std::vector<std::vector<std::vector<PA::Polynomial<Frac> > > > array_TPolynomials;
        std::vector<std::vector<Frac> > array_invWDiagonals;
        std::vector<std::vector<std::vector<Frac> > > array_C;
        PA::Polynomial<Frac> pol_20G;
        PA::Polynomial<Frac> pol_Uh;
        PA::Polynomial<Frac> pol_Qh;
    };
    TestSetup getFromFile(const std::string& inS){
        TestSetup res;
        int succ=getFromString(inS,res.array_TPolynomials,res.array_invWDiagonals);
        std::cout<<"Identities (hom, mult, orth, weights) = "<<succ<<"\n";
        succ=getFromString(inS, res.array_C, res.array_invWDiagonals);
        std::cout<<"Witnesses and trace sum (psd, trace) = "<<succ<<"\n";
        res.pol_20G=create20G(res.array_C,res.array_TPolynomials,res.array_invWDiagonals);
        res.pol_Uh=GL_C.paper_Uh;
        res.pol_Qh=GL_C.paper_Qh;
        return res;
    }
    template<typename TTT>
    struct Box{
        std::vector<TTT> ls;
        std::vector<TTT> hs;
        int inequalityCode;
        Box();
        int operator<(const Box<TTT>& ) const;
    };
    template<typename TTT>
    Box<TTT>::Box(){
        inequalityCode=0;
        ls.resize(3);hs.resize(3);
        ls[0]=-1;hs[0]=1;ls[1]=-1;hs[1]=1;ls[2]=-1;hs[2]=1;
    }
    template<typename TTT>
    int Box<TTT>::operator<(const Box<TTT>& oth)const{
        long i=0;
        TTT sz=hs[0]-ls[0];
        TTT oth_sz=oth.hs[0]-oth.ls[0];
        if(sz<oth_sz){return 1;}
        if(oth_sz<sz){return 0;}
        TTT distFromZero=ls[0]+ls[1]+ls[2];
        TTT othDistFromZero=oth.ls[0]+oth.ls[1]+oth.ls[2];
        if(distFromZero<othDistFromZero){return 1;}
        if(othDistFromZero<distFromZero){return 0;}
        while(i<3){
            if(ls[i]<oth.ls[i]){return 1;}
            if(oth.ls[i]<ls[i]){return 0;}
            if(hs[i]<oth.hs[i]){return 1;}
            if(oth.hs[i]<hs[i]){return 0;}
            ++i;
        }
        return 0;
    }
    template<typename TTT>
    void removeTopThenSplitAndAddToStack(std::stack<Box<TTT> >& _st){
        if(_st.size()==0){return;}
        Box<TTT> b=_st.top();
        _st.pop();
        TTT half(1,2);
        Box<TTT> tmp;
        tmp.ls[0]=b.ls[0];
        tmp.hs[0]=(b.hs[0]+b.ls[0])* half;
        tmp.ls[1]=b.ls[1];
        tmp.hs[1]=(b.hs[1]+b.ls[1])* half;
        tmp.ls[2]=b.ls[2];
        tmp.hs[2]=(b.hs[2]+b.ls[2])* half;
        _st.push(tmp);
        tmp.ls[0]=(b.hs[0]+b.ls[0])*half;
        tmp.hs[0]=b.hs[0];
       /* no need to do these instructions because they stayed from the old code:
        tmp.ls[1]=b.ls[1];
        tmp.hs[1]=(b.hs[1]+b.ls[1])* half;
        tmp.ls[2]=b.ls[2];
        tmp.hs[2]=(b.hs[2]+b.ls[2])* half;*/
        _st.push(tmp);
        tmp.ls[1]=(b.ls[1]+b.hs[1])*half;
        tmp.hs[1]=b.hs[1];
        // happy with previous [0] and [2]
        _st.push(tmp);
        tmp.ls[2]=(b.ls[2]+b.hs[2])*half;
        tmp.hs[2]=b.hs[2];
        // happy with previous [0] and [1]
        _st.push(tmp);
        tmp.ls[1]=b.ls[1];
        tmp.hs[1]=(b.hs[1]+b.ls[1])* half;
        // happy with previous [0] and [2]
        _st.push(tmp);
        tmp.ls[0]=b.ls[0];
        tmp.hs[0]=(b.hs[0]+b.ls[0])* half;
        // happy with previous [1] and [2]
        _st.push(tmp);
        tmp.ls[1]=(b.ls[1]+b.hs[1])*half;
        tmp.hs[1]=b.hs[1];
        //happy with previous [0] and [2]
        _st.push(tmp);
        tmp.ls[2]=b.ls[2];
        tmp.hs[2]=(b.ls[2]+b.hs[2])*half;
        //happy with previous [0] and [1]
        _st.push(tmp);
    }
}

#endif
