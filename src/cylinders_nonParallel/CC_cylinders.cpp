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


#ifndef __INCL_CC_CPP
#define __INCL_CC_CPP
namespace CC{
std::string helpString(){
    std::string res;
    res+="Syntax: ./cyl -i cyl_inputData.txt -o cyl_outputData\n";
    return res;
}
Frac power(Frac x, long a){
    Frac res=1;
    for(long i=0;i<a;++i){res*=x;}
    return res;
}
Frac findConservativeMax(long a, long b, long c, Frac coeff, const CPC::Box<Frac>& box){
    Frac xV[2], yV[2],zV[2];
    xV[0]=power(box.ls[0],a);xV[1]=power(box.hs[0],a);
    yV[0]=power(box.ls[1],b);yV[1]=power(box.hs[1],b);
    zV[0]=power(box.ls[2],c);zV[1]=power(box.hs[2],c);
    Frac extreme1=xV[0]*yV[0]*zV[0];
    Frac extreme2=extreme1;
    Frac tmp;
    for(long i=0;i<2;++i){
        for(long j=0;j<2;++j){
            for(long k=0;k<2;++k){
                tmp=xV[i]*yV[j]*zV[k];
                if(tmp<extreme1){extreme1=tmp;}
                if(extreme2<tmp){extreme2=tmp;}
            }
        }
    }
    extreme1*=coeff;extreme2*=coeff;
    if(extreme1<extreme2){return extreme2;}
    return extreme1;
}
Frac findConservativeBound(const std::vector<std::vector<std::vector<Frac> > > &mP,
                             const CPC::Box<Frac>& b,
                             int boundTypePm1){
    Frac res=0;
    long maxDegree=mP.size();
    for(long i=0;i<maxDegree;++i){
        for(long j=0;j<maxDegree;++j){
            for(long k=0;k<maxDegree;++k){
                if(!(mP[i][j][k]==0)){
                    if(boundTypePm1<0){
                        Frac mpijk=mP[i][j][k];
                        mpijk.setNum(-mpijk.getNum());
                        res-=findConservativeMax(i,j,k,mpijk,b);
                    }
                    else{
                        res+=findConservativeMax(i,j,k,mP[i][j][k],b);
                    }
                }
            }
        }
    }
    return res;
}
int whichInequalityCanBeGuaranteedOnBox(const std::vector<std::vector<std::vector<Frac> > > &mP,
                                        const CPC::Box<Frac>& b){
    // returns 1 if mP>=0 is guaranteed on the box
    // returns -1 if mP<0 is guaranteed on the box
    // returns 0 if the algorithm can guarantee neither mP>=0 nor mP<0
    Frac conservativeUpperBound=findConservativeBound(mP,b,1);
    Frac conservativeLowerBound=findConservativeBound(mP,b,-1);
    if(conservativeLowerBound.getNum()>=0){return 1;}
    if(conservativeUpperBound.getNum()<0){return -1;}
    return 0;
}
int succeedOnOneBox(const std::vector<std::vector<std::vector<Frac> > > &mGh,
                    const std::vector<std::vector<std::vector<Frac> > > &mUh,
                    const std::vector<std::vector<std::vector<Frac> > > &mQh,
                    const CPC::Box<Frac>& b){
    if(whichInequalityCanBeGuaranteedOnBox(mUh,b)==-1){return 1;}
    if(whichInequalityCanBeGuaranteedOnBox(mQh,b)==-1){return 2;}
    if(whichInequalityCanBeGuaranteedOnBox(mGh,b)== 1){return 3;}
    return 0;
}
struct TestResult{
public:
    ssm::set<CPC::Box<Frac> > allBoxes;
    long counterPartition;
    int success;
};
TestResult runTest(const PA::Polynomial<Frac>& Gh,
                   const PA::Polynomial<Frac>& Uh,
                   const PA::Polynomial<Frac>& Qh){
    TestResult res;
    res.success=0;
    res.counterPartition=0;
    std::vector<std::vector<std::vector<Frac> > > mGh=PA::polMatrix(Gh);
    std::vector<std::vector<std::vector<Frac> > > mUh=PA::polMatrix(Uh);
    std::vector<std::vector<std::vector<Frac> > > mQh=PA::polMatrix(Qh);
    CPC::Box<Frac> bigBox,currentBox;
    int successInd;
    std::stack<CPC::Box<Frac> > remaining;
    remaining.push(bigBox);
    CPC::removeTopThenSplitAndAddToStack(remaining);
    long giveUpCounter=5000000;
    while((remaining.size()>0)&&(giveUpCounter>0)){
        --giveUpCounter;currentBox=remaining.top();
        if(giveUpCounter%10000==0){
            std::cout<<"Give-up-counter is "<<
            giveUpCounter<<" remaining boxes = "
            <<remaining.size()<<" certified so far = "<<res.counterPartition
            <<" overflow flag = "<<FA::GL_OVERFLOW<<"\n";
            if(FA::GL_OVERFLOW){std::cout<<"OVERFLOW: aborting face.\n";return res;}
        }
        successInd=succeedOnOneBox(mGh,mUh,mQh,currentBox);
        if(successInd!=0){ 
            remaining.pop();
            ++res.counterPartition;
            currentBox.inequalityCode=successInd;
            res.allBoxes.insert(currentBox);
        }
        else{
            CPC::removeTopThenSplitAndAddToStack(remaining);
        }
    }
    if(remaining.size()==0){
        res.success=1;
    }
    return res;
}
std::string initialCommentBoxes(){
    std::string res;
    res+="This file contains the partition of the three-dimensional cube\n";
    res+="centered at origin. For each box in the partition, we write \n";
    res+="1: If we can guarantee that U_h<0\n";
    res+="2: If we can guarantee that Q_h<0\n";
    res+="3: If we can guarantee that G>=0\n";
    return res;
}
std::string printSetup(const CPC::TestSetup &ts){
    std::string res;
    for(long l=0;l<4;++l){
        long ell=l+1;
        int matrixStarted;
        res+="\\section{Term $\\ell="+std::to_string(ell)+"$}\n";
        res+="\\subsection{Matrix $D_"+std::to_string(ell)+"$}\n";
        res+="We will write the matrix $D^{-1}$ instead of $D$. The matrices $D$ and $D^{-1}$ ";
        res+="are diagonal, so one is easy to construct from the other. The matrix $D^{-1}$ is ";
        res+="more convenient because the denominators of diagonal terms are powers of $2$. Therefore, ";
        res+="their binary expressions are finite and there is no rounding in computer storage.";
        res+="\\[D_"+std::to_string(ell)+"^{-1}=\\left[\\begin{array}{";
        for(long i=0;i<2*ell+1;++i){res+="c";}
        res+="}\n";
        matrixStarted=0;
        for(long i=0;i<2*ell+1;++i){
            if(matrixStarted){res+="\\\\\n";}
            matrixStarted=1;
            int rowStarted=0;
            for(long j=0;j<2*ell+1;++j){
                if(rowStarted){
                    res+=" & ";
                }
                rowStarted=1;
                if(i==j){
                    res+=PA::printDiadicRational(ts.array_invWDiagonals[l][i]);
                }
                else{
                    res+=std::to_string(0);
                }
            }
        }
        res+="\n\\end{array}\\right].\\]\n";
        res+="\\subsection{Matrix $C_"+std::to_string(ell)+"$}\n";
        res+="\\[C_"+std::to_string(ell)+"=\\left[\\begin{array}{";
        for(long i=0;i<2*ell+1;++i){res+="c";}
        res+="}\n";
        matrixStarted=0;
        for(long i=0;i<2*ell+1;++i){
            if(matrixStarted){res+="\\\\\n";}
            matrixStarted=1;
            int rowStarted=0;
            for(long j=0;j<2*ell+1;++j){
                if(rowStarted){
                    res+=" & ";
                }
                rowStarted=1;
                res+=PA::printDiadicRational(ts.array_C[l][i][j]);
            }
        }
        res+="\n\\end{array}\\right].\\]\n";
        res+="\\subsection{Matrix $T_{h}^"+std::to_string(ell)+"(q_0,q_1,q_2,q_3)$}\n";
        for(long i=0;i<2*ell+1;++i){
            matrixStarted=1;
            for(long j=0;j<2*ell+1;++j){
                res+="\\begin{align*}(T_{h}^"+std::to_string(ell)+")_{"+std::to_string(i+1)+","+std::to_string(j+1)+"}&=";
                res+=(ts.array_TPolynomials[l][i][j]).toString();
                res+="\\end{align*}\n";
            }
        }
    }
    res+="\\section{Polynomial $U_h$}\n";
    res+="\\begin{align*}U_h(q_0,q_1,q_2,q_3)&="+ts.pol_Uh.toString()+"\\end{align*}\n";
    res+="\\section{Polynomial $Q_h$}\n";
    res+="\\begin{align*}Q_h(q_0,q_1,q_2,q_3)&="+ts.pol_Qh.toString()+"\\end{align*}\n";
    res+="\\section{Polynomial $G$}\n";
    res+="We store $20G$ instead of $G$ in the memory because $20G$ can have all of its coefficients stored with full precision, without rounding. Recall that computers work in binary and for them, $\\frac15$ is an infinite expression that must be rounded. Computers can only tolerate denominators that are powers of two.\n";
    res+="\\begin{align*}20G(q_0,q_1,q_2,q_3)&="+ts.pol_20G.toString()+"\\end{align*}\n";
    return res;
}
std::string polXYZSummary(const std::string& polName,const PA::Polynomial<Frac>& pol){
    std::string res;
    ssm::set<std::string> avr= PA::allVariables(pol);
    std::map<std::string,PA::Polynomial<Frac> > substMap;
    std::string vNames="xyz";
    for(long i=0;i<avr.size();++i){
        std::string tmpName;
        tmpName+=vNames[i];
        PA::Monomial<Frac> m(1,tmpName,1);
        PA::Polynomial<Frac> p;
        p+=m;
        substMap[avr[i]]=p;
    }
    PA::Polynomial<Frac> pxyz=PA::substitute(pol,substMap);
    res+="\\subsection{Polynomial $"+polName+"(x,y,z)$}\n";
    res+="\\begin{align*}"+polName+"(x,y,z)&="+pxyz.toString()+"\\end{align*}";
    return res;
}
std::string mainWork(const std::string& _srcTxt, const std::string& scDest){
    CPC::TestSetup ts=CPC::getFromFile(_srcTxt);
    if(ts.success==0){return "Failed.";}
    PA::Polynomial<Frac> p20G3;
    PA::Polynomial<Frac> pUh3;
    PA::Polynomial<Frac> pQh3;
    TestResult res1T;
    int success=1;
    std::string res,resB;
    long totalCounter=0;
    res+=printSetup(ts);
    resB+=initialCommentBoxes();
    for(long i=0;i<4;++i){
        resB+="\n\n\n*** Substitution q_"+std::to_string(i)+"=1 ***\n\n";
        res+="\\section{Substitution $q_"+std::to_string(i)+"=1$"+"}\n";
        std::cout<< ("\n\n\n*** Substitution q_"+std::to_string(i)+"=1 ***\n\n");
        std::map<std::string,Frac> vEval;
        vEval["q_"+std::to_string(i)]=1;
        p20G3=ts.pol_20G.evaluate(vEval);
        pUh3=ts.pol_Uh.evaluate(vEval);
        pQh3=ts.pol_Qh.evaluate(vEval);
        res+=polXYZSummary("U_h",pUh3);
        res+=polXYZSummary("Q_h",pQh3);
        res+=polXYZSummary("20G",p20G3);
        res1T=runTest(p20G3,pUh3,pQh3);
        for(long iBox=0;iBox<res1T.allBoxes.size();++iBox){
            CPC::Box<Frac> currentBox=res1T.allBoxes[iBox];
            resB+="[";
            for(long iCoord=0;iCoord<3;++iCoord){
                resB+="[["+PA::printDiadicRational(currentBox.ls[iCoord])+"][";
                resB+=PA::printDiadicRational(currentBox.hs[iCoord])+"]]";
            }
            resB+="["+std::to_string(currentBox.inequalityCode)+"]";
            resB+="]\n";
        }
        
        resB+="Total number of boxes in the partition = "+std::to_string(res1T.counterPartition)+"\n\n";
        totalCounter+=res1T.counterPartition;
        success*=res1T.success;
    }
    IOF::toFile(scDest+".tex",res);
    IOF::toFile(scDest+".txt",resB);
    return "Success = "+std::to_string(success)
            +"\nTotal number of certified boxes: "
            +std::to_string(totalCounter);;
}
}

#endif
