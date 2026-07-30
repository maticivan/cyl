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
void runTestP(long threadId, long* nThreads,
              fst::stack<CPC::Box<Frac> >* sBoxes,
              long* numStacksForOneEll,
              fst::stack<CPC::Box<Frac> >* a_r_B,
              long* a_forTotalCounters,
              int* a_forSuccess,
              std::vector<std::vector<std::vector<Frac> > >* a_mGh,
              std::vector<std::vector<std::vector<Frac> > >* a_mUh,
              std::vector<std::vector<std::vector<Frac> > >* a_mQh){
    long ellShiftStacks=0;
    long stackToWork;
    long writingPosition=threadId;
    int& ind_success=*a_forSuccess;
    long giveUpCounter;
    int successInd;
    CPC::Box<Frac> currentBox;
    for(long ell=0;ell<4;++ell){
        fst::stack<CPC::Box<Frac> >& wB= a_r_B[writingPosition];
        long& totalCounter=a_forTotalCounters[writingPosition];
        stackToWork=threadId;
        while(stackToWork<*numStacksForOneEll){
            giveUpCounter=5000000;
            fst::stack<CPC::Box<Frac> >& remaining=sBoxes[stackToWork+ellShiftStacks];
            while((remaining.size()>0)&&(giveUpCounter>0)){
                    --giveUpCounter;currentBox=remaining.top();
                    successInd=succeedOnOneBox(a_mGh[ell],a_mUh[ell],a_mQh[ell],currentBox);
                    if(successInd!=0){
                        remaining.pop();
                        ++totalCounter;
                        currentBox.inequalityCode=successInd;
                        wB.push(currentBox);
                    }
                    else{
                        CPC::removeTopThenSplitAndAddToStack(remaining);
                    }
                }
                if(remaining.size()!=0){
                    ind_success=0;
                }
            stackToWork+=*nThreads;
        }
        
        ellShiftStacks+=*numStacksForOneEll;
        writingPosition+=*nThreads;
    }
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
void splitUntilAcceptable(fst::stack<CPC::Box<Frac> >& allBoxes,
                          const CPC::Box<Frac>& boxToSplit,
                          Frac targetSize){
    if( !(targetSize<boxToSplit.hs[0]-boxToSplit.ls[0])){
        allBoxes.push(boxToSplit);
        return;
    }
    fst::stack<CPC::Box<Frac> > oneSplit(8);
    oneSplit.push(boxToSplit);
    CPC::removeTopThenSplitAndAddToStack(oneSplit);
    while(!oneSplit.empty()){
        splitUntilAcceptable(allBoxes,oneSplit.top(),targetSize);
        oneSplit.pop();
    }
}
std::pair<fst::stack<CPC::Box<Frac> >*,long> splitIntoBoxes(int exponentOfDenominator,int numCopies){
    // .second of the return value will be the number of stacks in one of the numCopies
    // .first will be an array of stacks of length (.second)*numCopies
    std::pair<fst::stack<CPC::Box<Frac> >*,long> res;
    CPC::Box<Frac> bigBox;
    Frac targetSize(1);Frac oneHalf(1,2);
    long numInOneCopy=1;
    for(int j=0;j<exponentOfDenominator;++j){
        targetSize*=oneHalf;
        numInOneCopy*=8;
    }
    numInOneCopy*=8;
    res.second=numInOneCopy;
    long total=numInOneCopy*numCopies;
    res.first=new fst::stack<CPC::Box<Frac> >[total];
    fst::stack<CPC::Box<Frac> > allBoxes(numInOneCopy);
    splitUntilAcceptable(allBoxes,bigBox,targetSize);
    std::vector<long> permutation=RNDF::genMTwistPermutation(numInOneCopy);
    long i=0;
    while(i<numInOneCopy){
        long j=permutation[i];
        while(j<total){
            ((res.first)[j]).push(allBoxes.top());
            j+=numInOneCopy;
        }
        allBoxes.pop();
        ++i;
    }
    return res;
}
std::string mainWork(const std::string& _srcTxt, const std::string& scDest){
    CPC::TestSetup ts=CPC::getFromFile(_srcTxt);
    if(ts.success==0){return "Failed.";}
    long nThreads=std::thread::hardware_concurrency();
    if(nThreads==0){nThreads=1;}
    std::cout<<"Partitioning the box and verifying the inequalities. This may take a while.\n";
    if(nThreads==1){
        std::cout<<"Your CPU has only one core. This may take up to 20 minutes.\n";
    }
    else{
        std::cout<<"Your CPU has multiple cores. I'll use them.\nNumber of threads that ";
        std::cout<<"are available to me is ";
        std::cout<<nThreads<<"\n";
    }
    int success=1;
    std::string res,resB;
    long totalCounter=0;
    res+=printSetup(ts);
    resB+=initialCommentBoxes();
    std::vector<std::vector<std::vector<Frac> > >* mGh=new std::vector<std::vector<std::vector<Frac> > >[4];
    std::vector<std::vector<std::vector<Frac> > >* mUh=new std::vector<std::vector<std::vector<Frac> > >[4];
    std::vector<std::vector<std::vector<Frac> > >* mQh=new std::vector<std::vector<std::vector<Frac> > >[4];
    PA::Polynomial<Frac>* a_pUh3=new PA::Polynomial<Frac>[4];
    PA::Polynomial<Frac>* a_pQh3=new PA::Polynomial<Frac>[4];
    PA::Polynomial<Frac>* a_p20G3=new PA::Polynomial<Frac>[4];
    for(long i=0;i<4;++i){
        std::map<std::string,Frac> vEval;
        vEval["q_"+std::to_string(i)]=1;
        a_pUh3[i]=ts.pol_Uh.evaluate(vEval);
        a_pQh3[i]=ts.pol_Qh.evaluate(vEval);
        a_p20G3[i]=ts.pol_20G.evaluate(vEval);
        
        
        mGh[i]=PA::polMatrix(a_p20G3[i]);
        mUh[i]=PA::polMatrix(a_pUh3[i]);
        mQh[i]=PA::polMatrix(a_pQh3[i]);
    }
    //presplitting to side  1/4 is lossless: an instrumented run of the adaptive algorithm showed
    //the largest box ever certified has side 1/4 on every face, so the adaptive tree contains the
    //uniform depth-3 tree.
    std::pair<fst::stack<CPC::Box<Frac> >*,long> bSplit=splitIntoBoxes(2,4);
    std::thread* allTh=new std::thread[nThreads];
    fst::stack<CPC::Box<Frac> >* a_Res_B=new fst::stack<CPC::Box<Frac> >[4*nThreads];
    for(long i=0;i<4*nThreads;++i){
        a_Res_B[i].setFastStorageCapacity(200000);
    }
    long* aTotalCounters=new long[4*nThreads];
    for(long i=0;i<4*nThreads;++i){aTotalCounters[i]=0;}
    int* aSuccesses=new int[nThreads];
    for(long i=0;i<nThreads;++i){aSuccesses[i]=1;}
    for(long i=0;i<nThreads;++i){
        allTh[i]=std::thread(runTestP,
                             i,&nThreads,
                             bSplit.first,
                             &(bSplit.second),
                             a_Res_B,
                             aTotalCounters,aSuccesses+i,
                             mGh, mUh, mQh
                             );
    }
    for(long i=0;i<nThreads;++i){
        allTh[i].join();
        for(long ell=0;ell<4;++ell){
            totalCounter+=aTotalCounters[i+ell*nThreads];
        }
        
        success*=aSuccesses[i];
    }
    long base=0;
    for(long ell=0;ell<4;++ell){
        resB+="\n\n\n*** Substitution q_"+std::to_string(ell)+"=1 ***\n\n";
        res+="\\section{Substitution $q_"+std::to_string(ell)+"=1$"+"}\n";
        res+=polXYZSummary("U_h",a_pUh3[ell]);
        res+=polXYZSummary("Q_h",a_pQh3[ell]);
        res+=polXYZSummary("20G",a_p20G3[ell]);
        long counterPartition=0;
        ssm::set<CPC::Box<Frac> > sAll;
        for(long i=0;i<nThreads;++i){
            while(!((a_Res_B[base+i]).empty())){
                sAll.insert((a_Res_B[base+i]).top());
                (a_Res_B[base+i]).pop();
            }
            counterPartition+=aTotalCounters[base+i];
        }
        for(long i=0;i<sAll.size();++i){
            CPC::Box<Frac> currentBox=sAll[i];
            resB+="[";
            for(long i=0;i<3;++i){
                resB+="[["+PA::printDiadicRational(currentBox.ls[i])+"][";
                resB+=PA::printDiadicRational(currentBox.hs[i])+"]]";
            }
            resB+="["+std::to_string(currentBox.inequalityCode)+"]";
            resB+="]\n";
        }
        resB+="Total number of boxes in the partition = "+std::to_string(counterPartition)+"\n\n";
        base+=nThreads;
    }
    delete[] aSuccesses;
    delete[] aTotalCounters;
    delete[] a_Res_B;
    delete[] allTh;
    delete[] bSplit.first;
    
    delete[] mGh;
    delete[] mUh;
    delete[] mQh;
    IOF::toFile(scDest+".tex",res);
    IOF::toFile(scDest+".txt",resB);
    return "Success = "+std::to_string(success)
            +"\nTotal number of certified boxes: "
            +std::to_string(totalCounter);
}
}

#endif
