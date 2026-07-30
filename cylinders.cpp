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

#include "src/fundamental/SHF_standardHeadersMinimal.cpp"
#include <thread>
#include "src/fundamental/ssm.cpp"
#include "src/fundamental/fst.cpp" 
#include "src/fundamental/MFRF_massiveFindReplace.cpp"
#include "src/fundamental/BF_bijFunction.cpp"
#include "src/fundamental/TMF_timer.cpp"
#include "src/fundamental/SF_strings.cpp"
#include "src/fundamental/GF_globalInitsMinimal.cpp"
#include "src/fundamental/RNDF_randomCodesAndPermutations.cpp"
#include "src/fundamental/IOF_io.cpp"
#include "src/tools/CLAT_commandLineArguments.cpp"
#include "src/algebra/FA_fractions.cpp"
#include "src/algebra/PA_polynomials.h"
#include "src/algebra/PA_polynomials.cpp"
#include "src/cylinders/CPC_cPolynomials.cpp"
#include "src/cylinders/CC_cylinders.cpp"
 
int main(int argc, char **argv){
    TMF::Timer tm;
    tm.start();
    CLAT::CommandLineArgs inputArgs=CLAT::getArguments(CLAT::argsInStrings(argc,argv));
    std::string sourceTxt;
    if(inputArgs.inputFile=="notFound"){
        inputArgs.helpNeeded=1;
    }
    else{
        sourceTxt=IOF::fileToString(inputArgs.inputFile);
        if(sourceTxt=="fileNotFound"){
            std::cout<<"Input file not found\n";
            return 0;
        }
    }
    if(inputArgs.helpNeeded){
        std::cout<<CC::helpString()<<"\n";
        return 0;
    }
    std::cout<<CC::mainWork(sourceTxt,inputArgs.outputFile)<<"\n";
    if((FA::GL_OVERFLOW)||(FA::GL_DIVBY0)){
        std::cout<<"WARNING: The program had a division by zero or an overflow in calculation. Results are not trustworthy. \n";
    }
    tm.end();
    std::cout<<"Total time: "<<tm.getTimeSec()<<" seconds.\n";
    return 0;
}

