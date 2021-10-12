#ifndef CAVIARMODEL_H
#define CAVIARMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

//#ifndef CAVIAR_POSTCAL_H
//#define CAVIAR_POSTCAL_H
#include "caviar_PostCal.h"
//#endif

using namespace std;
using namespace arma;

 
class CaviarModel{
	
public:
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
	double * sigma;
        double * stat;
        char * pcausalSet;
        int * rank;
        bool histFlag;
	postCal * post;
	string * snpNames;
	string ldFile;
        string zFile;
        string outputFileName;
        string geneMapFile;	

	CaviarModel(string ldFile, string zFile, string outputFileName, int totalCausalSNP, double NCP, double rho, bool histFlag, double gamma=0.01) {
//         cout<<"Perform caviar finemapping step"<<endl;

		int tmpSize = 0;
		this->histFlag = histFlag;
		this->NCP = NCP;
		this->rho = rho;
		this->gamma = gamma;
		this->ldFile = ldFile;
		this->zFile  = zFile;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
  //              cout<<"HERE1"<<endl;
		fileSize(ldFile, tmpSize);
    //            cout<<"HERE2"<<endl;
	        snpCount   = (int)sqrt(tmpSize);
         	sigma      = new double[snpCount * snpCount];
		stat       = new double[snpCount];
		pcausalSet = new char[snpCount];
		rank       = new int[snpCount];
		snpNames   = new string [snpCount];
      //          cout<<"HERE3"<<endl;
		importData(ldFile, sigma);
        //        cout<<"HERE4"<<endl;
		importDataFirstColumn(zFile, snpNames);
          //      cout<<"HERE5"<<endl;
		importDataSecondColumn(zFile, stat);
            //    cout<<"HERE6"<<endl;
		makeSigmaPositiveSemiDefinite(sigma, snpCount);
              //  cout<<"HERE7"<<endl;
		for (int i = 0; i < snpCount; i++){
			if (abs(stat[i]) > NCP) NCP = abs(stat[i]);
		}
//               cout<<"OK1"<<endl;
		post = new postCal(sigma, stat, snpCount, totalCausalSNP, snpNames, gamma);
	}
	void run() {
        	post->findOptimalSetGreedy(stat, NCP, pcausalSet, rank, rho, outputFileName);
	}
	void finishUp() {
		ofstream outputFile;
                string outFileNameSet = string(outputFileName)+"_set";
                outputFile.open(outFileNameSet.c_str());
                for(int i = 0; i < snpCount; i++) {
                        if(pcausalSet[i] == '1')
                                outputFile << snpNames[i] << endl;
                }
                post->printPost2File(string(outputFileName)+"_post");
                //output the histogram data to file
                if(histFlag)
                	post->printHist2File(string(outputFileName)+"_hist");
	}

	void printLogData() {
		//print likelihood
		//print all possible configuration from the p-causal set
		post->computeALLCausalSetConfiguration(stat, NCP, pcausalSet,outputFileName+".log");
	}

	~CaviarModel() {
	}

};
 
#endif
