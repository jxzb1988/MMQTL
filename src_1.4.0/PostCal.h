#ifndef POSTCAL_H
#define POSTCAL_H

#include <iostream>
#include <fstream>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>
#ifndef _UTIL_H
#define _UTIL_H

#include "Util.h"
#endif





using namespace std;
using namespace arma;

void printGSLPrint(mat A, int row, int col);
 
class PostCal{

private:
	double gamma;		// the probability of SNP being causal
	double * postValues;	//the posterior value for each SNP being causal
//	double * sigma;		//the LD matrix
	double * histValues;	//the probability of the number of causal SNPs, we make the histogram of the causal SNPs
	int snpCount;		//total number of variants (SNP) in a locus
	int maxCausalSNP;	//maximum number of causal variants to consider in a locus
	double sigmaDet;	//determinie of matrix

	double totalLikeLihood; //Compute the total likelihood of all causal status (by likelihood we use prior)

	mat sigmaMatrix;
	mat invSigmaMatrix;
	mat statMatrix;
        mat statMatrixtTran;
	double baseValue;			//base value used for calculation for overflow	
	string * snpNames;
        int nthread;
        int number;
        mat X;
        mat stat;
      //  double inv_ve=3.0;
public:
      //  int nthread;
	PostCal(double * stat, int snpCount, int maxCausalSNP, string * snpNames, double gamma,int nthread,int number, mat &X) {
                cout<<"Initiate PostCal"<<endl;
		baseValue = 0;
		this->gamma = gamma;
                this->nthread=nthread;
		this->snpNames = snpNames;
		this-> snpCount = snpCount;
		this-> maxCausalSNP = maxCausalSNP;
                //this->sigma = new double[snpCount * snpCount];
		this-> postValues = new double [snpCount];
		this-> histValues = new double [maxCausalSNP+1];             
                this-> stat       = mat(number,1,fill::zeros); 
                this-> number=number;
	        this-> X     =X;
		statMatrix                 = mat (number, 1,fill::zeros);
		statMatrixtTran            = mat (1, number,fill::zeros);
//		sigmaMatrix         	   = mat (snpCount, snpCount,fill::zeros);
             //   cout<<"Come into PostCal, and it is OK1"<<endl;
//		for(int i = 0; i < snpCount*snpCount; i++)
//			this->sigma[i] = sigma[i];
             //   cout<<"Come into PostCal, and it is OK1_1"<<endl;                
		for(int i = 0; i < snpCount; i++)
                        this->postValues[i] = 0;
             //   cout<<"Come into PostCal, and it is OK1_2"<<endl;
		for(int i= 0; i <= maxCausalSNP;i++)
			this->histValues[i] = 0;
             //   cout<<"Come into PostCal, and it is OK1_3"<<endl;
		for(int i = 0; i < number; i++) {
                	statMatrix(i,0) = stat[i];
        	        statMatrixtTran(0,i) = stat[i];
               //         cout<<"i is: "<<i<<" What's Wrong"<<endl;
               //         cout<<"Value is: "<<stat[i]<<endl;
                        this->stat(i,0)=stat[i];
	        }
	//	cout<<"Come into PostCal, and it is OK2"<<endl;
//		for(int i = 0; i < snpCount; i++) {
  //              	for (int j = 0; j < snpCount; j++)
    //                   		sigmaMatrix(i,j) = sigma[i*snpCount+j];
      // 		}
         //       cout<<"Come into PostCal, and it is OK3"<<endl;
//		invSigmaMatrix = inv(sigmaMatrix);
//		sigmaDet       = det(sigmaMatrix);
         //       cout<<"Come into PostCal, and it is OK4"<<endl;
	
	}
        ~PostCal() {
		delete [] histValues;
		delete [] postValues;
  //              delete [] sigma;
	}

        double likelihood(int * configure) ;
	int nextBinary(int * data, int size) ;
        int decomp(vector<int> &str, int *data, int num);
        string convert_symbol(int  *data, int num,vector<int> &output);

	double computeTotalLikelihood(double * stat,map<int,double>& Weight,int nthread) ;	
	double findOptimalSetGreedy(double * stat, char * configure, int *rank,  double inputRho, map<int,double>& Weight, int nthread);
	string convertConfig2String(int * config, int size);
	void printHist2File(string fileName) {
		exportVector2File(fileName, histValues, maxCausalSNP+1);
	}
	void printPost2File(string fileName) {
		double total_post = 0;
                ofstream outfile(fileName.c_str(), ios::out);	
		for(int i = 0; i < snpCount; i++)
                	total_post += postValues[i];
		for(int i = 0; i < snpCount; i++) {
			outfile << snpNames[i] << "\t" << postValues[i]/total_post << "\t" << postValues[i]/totalLikeLihood << endl;
		}
	}

};
 
#endif
