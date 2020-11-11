#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <armadillo>
#include <vector>
using namespace std;
using namespace arma;
int Readbim(string file);
int Readfam(string file);
void Readbim(string bimfile, map<int, vector<int> > &order,vector<string> &variant);
int Readbim(string bimfile, string output_file);
void Read_anno(string gene, string anno, int &chr, int &start, int &end);
void Match_ind(map<string, int> &ind_geno, vector<string> &ind_phe, vector<int> &index, map<string, int> &index2, vector<bool> &INDinc );
bool Match_cis(vector<int> &Order, int chr, int start, int end, int dis);
void Readfam(string file,map<string, int> &ind, vector<string> &individuals);
int Prepare_input(string input,          string geno_file,          string gene, string anno, vector<string> &geno_file_, vector<string> &allele_file_, vector<string> &y, map<string, bool> & pos, map<int,bool> &data_index, int colocalization, bool finemap,int cis_length, int cov_length, vector<string> &cov, int LD_index);

//int Prepare_input(string input_1_causal, string geno_file_1_causal, string gene, string anno, vector<string> &input_1,  vector<string> &allele_file_, vector<string> &y, map<string, bool> &pos,  map<int, bool> &data_index,int colocal,        bool finemap,int cis_length, int cov_length, vector<string> &covariate_1, int LD_index);








int Read_geno(string geno, vector<string> &geno_,vector<vector<string> > &Ind,vector<vector<double> > &Phe ,int chr, int start, int end, map<string, bool> &pos, string gene, map<int, bool> &data_index,int cis_length,int  cov_length, vector<string> &cov, int LD_index);
int Read_phe(string input, vector<vector<string> > &Ind,vector<vector<double> > &Phe, string gene, map<int, bool> data_index);

void outputcov(string cov,vector<bool> &INDinc);
