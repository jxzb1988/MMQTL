#ifndef _CONDITIONAL_ANALYSIS_H
#define _CONDITIONAL_ANALYSIS_H

#include<iostream>
#include<fstream>
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <regex>
#include <iomanip>


using namespace std;
using namespace arma;

int ReadInputFile(string input, vector<string> &output);
int ReadInputFile(string input, map<string, string> &output);
int ReadInputFile(string input, vector<string> &output, vector<string> &allele_file);

int ReadInputFile(string input, map<string, double> &output);

void remove_file (vector<string> &file);
void remove_file (vector<string> &file, string format);
void remove_file (string file);
void read_table(string data, vector< vector<double> > &test,  vector<string> &ind, vector<string> &variant);
void read_table(string data, vector< vector<double> > &test,  vector<string> &ind, vector<string> &variant, map<string, int> &VAR);
int match(vector<string> &data, string target);
double calculate_VIF(mat &input,vector<string> &variant,string target,vector<string> &detected);
void readgemma(string input, vector < vector<double> > &out_gemma, vector<string> & variant);

int calculate_Mvalue(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > > &beta_,  vector<vector<vector<double> > > &sd_, int index, string output, vector<vector<double> >  &variant, map<string,bool> &pos, map<string,bool> &pos1, map<int,bool> &data_index);

int create_bimbam(string file,string str);

double Han_adjust(int nstudy, mat &matrix, double Z_score);
void readgemma(string input, vector < vector<double> > &out_gemma, vector<string> & variant, vector<vector<double> > & z_score, vector<vector<double> > &beta, vector<vector<double> > &sd, int index_mixed);
void brent_test();
int  updatefam(vector< vector<string> > &input , mat &value,string output);
void copybim(string input, string output);
void copybim(string input, string output, vector<int> &index);
void copyfam(string input, string output);
void  outputbed(mat & _snp_2, mat & _snp_1, string output);
void readbim(string file,map<string, int> &order);
int readfam(string file);
int readfam(string file,vector<double> &pheno, vector< vector<string> > &Fam_infor);
int readfam(string file,vector<double> &pheno, vector< vector<string> > &Fam_infor,int index);
int readbim(string file);
void readbim(string file,vector<string>&variant);
void extractvariant(string input, string output, double r2_cutoff, string target);
int extractvariant(string input, string output, map<string,bool> &target1, string target);
void copybed(string input, string output);
int conditional_analysis(string input,string gene, int totalCausalSNP, float rho, bool histFlag, double gamma, string weight, int nthread, string covariate, vector<string> &grm_file, string outputFile, string geno_file);

int brent_method( vector<double>& A1,  vector<double>& A2, vector<double>& A3);


//int conditional_analysis(string file, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,string covariate, string grm_file,string outputFile,string geno_file,vector<vector<double> > &z_score, vector<vector<double> > &beta, vector<vector<double> > &sd, int index, vector<string>  &variant, vector<string> &removal, int I, int colocalization, vector<double> &inv_sample);

int conditional_analysis(string file, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,string covariate, string grm_file,string outputFile,string geno_file,vector<vector<double> > &z_score, vector<vector<double> > &beta, vector<vector<double> > &sd, int index, vector<string>  &variant, vector<string> &removal, int I, int colocalization, double &Variance);



//    conditional_analysis(file_[i],    gene,            totalCausalSNP,      rho,      histFlag,        gamma,       weight,     nthread,       covariate,         grm_file_[i],   OutputFile,       geno_file_[i],                     z_score,                          beta,                          sd,         1);

void read_allele(string file, map<string,string> &alelle);

void prepare_gwas_summary(string input, string output,map <string, vector<double>>  &Res_fixed, map<string, bool> &index,vector<string> &target, map<string, bool> &MAP);
void prepare_gwas_summary(map <string, vector<double>>  &Res_fixed, string output, map<string, bool> &index, vector<string> &target);
void prepare_ld(string input,string output, map<string,bool> & index);

void read_cis_summary(string file,vector<vector<double> > &z_score,vector<vector<double> > &beta,vector<vector<double> > &sd,vector<string> &variant, map<string, string> &allele);

int MM_conditional_analysis(vector<string> &input,string gene, int totalCausalSNP, float rho, bool histFlag, double gamma, string weight, int nthread, vector<string> &covariate_, vector<string> &grm_file, string outputFile, vector<string> &geno_file,string meta_mode, vector<string> &allele_file, map<string, bool> & pos, int colocalization, map<string, string> &GWAS_file, map<string, double > & MAF, map<int,bool> &data_index, int nested, bool finemap, double p_cutoff,bool cis_summary, string cor_file, vector<string> &causal, int eQTL_number, bool Han);

//int MM_conditional_analysis(vector<string> &input,string gene, int totalCausalSNP, float rho, bool histFlag, double gamma, string weight, int nthread, string covariate, vector<string> &grm_file, string outputFile, vector<string> &geno_file,string meta_mode, vector<string> &allele_file, map<string, bool> & pos, int colocalization, map<string, string> &GWAS_file, map<string, double > & MAF, map<int,bool> &data_index, int nested, bool finemap, double p_cutoff, bool cis_summary);


//int conditional_analysis(string input,string gene, int totalCausalSNP, float rho, bool histFlag);

//void localization(vector<vector<string> > &variant_, mat &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF);


int localization(vector<vector<string> > &variant_, map <string, vector<double>>  &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF);


int localization(vector<vector<string> > &variant_, map <string, vector<double>>  &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF, string colocalization_method,string gene, bool finemap);

double MetaAnalysis(vector<vector<vector<double> > > &z_score, vector<vector<vector<double> > > &beta_, vector<vector<double> > &sd, int x, mat & Res_fixed, mat &Res_random, string meta_mode, string output, vector<string>  &variant, map<string, bool> &pos,map<string, string> &GWAS_file, map<string, double > &MAF,string gene,vector<double> &test1_1,vector<double> &test2_1, vector<double> &test3_1, map<string, string> &ALLELE, bool write_cor);

double MetaAnalysis(vector<vector<vector<double> > > &z_score, vector<vector<vector<double> > > &beta_, vector<vector<double> > &sd, int x, mat & Res_fixed, mat &Res_random, string meta_mode, string output, vector<string>  &variant, map<string, bool> &pos,map<string, string> &GWAS_file, map<string, double > &MAF,string gene, bool INDEX, bool finemap, map<string, string> &ALLELE, string cor_file, bool write_cor);


void read_cor_file(mat &Cov_, string cor_file, int n);

void ConstructCov(vector<vector<vector<double> > > &z_score, vector<vector<vector<double> > > &beta_, vector<vector<double> > &sd, int x);

double ConstructCovRandom(vector<vector<vector<double> > > &z_score, vector<vector<vector<double> > > &beta_, vector<vector<double> > &sd, int x, mat & Res_fixed, mat &Res_random , string output, vector<vector<string> >  &variant, map<string, bool> &pos, map<string, bool> &pos1, map <string, vector<double>> &Res, map<string, double> &MAF, vector<double> &VARIANCE, map<int,bool> &data_index, vector<double> &test1_1,vector<double> &test2_1, vector<double> &test3_1, map<string, string> & ALLELE, string cor_file, bool write_cor);

#endif
