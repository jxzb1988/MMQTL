#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h> 
#include <armadillo>
#include <stdlib.h>
#include <stdio.h>
#ifndef _UTIL_H
#define _UTIL_H
#include "Util.h"
#endif

#ifndef _POSTCAL_H
#define _POSTCAL_H
#include "PostCal.h"
#endif

#ifndef _TOPKSNP_H
#define _TOPKSNP_H
#include "TopKSNP.h"
#endif

#ifndef _MEQTLPOLYGMODEL_H
#define _MEQTLPOLYGMODEL_H
#include "MeQTLPolyGModel.h"
#endif




#ifndef _PREPARE_INPUT_H
#define _PREPARE_INPUT_H

#include "prepare_input.h"
#endif


#ifndef _CAVIARMODEL_H
#define _CAVIARMODEL_H

#include "CaviarModel.h"
#endif




#include "conditional_function.h"

#include <sys/stat.h>

using namespace std;




void output_version(){

    cout<<"multivariate multiple eQTL, mmQTL"<<endl;
    cout<<"version 1.3.0a"<<endl;
    cout<<"Please report bugs to Biao Zeng <biao.zeng@mssm.edu> or gabriel.hoffman@mssm.edu or panagiotis.roussos@mssm.edu"<<endl;
}

void PrintLicense(void) {
  cout << endl;
  cout <<"Copyright 2020-present, Biao Zeng"<<endl;
  cout << "The Software Is Distributed Under GNU General Public "
       << "License, and  A copy of the GNU General Public License is attached along with this program." << endl;

}


int main( int argc, char *argv[]  ){
        
        
        if (argc <= 1) {
                cout<<"Error, parameter(s) needed"<<endl;
                return EXIT_SUCCESS;
        }
        if (argc==2 && argv[1][0] == '-' && argv[1][1] == 'c') {
                
                return EXIT_SUCCESS;
        }
        int totalCausalSNP = 1;
	double gamma = 0.01;
	float rho = 0.95;
	bool histFlag = false;
        bool primary_only = false;
        bool Han=false;
        int eQTL_number = -1;
	int oc = 0;
        int MM_mode =1;	
	string yFile  = "";
	string outputFileName = "test";
	string geneMapFile = "";	
        string weight = "";
        string covariate_file = "";
        string covariate  ="";
        string grm_test ="";
        string gene="PolyQTL_output";
        vector<string> grm_file;
        string Grm_file;
        string X_file="";
        string geno_file="";    
        string meta_mode="fixed";
        string ldFile="";    
        int number;
        int nthread=1;
        int mode=-1;
        int colocal=0;
        double NCP=6.0;
        bool bfile=true;
        string gwas_summary="";
        string gwas_file   ="";
        string input="";
        string anno="";
        map<string, double> MAF;
        string MAF_file="";
        string gene_file="";
        string col_method="eCAVIAR";
        int nested=-1;
        int LD_index=1;
        bool finemap=false;
        int cis_length=1000000;
        int cov_length=3000000;
        double p_cutoff=0.000001;
        bool cis_summary=false;
        string cor_file="";
        bool test_only=false;
        bool joint=false;
        string str="";
        int minimum_tissue=3;
        bool mode_silence=false;
        for(int i = 1; i < argc; i++) 
          {
                if (strcmp(argv[i], "-number")==0 || strcmp(argv[i], "--number")==0 || strcmp(argv[i], "-n")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
  //                      str.assign(argv[i]);
                        number= atoi(argv[i]);
                }
               if (strcmp(argv[i], "-minimum_tissue")==0 || strcmp(argv[i], "--minimum_tissue")==0 ) {
                     if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                     minimum_tissue= atoi(argv[i]);
                }
                else if (strcmp(argv[i], "-silence")==0) {
                        mode_silence=true;
                }
                
                else if (strcmp(argv[i], "-license")==0  || strcmp(argv[i], "--license")==0) {
                       PrintLicense();
                       return EXIT_SUCCESS;
                }

               else if(strcmp(argv[i], "-output")==0 || strcmp(argv[i], "--output")==0 || strcmp(argv[i], "-o")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        outputFileName =str;
                }
      
              

                else if(strcmp(argv[i], "-gene")==0 || strcmp(argv[i], "--gene")==0 || strcmp(argv[i], "-T")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        gene =str;
                }

                else if(strcmp(argv[i], "-input_phe")==0 || strcmp(argv[i], "--input_phe")==0 || strcmp(argv[i], "-P")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        input =str;
                }

               else if(strcmp(argv[i], "-cor_file")==0 || strcmp(argv[i], "--cor_file")==0 || strcmp(argv[i], "-E")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        cor_file =str;
                }

               else if(strcmp(argv[i], "-nested")==0 || strcmp(argv[i], "--nested")==0 || strcmp(argv[i], "-N")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
//                        str.assign(argv[i]);
                        nested =atoi(argv[i]);
                }

               

               else if(strcmp(argv[i], "-LD_index")==0 || strcmp(argv[i], "--LD_index")==0 || strcmp(argv[i], "-i")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
           //             str.clear();
             //           str.assign(argv[i]);
                        LD_index =atoi(argv[i]);
                }
 
               else if(strcmp(argv[i], "-gene_file")==0 || strcmp(argv[i], "--gene_file")==0 || strcmp(argv[i], "-L")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        gene_file =str;
                }

               else if(strcmp(argv[i], "-colocal")==0 || strcmp(argv[i], "--colocal")==0 || strcmp(argv[i], "-I")==0) {
//                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        colocal =1;
                }


                else if(strcmp(argv[i], "-finemap")==0 || strcmp(argv[i], "--finemap")==0 || strcmp(argv[i], "-F")==0) {
  //                      if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        finemap =true;
                }

                else if(strcmp(argv[i], "-cis_summary")==0 || strcmp(argv[i], "--cis_summary")==0 || strcmp(argv[i], "-D")==0) {
    //                    if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        cis_summary =true;
                }

                else if(strcmp(argv[i], "-primary_only")==0 || strcmp(argv[i], "--primary_only")==0) {
                                                primary_only =true;
                }
    

                else if(strcmp(argv[i], "-cis_length")==0 || strcmp(argv[i], "--cis_length")==0 || strcmp(argv[i], "-V")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        cis_length =atoi(argv[i]);
                }
                else if(strcmp(argv[i], "-eQTL_number")==0 || strcmp(argv[i], "--eQTL_number")==0 ) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        eQTL_number =atoi(argv[i]);
                }
               
 
                else if(strcmp(argv[i], "-version")==0 || strcmp(argv[i], "--version")==0 || strcmp(argv[i], "-v")==0) {
                         output_version();
                         exit(0);
                }
               
                else if(strcmp(argv[i], "-help")==0 || strcmp(argv[i], "--h")==0 || strcmp(argv[i], "-h")==0) {
                 
                                cout << "Options: " << endl;
                                cout << "-h, --help                     show this help message and exit " << endl;
                                cout << "-I,                     Perform colocalization analysis " << endl;
                //                cout << "-b,                      " << endl;
                                cout << "-o OUTFILE, --out=OUTFILE      specify the output file" << endl;
//                                cout << "-l LDFILE, --ld_file=LDFILE    the ld input file" << endl;
                                cout << "-p yFILE, --y_file=yFILE       phenotype" << endl;
  //                              cout << "-r RHO, --rho-prob=RHO         set $pho$ probability (default 0.95)" << endl;
//                                cout << "-g GAMMA, --gamma              set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
                                cout << "-c causal                      set the maximum number of causal SNPs" << endl;
                                cout << "-C covariate, --covariate      set the covariate matrix "<<endl;
                                cout << "-a gene annotation file"<<endl;
                                cout << "-geno_file FILE, --geno_file FILE	the file containing the path of plink-format genotype file " <<endl;
//                                cout << "-T output folder, --target      gene name "<<endl;
//                                cout << "-M 0/1,       run in meta manner or not "<<endl;
                                cout << "-grm_file, --grm_file files containing GRM file paths "<<endl;
                                cout << "-A Meta analysis pattern, either fixed or random "<<endl;
//                                cout << "-f 1                           to out the probaility of different number of causal SNP" << endl;
//                                cout << "-w Weight file, --weight       set the biological annotation to use" << endl;
                                cout << "-V cis windown length       " << endl;
                                cout << "-Han apply Han's method to adjust results of random-effect model"<<endl;
                                cout << "-B window to estimate covariance,        " << endl;
                                cout << "-S GWAS summary statistic results, --weight       Summary statistics to perform colocalization analysis" << endl;
  //                              cout << "-n Number of samples, --number       set the biological annotation to use" << endl;
    //                            cout << "-x genotype information,       genotype file for explored variants" << endl;
//                                cout << "-t threads to use, --nthread       set the threads to use" << endl;
                                exit(0);
                }
               else if(strcmp(argv[i], "-cov_length")==0 || strcmp(argv[i], "--cov_length")==0 || strcmp(argv[i], "-B")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        cov_length =atoi(argv[i]);
                }
                else if(strcmp(argv[i], "-Han")==0 || strcmp(argv[i], "--Han")==0 ) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {Han=true;continue;}
                //        ++i;
//                        Han=true;
                }
 
               else if(strcmp(argv[i], "-gwas_summary")==0 || strcmp(argv[i], "--gwas_summary")==0 || strcmp(argv[i], "-S")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        gwas_summary =str;
                }
 
               else if(strcmp(argv[i], "-gwas_file")==0 || strcmp(argv[i], "--gwas_file")==0 || strcmp(argv[i], "-s")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        gwas_file =str;
                }
 
               else if(strcmp(argv[i], "-yFile")==0 || strcmp(argv[i], "--yFile")==0 || strcmp(argv[i], "-p")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        yFile =str;
                }


               else if(strcmp(argv[i], "-rho")==0 || strcmp(argv[i], "--rho")==0 || strcmp(argv[i], "-r")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
  //                      str.assign(argv[i]);
                        rho =atof(argv[i]);
                }

               else if(strcmp(argv[i], "-xfile")==0 || strcmp(argv[i], "--xfile")==0 || strcmp(argv[i], "-x")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        X_file =str;
                }

               else if(strcmp(argv[i], "-p_cutoff")==0 || strcmp(argv[i], "--p_cutoff")==0 || strcmp(argv[i], "-X")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
  //                      str.assign(argv[i]);
                        p_cutoff =atof(argv[i]);
                }

              else if(strcmp(argv[i], "-meta_mode")==0 || strcmp(argv[i], "--meta_mode")==0 || strcmp(argv[i], "-A")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        meta_mode =str;
                }

              else if(strcmp(argv[i], "-anno")==0 || strcmp(argv[i], "--anno")==0 || strcmp(argv[i], "-a")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        anno =str;
                }

             else if(strcmp(argv[i], "-totalCausalSNP")==0 || strcmp(argv[i], "--totalCausalSNP")==0 || strcmp(argv[i], "-c")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
//                        str.assign(argv[i]);
                        totalCausalSNP =atoi(argv[i]);
                }

             else if(strcmp(argv[i], "-gamma")==0 || strcmp(argv[i], "--gamma")==0 || strcmp(argv[i], "-g")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
  //                      str.assign(argv[i]);
                        gamma =atof(argv[i]);
                }

             else if(strcmp(argv[i], "-weight")==0 || strcmp(argv[i], "--weight")==0 || strcmp(argv[i], "-w")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        weight =str;
                }
 
              else if(strcmp(argv[i], "-nthread")==0 || strcmp(argv[i], "--nthread")==0 || strcmp(argv[i], "-t")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
//                        str.clear();
//                        str.assign(argv[i]);
                        nthread =atoi(argv[i]);
                }

              else if(strcmp(argv[i], "-MM_mode")==0 || strcmp(argv[i], "--MM_mode")==0 || strcmp(argv[i], "-M")==0) {
//                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        MM_mode =1;
                }
 
              else if(strcmp(argv[i], "-MAF_file")==0 || strcmp(argv[i], "--MAF_file")==0 || strcmp(argv[i], "-m")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        MAF_file =str;
                }


              else if(strcmp(argv[i], "-histFlag")==0 || strcmp(argv[i], "--histFlag")==0 || strcmp(argv[i], "-f")==0) {
//                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        histFlag =true;
                }

               else if(strcmp(argv[i], "-test_only")==0 || strcmp(argv[i], "--test_only")==0 ) {
                        test_only=true;
                }
              else if(strcmp(argv[i], "-joint")==0 || strcmp(argv[i], "--joint")==0 ) {
                        joint=true;
                }
              else if(strcmp(argv[i], "-bfile")==0 || strcmp(argv[i], "--bfile")==0 || strcmp(argv[i], "-b")==0) {
//                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        bfile =true;
                }

              else if(strcmp(argv[i], "-grm_file")==0 || strcmp(argv[i], "--grm_file")==0 || strcmp(argv[i], "-G")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        grm_file.push_back(str);
                }
 
              else if(strcmp(argv[i], "-grm_file_single")==0 || strcmp(argv[i], "--grm_file_single")==0 || strcmp(argv[i], "-R")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        Grm_file=str;
                }


              else if(strcmp(argv[i], "-ldFile")==0 || strcmp(argv[i], "--ldFile")==0 || strcmp(argv[i], "-l")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        ldFile=str;
                }

              else if(strcmp(argv[i], "-covariate")==0 || strcmp(argv[i], "--covariate")==0 || strcmp(argv[i], "-C")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                        covariate_file=str;
                }


              else if(strcmp(argv[i], "-geno_file")==0 || strcmp(argv[i], "--geno_file")==0 || strcmp(argv[i], "-Z")==0) {
                        if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
                        ++i;
                        str.clear();
                        str.assign(argv[i]);
                       geno_file=str;
                }
          }
//        cout<<"cis_length is: "<<cis_length<<endl;
/*
	while ((oc = getopt(argc, argv, "vhMDYbFIl:t:o:E:X:x:i:Z:C:B:V:z:p:L:s:P:n:w:g:r:R:N:c:G:w:f:m:A:a:T:S:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 1.0:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
                                cout << "-I,                     Perform colocalization analysis " << endl;
                                cout << "-b,                      " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-p yFILE, --y_file=yFILE	phenotype" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability (default 0.95)" << endl;
				cout << "-g GAMMA, --gamma		set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
                                cout << "-C covariate, --covariate      set the covariate matrix "<<endl;
                                cout << "-a gene annotation file"<<endl;
                                cout << "-T output folder, --target      gene name "<<endl;
                                cout << "-M 0/1,       run in meta manner or not "<<endl;
                                cout << "-G genetic relatedness matrix, --GRM  set the genetic relatedness matrix "<<endl;
                                cout << "-A Meta analysis pattern, either fixed or random "<<endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
                                cout << "-w Weight file, --weight       set the biological annotation to use" << endl;
                                cout << "-V cis windown length       " << endl;
                                cout << "-B window to estimate covariance,        " << endl;
                                cout << "-S GWAS summary statistic results, --weight       Summary statistics to perform colocalization analysis" << endl;
                                cout << "-n Number of samples, --number       set the biological annotation to use" << endl;
                                cout << "-x genotype information,       genotype file for explored variants" << endl;
                                cout << "-t threads to use, --nthread       set the threads to use" << endl;
				exit(0);
                        case 'n':
                                number = atoi(optarg);
                                break;
			case 'o':
				outputFileName = string(optarg);
				break;
                        case 'T':
                                gene           =string(optarg);
                                break;
                        case 'P':
                                input          = string(optarg);
                                break;
                        case 'E':
                                cor_file          = string(optarg);
                                break;  
                        case 'N':
                                nested = atoi(optarg);
                                break;
                        case 'L':
                                gene_file = string(optarg);
                                break;
                        case 'I':
                                colocal =1;
                                break;
                        case 'i':
                                LD_index=atoi(optarg);
                                break;
                        case 'F':
                                finemap =true;
                                break;
                        case 'D':
                                cis_summary =true;
                                break;
                        case 'V':
                                 cis_length    = atoi(optarg);
                                 break;
                        case 'B':
                                 cov_length    = atoi(optarg);
                                 break;
                        case 'S':
                                gwas_summary = string(optarg);
                                break;
		        case 's':
                                gwas_file    = string(optarg);
                                break;
			case 'p':
				yFile = string(optarg);
				break;
			case 'r':
				rho = atof(optarg);
				break;
                        case 'x':
                                X_file = string(optarg);
                                break;
                        case 'X':
                                p_cutoff = atof(optarg);
                                break;
                        case 'A':
                                meta_mode = string(optarg);
                                break;
                        case 'a':
                                anno = string(optarg);
                                break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'g':
				gamma = atof(optarg);
				break;
                        case 'w':
                                weight=string(optarg);
                                break;
                        case 't':
                                nthread=atoi(optarg);
                                break;
                        case 'M':
                                MM_mode=1;
                                break;
                        case 'm':
                                MAF_file=string(optarg);
                                break;
                        case 'Y':
                                test_only=true;
                                break;
			case 'f':
                                histFlag = true;
                                break;
                        case 'b':
                                bfile = true;
                                break;
                        case 'G':
                                grm_test     =string(optarg);
                                grm_file.push_back(grm_test);
                                break;
                        case 'R':
                                Grm_file     =string(optarg);
                                break;
                        case 'l':
                                ldFile     =string(optarg);
                                break;
                        
                        case 'C':
                                covariate=string(optarg);
                                break;
                        case 'Z':
                                geno_file=string(optarg);
                                break;
			case ':':
			case '?':
			default:
				cout << "Strange" << endl;
				break;
		}
	}
*/
        if(primary_only)
         {
           eQTL_number=0;
         }
        LD_index--;
        if(gene_file !="")
         {
           vector<string> gene_list;
           ReadInputFile(gene_file,gene_list);
           for(int i=0;i<gene_list.size();i++)
            {
           //    string gene = gene_list[i];   //Changed by Biao Zeng, 10/08/2020, 06:16 am
               gene.clear();
               gene = gene_list[i];
               if(gene=="")
         {
           cout<<"Gene ID is neded to be provided, EXIT"<<endl;
           exit(EXIT_FAILURE);
         }
        ifstream check_dir2("output/");
        if (!check_dir2) {
                mkdir("output", S_IRWXU|S_IRGRP|S_IROTH);
        }
        outputFileName =   gene+"/"+outputFileName;
        ifstream check_dir3(gene.c_str());
        if (!check_dir3) {
                mkdir(gene.c_str(), S_IRWXU|S_IRGRP|S_IROTH);
        }
        map<string,string> gwas_file_;
        if(colocal==1)
         {
           ReadInputFile(gwas_file,gwas_file_);
         }
        if(MM_mode!=0)
          {
            vector<string> input_;
            vector<string> grm_file_;
            vector<string> geno_file_;
            vector<string> Geno_file_;
            vector<string> allele_file_;
            vector<string> covariate_file_;
            vector<string> y;
            vector<string> causal;
            map<int, bool> data_index;
            map<string, bool> pos;
           if(!bfile)
            {
              ReadInputFile(input, input_, allele_file_);
              ReadInputFile(Grm_file, grm_file_);
              ReadInputFile(geno_file, geno_file_);

            } else
            {
              if(!cis_summary)
               {
             if(!test_only)
              {
              int state =Prepare_input(input, geno_file, gene, anno, input_, allele_file_,y,pos,data_index,colocal, finemap, cis_length, cov_length,covariate_file_, LD_index);
              if(state!=0)
               {
                 cout<<"Please have a check on input genotype or phenotype"<<endl;
                 return -1;
               }
               } else
               {
                 ReadInputFile(input, input_);
                 if(input_.size()<=0)
                  {
                    return -1;
                  }
                 for(int i=0;i<input_.size();i++)
                  {
                    allele_file_.push_back(string("NA"));
                    y.push_back(input_[i]);
                    data_index.insert(pair<int, bool>(i, true));
                  }
               }


              }
              if(Grm_file!="")
               {
                 ReadInputFile(Grm_file, grm_file_);
               } else
               {
                 {
                   grm_file_.push_back(string("NA"));
                 }
               }
            }
                  for(int i=0;i<input_.size();i++)
                   {
                     Geno_file_.push_back(string("NA"));
//                     covariate_.push_back(string("NA"));
                   }

             if(covariate_file!="")
               {

                 ReadInputFile(covariate_file, covariate_file_);
               } else
               {
                 for(int i=0;i<input_.size();i++)
                 {
                   covariate_file_.push_back(string("NA"));
                 }
               }



            if(grm_file_.size()==0)
             {
               for(int i=0;i<input_.size();i++)
                {
                  grm_file_.push_back(string("NA"));
                  
                }
             }

            if(covariate_file_.size()==0)
             {
               for(int i=0;i<input_.size();i++)
                {
                  covariate_file_.push_back(string("NA"));
                }
             }
          
//           cout<<"Output covariate file"<<endl;
//           for(int i =0;i<input_.size();i++)
//            {
//              cout<<covariate_file_[i]<<endl;
//            }

             if(!test_only)
              {
             
              
               MM_conditional_analysis(y,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate_file_, grm_file_,outputFileName, Geno_file_, meta_mode, allele_file_, pos, colocal, gwas_file_, MAF, data_index, nested, finemap,p_cutoff, cis_summary,cor_file, causal, eQTL_number, minimum_tissue, Han);


              } else
              {
               cout<<"The analysis is just for test"<<endl;
//               causal.push_back(string("rs11899598"));
//               causal.push_back(string("rs6747985"));
               causal.push_back(string("rs832190"));
               causal.push_back(string("rs72875180"));
               causal.push_back(string("rs73120888"));
              }
               if(cis_summary)
                {
                  return 0;
                }
            vector<string> geno_file_1;
            ReadInputFile(geno_file, geno_file_1);
           if(causal.size()>1)
            { 
            for(int J=0;J<causal.size();J++)
             {
//                cout<<"Analysis for peak variant"<<J<<endl;
                map<string, bool> Causal;
                vector<string> covariate_1;
                for(int K=0;K<causal.size();K++)
                 {
                   if(J==K)
                    {
                      Causal.insert(pair<string, bool>(causal[J],true));
                    } else
                    {
                       Causal.insert(pair<string, bool>(causal[J],false));
                    }
                 }
                vector<string> geno_file_1;
                vector<string> phe_file_1;
                vector<string> grm_file_1;    
                vector<string> allele_file_1; 
                for(int L=0;L<Geno_file_.size();L++)
                 {
                   string Output=outputFileName+string("_geno_causal")+to_string(L);
                   int result=extractvariant(Geno_file_[L], Output, Causal, causal[J]);
                   if(result==0)
                    {
                      covariate_1.push_back(Output+string("_cov"));                      
                      geno_file_1.push_back(Output);
                      phe_file_1.push_back(input_[L]); 
                      grm_file_1.push_back(grm_file_[L]);
                      allele_file_1.push_back(allele_file_[L]);
                    } else if(result==-1)
                    {
                      covariate_1.push_back(string("NA"));
                      geno_file_1.push_back(Output);
                      phe_file_1.push_back(input_[L]);
                      grm_file_1.push_back(grm_file_[L]);
                      allele_file_1.push_back(allele_file_[L]);
                    } 

                    
                 }
                string geno_file_1_causal=outputFileName+string("_geno_file_causal_")+to_string(J);
                
//                string A;
                fstream  FHOU(geno_file_1_causal.c_str(),ios::out);
                if(!FHOU)
                 {
                   cout<<"Error for input or output"<<endl;
                   return (1);
                 }
//                for(map<string,bool>::iterator i=included.begin();i!=included.end();i++)
                for(int i=0;i<geno_file_1.size();i++)
                 {
                      FHOU<<geno_file_1[i]<<endl;
                 }
               FHOU.clear();
               FHOU.clear();

               string input_1_causal=outputFileName+string("_phe_file_causal_")+to_string(J);


                fstream  FHOU1(input_1_causal.c_str(),ios::out);
                if(!FHOU1)
                 {
                   cout<<"Error for input or output"<<endl;
                   return (1);
                 }
//                for(map<string,bool>::iterator i=included.begin();i!=included.end();i++)
                for(int i=0;i<phe_file_1.size();i++)
                 {
                      FHOU1<<phe_file_1[i]<<endl;
                 }
               FHOU1.clear();


               vector<string> causal_1; 
               vector<string> input_1;
               vector<string> y_1;
               for(int i=0;i<y.size();i++)
                {
                  y_1.push_back(y[i]);
                }
               
                int LD_index_test=-1;
                colocal=-1;
                int state =Prepare_input(input_1_causal, geno_file_1_causal, gene, anno, input_1, allele_file_,y,pos,data_index,colocal, false, cis_length, cov_length, covariate_1, LD_index); 
                                      
                


//Prepare_input(string input_1_causal, string geno_file_1_causal, string gene, string anno, string input_1, vector<string> &allele_file_,vector<string> &y,map<string, bool> &pos, map<int, bool> &data_index,int colocal, bool finemap, int cis_length, int cov_length, vector<string> covariate_1, int LD_index);














 
               string cor_file_1=outputFileName+string("_peak_1_statistical_signal_summary_correlation");
               
                MM_conditional_analysis(y_1,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate_1, grm_file_1,outputFileName, Geno_file_, meta_mode, allele_file_1, pos, colocal, gwas_file_, MAF, data_index, nested, finemap,p_cutoff, cis_summary,cor_file_1, causal_1, eQTL_number, minimum_tissue, Han);



              }
             }
            for(int i=0;i<y.size();i++)
             {
               string file1=gene+string("/affected_allele_tissue_")+to_string(i);
               char c[file1.size() + 1];
               strcpy(c, file1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
             }
            for(int i=0;i<y.size();i++)
             {
               if(!(y[i] =="NA"))
                {
                  char c[y[i].size() + 1];
                  strcpy(c, y[i].c_str());
                  if(remove(c))
                   {
                     cout<<"Successfully delete"<<endl;
                   }
                }
             }
//   output_of_GEMMA_7_ENSG00000034053.assoc.txt



            for(int i=0;i<y.size();i++)
             {
               string file1=string("output")+string("/output_of_GEMMA_")+to_string(i)+gene+string(".assoc.txt");
               char c[file1.size() + 1];
               strcpy(c, file1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
             }



            string ldfile=gene+string("/ldFile");
            char c1[ldfile.size() + 1];
               strcpy(c1, ldfile.c_str());
               if(remove(c1))
                {
                  cout<<"Successfully delete"<<endl;
                }
            
            string ldfile2=gene+string("/ldFile2");
            char c2[ldfile2.size() + 1];
               strcpy(c2, ldfile2.c_str());
               if(remove(c2))
                {
                  cout<<"Successfully delete"<<endl;
                }
              
            }
        }
       return 0;
     }
//      cout<<"Run analysis for gene "<<gene<<endl;
        if(gene=="")
         {
           cout<<"Gene ID is neded to be provided, EXIT"<<endl;
           exit(EXIT_FAILURE);
         }

        ifstream check_dir2("output/");
        if (!check_dir2) {
                mkdir("output", S_IRWXU|S_IRGRP|S_IROTH);
        }

        outputFileName =   gene+"/"+outputFileName;

        ifstream check_dir3(gene.c_str());
        if (!check_dir3) {
                mkdir(gene.c_str(), S_IRWXU|S_IRGRP|S_IROTH);
        }
       

        map<string,string> gwas_file_;
        if(colocal==1)
         {
           ReadInputFile(gwas_file,gwas_file_);
          if(col_method =="coloc")
           {
             ReadInputFile(MAF_file, MAF);
           }
         }
        if(MM_mode!=0)
          {
            vector<string> input_;
            vector<string> grm_file_;
            vector<string> geno_file_;
            vector<string> Geno_file_;
            vector<string> allele_file_;
            vector<string> covariate_file_;
            vector<string> y;
            map<string, bool> pos;
            map<int, bool> data_index;
            map<int, bool> data_index1;
           if(!bfile)
            {
              ReadInputFile(input, input_, allele_file_); 
              ReadInputFile(Grm_file, grm_file_);
              ReadInputFile(geno_file, geno_file_);          
              
            } else
            {
//             cout<<"Come to prepare input file"<<endl;
//             cout<<"cis_summary is: "<<cis_summary<<endl;
//             if(!test_only)
//              {
              if(!cis_summary)
               {
                 int state =Prepare_input(input, geno_file, gene, anno, input_, allele_file_,y,pos,data_index,colocal, finemap, cis_length, cov_length,covariate_file_, LD_index);
                 if(state!=0)
                  {
                    cout<<"Please have a check on input genotype or phenotype"<<endl;
                    return -1;
                  }
               } else
               {
//                 cout<<"Perform meta-analysis on eQTL summary result"<<endl;
                // cout<<"Find summary result place"<<endl;
                 ReadInputFile(input, input_);
                // cout<<"Summary result file finding is over"<<endl;
                // cout<<"Number of summary results is: "<<input_.size()<<endl;
                 if(input_.size()<=0)
                  {
                    return -1;
                  }
                 for(int i=0;i<input_.size();i++)
                  {
                    allele_file_.push_back(string("NA"));
                    y.push_back(input_[i]);
                    data_index.insert(pair<int, bool>(i, true));
                    data_index1.insert(pair<int, bool>(i, true));
                  }
               }
//              }
//              cout<<"Prepare input is over"<<endl;
              if(Grm_file=="")
               {
                for(int i=0;i<input_.size();i++)
                 {
                   grm_file_.push_back(string("NA"));
                 }
               } else
               {
                 vector<string> grm_file_tem;
                 ReadInputFile(Grm_file, grm_file_tem);
                 for(int i=0;i<grm_file_tem.size();i++)
                  {
                    map<int,bool>::iterator ITE;
                    ITE=data_index.find(i);
                    if(ITE!=data_index.end())
                     {
                       if(ITE->second)
                        {
                          grm_file_.push_back(grm_file_tem[i]);
                        }
                     }
                  }
               }
//               cout<<"Prepare GRM file is over"<<endl;
              }
              for(int i=0;i<input_.size();i++)
               {
                 Geno_file_.push_back(string("NA"));
//                 covariate_.push_back(string("NA"));
               }

               if(covariate_file!="")
               {

                 ReadInputFile(covariate_file, covariate_file_);
               } else
               {
                 for(int i=0;i<input_.size();i++)
                 {
                   covariate_file_.push_back(string("NA"));
                 }
               }
 


              if(grm_file_.size()==0)
               {
                 for(int i=0;i<input_.size();i++)
                  {
                    grm_file_.push_back(string("NA"));
                  } 
               }
 
              if(covariate_file_.size()==0)
               {
                  for(int i=0;i<input_.size();i++)
                   {
                     covariate_file_.push_back(string("NA"));
                   }
               }
//              cout<<"Output covariate file"<<endl;
//              for(int i =0;i<covariate_file_.size();i++)
//               {
//                 cout<<covariate_file_[i]<<endl;
//               }

              if(input_.size()==0)
               {
                 if(!test_only)
                  {
                    exit(EXIT_FAILURE);
                  }
               }
              vector<string> causal;
//             cout<<"Start to the main function to perform meta-eQTL"<<endl;
            if(!test_only)
             {
               
              MM_conditional_analysis(y,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate_file_, grm_file_,outputFileName, Geno_file_, meta_mode, allele_file_, pos, colocal, gwas_file_, MAF, data_index,nested, finemap,p_cutoff, cis_summary, cor_file, causal, eQTL_number, minimum_tissue, Han);
             } else
             {
//               cout<<"This is just for test"<<endl;
//               causal.push_back(string("rs11899598"));
//               causal.push_back(string("rs6747985"));
               causal.push_back(string("rs832190"));
               causal.push_back(string("rs72875180"));
               causal.push_back(string("rs73120888"));
             }
          //  cout<<"meta eQTL detection is over"<<endl;
            if(cis_summary)
            {
               return 0;
            } 
             
//           if(eQTL_number>0 && causal.size()>=eQTL_number)
//            {
//              return 0;
//            } 

           if(causal.size()>1 && joint)
            {
//            cout<<"Yes, there are multiple eQTL found"<<endl;
            vector<string> geno_file_1;
            ReadInputFile(geno_file, geno_file_1);

//            for(int x=0;x<geno_file_1.size();x++)
//             {
//               cout<<"Genotype file is "<<geno_file_1[x]<<endl;
//             }
            for(int J=0;J<causal.size();J++)
             {
  //              cout<<"Perform joint analysis for peak variant "<<causal[J]<<endl;
                map<string, bool> Causal;
                vector<string> covariate_1;
                for(int K=0;K<causal.size();K++)
                 {
                   if(J==K)
                    {
                      Causal.insert(pair<string, bool>(causal[K],true));
                    } else
                    {
                       Causal.insert(pair<string, bool>(causal[K],false));
                    }
                 }
//                cout<<"size of Causal is: "<<Causal.size()<<endl;
                vector<string> geno_file_2;
                vector<string> phe_file_1;
                vector<string> grm_file_1;
                vector<string> allele_file_1;
                vector<string> input_2; 
                ReadInputFile(input, input_2);
                for(int L=0;L<geno_file_1.size();L++)
                 {
//                   cout<<"Start to prepare genotype for data"<<L<<endl;
                   string Output=outputFileName+string("_causal_")+to_string(J)+string("_data_")+to_string(L);
//                   cout<<"Start to extract genotype, and store in plink binary format file"<<endl;
//                   cout<<"Input is: "<<geno_file_1[L]<<endl;
//                   cout<<"Output is: "<<Output<<endl;
                   int result=extractvariant(geno_file_1[L], Output, Causal, causal[J]);
//                   cout<<"result is "<<result<<endl;
                   if(result==0)
                    {
                      covariate_1.push_back(Output+string(".bed_cov"));
                      geno_file_2.push_back(Output);
                      phe_file_1.push_back(input_2[L]);
                      grm_file_1.push_back(grm_file_[L]);
                      allele_file_1.push_back(allele_file_[L]);
                    } else if(result==-1)
                    {
                      covariate_1.push_back(string("NA"));
                      geno_file_2.push_back(Output);
                      phe_file_1.push_back(input_2[L]);
                      grm_file_1.push_back(grm_file_[L]);
                      allele_file_1.push_back(allele_file_[L]);
                    }


                 }
//                cout<<"Create genotype list file"<<endl;
                string geno_file_1_causal=outputFileName+string("_geno_file_causal_")+to_string(J);

                fstream  FHOU(geno_file_1_causal.c_str(),ios::out);
                if(!FHOU)
                 {
                   cout<<"Error for input or output"<<endl;
                   return (1);
                 }
                for(int i=0;i<geno_file_2.size();i++)
                 {
                      FHOU<<geno_file_2[i]<<endl;
                 }
               FHOU.clear();
               FHOU.close();
               string input_1_causal=outputFileName+string("_phe_file_causal_")+to_string(J);


                fstream  FHOU1(input_1_causal.c_str(),ios::out);
                if(!FHOU1)
                 {
                   cout<<"Error for input or output"<<endl;
                   return (1);
                 }
                for(int i=0;i<phe_file_1.size();i++)
                 {
                      FHOU1<<phe_file_1[i]<<endl;
                 }
               FHOU1.clear();


               vector<string> causal_1;
               vector<string> input_1;
               vector<string> y_1;
//               for(int i=0;i<y.size();i++)
//                {
//                  y_1.push_back(y[i]);
//                }
//               cout<<"Perform backward selection"<<endl;
//               cout<<"input_1_causal is: "<<input_1_causal<<endl;
//               cout<<"geno_file_1_causal is: "<<geno_file_1_causal<<endl;
               int colocal_1=-1;
//                return 0;
//                 
               map<int, bool> data_index2;   
               for(map<int,bool>::iterator it=data_index1.begin();it!=data_index1.end();it++ )
                 {
                   data_index2.insert(pair<int, bool>(it->first,it->second));
                 }
                int state =Prepare_input(input_1_causal, geno_file_1_causal, gene, anno, input_1, allele_file_1,y_1,pos,data_index2,colocal_1, false, cis_length, cov_length, covariate_1, LD_index);
               vector<string> covariate_2;
               for(int I=0;I<covariate_1.size();I++)
                {
                  if(covariate_1[I] =="NA")
                   {
                     covariate_2.push_back(string("NA"));
                   } else
                   {
                     string test=covariate_1[I]+string("_adjust");
                     covariate_2.push_back(test);
                   }
                }
//               for(int I=0;I<Geno_file_.size();I++)
//                {
//                  Geno_file_[I]="NA";
//                }
               string cor_file_1=outputFileName+string("_peak_1_statistical_signal_summary_correlation");
               string outputFileName1=outputFileName+string("_causal")+to_string(J);
//               cout<<"Preparing input file is over"<<endl;
               MM_conditional_analysis(y_1,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate_2, grm_file_1,outputFileName1, Geno_file_, meta_mode, allele_file_1, pos, -1, gwas_file_, MAF, data_index2, nested, false,p_cutoff, cis_summary,cor_file_1, causal_1, eQTL_number, minimum_tissue, Han);
//               cout<<"backward selection for variant "<<J<<" is over"<<endl;
               remove_file(covariate_2);
               remove_file(geno_file_2,"bfile");
               remove_file(covariate_1);
               remove_file(y_1);
//               remove_file(cor_file_1);
               remove_file(geno_file_1_causal);
               remove_file(input_1_causal);

              }
             }





 
            for(int i=0;i<y.size();i++)
             {
               string file1=gene+string("/affected_allele_tissue_")+to_string(i);
               char c[file1.size() + 1];
               strcpy(c, file1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
             }
            
            for(int i=0;i<y.size();i++)
             {
               if(!(y[i] =="NA"))
                {         
                  char c[y[i].size() + 1];
                  strcpy(c, y[i].c_str());
                  if(remove(c))
                   {
                     cout<<"Successfully delete"<<endl;
                   }
                }
             }



            for(int i=0;i<y.size();i++)
             {
               string file1=string("output")+string("/output_of_GEMMA_")+to_string(i)+string("_")+gene+string(".assoc.txt");
               char c[file1.size() + 1];
               strcpy(c, file1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
             }

            string ldfile=gene+string("/ldFile");
            char c1[ldfile.size() + 1];
               strcpy(c1, ldfile.c_str());
               if(remove(c1))
                {
                  cout<<"Successfully delete"<<endl;
                }

            string ldfile2=gene+string("/ldFile2");
            char c2[ldfile2.size() + 1];
               strcpy(c2, ldfile2.c_str());
               if(remove(c2))
                {
                  cout<<"Successfully delete"<<endl;
                }
            string probe0=gene+string("/probe0");
            char c3[probe0.size() + 1];
               strcpy(c3, probe0.c_str());
               if(remove(c3))
                {
                  cout<<"Successfully delete"<<endl;
                }
 
  
            string probe0_phe=gene+string("/probe0_phe");
            char c4[probe0_phe.size() + 1];
               strcpy(c4, probe0_phe.c_str());
               if(remove(c4))
                {
                  cout<<"Successfully delete"<<endl;
                }


            return 0;
          }
        if(input!="")
         {
	   ifstream check_dir(gene);
           if (!check_dir) {
                mkdir(gene.c_str(), S_IRWXU|S_IRGRP|S_IROTH);
             }
           conditional_analysis(input,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,outputFileName, geno_file);                     
           return 0;
         } else
         {
           MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
	   MeQTLPoly.run();
	   MeQTLPoly.finishUp();		
	   return 0;
         }
        if(colocal==1 && gwas_summary!="")
         {
            MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName + "_1", totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
            CaviarModel gwasModel(ldFile, gwas_summary, outputFileName + "_2", totalCausalSNP, NCP, rho, histFlag);
            gwasModel.run();
            gwasModel.finishUp();

            int snpCount  = gwasModel.snpCount;
            string *snpNames  = new string [snpCount];
            double *stat1    = new double [snpCount];
            double *stat2     = new double [snpCount];
            importDataFirstColumn(outputFileName+"_1"+"_post", snpNames, 1);
            importDataNthColumn(outputFileName+"_1"+"_post", stat1, 3, 1);
            importDataNthColumn(outputFileName+"_2"+"_post", stat2, 3, 1);
            ofstream outfile( (outputFileName+"_col").c_str(), ios::out);
            double sumCLPP = 0;
            for(int i = 0; i < snpCount; i++) 
             {
                sumCLPP = sumCLPP + stat1[i] * stat2[i];
             }
            outfile << "SNP_ID\tProb_in_pCausalSet\tCLPP" << endl;
            for(int i = 0; i < snpCount; i++) {
                outfile << snpNames[i] << "\t" << (stat1[i] * stat2[i])/sumCLPP << "\t" << stat1[i] * stat2[i] << endl;
             }
            outfile.close();
         }

  }
