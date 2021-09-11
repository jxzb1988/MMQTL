#ifndef GMODEL_H
#define GMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>
#include <sstream>
#ifndef _CAVIAR_POSTCAL_H
#define _CAVIAR_POSTCAL_H


#include "PostCal.h"
#endif

using namespace std;
using namespace arma;

class MeQTLPolyGModel{
	
public:
	double rho;
	double gamma;
	int snpCount;
//        int Number_sample;
	int totalCausalSNP;
//	double * sigma;
        double * stat;
        double * stat_snp;
        char * configure;
        int * rank;
        bool histFlag;
	PostCal * post;
	string * snpNames;
//	string ldFile;
        string yFile;
        string outputFileName;
        string geneMapFile;	
        map<int,double> Weight;
        string weight;
        int nthread;
        int number;
        string covariate;
        vector<string> grm_file;
	MeQTLPolyGModel(string yFile, string outputFileName, int totalCausalSNP,double rho,bool histFlag,double gamma,string  weight, int nthread,string covariate, vector<string> grm_file,string X_file )           {
             //   cout<<"Get into CaviarModel"<<endl;
		int tmpSize = 0;
		this->histFlag = histFlag;
		this->rho = rho;
		this->gamma = gamma;
	//	this->ldFile = ldFile;
		this->yFile  = yFile;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
                this->weight = weight;
                this->nthread= nthread;
                this->covariate=covariate;
                this->grm_file      =grm_file;
                this->number        =0;
                this->snpCount      =0;
                vector<mat> grm;
                mat cov;
              
		calculate_variant_number(X_file, snpCount);
            //    cout<<"Get explored variant number is over"<<endl;
             //   cout<<"snpCount is: "<<snpCount<<endl;
                calculate_sample_size(yFile, number);
            //    cout<<"Get sample size is over: "<<endl;
             //   cout<<"Sample size is: "<<number<<endl;
	      //  snpCount = (int)sqrt(tmpSize);
              //  snpCount = (int) tmpSize;
                mat X=mat(number, snpCount, fill::zeros);
            //    cout<<"snpCount is: "<<snpCount<<endl;
            //    cout<<"number is: "<<number<<endl;
                //cout<<"X is: "<<X<<endl;
                import_geno(X_file,X);
              //  cout<<"Extracting genotype is over"<<endl;
         //	sigma     = new double[snpCount * snpCount];
		stat      = new double[number];
                stat_snp  = new double[snpCount];
		configure = new char[snpCount];
		rank      = new int[snpCount];
		snpNames  = new string [snpCount];
          //      cout<<"Come to extract ldFile"<<endl;
	  //	importData(ldFile, sigma);
             //   cout<<"Come to makeSigmaPositiveSemiDefinite"<<endl;
	//	makeSigmaPositiveSemiDefinite(sigma, snpCount);
//                cout<<"Come to importDataFirstRow"<<endl;
		importDataFirstRow(X_file, snpNames);
	//	importDataSecondColumn(yFile, stat);
//	        cout<<"Come to importDataSecondColumn"<<endl;
                importDataFirstColumn(yFile, stat);
                double sum=0;
                for(int i=0;i<number;i++)
                 {
                   sum+=stat[i];
                 }
                sum/=number;
                for(int i=0;i<number;i++)
                 {
                   stat[i]-=sum;
                 }
              //  cout<<"extracting is over"<<endl;
                int peak_index=0;
                double * stat_for_peak;
                stat_for_peak      = new double[number]; 
                for(int i=0;i<number;i++)
                 {
                   stat_for_peak[i]=stat[i];
                 }
  //              cout<<"Get peak signal"<<endl;
                int test=calculate_stat(stat_for_peak,X,peak_index);
  //              cout<<"test is: "<<test<<endl;
  //              cout<<"Peak variant is detected"<<endl; 
  //              cout<<"peak_index is: "<<peak_index<<endl;
                delete stat_for_peak;
               // calculate_stat(stat,X);
                if(grm_file.size()!=0)
                 {
                   importgrm(grm_file, grm, number);
            //       peak_index=9;
               //    cout<<"Output genotype"<<endl;
                //   cout<<"***********************************"<<endl;
                 //  cout<<"X is: "<<X<<endl; 
                 //  cout<<"***********************************"<<endl;
                   mat peak_geno=X.col(peak_index);
                   if(covariate!="")
                     {
                       importcov(covariate,cov);
                       remove_poly(stat,grm,X,peak_geno);
                     } else
                     {
                       remove_poly(stat,grm,X,peak_geno);
                     }
                 }
         //       for(int i=0;i<number;i++)
         //        {
         //          stat_for_peak[i]=stat[i];
         //        }
         //       cout<<"Get peak signal after removing polygenic effect"<<endl;
         //       calculate_stat(stat_for_peak,X,peak_index);

    //            cout<<"Remove polygenetic influence is finished"<<endl;
             //   int peak_index=calculate_stat(stat,X);
             //   if(covariate!="")
             //    {
             //      post = new PostCal(&sigma, &stat_snp, snpCount, totalCausalSNP, &snpNames, gamma,nthread,cov); 
             //    } 
             //   else
             //    {
		   post = new PostCal(stat, snpCount, totalCausalSNP, snpNames, gamma, nthread,number, X);
             //    }
	}
        int calculate_stat(double *stat, mat &X, int &peak)
         {
        //    cout<<"Come to get the peak signal"<<endl;
        //    cout<<"number is: "<<number<<endl;
            mat test=mat(number,1,fill::zeros);
            int row_num=X.n_rows;
              int col_num=X.n_cols;
      //        cout<<"row_num is: "<<row_num<<endl;
      //        cout<<"col_num is: "<<col_num<<endl;
            //int len=stat.size();
            //  cout<<"len is: "<<len<<endl;
      //        cout<<"peak is: "<<peak<<endl;
        //    cout<<"Have a test"<<endl;
            mat stat_test=mat(number,1,fill::zeros);
        //    cout<<"Initiation is OK"<<endl;
            for(int i=0;i<number;i++)
             {
               stat_test(i,0)=stat[i];
             }
            mat beta =mat(snpCount,1,fill::zeros);
            for(int i=0;i<snpCount;i++)
             {
               vec X_test=X.col(i);
               mat X_chosen=mat(number,1,fill::zeros);
               for(int j=0;j<number;j++)
                {
                  X_chosen(j,0)=X_test[j];
                }
               mat XtX=trans(X_chosen)*X_chosen;
               mat inv_XtX=pinv(XtX);
               for(int x=0;x<number;x++)
                {
                }
               mat xy=trans(X_chosen)*stat_test;
               mat res_test=inv_XtX*trans(X_chosen)*stat_test;
               mat resi=stat_test-X_chosen*res_test;
               double  phe_variance=stddev(resi.col(0));
               colvec residuals = stat_test - X_chosen * res_test;
               double s2=0;
               for(int i_x=0;i_x<number;i_x++)
                {
                  s2+=residuals(i_x)*residuals(i_x);
                }
               s2/=number;
               colvec sderr = sqrt(s2 *  diagvec(pinv(trans(X_chosen)*X_chosen))); 
               beta(i,0)=res_test(0,0)/sderr(0);
             }
            double peak_signal=0;
            double *beta_test=new double[snpCount];
            for(int i=0;i<snpCount;i++)
             {
               beta_test[i]=beta(i,0);
               if(abs(beta_test[i])>peak_signal)
                {
                  peak=i;
                  peak_signal=abs(beta_test[i]);
                }
             }
         //   cout<<"peak_signal is: "<<peak_signal<<endl;
         //   cout<<"peak is: "<<peak<<endl;
         //   cout<<"Detecting peak is over"<<endl;
         //   return 0;
         //     cout<<"X is; "<<X<<endl;
         //     cout<<"peak is: "<<peak<<endl;
              row_num=X.n_rows;
              col_num=X.n_cols;
         //     cout<<"row_num is: "<<row_num<<endl;
         //     cout<<"col_num is: "<<col_num<<endl;
              delete beta_test;
              return 0;
         }
        void import_geno(string & X_file, mat & X)
         {
           ifstream fin(X_file.c_str(), std::ifstream::in);
           vector<double> words;
           string word;
           string line;
           int x=0;
           int y=0;
           string::size_type sz;
           getline(fin, line);
           stringstream ss(line);
           y=0;
           string str;
           int num_test=y;
           vector<int> row;
           string test;
           while(fin && getline(fin, line))
            {
              stringstream ss(line);
              y=0;
              while(ss && ss >> word)
               {
                 if(word.compare("NA")==0||word.compare("na")==0)
                  {
                    row.push_back(x*num_test+y);
                    X(x,y)=-999;
                    y++;
                    continue;
                                    }
                 X(x,y)=atof(word.c_str());
                 y++;
               }
              x++;
            }
          for(int i=0;i<y;i++)
            {
              double mean=0;
              int number=0;
              for(int j=0;j<x;j++)
               {
                 if(X(j,i)!=-999)
                  {
                    mean+=X(j,i);
                    number++;
                  }
               }
              mean/=number;
              for(int j=0;j<x;j++)
               {
		 if(X(j,i)==-999)
                  {
                    X(j,i)=0;
                    continue;
                  }
    
                 X(j,i)-=mean;
               //  mat ran=randu<mat>(1,1);
               //  X(j,i)=0.05*ran(0,0);
               }
	      num_test=y;
              if(row.size()>0)
               {
       //          for(int x_t=0;x_t<row.size();x_t++)
       //           {
       //             int j_t=int(row[x_t]/num_test);
       //             int i_t=int(row[x_t]%num_test);
       //             if(i_t==i)
       //              {
       //                X(j_t,i_t)=0;
       //              }
        //          }
               }
            }
     
   //return 0;
  }
         
        void importgrm(vector<string> &grm_file, vector<mat> &grm, int N)
         {
           int index = 0;
           for(int m=0;m<grm_file.size();m++)
            {
              mat  test=mat(N,N,fill::zeros);
              string file=grm_file[m];
              ifstream fin(file.c_str(), std::ifstream::in);
              vector<double> words;
              string line;
              int x=0;
              int y=0;
              string::size_type sz;
              while(fin && getline(fin, line))
               {
                 double word;
                 stringstream ss(line);
                 y=0;
                 while(ss && ss >> word)
                  {
           //         cout << word << "\t";
                    words.push_back(word);
                 //   cout<<"here"<<endl;
                 //   cout<<"x is: "<<x<<", y is: "<<y<<endl;
                    test(x,y)=word;
                 //   cout<<"here1"<<endl;
                    y++;
                 }
               // cout << "(newline)\n";
                x++;
               }
            //  cout<<"mat test is: "<<test<<endl;
              grm.push_back(test);
            }
         //  return(1);
         }
         int importcov(string &cov_file, mat &cov)
         {
         //  cout<<"file is: "<<file<<endl;
           int index = 0;
              ifstream fin(cov_file.c_str(), std::ifstream::in);
              vector<double> words;
              string line;
              int x=0;
              int y=0;
              string::size_type sz;
              while(fin && getline(fin, line))
               {
                 double word;
                 stringstream ss(line);
                 y=0;
                 while(ss && ss >> word)
                  {
   //                 cout << word << "\t";
                    words.push_back(word);
   //                 cout<<"here"<<endl;
   //                 cout<<"x is: "<<x<<", y is: "<<y<<endl;
                    cov(x,y)=word;
   //                 cout<<"here1"<<endl;
                    y++;
                  }
   //             cout << "(newline)\n";
                x++;
               }
           return(1);
         }
       void importDataFirstRow(string fileName, string * list)
        {
          int index = 0;
          string data = "";
          string line = "";
          string word = "";
          ifstream fin(fileName.c_str(), std::ifstream::in);
          if(getline(fin, line) )
           {
             istringstream iss(line);
             while(iss && iss >> word)
              {
                list[index] = word;
                index++;
              }
           }
   //     cout << "FINISH" << endl;
        fin.close();
        }
       void remove_poly(double *stat_ori, vector<mat> &grm , mat &X_all, mat &X)
         {
           int r=grm.size()+1;
           mat Var=mat(r,1,fill::zeros);
           mat stat =mat(number,1,fill::zeros);
           for(int i=0;i<number;i++)
            { 
              stat(i,0)=stat_ori[i];
            }
           for(int x=0;x<X.n_rows;x++)
            {
              for(int y=0;y<X.n_cols;y++)
               {
               }
            }  
           for(int i=0;i<r;i++)
            {
              Var(i,0)=1.0/r;
            }
           int N=stat.n_rows;
           mat A_test =mat(N,N,fill::eye);
           grm.push_back(A_test);
           for(int X_t=0;X_t<grm.size();X_t++)
            {
              mat test=grm[X_t];
              int test_row=test.n_cols;
              int test_col=test.n_rows; 
            }
           vector<mat> A=grm;
           mat AI = mat(r,r,fill::zeros);
           mat S  = mat(r,r,fill::zeros);
           mat s  = mat(r,1,fill::zeros);
           double l_dif =10.0;
           int it =0;
           mat y=stat;
           mat var1=var(y);
           Var =var1(0,0) * Var;
           mat V=mat(N,N,fill::zeros);
           for(int i=0;i<r;i++)
            {
              V =V+Var(i,0)*A[i];
            }
           mat Vinv =pinv(V);
           mat VinvX=Vinv*X;
           mat pinv_Xt_Vinv_X=pinv(trans(X)*Vinv*X);
           mat Xt_Vinv=trans(X)*Vinv;
           mat P =Vinv-Vinv*X*pinv(trans(X)*Vinv*X)*trans(X)*Vinv;
           for(int i=0;i<r;i++)
            {
              mat test5=(Var(i,0)*Var(i,0)*trans(y)*P*A[i]*P*y+sum(diagvec(Var(i,0)*mat(N,N,fill::eye)-Var(i,0)*Var(i,0)*P*A[i])))/N;
              Var(i,0)=test5(0,0);
            }
           V=mat(N,N,fill::zeros);
           for(int i=0;i<r;i++)
            {
              V =V + A[i] * Var(i,0);
            }
           Vinv =pinv(V);
           VinvX=Vinv*X;
        //   cout<<"6 pinv"<<endl;
           pinv_Xt_Vinv_X=pinv(trans(X)*Vinv*X);
           Xt_Vinv=trans(X)*Vinv;
           P=Vinv-Vinv*X*pinv(trans(X)*Vinv*X)*trans(X)*Vinv;
           double value, sign;
           log_det(value,sign,V);
           double value1=value;
           log_det(value,sign,trans(X)*Vinv*X);
           double value2=value;
           mat test_ypy=trans(y)*P*y;
           double test4=-0.5*(value1+value2+test_ypy(0,0));
           double logL=test4;
           it=0;
           l_dif=10.0;
           while ( it < 100 & ( abs(l_dif) >= 0.0001))
            {
          //    cout<<"l_dif is: "<<l_dif<<endl;
          //    cout<<"Var is: "<<Var<<endl;
              if(abs(l_dif) >= 10^-4)
               {
            //     cout<<"Yes1"<<endl;
               }
              if(abs(l_dif)<10^-2)
               {
              //   cout<<"Yes2"<<endl;
               }
              if(-1.0*l_dif >0)
               {
                // cout<<"Yes3"<<endl;
               }
            //  cout<<"Come into REML step"<<endl;
              it = it + 1;
              mat AI =mat(r,r,fill::zeros);
              for (int i=0;i<=r-1;i++)
               {
                 for (int ii=0;ii<=r-1;ii++)
                  {
                    if ( i == r-1 && ii == r-1 )
                     {
                       mat res_test=trans(y) * P * P * P * y;
                       AI(r-1,r-1) = res_test(0,0);
                     }
                    else if ( i == r-1 )
                     {
                       mat res_test=trans(y) * P * P * A[ii] * P * y;
                       AI(r-1,ii) = res_test(0,0);
                     }
                    else if ( ii == r-1 )
                     {
                       mat res_test=trans(y) * P * A[i] * P * P * y;
                       AI(i,r-1) = res_test(0,0);
                     }
                    else
                     {
                       mat res_test=trans(y) * P * A[i] * P * A[ii] * P * y;
                       AI(i,ii) = res_test(0,0);
                     }
                 }
               }
           //   cout<<"AI is: "<<AI<<endl;
           //   cout<<"It is OK1"<<endl;
              AI = 0.5*AI;
              for (int i=0;i<=r-1;i++)
               {
                 if ( i == r-1 )
                  {
                    double Sum=sum(diagvec(( P )));
                    mat res_test=trans(y) * P * P * y;
             //       cout<<"i is: "<<i<<" Sum1 is: "<<Sum<<"  and Sum2 is: "<<res_test<<endl;
                    s(i,0) = sum(diagvec(( P ))) - res_test(0,0);
             //       cout<<"s is: "<<s<<endl;
                  } else
                  {
                    double Sum=sum(diagvec(( P *A[i])));
                    //mat res_test=trans(y) * P * P * y;
                    //cout<<"Sum1 is: "<<Sum<<"  and Sum2 is: "<<res_test<<endl;
                    mat res_test=trans(y) * P * A[i] * P * y;
               //     cout<<"i is: "<<i<<" Sum1 is: "<<Sum<<"  and Sum2 is: "<<res_test<<endl;
                    s(i,0) = sum(diagvec(( P * A[i] ))) - res_test(0,0);
                //    cout<<"s is: "<<s<<endl;
                  }
               }
              s = -0.5*s;
           //   cout<<"Output s matrix"<<endl;
              for(int i_x=0;i_x<s.n_rows;i_x++)
               {
                 for(int j_x=0;j_x<s.n_cols;j_x++)
                  {
             //       cout<<" "<<s(i_x,j_x);
                  }
               //  cout<<endl;
               }
         //     cout<<endl;
         //     cout<<"s is: "<<s<<endl;
         //     cout<<"It is OK2"<<endl;
              if ( l_dif > 1 )
                {
           //       cout<<"l_dif is: "<<l_dif<<endl;
           //       cout<<"8_1 pinv"<<endl;
           //       cout<<"Output AI"<<endl;
                  for(int m=0;m<AI.n_rows;m++)
                   {
                     for(int n=0;n<AI.n_cols;n++)
                      {
             //           cout<<" "<<AI(m,n);
                      }
               //      cout<<endl;
                   }
                //  cout<<endl;
                  Var = Var + 0.316*(pinv(AI) * s);
                }else 
                {
             //     cout<<"l_dif is: "<<l_dif<<endl;
             //     cout<<"8_2 pinv"<<endl;
             //     cout<<"Output AI"<<endl;
                  for(int m=0;m<AI.n_rows;m++)
                   {
                     for(int n=0;n<AI.n_cols;n++)
                      {
               //         cout<<" "<<AI(m,n);
                      }
               //      cout<<endl;
                   }
               //   cout<<endl;

                  Var = Var + pinv(AI) * s;
                }
       //       cout<<"Var is: "<<Var<<endl;
       //       cout<<"It is OK2_1"<<endl;
//              V = mat(N,N,fill::zeros);
              V = 0.001*mat(N,N,fill::eye);
       //       cout<<"It is OK2_2"<<endl;
              for (int i=0;i<=r-1; i++ )
               {
                 V = V + A[i] * Var(i,0);
               }
       //       cout<<"It is OK2_3"<<endl;
       //       cout<<"V is: "<<endl;
       //       for(int i_x=0 ;i_x<V.n_rows;i_x++)
       //        {
       //          for(int j_x=0; j_x<V.n_cols;j_x++)
       //           {
       //             cout<<" "<<V(i_x,j_x);
       //           }
       //          cout<<endl;
       //        }
       //       cout<<endl;
       //       cout<<"Inverse covariance matrix"<<endl;
              for(int i_x=0 ;i_x<V.n_rows;i_x++)
               {
                 for(int j_x=0; j_x<V.n_cols;j_x++)
                  {
       //             cout<<" "<<V(i_x,j_x);
                  }
       //          cout<<endl; 
               }
       //       cout<<endl;
              cout<<"Perform inv of V"<<endl;
              Vinv = inv(V);
       //       cout<<"Inverse is over"<<endl;
       //       cout<<"It is OK3"<<endl;
       //       cout<<"9 pinv"<<endl;
              P = Vinv - Vinv * X * pinv( trans(X) * Vinv * X ) * trans(X) * Vinv;
              log_det(value,sign,V);
              value1=value;
              log_det(value,sign,trans(X)*Vinv*X);
              value2=value;
              test_ypy=trans(y) * P * y;
              double  res_logL=value1 + value2 + test_ypy(0,0);
              
              double new_logL = -0.5 * res_logL;
       //       cout<<"test4 is: "<<new_logL<<endl;
              l_dif = new_logL - logL;
       //       cout<<"l_dif is: "<<l_dif<<endl;
              logL = new_logL;
       //       cout<<"It is OK4"<<endl;
            }
       //    cout<<"Here is OK12"<<endl;
       //    cout<<"N is: "<<N<<endl;
           mat  Cov=mat(N, N, fill::zeros);
           if(Var(0,0)>0)
            {
              
            
           for(int i=0;i<r;i++)
            {
         //     cout<<"Var is: "<<Var(i,0)<<endl;
              Cov+=Cov+grm[i]*Var(i,0);
            }
         //   cout<<"Here is OK13"<<endl;
         //  cout<<"10 pinv"<<endl;
           mat invCov=inv(Cov);
           cx_mat  square_root_invCov=sqrtmat(invCov);
           stat=real(square_root_invCov)*stat;
           X_all   =real(square_root_invCov)*X_all; 
           }
//           cout<<"Transporming is over"<<endl;
  //         cout<<"Transformed phenotype is: "<<endl;
    //       cout<<stat<<endl;
      //     cout<<"Transformed genotype is: "<<endl;
        //   cout<<X_all<<endl;
        //   cout<<"Here is OK14"<<endl;
        //   cout<<"New phenotype is: "<<stat<<endl;
           for(int i=0;i<stat.n_rows;i++)
            {
        //      cout<<"Transformed phenotype is: "<<stat(i,0)<<endl;
              stat_ori[i]=stat(i,0);
            } 
         //  for(int i=0;i<stat.n_rows;i++)
         //   {
         //     stat_ori[i]=stat(i,0);
         //   }
         }
        void extract_weight(string weight, map<int,double>& Weight)
          {
            string line = "";
            string fileName=weight;
            ifstream fin(fileName.c_str(), std::ifstream::in);
            int index=-1;
            while (getline(fin, line))
             {
               string snp;
               index++;
               double weight;
               istringstream iss(line);
               iss >> snp;
               iss >> weight;
               weight = 1 / (1 + exp(7.2 - weight*3));
               Weight.insert(map<int, double>::value_type(index, weight));
             } 
          }
        void calculate_variant_number(string fileName, int & index)
          {
          //  string fileName="test.geno";
           // int index = 0;
            string data = "";
            string line = "";
            string word = "";
            ifstream fin(fileName.c_str(), std::ifstream::in);
            if(getline(fin, line) )
             {
               istringstream iss(line);
               while(iss && iss >> word)
                {
                  index++;
                }
             }
     //   cout << "Number of variants is: " <<index<<endl;
     //   cout << "FINISH" << endl;
            fin.close();
          }
        void calculate_sample_size(string fileName, int & index)
          {
            string data = "";
            string line = "";
            string word = "";
            ifstream fin(fileName.c_str(), std::ifstream::in);
            while(getline(fin, line) )
             {
               istringstream iss(line);
              // while(iss && iss >> word)
              //  {
               index++;
         //      cout<<"Index is: "<<index<<endl;
              //  }
             }
            fin.close();
         //   cout<<"Index is: "<<index<<endl;
          }
	void run() {
        	post->findOptimalSetGreedy(stat, configure, rank, rho, Weight, nthread);
	}
	void finishUp() {
                string outFileNameSet = string(outputFileName)+"_set";
                ofstream outputFile;
                outputFile.open(outFileNameSet.c_str(),ios_base::out);
                for(int i = 0; i < snpCount; i++) {
                        if(configure[i] == '1')
                                outputFile << snpNames[i] << endl;
                }
                post->printPost2File(string(outputFileName)+"_post");
                if(histFlag)
                	post->printHist2File(string(outputFileName)+"_hist");
	}
        ~MeQTLPolyGModel() {
	}

};
 
#endif
