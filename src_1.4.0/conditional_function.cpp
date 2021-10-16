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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <regex>
#include "gemma_param.h"
#include "gemma.h"
#include <iomanip>
#ifndef _MEQTLPOLYGMODEL_H
#define _MEQTLPOLYGMODEL_H
#include "MeQTLPolyGModel.h"
#endif


#ifndef _CAVIARMODEL_H
#define _CAVIARMODEL_H

#include "CaviarModel.h"
#endif


#ifndef _UTIL_H
#define _UTIL_H

#include "Util.h"
#endif



#include "conditional_function.h"
#include <boost/math/tools/minima.hpp>
#include <boost/math/distributions/beta.hpp>


using boost::math::tools::brent_find_minima;

using namespace std;
using namespace arma;


struct func_dist_square
 {
   func_dist_square(const vector<double>& a1,const vector<double>& a2,const vector<double>& a3):A1(a1),A2(a2),A3(a3){};
   double operator()(const double& x)
    {
         double d1 = 0.0;
         double d2 = 0.0;
         for (int b = 0; b < A1.size(); b++)
          {
            d1 += log(A1[b] + x);
          }
         for (int b = 0; b < A2.size(); b++)
          {
            d2 += A2[b] / (A3[b] + x);
          }


          return 0.5 * (A1.size() * log(6.283185307179586) + d1 + d2);
       }
     private:
     vector<double> A1;
     vector<double> A2;
     vector<double> A3;

 };



int brent_method( vector<double>& A1,  vector<double>& A2, vector<double>& A3, double& start, double& end, double& result1, double& result2)
 {
   const int double_bits = numeric_limits<double>::digits;
   pair<double, double> r = brent_find_minima(func_dist_square(A1,A2,A3), start,end, double_bits);
   streamsize precision_1 = cout.precision(numeric_limits<double>::digits10);
//   cout << "x at minimum = " << r.first
//       << ", f(" << r.first << ") = " << r.second << endl;
   vector<double> a;
   double a1=1.0;
   double a2=2.;
   a.push_back(a1);
   a.push_back(a2);
   try
    {
      pair<double, double> r2 = brent_find_minima(func_dist_square(A1,A2,A3), start, end, double_bits);
//      cout << "x at minimum = " << r2.first       << ", f(" << r2.first << ") = " << r2.second << endl;
      result1 = r2.first;
      result2 = r2.second;
      return 0;
    }
   catch(int e)
    {
      cout << "An exception occurred when apply brent method to find minimum value. Exception Nr. " << e <<endl;
      exit(0);
    }

 }



void read_allele(string file,map<string, string> &allele)
 {
//   cout<<"Come into read_allele function"<<endl;
   ifstream FHIN1(file.c_str(),ifstream::in);
   string a;
   int b=0;
   while(FHIN1 &&getline(FHIN1,a))
     {
          stringstream ss(a);
          string word;
          ss>>word;
          string word1;
          ss>>word1;
          allele.insert(pair<string,string>(word,word1));
      }

    FHIN1.clear();
    FHIN1.close();
//   cout<<"Allele infor reading is over"<<endl;
 }




//Read phe+genotype table
void read_table(string data, vector< vector<double> > &test,  vector<string> &ind, vector<string> &variant, map<string, int> & VAR)
  {
    ifstream FHIN1(data.c_str(),ifstream::in);
    string a;
    int b=0;
    while(FHIN1 &&getline(FHIN1,a))
     {
       b++;
       if(b==1)
        {
          int INDEX=-1;
          stringstream ss(a);
          string word;
          ss>>word;
          while(ss && ss>>word)
           {
             INDEX++;
             VAR.insert(pair<string, int>(word, INDEX));
             variant.push_back(word);
           }
         } else
         {
           double c=0.0;
           stringstream ss(a);
           string waste;
           ss>>waste;
           ind.push_back(waste);
           vector<double> T;
           while(ss && ss>>c)
            {
                 T.push_back(c);
            }
           test.push_back(T);
         }
      }
    FHIN1.clear();
    FHIN1.close();
  }





void read_table(string data, vector< vector<double> > &test,  vector<string> &ind, vector<string> &variant)
  {
    ifstream FHIN1(data.c_str(),ifstream::in);
    string a;
    int b=0;
    while(FHIN1 &&getline(FHIN1,a))
     {
       b++;
       if(b==1)
        {
//          int INDEX=-1;
          stringstream ss(a);
          string word;
          ss>>word;
          while(ss && ss>>word)
           {
//             INDEX++;
//             VAR.insert(pair<string, int>(word, INDEX));
             variant.push_back(word);
           }
         } else
         {
           double c=0.0;
           stringstream ss(a);
           string waste;
           ss>>waste;
           ind.push_back(waste);
           vector<double> T;
           while(ss && ss>>c)
            {
                 T.push_back(c);
            }
           test.push_back(T);
         }
      }
    FHIN1.clear();
    FHIN1.close();
  }





//Find the index of a specific string in a string vector
int match(vector<string> &data, string target)
 {
   
   int j=-1;
   int i=-1;
   vector<string>::iterator Col;

   for(Col=data.begin();Col!=data.end();Col++)
    {
      j++;
      if((*Col)==target)
       {
         i=j;
         break;
       }
    }
   return i;
 }

//Calculate the variance inflation factor to determine whether include the explored variant or not
double calculate_VIF(mat &input,vector<string> &variant,string target,vector<string> &detected)
 {
//   cout<<"Come to function calculate_VIF"<<endl;
   int res=match(variant,target);
   mat X_chosen=mat(input.n_rows,detected.size(),fill::zeros);
   vector<int> DELETE;
   for(int i=0;i<detected.size();i++)
    {
      int ind=match(variant, detected[i]);
      if(ind>0)
       { 
         X_chosen.col(i)=input.col(ind);
       } else
       {
         DELETE.push_back(i);
       }
    }
   if(DELETE.size()==detected.size())
    {
      return 1.0;
    }
   if(DELETE.size()>0)
    {
      uvec Indice = conv_to< uvec >::from(DELETE);
      X_chosen.shed_cols(Indice);
    }
//   cout<<"Genotype get, and start to calculate VIF"<<endl;
//   cout<<"genotype y "<<endl;
//   cout<<input.col(res);
//   cout<<"genotype x "<<endl;
//   cout<<X_chosen<<endl;
   mat stat_test=input.col(res);
   mat XtX=trans(X_chosen)*X_chosen;
   mat inv_XtX=pinv(XtX);

   mat res_test=inv_XtX*trans(X_chosen)*stat_test;
   mat resi=stat_test-X_chosen*res_test;
   double  phe_variance=stddev(resi.col(0));
   colvec residuals = stat_test - X_chosen * res_test;
   double s2=0;
   double s=0;
   for(int i_x=0;i_x<input.n_rows;i_x++)
    {
      s2+=residuals(i_x)*residuals(i_x);
      s+= stat_test(i_x,0)*stat_test(i_x,0);
    }
//   double r2=1-s2/s;
    
   double r2=s2/s;
   if(r2<0.00001)
    {
      r2=0.00001;
    }
   double vif=1/r2;
   return vif;
 }




//read GEMMA result output
void readgemma(string input, vector < vector<double> > &out_gemma, vector<string> & variant)
 {
   ifstream FHIN(input.c_str(),ios::in);
   if(!FHIN)
    {
      return ;
    }
   string A;
   getline(FHIN,A);
   while(FHIN && getline(FHIN, A))
    {
         stringstream ss(A);
         vector<string> X;
         string B;
         vector<double> Y;
         ss>>B;
         ss>>B;
         variant.push_back(B);
         for(int i=3;i<8;i++)
          {
            ss>>B;
          }
         double test;
         double test1;
         double test2;
         ss>>test;
         Y.push_back(test);
         ss>>test1;
         ss>>test2;
         ss>>test2;
         Y.push_back(test2);
         if(test1<=0.000000001)
          {
            test=0;
          }else
          {
            test=test/test1;
          }
         if(test<0)
          {
            test=-test;
          }
         Y.push_back(test);
         out_gemma.push_back(Y);
  //       cout<<"Statistical evidence"<<Y[0]<<" "<<Y[1]<<" "<<Y[2]<<endl;
    }
   FHIN.clear();
   FHIN.close();

 }

void readgemma(string input, vector < vector<double> > &out_gemma, vector<string> & variant, vector <vector<double> > &z_score, vector <vector<double> > &beta,vector <vector<double> > &sd, int index_mixed)
 {
   ifstream FHIN(input.c_str(),ios::in);
   if(!FHIN)
    {
      return ;
    }
   string A;
   vector<double> Z_score;
   vector<double> Beta;
   vector<double> Sd;
   getline(FHIN,A);
//   cout<<"index_mixed is: "<<index_mixed<<endl;
   while(FHIN && getline(FHIN, A))
    {
     if(index_mixed==1)
     {
      stringstream ss(A);
      vector<string> X;
      string B;
      vector<double> Y;
      ss>>B;
      ss>>B;
      variant.push_back(B);
      for(int i=3;i<8;i++)
       {
         ss>>B;
       }
      double test;
      double test1;
      double test2;
     //Changed by Biao Zeng, 04/09/2021, 03:03 pm
      double test3;
      ss>>test;
      Y.push_back(test);
      ss>>test1;
      Sd.push_back(test1);
      ss>>test2;
      ss>>test2;
      Beta.push_back(test);
      
      Y.push_back(test2);
      test=test/test1;
      //Changed by Biao Zeng, 04/09/2021, 03:03 pm
      test3=test; 
      if(test<0)
       {
         test=-test;
       }

      Y.push_back(test);
      
      out_gemma.push_back(Y);
      //Changed by Biao Zeng, 04/09/2021, 03:03 pm
//      Z_score.push_back(test);
      Z_score.push_back(test3);
     } else
     {
       stringstream ss(A);
      vector<string> X;
      string B;
      vector<double> Y;
      ss>>B;
      ss>>B;
      variant.push_back(B);
//      cout<<B<<" ";
      for(int i=3;i<9;i++)
       {
         ss>>B;
       }
      double test;
      double test1;
      double test2;
      double test3;
      ss>>test;
      Y.push_back(test);
      ss>>test1;
      Sd.push_back(test1);
      ss>>test2;
//      ss>>test2;
      Beta.push_back(test);

      Y.push_back(test1);
      test=test/test1;
      test3=test;
      if(test<0)
       {
         test=-test;
       }

      Y.push_back(test);

      out_gemma.push_back(Y);
  //    cout<<"Statistical evidence"<<Y[0]<<" "<<Y[1]<<" "<<Y[2]<<endl;
      Z_score.push_back(test3);
     }
    }
   z_score.push_back(Z_score);
   beta.push_back(Beta);
   sd.push_back(Sd);
   FHIN.clear();
   FHIN.close();

 }

int updatefile(vector< vector<string> > &input , mat &value,string output)
 {
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
      return (1);
    }
   
   if(value.n_rows!=input.size()-1)
    {
         return (1);
    }
   FHOU<<input[0][0];
   for(int j=1;j<input[0].size();j++)
    {
      FHOU<<"	"<<input[0][j];    
    }
   FHOU<<endl;
   for(int i=1;i<input.size();i++)
    {
      FHOU<<input[i][0];
      FHOU<<"	"<<value(i-1,0); 
      for(int j=2;j<input[i].size();j++)
       {
         FHOU<<"	"<<input[i][j];
       }
      FHOU<<endl;
    }
   FHOU.clear();
   FHOU.close();
   return (0);    
 }





//Attached new residual phenotype, removing effect of detected SNPs, to plink fam file.
int  updatefam(vector< vector<string> > &input , mat &value,string output)
 {
//   cout<<"Come into updatefam fucntion"<<endl;
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }
   if(value.n_rows!=input.size())
    {
      return (1);
    }
   for(int i=0;i<input.size();i++)
    {
      FHOU<<input[i][0];
      for(int j=1;j<input[i].size();j++)
       {
         FHOU<<"        "<<input[i][j];
       }
      FHOU<<"   "<<value(i,0)<<endl;
    }
   FHOU.clear();
   FHOU.clear();
   return (0);
 }



//Make a copy of plink bim file
void copybim(string input, string output)
 {
   ifstream FHIN(input.c_str(),ios::in);
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHIN || !FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return;
    }
   while(FHIN && getline(FHIN,A))
    {
      FHOU<<A;
      FHOU<<endl;
    }
   FHIN.clear();
   FHOU.clear();
   FHIN.close();
   FHOU.close();
 }

//Extract a specific set of variants, and output a new bim file
void copybim(string input, string output, vector<int> &index)
 {
   ifstream FHIN(input.c_str(),ios::in);
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHIN || !FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return;
    }
   int val=-1;
   while(FHIN && getline(FHIN,A))
    {
      val++;
      if(find(index.begin(), index.end(), val) != index.end())
       {       
         FHOU<<A;
         FHOU<<endl;
       }
    }
   FHIN.clear();
   FHOU.clear();
   FHIN.close();
   FHOU.close();
 }


//Make a copy of fam file
void copyfam(string input, string output)
 {
   ifstream FHIN(input.c_str(),ios::in);
   string A;
   fstream FHOU(output.c_str(),ios::out);
   if(!FHIN || !FHOU)
    {
      cout<<"Error for input or output"<<endl;
    }
   while(FHIN && getline(FHIN,A))
    {
      stringstream ss(A);
      vector<string> B;
      string Test;
      while(ss && ss>>Test)
       {
         FHOU<<Test<<"	"; 
       }
      FHOU<<endl;
    }
   FHIN.clear();
   FHIN.close();
   FHOU.clear();
   FHOU.close();
 }

//Make a copy of bed file
void  outputbed(mat & _snp_2, mat & _snp_1, string output)
 {
   int i=0, pos=0, j=0;
   string OutBedFile=output;
//   cout<<"Come into outputbed to make a bed file"<<endl;
//   cout<<"output file is: "<<output<<endl;
   fstream OutBed(OutBedFile.c_str(), ios::out|ios::binary);
   if(!OutBed) throw("Error: can not open the file ["+OutBedFile+"] to write.");
//   cout<<"Writing genotypes to PLINK BED file ["+OutBedFile+"] ..."<<endl;
//   cout<<"Stpe1 ok"<<endl;
   bitset<8> b;
   char ch[1];
   b.reset();
   b.set(2);  b.set(3);  b.set(5);  b.set(6);
   ch[0] = (char)b.to_ulong();
   OutBed.write(ch,1);
   b.reset();
   b.set(0);  b.set(1);  b.set(3);  b.set(4);
   ch[0] = (char)b.to_ulong();
   OutBed.write(ch,1);
   b.reset();
   b.set(0);
   ch[0] = (char)b.to_ulong();
//   cout<<"Step2 ok"<<endl;
   OutBed.write(ch,1);
   int n_sid=_snp_2.n_rows;
   int n_ind=_snp_2.n_cols;
//   cout<<"Number of individuals is: "<<n_ind<<endl;
//   cout<<"Number of variants is: "<<n_sid;
//   cout<<"Right now, it is OK"<<endl;
   for(i=0; i<n_sid; i++){
      pos=0;
      b.reset();
      for(j=0; j<n_ind; j++){
        b[pos++]=(!_snp_2(i,j));
        b[pos++]=(!_snp_1(i,j));
        if(pos>7 || j==n_ind-1){
          ch[0]=(char)b.to_ulong();
          OutBed.write(ch,1);
          pos=0;
          b.reset();
        }
      }
    }
   OutBed.close();
//   cout<<"Function outputbed over"<<endl;
 }

//Extract a specific set of variants from a bim file 
void readbim(string file,map<string, int> &order)
 {

   int ibuf=-1;
   string A;
   string cbuf="0";
   double dbuf=0.0;
   string str_buf;
   ifstream Bim(file.c_str());
   if(!Bim) throw("Error: can not open the file ["+file+"] to read.");
//   cout<<"Reading PLINK BIM file from ["+file+"]."<<endl;
   while(Bim && getline(Bim, A))
    {
      ibuf++;
      stringstream ss(A);
      string waste;
      ss>>waste;
      ss>>waste;
//      cout<<"waste is: "<<waste<<", and ibuf is: "<<ibuf<<endl;
      order.insert(pair<string,int>(waste,ibuf));
    }
//   cout<<"Read bim file is over"<<endl;
   Bim.close();
 }

//Read a bim file
void readbim(string file,vector<string>&variant)
 {

   int ibuf=-1;
   string A;
   string cbuf="0";
   double dbuf=0.0;
   string str_buf;
   ifstream Bim(file.c_str());
   if(!Bim) throw("Error: can not open the file ["+file+"] to read.");

   while(Bim && getline(Bim, A))
    {
      ibuf++;
      stringstream ss(A);
      string waste;
      ss>>waste;
      ss>>waste;
      variant.push_back(waste);
    }
   Bim.close();
 }


//Read a fam file, and calcuate the individual number 
int readfam(string file)
 {
   ifstream Fam(file.c_str());
   string A;
   if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
   int i=0;
   while(Fam && getline(Fam, A))
    {
      i++;
    }
   Fam.clear();
   Fam.close();
   return i;
 }

//Read a fam file and extract phenotype 




int readfam(string file,vector<double> &pheno, vector< vector<string> > &Fam_infor,int index)
 {
   ifstream Fam(file.c_str());
   string A;
   if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
//   cout<<"Reading PLINK FAM file from ["+file+"]."<<endl;
   while(Fam && getline(Fam, A))
    {
      stringstream ss(A);
      vector<string> T;
      string TT;
  //    for(int j=0;j<5;j++)
      while(ss && ss>>TT)
       {
                  //    ss>>TT;
         T.push_back(TT);
       }
      Fam_infor.push_back(T);
      double TTT;
      ss>>TTT;
      pheno.push_back(TTT);
    }
   Fam.clear();
   Fam.close();
 }
                                                                     

int readfam(string file,vector<double> &pheno, vector< vector<string> > &Fam_infor)
 {
   ifstream Fam(file.c_str());
   string A;
   if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
//   cout<<"Reading PLINK FAM file from ["+file+"]."<<endl;
   while(Fam && getline(Fam, A))
    {
      stringstream ss(A);
      vector<string> T;
      string TT;
      for(int j=0;j<5;j++)
  //    while(ss && ss>>TT)
       {
         ss>>TT;
         T.push_back(TT);
       }
      Fam_infor.push_back(T);
      double TTT;
      ss>>TTT;
      pheno.push_back(TTT);
    }
   Fam.clear();
   Fam.close();
 }

//Read a bim file, and calculate the number of variants
int readbim(string file)
 {
   int ibuf=0;
   string A;
   string cbuf="0";
   double dbuf=0.0;
   string str_buf;
   ifstream Bim(file.c_str());
   if(!Bim) throw("Error: can not open the file ["+file+"] to read.");
   while(Bim && getline(Bim, A))
    {
      ibuf++;
    }
   Bim.close();
   return ibuf;
 }

void remove_file (vector<string> &file)
 {
   for(int i=0;i<file.size();i++)
             {
               cout<<"Removed file is: "<<i<<", "<<file[i]<<endl;
               if(file[i]=="NA")
                {
                  continue;
                }
               string file1=file[i];
               char c[file1.size() + 1];
               strcpy(c, file1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
             }

 }

void remove_file (vector<string> &file, string format)
 {
   if(format =="bfile")
    {
   for(int i=0;i<file.size();i++)
             {
               cout<<"Removed file is: "<<i<<", "<<file[i]<<endl;
                if(file[i]=="NA")
                {
                  continue;
                }
               string bed=file[i]+string(".bed");
               char c1[bed.size() + 1];
               strcpy(c1, bed.c_str());
               if(remove(c1))
                {
                  cout<<"Successfully delete"<<endl;
                }
               string bim=file[i]+string(".bim");
               char c2[bim.size() + 1];
               strcpy(c2, bim.c_str());
               if(remove(c2))
                {
                  cout<<"Successfully delete"<<endl;
                }
              
               string fam=file[i]+string(".fam");
               char c3[fam.size() + 1];
               strcpy(c3, fam.c_str());
               if(remove(c3))
                {
                  cout<<"Successfully delete"<<endl;
                }
             }
     }

 }

void remove_file (string file)
 {
    cout<<"Removed file "<<file<<endl;
              char c[file.size() + 1];
               strcpy(c, file.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }

 }
int extractvariant(string input, string output, map<string,bool> &target1, string target)
 {
 //   cout<<"Come into extractvariant function to extract genotype for tartet variants"<<endl; 
    int i=0, j=0, k=0;
    input=input+string(".bed");
    output=output+string(".bed");
    string strfam="fam";
    string strbim="bim";
    string famfile=input;
    string bimfile=input;
    famfile.replace(famfile.end()-3,famfile.end(),strfam);
    bimfile.replace(bimfile.end()-3,bimfile.end(),strbim);
    string famout=output;
    string bimout=output;
    famout.replace(famout.end()-3,famout.end(),strfam);
    bimout.replace(bimout.end()-3,bimout.end(),strbim);
    int nsnp=readbim(bimfile);
    int nind=readfam(famfile);
    map<string,int> order;
    readbim(bimfile,order);

//    for(map<string,int>::iterator Ite=order.begin();Ite!=order.end();Ite++)
//     {
//       cout<<"Variant is: "<<Ite->first<<", and order is: "<<Ite->second<<endl;
//     }

    map<int, bool> INDEX;
    int index1=-1;
    for(map<string,bool>::iterator Ite=target1.begin();Ite!=target1.end();Ite++)
     {
//       bool index1=false;
//       cout<<"tested variant is: "<<Ite->first<<endl;
       if(order.find(Ite->first)!=order.end())
        {
          if(Ite->first ==target)
           {
  //           cout<<"Peak variant found"<<Ite->first<<endl;
    //         cout<<"Order of this variant is: "<<order[Ite->first]<<endl;
             index1=order[Ite->first];
             INDEX.insert(pair<int,bool>(order[Ite->first],true));
           } else
           {
      //       cout<<"Covariate variant found"<<endl;
             INDEX.insert(pair<int,bool>(order[Ite->first],false));
           }
        }
       
     }
//    cout<<"index1 is: "<<index1<<endl;
    if(index1==-1)
     {
       return (-2);
     }   
     
    mat _snp_2 =mat(INDEX.size(),nind, fill::zeros);
    mat _snp_1 =mat(INDEX.size(),nind, fill::zeros);
    char ch[1];
    bitset<8> b;
//    cout<<"Size of INDEX is: "<<INDEX.size()<<endl;
//    cout<<"input is: "<<input<<endl;
    fstream BIT(input.c_str(), ios::in|ios::binary);
//    cout<<"Start to read genotype"<<endl;
    if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
    for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
    int snp_indx=0, indi_indx=0;
    map<int, bool> INDEX1;
//    cout<<"I am here1"<<endl;
//    cout<<"nsnp is: "<<nsnp<<endl;
    for(j=0, snp_indx=0; j<nsnp; j++)
     {
       if(INDEX.find(j)==INDEX.end())
        {
//          cout<<"Not variant I interest"<<endl;
//         for(int i=0;i<nind;)
         for(int i=0;i<nind;i+=4)
          {
          BIT.read(ch,1);
/*
          b=ch[0];
          k=0;
          while(k < 7 && i < nind)
           {
             k++;
             k++;
             indi_indx++;
             i++;
           }
*/
          }
//          cout<<"Is wrong?"<<endl;
          continue;
        }
      if(j==index1)
       {
         INDEX1.insert(pair<int,bool> (snp_indx,true));
       } else
       {
         INDEX1.insert(pair<int,bool> (snp_indx,false));
       }
       for(i=0, indi_indx=0; i<nind;)
        {
          BIT.read(ch,1);
//          if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
          b=ch[0];
          k=0;
          while(k < 7 && i < nind)
           {
             _snp_2(snp_indx,indi_indx)=(!b[k++]);
             _snp_1(snp_indx,indi_indx)=(!b[k++]);
             indi_indx++;
             i++;
           }
        }
       snp_indx++;
     }
    mat dosage=mat(INDEX.size(),nind, fill::zeros);
    for(int i=0;i<_snp_2.n_rows;i++)
     {
       for(int j=0;j<_snp_2.n_cols;j++)
        {
          if(_snp_2(i,j)==0 && _snp_1(i,j)==1)
           {
             dosage(i,j)=-9;
           } else
           {
             dosage(i,j)=_snp_2(i,j)+_snp_1(i,j);
           }
        }
     }
    for(int i=0;i<dosage.n_rows;i++)
     {
       double Mean=0;
       int num=0;
       for(int j=0;j<dosage.n_cols;j++)
        {
          if(dosage(i,j)!=-9)
           {
             Mean+=dosage(i,j);
             num++;
           }
        }
       Mean/=num;
       for(int j=0;j<dosage.n_cols;j++)
        {
          if(dosage(i,j)==-9)
           {
             dosage(i,j)=Mean;
           }
        }
     }
 //   map<string,int> order;
    map<int,string> order1;
 //   readbim(bimfile,order);
    map<string,int>::iterator Ite;
    Ite=order.find(target);
    if(Ite==order.end())
     {
       return -2;
     }
    map<string,int>::iterator Ite1;
    for(Ite1=order.begin();Ite1!=order.end();Ite1++)
     {
       order1.insert(pair<int, string>(Ite1->second, Ite1->first));
     }
    int T_index=Ite->second;
    vector<int> index;
    vector< vector<int> > _Snp_2;
    vector< vector<int> > _Snp_1;
    vector< vector<int> > _Snp_2_cov;
    vector< vector<int> > _Snp_1_cov;
    vector <int> C;
    bool cov_exists=false;
    for(int x=0;x<_snp_2.n_rows;x++)
     {
    //   mat r2=cor(dosage.row(x),dosage.row(T_index));
//       if(r2(0,0)*r2(0,0)>=r2_cutoff)
       if(INDEX1[x])
        {
      //    C.push_back(x);
          vector<int> Test1;
          vector<int> Test2;
          for(int y=0;y<_snp_2.n_cols;y++)
           {
             Test1.push_back(_snp_1(x,y));
             Test2.push_back(_snp_2(x,y));
           }
          _Snp_2.push_back(Test2);
          _Snp_1.push_back(Test1);
        } else 
        {
 //         if(target1.find(order1[x])==target1.end())
 //          { 
 //            continue;
 //          }
          cov_exists=true;
    //      C.push_back(x);
          vector<int> Test1;
          vector<int> Test2;
          for(int y=0;y<_snp_2.n_cols;y++)
           {
             Test1.push_back(_snp_1(x,y));
             Test2.push_back(_snp_2(x,y));
           }
          _Snp_2_cov.push_back(Test2);
          _Snp_1_cov.push_back(Test1);
        }
     }
//    if(!cov_exists)
//     {
//       return -1;
//     }
   int Ret=-1;
   if(INDEX1.size()>1)
   {
     Ret=0;
    string output_cov=output+string("_cov");
  
    
    fstream  FHOU(output_cov.c_str(),ios::out);
    if(!FHOU)
    {
      return (1);
    }
//    cout<<"dosage is: "<<endl;
//    cout<<dosage<<endl;
//   FHOU<<input[0][0];
//   for(int i=0;i<C.size();i++)
   for(int i=0;i<dosage.n_cols;i++)
    {
//      if(INDEX1[i])
//       {
//         continue;
//       }
      FHOU<<"1";
      for(int j=0;j<dosage.n_rows;j++)
       {
         if(INDEX1[j])
          {
            continue;
          }
         FHOU<<"	"<<dosage(j,i);
       }
      FHOU<<endl;
    }
   FHOU.clear();
   FHOU.close(); 
   }
//   return (0);   
 
    mat _snp_2_target=mat(_Snp_2.size(),_Snp_2[0].size(),fill::zeros);
    mat _snp_1_target=mat(_Snp_1.size(),_Snp_1[0].size(),fill::zeros);
    for(int i=0;i<_Snp_1.size();i++)
     {
       for(int j=0;j<_Snp_1[0].size();j++)
        {
          _snp_1_target(i,j)=_Snp_1[i][j];
          _snp_2_target(i,j)=_Snp_2[i][j];
        }
     }
    outputbed(_snp_2_target,_snp_1_target,output);
    copyfam(famfile,famout);
//    cout<<"index1 is: "<<index1<<endl;
    C.push_back(index1);
    copybim(bimfile,bimout,C);
    BIT.clear();
    BIT.close();
    return Ret;
 }


//Given a r2 cutoff, extract the variants locating in high LD with a specific varaint
void extractvariant(string input, string output, double r2_cutoff, string target)
 {
    int i=0, j=0, k=0;
    string strfam="fam";
    string strbim="bim";
    string famfile=input;
    string bimfile=input;
    famfile.replace(famfile.end()-3,famfile.end(),strfam);
    bimfile.replace(bimfile.end()-3,bimfile.end(),strbim);
    string famout=output;
    string bimout=output;
    famout.replace(famout.end()-3,famout.end(),strfam);
    bimout.replace(bimout.end()-3,bimout.end(),strbim);
    int nsnp=readbim(bimfile);
    int nind=readfam(famfile);
    mat _snp_2 =mat(nsnp,nind, fill::zeros);
    mat _snp_1 =mat(nsnp,nind, fill::zeros);
    char ch[1];
    bitset<8> b;
    fstream BIT(input.c_str(), ios::in|ios::binary);
    if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
    for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
    int snp_indx=0, indi_indx=0;
    for(j=0, snp_indx=0; j<nsnp; j++)
     {
       for(i=0, indi_indx=0; i<nind;)
        {
          BIT.read(ch,1);
          if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
          b=ch[0];
          k=0;
          while(k < 7 && i < nind)
           {
             _snp_2(snp_indx,indi_indx)=(!b[k++]);
             _snp_1(snp_indx,indi_indx)=(!b[k++]);
             indi_indx++;
             i++;
           }
        }
       snp_indx++;
     }
    mat dosage=mat(nsnp,nind, fill::zeros);
    for(int i=0;i<_snp_2.n_rows;i++)
     {
       for(int j=0;j<_snp_2.n_cols;j++)
        {
          if(_snp_2(i,j)==0 && _snp_1(i,j)==1)
           {
             dosage(i,j)=-9;
           } else
           {
             dosage(i,j)=_snp_2(i,j)+_snp_1(i,j);
           }
        }
     }
    for(int i=0;i<dosage.n_rows;i++)
     {
       double Mean=0;
       int num=0;
       for(int j=0;j<dosage.n_cols;j++)
        {
          if(dosage(i,j)!=-9)
           {
             Mean+=dosage(i,j);
             num++;
           }
        }
       Mean/=num;
       for(int j=0;j<dosage.n_cols;j++)
        {
          if(dosage(i,j)==-9)
           {
             dosage(i,j)=Mean;
           }
        }
     }

    map<string,int> order;
    readbim(bimfile,order);
    map<string,int>::iterator Ite; 
    Ite=order.find(target);
    if(Ite==order.end())
     {
       return ;
     }
    int T_index=Ite->second;
    vector<int> index;
    vector< vector<int> > _Snp_2;
    vector< vector<int> > _Snp_1;
    vector <int> C;
    for(int x=0;x<_snp_2.n_rows;x++)
     {
       mat r2=cor(dosage.row(x),dosage.row(T_index));
       if(r2(0,0)*r2(0,0)>=r2_cutoff)
        {
          C.push_back(x);
          vector<int> Test1;
          vector<int> Test2;
          for(int y=0;y<_snp_2.n_cols;y++)
           {
             Test1.push_back(_snp_1(x,y));
             Test2.push_back(_snp_2(x,y));
           }
          _Snp_2.push_back(Test2);
          _Snp_1.push_back(Test1);
        }
     }

    mat _snp_2_target=mat(_Snp_2.size(),_Snp_2[0].size(),fill::zeros);
    mat _snp_1_target=mat(_Snp_1.size(),_Snp_1[0].size(),fill::zeros);
    for(int i=0;i<_Snp_1.size();i++)
     {
       for(int j=0;j<_Snp_1[0].size();j++)
        {
          _snp_1_target(i,j)=_Snp_1[i][j];
          _snp_2_target(i,j)=_Snp_2[i][j];
        }
     }
    outputbed(_snp_2_target,_snp_1_target,output);
    copyfam(famfile,famout);
    copybim(bimfile,bimout,C);
    BIT.clear();
    BIT.close();
 }

//Make a copy of plink bed file
void copybed(string input, string output)
 {
 //   cout<<"Come into copybed function"<<endl;
    int i=0, j=0, k=0;
    string strfam="fam";
    string strbim="bim";
    string famfile=input;
    string bimfile=input;

    famfile.replace(famfile.end()-3,famfile.end(),strfam);
    bimfile.replace(bimfile.end()-3,bimfile.end(),strbim);
    int nsnp=readbim(bimfile);
    int nind=readfam(famfile);
//    cout<<"Number of variants is: "<<nsnp<<endl;
//    cout<<"Number of individuals is: "<<nind<<endl;
    mat _snp_2 =mat(nsnp,nind, fill::zeros);
    mat _snp_1 =mat(nsnp,nind, fill::zeros);
    char ch[1];
    bitset<8> b;
    fstream BIT(input.c_str(), ios::in|ios::binary);
    if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
 //   cout<<"Reading PLINK BED file from ["+input+"] in SNP-major format ..."<<endl;
    for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
//    cout<<"Here1 is OK"<<endl;
    int snp_indx=0, indi_indx=0;
    for(j=0, snp_indx=0; j<nsnp; j++)
     {
    //   cout<<"Here2 is OK"<<endl;
       for(i=0, indi_indx=0; i<nind;)
        {
          BIT.read(ch,1);
          if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
          b=ch[0];
          k=0;
          while(k < 7 && i < nind)
           {
             _snp_2(snp_indx,indi_indx)=(!b[k++]);
             _snp_1(snp_indx,indi_indx)=(!b[k++]);
             indi_indx++;
             i++;
           }
        }
       snp_indx++;
     }
//    cout<<"Output bed file"<<endl;
    outputbed(_snp_2,_snp_1,output);
    BIT.clear();
    BIT.close();
//    cout<<"Function copybed is over"<<endl;
 }

void prepare_geno_phe(string input, string X_file, string YFile)
 {
//   cout<<"Come into function to prepare genotype and phenotype"<<endl;
   int i,j,k;
   string bedfile=input+".bed";
   string famfile=input+".fam";
   string bimfile=input+".bim";
//   cout<<"bedfile is: "<<bedfile<<endl;
//   cout<<"famfile is: "<<famfile<<endl;
//   cout<<"bimfile is: "<<bimfile<<endl;
   int nsnp=readbim(bimfile);
   int nind=readfam(famfile);
   vector<string> variant;
   readbim(bimfile,variant);
   mat _snp_2 =mat(nsnp,nind, fill::zeros);
   mat _snp_1 =mat(nsnp,nind, fill::zeros);
//   cout<<"nsnp is: "<<nsnp<<endl;
//   cout<<"nind is: "<<nind<<endl;
   char ch[1];
   bitset<8> b;
   fstream BIT(bedfile.c_str(), ios::in|ios::binary);
   if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
   for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
   int snp_indx=0, indi_indx=0;
   for(j=0, snp_indx=0; j<nsnp; j++)
    {
      for(i=0, indi_indx=0; i<nind;)
       {
         BIT.read(ch,1);
         if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
         b=ch[0];
         k=0;
         while(k < 7 && i < nind)
          {
            _snp_2(snp_indx,indi_indx)=(!b[k++]);
            _snp_1(snp_indx,indi_indx)=(!b[k++]);
            indi_indx++;
            i++;
          }
       }
      snp_indx++;
    }
//   cout<<"Reading PLINK BED file is over"<<endl;
   mat dosage=mat(nsnp,nind, fill::zeros);
   for(int i=0;i<_snp_2.n_rows;i++)
    {
      for(int j=0;j<_snp_2.n_cols;j++)
       {
         if(_snp_2(i,j)==0 && _snp_1(i,j)==1)
          {
            dosage(i,j)=-9;
          } else
          {
            dosage(i,j)=_snp_2(i,j)+_snp_1(i,j);
          }
       }
    }
   for(int i=0;i<dosage.n_rows;i++)
    {
      double Mean=0;
      int num=0;
      for(int j=0;j<dosage.n_cols;j++)
       {
         if(dosage(i,j)!=-9)
          {
            Mean+=dosage(i,j);
            num++;
          }
       }
      Mean/=num;
      for(int j=0;j<dosage.n_cols;j++)
       {
         if(dosage(i,j)==-9)
          {
            dosage(i,j)=Mean;
          }
       }
    }
   ofstream outfile1(X_file.c_str(), ios::out);
   outfile1<<variant[0];
   for(int x=1;x<variant.size();x++)
    {
      outfile1<<" "<<variant[x];
    }
   outfile1<<endl;
   for(int i=0;i<dosage.n_cols;i++)
    {
      outfile1<<dosage(0,i);
      for(int j=1;j<dosage.n_rows;j++)
       {
         outfile1<<" "<<dosage(j,i);
       }
      outfile1<<endl;
    }
   outfile1.close();

   vector<double>  pheno;
   vector < vector<string> > fam_infor;
   readfam(famfile,pheno,fam_infor,1);
   ofstream outfile2(YFile.c_str(), ios::out);
   for (int i = 0; i < pheno.size(); i++)
     outfile2 << pheno[i] << endl;
   outfile2.close();

 }



int calculate_Mvalue(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > > &beta_,  vector<vector<vector<double> > > &sd_, int Index, string output, vector<vector<string> >  &variant_, map<string, bool> &pos, map<string, bool> &pos1, map<int,bool> &data_index)
 {

//   cout<<"Come into calculate_Mvalue function"<<endl;
   mat cov_z_score = mat(z_score_[0][Index].size(),z_score_.size(), fill::zeros);
   mat cov_beta    = mat(z_score_[0][Index].size(),z_score_.size(), fill::zeros);
   mat cov_z_score_copy = mat(z_score_[0][Index].size(),z_score_.size(), fill::zeros);
   mat cov_sd      = mat(z_score_[0][Index].size(),z_score_.size(), fill::zeros);
   mat cov_sd_copy      = mat(z_score_[0][Index].size(),z_score_.size(), fill::zeros);
//   cout<<"Extract z_score"<<endl;
   for(int i=0;i<z_score_.size();i++)
    {
      for(int j=0;j<z_score_[i][Index].size();j++)
       {
         cov_z_score(j,i)=z_score_[i][Index][j];
         cov_beta(j,i)   =beta_[i][Index][j];
         cov_z_score_copy(j,i)=z_score_[i][Index][j];
         cov_sd(j,i)=sd_[i][Index][j];
         cov_sd_copy(j,i)=sd_[i][Index][j];
       }
    }
   vector<int> indice;
//   cout<<"Remove significant signals"<<endl;
//   cout<<"Size of cov_z_score is: "<<cov_z_score.n_rows<<"  and"<<cov_z_score.n_cols<<endl;
   for(int i=0;i<cov_z_score.n_rows;i++)
    {
      int sym=0;

      for(int j=0;j<cov_z_score.n_cols;j++)
       {
//         cout<<i<<" "<<j<<endl;
         if(abs(cov_z_score(i,j))>2)
          {
             sym=1;
          }
       }
      if(sym==1)
       {
//         cout<<"I have to remove the variant "<<i<<endl;
         indice.push_back(i);
       }
    }
//   cout<<"Start to remove rows "<<indice.size()<<endl;
   uvec Indice = conv_to< uvec >::from(indice);
   cov_z_score_copy.shed_rows(Indice);
   cov_sd_copy.shed_rows(Indice);
//   cout<<"Variants with significant signals are removed"<<endl;
 
//   cout<<"Data used to contruct covariance: "<<endl;
//      cout<<cov_z_score_copy<<endl;
//      cout<<"Output data used to construct covariance is done"<<endl;

   mat Cov_1 =cov(cov_z_score_copy);
   
   int I_biao=cov_z_score.n_rows;
   mat MVALUE = mat(z_score_[0][Index].size(),z_score_.size(), fill::zeros);
   for(int i=0;i<cov_z_score.n_rows;i++)
    {
      if(!pos[variant_[0][i]] || !pos1[variant_[0][i]])
       {
         continue;
       }
      mat Cov = Cov_1;
      int na_index=0;
      vector<int> index;
      map<int, bool> Index2;
      for(int X=0;X<Cov.n_rows;X++)
       {
         for(int Y=0;Y<Cov.n_rows;Y++)
           {
             if(X==Y)
               {
                 Cov(X,Y)=cov_sd(i,X)*cov_sd(i,Y);
               } else
               {
                 Cov(X,Y)=Cov(X,Y)*cov_sd(i,X)*cov_sd(i,Y);
               }
           }
         if(cov_sd(i,X)<=0.0000001)
          {
            na_index=1;
            index.push_back(X);
            Index2.insert(pair<int, bool>(X, true));
          }else
          {
            Index2.insert(pair<int, bool>(X, false));
          }
       }
//      cout<<"Cov used for Mvalue calculation is: "<<endl;
//     cout<<Cov<<endl;
      uvec Indice_1 = conv_to< uvec >::from(index);
      if(index.size()>0 && Cov.n_rows-index.size()>=3)
       {
//         uvec Indice_1 = conv_to< uvec >::from(index);
         Cov.shed_cols(Indice_1);
         Cov.shed_rows(Indice_1);
       } else if(Cov.n_rows-index.size()<3)
       {
         for(int BIAO=0;BIAO<Cov.n_rows;BIAO++)
          {
            MVALUE(i,BIAO)=-9.0;
          } 
         continue;
       }
      mat inv_Cov = inv(Cov);
      mat Cov_new = mat(Cov.n_rows,Cov.n_rows,fill::zeros);
      int neg_index=0;
      for(int A=0;A<Cov.n_rows;A++)
       {
         Cov_new(A,A)=0.0;
         for(int C=0;C<Cov.n_rows;C++)
          {
            Cov_new(A,A)+=inv_Cov(A,C);
 
          }
         if(Cov_new(A,A)<=0.0)
          {
            neg_index=1; 
          }
         Cov_new(A,A)=1.0/Cov_new(A,A);
       }
      if(neg_index==1)
       {
         Cov_new = Cov;
       }
      mat Cov_beta= mat(Cov_1.n_rows,1,fill::zeros);
//      Cov_beta.shed_rows(Indice_1);
  
      for(int J=0;J<Cov_1.n_cols;J++)
       {
         Cov_beta(J,0) =  cov_beta(i,J);
       }
      int start=0;
      int number_tissue=Cov.n_cols;
  //    cout<<"Number of tissues is: "<<number_tissue<<endl;
      int end =exp2(number_tissue);
  //    cout<<"Value of end is: "<<end<<endl;
      int Interval=3;
      double arrayOfDouble5 [Interval][number_tissue];
      double arrayOfDouble4 [Interval][number_tissue];
      for(int B=0;B<Interval;B++)
       {
          for (int b5 = 0; b5 < number_tissue; b5++)
            {
              arrayOfDouble5[B][b5]=0.0;
              arrayOfDouble4[B][b5]=0.0;
            }
       }
      for(int B=0;B<Interval;B++)
       {
      for(int M=0;M<end;M++)
       {
  //       cout<<"M is: "<<M<<endl;
         int index[number_tissue];
         for(int I=0;I<number_tissue;I++)
          {
            index[I]=0;
          }
         int number_positive=0;
         int Ene =end;
         int End =M;
         for(int I=0;I<number_tissue;I++)
          {
            if(End%2)
             {
               index[I]=1;
               number_positive++;
             }
            End = End/2;
          }
         double Beta =1.0;
         double Alpha =1.0;
         double prior=log(gsl_sf_beta(number_positive+Alpha, number_tissue-number_positive+Beta)/gsl_sf_beta(Alpha,Beta));
         double d1 = 0.0;
         for (int b = 0; b < number_tissue; b++)
          {
            if (index[b]==0)
             {
               d1 += -0.5 * log(Cov_new(b,b)) - log(6.2831853071795862)/2.0 - 1.0/Cov_new(b,b) * Cov_beta(b,0) * Cov_beta(b,0) / 2.0;
             }
          }
         double d2 = 0.0;
         if (number_positive > 0)
          {
             double d3 = 0.0;
             double d4 = 0.0;
             double d5 = 0.0;
             double d6 = 0.0;
             for (int b1 = 0; b1 < number_tissue; b1++)
              {
                if (index[b1]==1)
                 {
                   d3 += 1.0/Cov_new(b1,b1);
                   d4 += 1.0/Cov_new(b1,b1) *  Cov_beta(b1,0);
                   d5 += 1.0/Cov_new(b1,b1)* Cov_beta(b1,0)* Cov_beta(b1,0);
                   d6 += log(1.0/Cov_new(b1,b1));
                 }
              }
            double d7 = d4 / d3;
            double d8 = d3;
            double d9 = 1.0 / (1.0 / d8 + 0.03*(B+1.0));
            double d10 = -(number_positive - 1) * log(6.2831853071795862)/2.0 + 0.5 * d6 - 0.5 * log(d3) - (d5 - d4 * d4 / d3) / 2.0;
            double d11 = 0.5 * log(d9) - log(6.2831853071795862)/2.0 - d9 * d7 * d7 / 2.0;
            d2 = d11 + d10;
          }
        d1= d1 + d2;
        for (int b5 = 0; b5 < number_tissue; b5++)
         {
           if (index[b5]==1)
            {
              arrayOfDouble5[B][b5] = arrayOfDouble5[B][b5] + exp(d1);
            } else
            {
              arrayOfDouble4[B][b5] = arrayOfDouble4[B][b5] + exp(d1);
            }
         }

        }
        }
   //   cout<<"Output variant"<<i<<endl;
//      for (int b2 = 0; b2 < number_tissue; b2++)
      int b2=0;
      int b1=0;
      while(b2<cov_z_score.n_cols && b1<cov_z_score.n_cols)
       {
         if(cov_sd(i,b1)>=0.0000001)
          {
            double d1=0.0; 
            double pos =0.0;
            double neg =0.0;
            for(int B=0;B<Interval;B++)
             {
//            cout<<"Variant_"<<i<<", and tissue_"<<number_tissue;
               pos +=arrayOfDouble5[B][b2];
               neg +=arrayOfDouble4[B][b2];
//            d1 += arrayOfDouble5[B][b2] / (arrayOfDouble4[B][b2] + arrayOfDouble5[B][b2]);
      //      cout<<", Mvalue is: "<<d1<<endl;
      //      MVALUE(i,b2)=d1;  
             }
         MVALUE(i,b1)=pos/(neg+pos);
//         cout<<"Variant_"<<i<<", and tissue_"<<b1;
  //       cout<<", Mvalue is: "<<MVALUE(i,b1)<<endl;
         b2++;
         b1++;
          } else
          {
            MVALUE(i,b1)=0;
    //        cout<<"Variant_"<<i<<", and tissue_"<<b1;
      //      cout<<", Mvalue is: "<<MVALUE(i,b1)<<endl;
            b1++;
          }
       }
       
   }
//   cout<<"Resulted Mvalue is: "<<endl;
//   cout<<MVALUE<<endl;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
//      cout<<"Error for input or output"<<endl;
      return (1);
    }
   FHOU<<"Variant";
   for(int j1=0,indx=0;j1<data_index.size();)
//   for(int j1=0;j1<MVALUE.n_cols;j1++)
    {
      if(!data_index[j1])
       {
         j1++;
         continue;
       }
      
      FHOU<<"\t"<<"Tissue"<<j1;
      j1++;
      indx++;
    }
   FHOU<<endl;
   for(int j1=0;j1<MVALUE.n_rows;j1++)
    {
      if(MVALUE(j1,0)==-9.0)
       {
         continue;
       }
//      FHOU<<"Variant"<<j1;
      if(!pos[variant_[0][j1]] || !pos1[variant_[0][j1]])
       {
         continue;
       }
      FHOU<<variant_[0][j1];
      for(int j2=0;j2<Cov_1.n_cols;)
       {
         FHOU<<"\t"<<MVALUE(j1,j2);
         j2++;
       }

      FHOU<<endl;
    }
//   FHOU<<MVALUE<<endl;
   FHOU.clear();
   FHOU.close();
   return 0;
 }












void read_cor_file(mat &Cov_, string cor_file, int n)
 {
   int ibuf=0;
   int ind1=-1;
   int ind2=-1;
//   cout<<"cor_file is: "<<cor_file<<endl;
//   cout<<"n is: "<<n<<endl;
   string a;
   string str_buf;
   ifstream Bim(cor_file.c_str());
   if(!Bim) throw("Error: can not open the file ["+cor_file+"] to read.");
//   cout<<"I am here in read_cor_file function"<<endl;
   while(Bim &&getline(Bim,a))
             {
               ind1++;
               ind2=-1;
               if(ind1<n)
                {
                   double c=0;
                   stringstream ss(a);

//                      vector<double> T;
                      while(ss && ss>>c)
                       {
                         cout<<" "<<c;
                         ind2++;
                         if(ind2<n)
                          {
                            Cov_(ind1,ind2)=c;
                          }
                       }

                }
               cout<<endl;
              }
   Bim.clear();
   Bim.close();
//   return ibuf;
 }



















double ConstructCovRandom(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > > &beta_,  vector<vector<vector<double> > > &sd_, int index, mat &Res_fixed, mat &Res_random, string meta_mode, string output, vector<vector<string> >  &variant_, map<string, bool> &pos, map<string, bool> &pos1, map <string, vector<double>> &Res, map<string, double > &MAF, vector<double> &VARIANCE, map<int,bool> &data_index, vector<double> &z_score, vector<double> &beta, vector<double> &sd, map<string, string> &ALLELE,string cor_file, bool write_cor, bool Han)
 {

   mat cov_z_score = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_beta    = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_z_score_copy = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_sd      = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_sd_copy      = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);


   for(int i=0;i<z_score_.size();i++)
    {  
      for(int j=0;j<z_score_[i][index].size();j++)
       {
         cov_z_score(j,i)=z_score_[i][index][j];
         cov_beta(j,i)   =beta_[i][index][j];
         cov_z_score_copy(j,i)=z_score_[i][index][j];
         cov_sd(j,i)=sd_[i][index][j];
         cov_sd_copy(j,i)=sd_[i][index][j];
       }
    }

//   cout<<"cov_z_score_copy is: "<<endl;
//   cout<<cov_z_score_copy<<endl;


   vector<int> indice;
   mat Cov_ =mat(cov_z_score_copy.n_cols,cov_z_score_copy.n_cols, fill::ones);
   if(cor_file=="")
    { 
   for(int i=0;i<cov_z_score.n_rows;i++)
    {
        
      int sym=0;
      for(int j=0;j<cov_z_score.n_cols;j++)
       {
  //       cout<<"Z_score is: "<<cov_z_score(i,j)<<endl;
  //       Changed by Biao Zeng, 11/08/2020, 11:38 PM.
//         if(abs(cov_z_score(i,j))>2 || abs(cov_z_score(i,j))<=0.00000001)
         if(abs(cov_z_score(i,j))>2 || abs(cov_z_score(i,j))<=0.00000001)
          {
             sym=1;
          }
       }
      if(sym==1)
       {
    //     cout<<"Will be removed"<<endl;
         indice.push_back(i);
       }
    }
//   cout<<"indice"<<endl;
   uvec Indice = conv_to< uvec >::from(indice);
   cov_z_score_copy.shed_rows(Indice);
   cov_sd_copy.shed_rows(Indice);
//   mat Cov_ =mat(cov_z_score_copy.n_cols,cov_z_score_copy.n_cols, fill::ones);
   }
   if(cor_file!="")
    {
//      cout<<"Go to extract summary correlation file"<<endl;
      read_cor_file(Cov_, cor_file, cov_z_score_copy.n_cols);
//      cout<<"Reading correlation file is done"<<endl;
    } else
    {
//      cout<<"Data used to contruct covariance: "<<endl;
//      cout<<cov_z_score_copy<<endl;
//      cout<<"Output data used to construct covariance is done"<<endl;
      Cov_ =cov(cov_z_score_copy);
    }

   if(write_cor)
    {
      string cor_output=output+string("_summary_correlation");
      fstream  FHOU(cor_output.c_str(),ios::out);
      for(int i=0;i<Cov_.n_rows;i++)
       {
         FHOU<<Cov_(i,0);
         for(int j=1;j<Cov_.n_cols;j++)
          {
            FHOU<<"\t"<<Cov_(i,j);
          }
         FHOU<<endl;
       }
      FHOU.clear();
      FHOU.close();
    }
   
//   cout<<"Cov_ is: "<<endl;
//   cout<<Cov_<<endl;

   int I_biao=cov_z_score.n_rows;
   int peak_index_fixed=-9999;
   double peak_beta_fixed=-9999;
   double peak_z_fixed   =0.0;
   int peak_index_random=-9999;
   double peak_beta_random=-9999;
   double peak_z_random   =0.0;
   mat Output_stat =mat(cov_z_score.n_rows, 3*(cov_z_score.n_cols+1)+1,fill::zeros);
 if(Cov_.n_rows>=3)
  {
   for(int i=0;i<cov_z_score.n_rows;i++)
    {
      mat Cov = Cov_;
      mat Cov1 = Cov_;
      int peak_index=-9999;
      vector<int> index;
      map<int,bool> INDEX;
      for(int X=0;X<Cov.n_rows;X++)
       {
         int na_index=0;
         for(int Y=0;Y<Cov.n_rows;Y++)
          {
             if(X==Y)
              {
                Cov(X,Y)=cov_sd(i,X)*cov_sd(i,Y);
                
              } else
              {
                //Changed by Biao Zeng, 11/09/2020, 12:12 AM
//                Cov(X,Y)=Cov(X,Y)*cov_sd(i,X)*cov_sd(i,Y);
                 if(Cov(X,Y)>1)
                  { 
                    Cov(X,Y)=0.9999*cov_sd(i,X)*cov_sd(i,Y);
                  } else
                  {
                    Cov(X,Y)=Cov(X,Y)*cov_sd(i,X)*cov_sd(i,Y);
                  }
              }
          }
         if(cov_sd(i,X)<=0.00000001)
          {
            na_index=1;
          }
         Output_stat(i,3*X+0) = cov_beta(i,X);
         Output_stat(i,3*X+1) = cov_sd(i,X);
         if(cov_sd(i,X)==0)
          {
            Output_stat(i,3*X+2) = 0;
          } else
          {
            Output_stat(i,3*X+2) = cov_beta(i,X)/cov_sd(i,X);
          }
         if(na_index==1)
          {
            index.push_back(X);
            INDEX.insert(pair<int,bool>(X,true));
          }
       }
      double meta_logP;
      uvec Indice_1 = conv_to< uvec >::from(index); 
      if(index.size()>0 && Cov.n_rows-index.size()>=3)
       {
         Cov.shed_cols(Indice_1); 
         Cov.shed_rows(Indice_1);     
         Cov1.shed_cols(Indice_1);
         Cov1.shed_rows(Indice_1);
       } else if(Cov.n_rows-index.size()<3 )
       {
//         cout<<"Explored variant is: "<<variant_[0][i];
//         cout<<", Number of tissues available is: "<<Cov.n_rows<<", the tissues to remove is: "<<index.size();
//         cout<<", and No enough tissue left"<<endl;
         Output_stat(i,3*Cov.n_rows+0) = 0.0;
         Output_stat(i,3*Cov.n_rows+1) = 0.0;
         Output_stat(i,3*Cov.n_rows+2) = 0.0;
         Output_stat(i,3*(Cov.n_rows+1)+0) = 0.0;
         vector<double> COL;
         COL.push_back(-9);
         COL.push_back(-9);
         z_score.push_back(0.0);
         beta.push_back(0.0);
         sd.push_back(0.0);
         continue;
       }
      vec eigval; //doubleMatrix1D1
      mat eigvec;  //doubleMatrix2D3
      eig_sym(eigval, eigvec, Cov) ;
       bool indicator=false;
      for(int indx=0;indx<eigval.n_elem;indx++)
       {
//         cout<<"Variant "<<i<<", Eigen value for PC"<<indx<<" is: "<<eigval(indx)<<endl;
         if(eigval(indx)<0)
          {
            indicator=true;
          }
       }
      if(indicator)
       {
         Output_stat(i,3*Cov.n_rows+0) = 0.0;
         Output_stat(i,3*Cov.n_rows+1) = 0.0;
         Output_stat(i,3*Cov.n_rows+2) = 0.0;
         Output_stat(i,3*(Cov.n_rows+1)+0) = 0.0;
         vector<double> COL;
         COL.push_back(-9);
         COL.push_back(-9);
         z_score.push_back(0.0);
         beta.push_back(0.0);
         sd.push_back(0.0);
         continue;
       }
      mat inv_Cov = inv(Cov);
      mat Cov_beta= mat(Cov.n_rows,1,fill::zeros);
      mat Cov_r= mat(Cov.n_rows,1,fill::zeros);
      for(int J=0,indx=0;J<Cov.n_cols+index.size();)
       {
         map<int,bool>::iterator ITE;
         ITE=INDEX.find(J);
         if(ITE!=INDEX.end())
          {
            J++;
            continue;
          }
         Cov_beta(indx,0) =  cov_beta(i,J);
         Cov_r(indx,0) =  cov_beta(i,J)/sqrt(VARIANCE[J]);
         J++;indx++;
       }
      mat One =mat(1,Cov_beta.n_rows,fill::ones);
      mat meta_beta=One*inv_Cov*Cov_beta*1.0/(One*inv_Cov*trans(One));
      mat meta_r   =One*inv_Cov*Cov_r*1.0/(One*inv_Cov*trans(One));
      mat Biao=One*inv_Cov*trans(One);
      double  meta_sd=sqrt(1/Biao(0,0));
      mat meta_z = meta_beta / meta_sd;
//      cout<<"Explored variant is: "<<i<<" "<<variant_[0][i]<<" meta_beta is: "<<meta_beta(0,0)<<", meta_sd is: "<<meta_sd<<", and meta_z is: "<<meta_z<<endl;
      vector<double> COL;
      
   //   double geno_variance =2*MAF[variant_[0][i]]*(1-MAF[variant_[0][i]]); //Changed by Biao Zeng 10/08/2020, 07:06 am
      COL.push_back(meta_r(0,0));
      COL.push_back(meta_z(0,0));
      Res.insert(pair<string, vector<double>>(variant_[0][i],COL));
      Output_stat(i,3*Cov_.n_rows+0) = meta_beta(0,0);
      Output_stat(i,3*Cov_.n_rows+1) = meta_sd;
      Output_stat(i,3*Cov_.n_rows+2) = meta_z(0,0);
      z_score.push_back(meta_z(0,0));
      beta.push_back(meta_beta(0,0));
      sd.push_back(meta_sd);
          if(meta_beta(0,0)<=0)
            {
               meta_logP = log10(2.0*normcdf(meta_z(0,0)));
            } else
            {
              meta_logP = log10(2.0*(1-normcdf(meta_z(0,0))));
            }
     if(meta_z(0,0)*meta_z(0,0)>peak_z_fixed*peak_z_fixed && pos[variant_[0][i]]) 
      {
        peak_index_fixed=i;
        peak_beta_fixed= meta_beta(0,0);
        peak_z_fixed   = meta_z(0,0);
      }
      int I_biao=Cov_beta.n_rows;
      cx_mat  square_root_invCov=sqrtmat(inv_Cov);
      mat Cov_beta_adjust=real(square_root_invCov)*Cov_beta;
      mat Cov_beta_adjust_sqr= trans(Cov_beta_adjust) * Cov_beta_adjust;
      mat doubleMatrix1D1 = eigval;
      int I = I_biao;
      mat doubleMatrix2D4 = mat(I,I,fill::eye);
      mat doubleMatrix2D5 = mat(I,I,fill::ones)/I;
      doubleMatrix2D4    = doubleMatrix2D4 - doubleMatrix2D5;
      mat doubleMatrix2D6 = mat(I,I,fill::eye);
      doubleMatrix2D6 = doubleMatrix2D6 + Cov;
      mat  doubleMatrix2D7 = doubleMatrix2D4 * trans(doubleMatrix2D6) * trans(doubleMatrix2D4);
      vec eigval2;
      mat eigvec2;
      eig_sym( eigval2, eigvec2, doubleMatrix2D7);
      vec eigval2_1 = eigval2-1;
      mat eigvec2_1 = eigvec2;
      vector<int> test_index;
      test_index.push_back(0);
      uvec Indices = conv_to< uvec >::from(test_index);
      eigvec2_1.shed_cols(Indices);
      mat doubleMatrix2D8 = eigvec2_1;
      int Test=I-1;
      mat doubleMatrix1D2 = mat(I-1,1,fill::zeros);
      for(int X=1;X<I;X++)
       {
         doubleMatrix1D2(X-1,0)=eigval2_1(X);
       }
      mat doubleMatrix1D3 = mat(I,1,fill::zeros);
      mat  doubleMatrix1D4 = trans(doubleMatrix2D8) * Cov_beta;
      mat  doubleMatrix1D5 = doubleMatrix1D4 % doubleMatrix1D4;
      mat doubleMatrix1D6 = doubleMatrix1D1;
      mat doubleMatrix1D7 = doubleMatrix1D2;
      double start=0.0;
      double end = 10000;
      vector<double> DoubleMatrix1D5;
      vector<double> DoubleMatrix1D6;
      vector<double> DoubleMatrix1D7;
      for(int j=0;j<doubleMatrix1D5.n_rows;j++)
       {
         DoubleMatrix1D5.push_back(doubleMatrix1D5(j,0));
       }
      for(int j=0;j<doubleMatrix1D6.n_rows;j++)
       {
         DoubleMatrix1D6.push_back(doubleMatrix1D6(j,0));
       }
      for(int j=0;j<doubleMatrix1D7.n_rows;j++)
       {
         DoubleMatrix1D7.push_back(doubleMatrix1D7(j,0));
       }
      double result1=-9999.0;
      double result2=-9999.0;
      brent_method(DoubleMatrix1D6, DoubleMatrix1D5, DoubleMatrix1D7,start, end, result1, result2);
      double d2 =-result2;
      double d4 = Cov_beta_adjust_sqr(0,0);
      double d3 = -0.5 * (I * log(6.2831853071795862) + sum(log(eigval)) + d4);
      double d5 = -2.0 * (d3 - d2);
      double statisticRandomEffects2_ = d5;
      //Changed ny Biao Zeng, 10/28/2020, 2:35 pm
//      if(statisticRandomEffects2_<0) 
       if(statisticRandomEffects2_<0.01)
       {
         statisticRandomEffects2_=0;
       }
      double pvalueRandomEffects2Asymptotic_ = 0.5 * gsl_cdf_chisq_Q (statisticRandomEffects2_, 1.0) + 0.5 * gsl_cdf_chisq_Q (statisticRandomEffects2_, 2.0);
      if(pvalueRandomEffects2Asymptotic_ >=0.000001)
       {
         pos1[variant_[0][i]]=false;
       }

   
      //Changed by Biao Zeng, 10/28/2020, 4:01 pm
      if(statisticRandomEffects2_<=9)
       { 
         statisticRandomEffects2_= meta_z(0,0) * meta_z(0,0);
       } else if(Han)
       {
         cout<<"Apply Han method to correct"<<endl;
         double p_Han_adjust=Han_adjust( Cov1.n_cols, Cov1,  sqrt(statisticRandomEffects2_));
         statisticRandomEffects2_= abs(gsl_cdf_ugaussian_Pinv (p_Han_adjust/2));
         statisticRandomEffects2_=statisticRandomEffects2_*statisticRandomEffects2_;
       }

      Output_stat(i,3*(Cov_.n_rows+1)+0) = sqrt(statisticRandomEffects2_);
      if(statisticRandomEffects2_*statisticRandomEffects2_>peak_z_random*peak_z_random  && pos[variant_[0][i]])
      {
        peak_index_random=i;
        peak_z_random   = statisticRandomEffects2_;
      }
   }
  if(peak_index_random<0)
   {
     return -1;
   }
  for(int j=0;j<cov_z_score.n_cols;j++)
   {
     Res_random(j,0)=peak_index_random;
     Res_random(j,1)=cov_beta(peak_index_random,j);
     Res_random(j,2)=peak_z_random;     
   }
  for(int j=0;j<cov_z_score.n_cols;j++)
   {
     Res_fixed(j,0)=peak_index_fixed;
     Res_fixed(j,1)=peak_beta_fixed;
     Res_fixed(j,2)=peak_z_fixed;
   }
 fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }
   FHOU<<"Variant"<<"\t"<<"Allele";
   for(int j=0,indx=0;j<data_index.size();)
    {
      if(!data_index[j])
       {
         j++;
         continue;
         
       }
      
      FHOU<<"\t"<<"beta_tissue_"<<j<<"\t"<<"sd_tissue_"<<j<<"\t"<<"z_tissue_"<<j;
      j++;
      indx++;
    }
   FHOU<<"\t"<<"fixed_beta"<<"\t"<<"fixed_sd"<<"\t"<<"fixed_z";
   FHOU<<"\t"<<"Random_Z"<<endl;
   for(int j1=0;j1<Output_stat.n_rows;j1++)
    {
      if(!pos[variant_[0][j1]])
       {
         continue;
       }
      string test = variant_[0][j1];
      FHOU<<test<<"\t"<<ALLELE[test];
      for(int j2=0;j2<Output_stat.n_cols;j2++)
       {
         FHOU<<"\t"<<Output_stat(j1,j2); 
       }
      FHOU<<endl;
    }
   FHOU.clear();
   FHOU.close();
  } else if(Cov_.n_rows==2)
  {
//    cout<<"There are only 2 data set available for analysis"<<endl;
    for(int i=0;i<cov_z_score.n_rows;i++)
    {
      mat Cov = Cov_;
      int peak_index=-9999;
      vector<int> index;
      map<int,bool> INDEX;
      for(int X=0;X<Cov.n_rows;X++)
       {
         int na_index=0;
         for(int Y=0;Y<Cov.n_rows;Y++)
          {
             if(X==Y)
              {
                Cov(X,Y)=cov_sd(i,X)*cov_sd(i,Y);

              } else
              {
                Cov(X,Y)=Cov(X,Y)*cov_sd(i,X)*cov_sd(i,Y);
              }
          }
         if(cov_sd(i,X)<=0.00000001)
          {
            na_index=1;
          }
         Output_stat(i,3*X+0) = cov_beta(i,X);
         Output_stat(i,3*X+1) = cov_sd(i,X);
         if(cov_sd(i,X)==0)
          {
            Output_stat(i,3*X+2) = 0;
          } else
          {
            Output_stat(i,3*X+2) = cov_beta(i,X)/cov_sd(i,X);
          }
         if(na_index==1)
          {
            index.push_back(X);
            INDEX.insert(pair<int,bool>(X,true));
          }
       }
      double meta_logP;
      uvec Indice_1 = conv_to< uvec >::from(index);
      if(index.size()>0 && Cov.n_rows-index.size()>=2)
       {
         Cov.shed_cols(Indice_1);
         Cov.shed_rows(Indice_1);
       } else if(Cov.n_rows-index.size()<2 )
       {
        Output_stat(i,3*Cov.n_rows+0) = 0.0;
         Output_stat(i,3*Cov.n_rows+1) = 0.0;
         Output_stat(i,3*Cov.n_rows+2) = 0.0;
         Output_stat(i,3*(Cov.n_rows+1)+0) = 0.0;
         vector<double> COL;
         COL.push_back(-9);
         COL.push_back(-9);
         z_score.push_back(0.0);
         beta.push_back(0.0);
         sd.push_back(0.0);
         continue;
       }
      vec eigval; //doubleMatrix1D1
      mat eigvec;  //doubleMatrix2D3
      eig_sym(eigval, eigvec, Cov) ;
       bool indicator=false;
      for(int indx=0;indx<eigval.n_elem;indx++)
       {
         if(eigval(indx)<0)
          {
            indicator=true;
          }
       }
      if(indicator)
       {
         Output_stat(i,3*Cov.n_rows+0) = 0.0;
         Output_stat(i,3*Cov.n_rows+1) = 0.0;
         Output_stat(i,3*Cov.n_rows+2) = 0.0;
         Output_stat(i,3*(Cov.n_rows+1)+0) = 0.0;
         vector<double> COL;
         COL.push_back(-9);
         COL.push_back(-9);
         z_score.push_back(0.0);
         beta.push_back(0.0);
         sd.push_back(0.0);
         continue;
       }
      mat inv_Cov = inv(Cov);
      mat Cov_beta= mat(Cov.n_rows,1,fill::zeros);
      mat Cov_r= mat(Cov.n_rows,1,fill::zeros);
      for(int J=0,indx=0;J<Cov.n_cols+index.size();)
       {
         map<int,bool>::iterator ITE;
         ITE=INDEX.find(J);
         if(ITE!=INDEX.end())
          {
            J++;
            continue;
          }
         Cov_beta(indx,0) =  cov_beta(i,J);
         Cov_r(indx,0) =  cov_beta(i,J)/sqrt(VARIANCE[J]);
         J++;indx++;
       }
      mat One =mat(1,Cov_beta.n_rows,fill::ones);
      mat meta_beta=One*inv_Cov*Cov_beta*1.0/(One*inv_Cov*trans(One));
      mat meta_r   =One*inv_Cov*Cov_r*1.0/(One*inv_Cov*trans(One));
      mat Biao=One*inv_Cov*trans(One);
      double  meta_sd=sqrt(1/Biao(0,0));
      mat meta_z = meta_beta / meta_sd;
      vector<double> COL;
      
   //   double geno_variance =2*MAF[variant_[0][i]]*(1-MAF[variant_[0][i]]); // Changed by Biao Zeng, 10/08/2020, 07:07 am
      COL.push_back(meta_r(0,0));
      COL.push_back(meta_z(0,0));
      Res.insert(pair<string, vector<double>>(variant_[0][i],COL));
      Output_stat(i,3*Cov_.n_rows+0) = meta_beta(0,0);
      Output_stat(i,3*Cov_.n_rows+1) = meta_sd;
      Output_stat(i,3*Cov_.n_rows+2) = meta_z(0,0);
      z_score.push_back(meta_z(0,0));
      beta.push_back(meta_beta(0,0));
sd.push_back(meta_sd);
          if(meta_beta(0,0)<=0)
            {
               meta_logP = log10(2.0*normcdf(meta_z(0,0)));
            } else
            {
              meta_logP = log10(2.0*(1-normcdf(meta_z(0,0))));
            }
     if(meta_z(0,0)*meta_z(0,0)>peak_z_fixed*peak_z_fixed && pos[variant_[0][i]])
      {
        peak_index_fixed=i;
        peak_beta_fixed= meta_beta(0,0);
        peak_z_fixed   = meta_z(0,0);
        peak_index_random=i;
        peak_beta_random= meta_beta(0,0);
        peak_z_random   = meta_z(0,0);
      }
  
      double pvalueRandomEffects2Asymptotic_ = 2.0*(1-normcdf(meta_z(0,0)));
      if(pvalueRandomEffects2Asymptotic_ >=0.000001)
       {
         pos1[variant_[0][i]]=false;
       }

      Output_stat(i,3*(Cov_.n_rows+1)+0) = meta_z(0,0);
//      if(statisticRandomEffects2_*statisticRandomEffects2_>peak_z_random*peak_z_random  && pos[variant_[0][i]])
//      {
//        peak_index_random=i;
//        peak_z_random   = statisticRandomEffects2_;
//      }
   }
  if(peak_index_fixed<0)
   {
     return -1;
   }
  for(int j=0;j<cov_z_score.n_cols;j++)
   {
     Res_random(j,0)=peak_index_fixed;
     Res_random(j,1)=peak_beta_fixed;;
     Res_random(j,2)=peak_z_random;
   }
  for(int j=0;j<cov_z_score.n_cols;j++)
   {
     Res_fixed(j,0)=peak_index_fixed;
     Res_fixed(j,1)=peak_beta_fixed;
     Res_fixed(j,2)=peak_z_random;
   }
 fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }
   FHOU<<"Variant"<<"\t"<<"Allele";
   for(int j=0,indx=0;j<data_index.size();)
    {
      if(!data_index[j])
       {
         j++;
         continue;

       }

      FHOU<<"\t"<<"beta_tissue_"<<j<<"\t"<<"sd_tissue_"<<j<<"\t"<<"z_tissue_"<<j;
      j++;
      indx++;
    }
   FHOU<<"\t"<<"fixed_beta"<<"\t"<<"fixed_sd"<<"\t"<<"fixed_z";
   FHOU<<"\t"<<"Random_Z"<<endl;
   for(int j1=0;j1<Output_stat.n_rows;j1++)
    {
      if(!pos[variant_[0][j1]])
       {
         continue;
       }
      FHOU<<variant_[0][j1]<<"\t"<<ALLELE[variant_[0][j1]];
      for(int j2=0;j2<Output_stat.n_cols;j2++)
       {
         FHOU<<"\t"<<Output_stat(j1,j2);
       }
      FHOU<<endl;
    }
   FHOU.clear();
   FHOU.close(); 
  } else if(Cov_.n_rows==1)
  {
//   cout<<"There is only one single data"<<endl;
   for(int i=0;i<cov_z_score.n_rows;i++)
    {
      mat Cov = Cov_;
      int peak_index=-9999;
      vector<int> index;
      map<int,bool> INDEX;
      double meta_logP;



      for(int X=0;X<Cov.n_rows;X++)
       {
         int na_index=0;
         for(int Y=0;Y<Cov.n_rows;Y++)
          {
             if(X==Y)
              {
                Cov(X,Y)=cov_sd(i,X)*cov_sd(i,Y);

              } else
              {
                Cov(X,Y)=Cov(X,Y)*cov_sd(i,X)*cov_sd(i,Y);
              }
          }
         if(cov_sd(i,X)<=0.00000001)
          {
            na_index=1;
          }
         Output_stat(i,3*X+0) = cov_beta(i,X);
         Output_stat(i,3*X+1) = cov_sd(i,X);
         if(cov_sd(i,X)==0)
          {
            Output_stat(i,3*X+2) = 0;
          } else
          {
            Output_stat(i,3*X+2) = cov_beta(i,X)/cov_sd(i,X);
          }
         if(na_index==1)
          {
            index.push_back(X);
            INDEX.insert(pair<int,bool>(X,true));
          }
       }


/*

      uvec Indice_1 = conv_to< uvec >::from(index);
      if(index.size()>0 && Cov.n_rows-index.size()>=2)
       {
         Cov.shed_cols(Indice_1);
         Cov.shed_rows(Indice_1);
       } else if(Cov.n_rows-index.size()<2 )
       {
        Output_stat(i,3*Cov.n_rows+0) = 0.0;
         Output_stat(i,3*Cov.n_rows+1) = 0.0;
         Output_stat(i,3*Cov.n_rows+2) = 0.0;
         Output_stat(i,3*(Cov.n_rows+1)+0) = 0.0;
         vector<double> COL;
         COL.push_back(-9);
         COL.push_back(-9);
         z_score.push_back(0.0);
         beta.push_back(0.0);
         sd.push_back(0.0);
         continue;
       }
      vec eigval; //doubleMatrix1D1
      mat eigvec;  //doubleMatrix2D3
      eig_sym(eigval, eigvec, Cov) ;
       bool indicator=false;
      for(int indx=0;indx<eigval.n_elem;indx++)
       {
         if(eigval(indx)<0)
          {
            indicator=true;
          }
       }
      if(indicator)
       {
         Output_stat(i,3*Cov.n_rows+0) = 0.0;
         Output_stat(i,3*Cov.n_rows+1) = 0.0;
         Output_stat(i,3*Cov.n_rows+2) = 0.0;
         Output_stat(i,3*(Cov.n_rows+1)+0) = 0.0;
         vector<double> COL;
         COL.push_back(-9);
         COL.push_back(-9);
         z_score.push_back(0.0);
         beta.push_back(0.0);
         sd.push_back(0.0);
         continue;
       }

      mat inv_Cov = inv(Cov);
      mat Cov_beta= mat(Cov.n_rows,1,fill::zeros);
      mat Cov_r= mat(Cov.n_rows,1,fill::zeros);
      for(int J=0,indx=0;J<Cov.n_cols+index.size();)
       {
         map<int,bool>::iterator ITE;
         ITE=INDEX.find(J);
         if(ITE!=INDEX.end())
          {
            J++;
            continue;
          }
         Cov_beta(indx,0) =  cov_beta(i,J);
         Cov_r(indx,0) =  cov_beta(i,J)/sqrt(VARIANCE[J]);
         J++;indx++;
       }
      mat One =mat(1,Cov_beta.n_rows,fill::ones);
      mat meta_beta=One*inv_Cov*Cov_beta*1.0/(One*inv_Cov*trans(One));

      mat meta_r   =One*inv_Cov*Cov_r*1.0/(One*inv_Cov*trans(One));
      mat Biao=One*inv_Cov*trans(One);
      double  meta_sd=sqrt(1/Biao(0,0));
      double meta_z = meta_beta / meta_sd;
*/
      vector<double> COL;
//      double geno_variance =2*MAF[variant_[0][i]]*(1-MAF[variant_[0][i]]);
      COL.push_back(1);
//      COL.push_back(meta_z);
//      Res.insert(pair<string, vector<double>>(variant_[0][i],COL));

      double meta_beta=cov_beta(i,0);
      double meta_sd = cov_sd(i,0);
      double meta_z  =0.0;
      if(meta_sd>0)
       {
         meta_z = meta_beta/meta_sd;
       }
      Output_stat(i,3*Cov_.n_rows+0) = meta_beta;
      Output_stat(i,3*Cov_.n_rows+1) = meta_sd;
      Output_stat(i,3*Cov_.n_rows+2) = meta_z;
//      vector<double> COL;
      

      COL.push_back(meta_z);
      Res.insert(pair<string, vector<double>>(variant_[0][i],COL));
      
      z_score.push_back(meta_z);
      beta.push_back(meta_beta);
      sd.push_back(meta_sd);
      if(meta_beta<=0)
            {
               meta_logP = log10(2.0*normcdf(meta_z));
            } else
            {
              meta_logP = log10(2.0*(1-normcdf(meta_z)));
            }
     if(meta_z*meta_z>peak_z_fixed*peak_z_fixed && pos[variant_[0][i]])
      {
        peak_index_fixed=i;
        peak_beta_fixed= meta_beta;
        peak_z_fixed   = meta_z;
        peak_index_random=i;
        peak_beta_random= meta_beta;
        peak_z_random   = meta_z;
      }

      double pvalueRandomEffects2Asymptotic_ = 2.0*(1-normcdf(meta_z));
      if(pvalueRandomEffects2Asymptotic_ >=0.000001)
       {
         pos1[variant_[0][i]]=false;
       }

      Output_stat(i,3*(Cov_.n_rows+1)+0) = meta_z;
//      if(statisticRandomEffects2_*statisticRandomEffects2_>peak_z_random*peak_z_random  && pos[variant_[0][i]])
//      {
//        peak_index_random=peak_index_fixed;
//        peak_z_random   = meta_z;
//      }
   }
  if(peak_index_fixed<0)
   {
     return -1;
   }
  for(int j=0;j<cov_z_score.n_cols;j++)
   {
     Res_random(j,0)=peak_index_fixed;
     Res_random(j,1)=peak_beta_fixed;
     Res_random(j,2)=peak_z_fixed;
   }
  for(int j=0;j<cov_z_score.n_cols;j++)
   {
     Res_fixed(j,0)=peak_index_fixed;
     Res_fixed(j,1)=peak_beta_fixed;
     Res_fixed(j,2)=peak_z_fixed;
   }
 fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }
   FHOU<<"Variant"<<"\t"<<"Allele";
   for(int j=0,indx=0;j<data_index.size();)
    {
      if(!data_index[j])
       {
         j++;
         continue;

       }

      FHOU<<"\t"<<"beta_tissue_"<<j<<"\t"<<"sd_tissue_"<<j<<"\t"<<"z_tissue_"<<j;
      j++;
      indx++;
    }
   FHOU<<"\t"<<"fixed_beta"<<"\t"<<"fixed_sd"<<"\t"<<"fixed_z";
   FHOU<<"\t"<<"Random_Z"<<endl;
   for(int j1=0;j1<Output_stat.n_rows;j1++)
    {
      if(!pos[variant_[0][j1]])
       {
         continue;
       }
      FHOU<<variant_[0][j1]<<"\t"<<ALLELE[variant_[0][j1]];
      for(int j2=0;j2<Output_stat.n_cols;j2++)
       {
         FHOU<<"\t"<<Output_stat(j1,j2);
       }
      FHOU<<endl;
    }
   FHOU.clear();
   FHOU.close();
  
  }
if(Cov_.n_rows>=3)
 {
 if(meta_mode =="random")
  {
    return 0.5 * gsl_cdf_chisq_Q (peak_z_random, 1.0) + 0.5 * gsl_cdf_chisq_Q (peak_z_random, 2.0);
  }else
  {
    return 2.0*(1-normcdf(abs(peak_z_fixed)));
  }
 } else if(Cov_.n_rows<=2)
 {
//    return 2.0*(1-normcdf(abs(peak_z_fixed)));
    return 2.0*(1-normcdf(abs(peak_z_fixed)));
 }
}


void ConstructCov(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > > &beta_,  vector<vector<vector<double> > > &sd_, int index)
 {

//   cout<<"Come into ConstructCov function"<<endl;


   mat cov_z_score = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_beta    = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_z_score_copy = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_sd      = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_sd_copy      = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
//   cout<<"Extract z_score"<<endl;
   for(int i=0;i<z_score_.size();i++)
    {  
      for(int j=0;j<z_score_[i][index].size();j++)
       {
         cov_z_score(j,i)=z_score_[i][index][j];
         cov_beta(j,i)   =beta_[i][index][j];
         cov_z_score_copy(j,i)=z_score_[i][index][j];
         cov_sd(j,i)=sd_[i][index][j];
         cov_sd_copy(j,i)=sd_[i][index][j];
       }
    }
   vector<int> indice; 
//   cout<<"Remove significant signals"<<endl;
//   cout<<"Size of cov_z_score is: "<<cov_z_score.n_rows<<"  and"<<cov_z_score.n_cols<<endl;
   for(int i=0;i<cov_z_score.n_rows;i++)
    {
      int sym=0;
      
      for(int j=0;j<cov_z_score.n_cols;j++)
       {
//         cout<<i<<" "<<j<<endl;
         if(abs(cov_z_score(i,j))>2)
          {
             sym=1;
          }
       }
      if(sym==1)
       {
//         cout<<"I have to remove the variant "<<i<<endl;
         indice.push_back(i);
       }
    }
//   cout<<"Start to remove rows "<<indice.size()<<endl;
   uvec Indice = conv_to< uvec >::from(indice);
   cov_z_score_copy.shed_rows(Indice);
   cov_sd_copy.shed_rows(Indice);
//   cout<<"Variants with significant signals are removed"<<endl;

//   cout<<"Data used to contruct covariance: "<<endl;
//      cout<<cov_z_score_copy<<endl;
//      cout<<"Output data used to construct covariance is done"<<endl;  
 
   mat Cov =cov(cov_z_score_copy);


   for(int i=0;i<cov_z_score.n_rows;i++)
    {
      for(int X=0;X<Cov.n_rows;X++)
       {
         for(int Y=0;Y<Cov.n_rows;Y++)
          {
             if(X==Y)
              {
                Cov(X,Y)=cov_sd(i,X)*cov_sd(i,Y);
                
              } else
              {
                Cov(X,Y)=Cov(X,Y)*cov_sd(i,X)*cov_sd(i,Y);
              }
          }
       }
      mat inv_Cov = inv(Cov);
      mat One =mat(1,cov_beta.n_cols,fill::ones);
      mat meta_beta=One*inv_Cov*trans(cov_beta);
      mat Biao=One*inv_Cov*trans(One);
      double  meta_sd=sqrt(1/Biao(0,0));
      mat meta_z = meta_beta / meta_sd;
//      cout<<"Summary result for peak signal "<<index<<endl;
          double meta_logP;
          if(meta_beta(0,0)<=0)
            {
               meta_logP = log10(2.0*normcdf(meta_z(0,0)));
            } else
            {
              meta_logP = log10(2.0*(1-normcdf(meta_z(0,0))));
            }
  //        cout<<"Variant "<<i<<"  is: "<<meta_z(0,0)<<" "<<meta_logP<<endl;      
   }
  }

void ReadSummary(string input, map<string, vector<double>> &summary)
 {

   ifstream FHIN(input.c_str(),ios::in);
   if(!FHIN)
    {
      return ;
    }
   string A;
   int b=-1;
   while(FHIN && getline(FHIN, A))
    {
      b++;
      if(b==0)
       {
         continue;
       }
         stringstream ss(A);
         vector<string> X;
         string B;
         string variant;
         vector<double> Y;
         ss>>variant;
         double beta;
         ss>>beta;
         double sd;
         ss>>sd;
         Y.push_back(beta);
         Y.push_back(sd);
         summary.insert(pair<string, vector<double>>(variant, Y));
    }
   FHIN.clear();
   FHIN.close();     
 }


//int localization(vector<vector<string> > &variant_, mat &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF)
int localization(vector<vector<string> > &variant_, map <string, vector<double>>  &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF)
 {
//   cout<<"Come to localization function to perform ecaviar-like analysis"<<endl;
   for(map<string,string>::iterator it=GWAS_file.begin();it!=GWAS_file.end();it++ )
    {
      map<string, vector<double>> summary;
      ReadSummary(it->second, summary);

      map<string, vector<double>> H4;
      for(int i=0;i<variant_[0].size();i++)
       {
         if(summary.find(variant_[0][i])!=summary.end())
          {
            double r1;
            double r2;
            double z1;
            double z2;
            vector<double> res;
            double prior1=0.0001;
            double prior2=0.0001;
            double prior12=0.000001;
//            double beta1=summary[variant_[0][i]][0];
//            z1=beta1/summary[variant_[0][i]][1];

            double beta1=summary[variant_[0][i]][0];
            z1=beta1/summary[variant_[0][i]][1];


            if(Res_fixed.find(variant_[0][i])==Res_fixed.end())
             {
               continue;
             }
            double beta2=Res_fixed[variant_[0][i]][0];
            z2=Res_fixed[variant_[0][i]][1];

//            double beta2=Res_fixed[variant_[0][i]][0];
//            z2=Res_fixed[variant_[0][i]][1];
//            if(beta2<0 && z2<0)
//             {
//               continue;
//             } 
            double var1 = (beta1/z1)*(beta1/z1);
            double var2 = (beta2/z2)*(beta2/z2);

            r1 = 0.15*0.15/(0.15*0.15+var1);
            r2 = 0.2* 0.2/(0.2*0.2+var2);
            
            double abf1=0.5*(log(1-r1)+r1*z1*z1);
            double abf2=0.5*(log(1-r2)+r2*z2*z2);
            double h0=0.0;
            double h1=abf1+log(prior1);
            double h2=abf2+log(prior2);
            double h3=abf1+abf2+log(prior1)+log(prior2);
            double h4=abf1+abf2+log(prior12);
            double sum=(exp(h0)+exp(h1)+exp(h2)+exp(h3)+exp(h4));
            h0=exp(h0)/sum;
            h1=exp(h1)/sum;
            h2=exp(h2)/sum;
            h3=exp(h3)/sum;
            h4=exp(h4)/sum; 
            res.push_back(h0);res.push_back(h1);res.push_back(h2);res.push_back(h3);res.push_back(h4);
            H4.insert(pair<string, vector<double>>(variant_[0][i],res));
          }
       }

      if(H4.size()==0)
       {
         continue;
       }
      string Output=output_prefix+string("_colocalized_with_")+it->first;
      fstream  FHOU(Output.c_str(),ios::out);
      if(!FHOU)
       {
         cout<<"Error for input or output"<<endl;
         return (1);
       }
     FHOU<<"variant	H0	H1	H2	H3	H4"<<endl;
     for(map<string,vector<double>>::iterator T=H4.begin();T!=H4.end();T++ )
      {
        FHOU<<T->first<<"	"<<H4[T->first][0]<<"	"<<H4[T->first][1]<<"	"<<H4[T->first][2]<<"	"<<H4[T->first][3]<<"	"<<H4[T->first][4]<<endl;
      }
     FHOU.clear();
     FHOU.close(); 
      
   }
  return 0;
 }



void prepare_gwas_summary(string input, string output,map <string, vector<double>>  &Res_fixed, map<string, bool> &index,vector<string> &target,map<string,bool> &MAP)
      {
         ifstream FHIN(input.c_str(),ios::in);
         if(!FHIN)
          {
            return ;
          }


//         cout<<"Output variants with LD information available"<<endl;
//         cout<<"*****************************************************"<<endl;

  //       for(map<string,bool>::iterator it=MAP.begin();it!=MAP.end();it++)
    //      {
      //      cout<<it->first<<endl;
        //  }

//         cout<<"********************************************************"<<endl;
         fstream  FHOU(output.c_str(),ios::out);
         if(!FHOU)
          {
            cout<<"Error for input or output"<<endl;
            return ;
          }
         string A;
         int b=-1;
         while(FHIN && getline(FHIN, A))
          {
            b++;
//            if(b==0)
//             {
//               continue;
//             }
            stringstream ss(A);
            vector<string> X;
            string B;
            string variant;
            vector<double> Y;
            ss>>variant;
//            cout<<"The explored in GWAS is "<<variant<<endl;
          
            if(Res_fixed.find(variant)==Res_fixed.end() || MAP.find(variant)==MAP.end())
             {
//               cout<<"There is no eQTL infor or LD"<<endl;
               index.insert(pair<string, bool>(variant, false));
               continue;
             }
//            cout<<"Variant "<<variant<<" can be explored"<<endl; 
            double beta;
            ss>>beta;
            double sd;
            ss>>sd;
            double z=beta/sd;
            FHOU<<variant<<"\t";
            FHOU<<z<<endl;
            index.insert(pair<string, bool>(variant, true));
            target.push_back(variant);
         }
        FHIN.clear();
        FHIN.close();
        FHOU.clear();
        FHOU.close();

      }


void prepare_gwas_summary(map <string, vector<double>>  &Res_fixed, string output, map<string, bool> &index, vector<string> &target)
      {
         fstream  FHOU(output.c_str(),ios::out);
         if(!FHOU)
          {
            cout<<"Error for input or output"<<endl;
            return ;
          }
         for(int i=0;i<target.size();i++)
          {
            if(Res_fixed.find(target[i])!=Res_fixed.end())
             {
               if(index.find(target[i])!=index.end())
                {
                  if(index[target[i]])
                   {
                     double z=Res_fixed[target[i]][1];
                     FHOU<<target[i]<<"\t"<<z<<endl;
                   }
                }
             }
          }
        FHOU.clear();
        FHOU.close();
      }


 void prepare_ld(string input,string output, map<string,bool> & index)
      {
/*
        cout<<"Come into function prepare_ld"<<endl;
        cout<<"input file is "<<input<<endl;
        cout<<"Output index"<<endl;
        for(map<string,bool>::iterator it=index.begin();it!=index.end();it++ )
         {
           if(it->second)
            {
              cout<<" "<<it->first;
            }
         }  
        cout<<endl;
*/
        ifstream FHIN(input.c_str(),ios::in);
         if(!FHIN)
          {
            return ;
          }
         fstream  FHOU(output.c_str(),ios::out);
         if(!FHOU)
          {
            cout<<"Error for input or output"<<endl;
            return ;
          }
         string A;
         int b=-1;
         map<int, bool> target;
         while(FHIN && getline(FHIN, A))
          {
            b++;
            if(b==0)
             {
               stringstream ss(A);
               string waste;
               ss>>waste;
               string variant;
               int indx=1;
               while(ss>>variant)
                {
                  if(index.find(variant)!=index.end())
                   {
                     target.insert(pair<int, bool>(indx,true));
                   } 
                  indx++;
                }
               continue;
             }
            if(target.find(b)==target.end())
             {
               continue;
             }
//            cout<<"Great, one variant for colocalization detected: "<<b<<endl;
            stringstream ss(A);
            string variant;
            ss>>variant;
//            cout<<"Variant is: "<<variant<<endl;
            double r;
            int b1=0;
            bool indx2=false;
            while(ss>>r)
             {
               b1++;
               if(target.find(b1)!=target.end())
                {
                  if(!indx2)
                   {
                     FHOU<<r;
      //               cout<<" "<<r;
                     indx2=true;

                   }else
                   {
                     FHOU<<" "<<r;
        //             cout<<" "<<r;
                   }
                }
             }
            FHOU<<endl;
        //    cout<<endl;
         }
        FHIN.clear();
        FHIN.close();
        FHOU.clear();
        FHOU.close();


      }



int localization(vector<vector<string> > &variant_, map <string, vector<double>>  &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF, string colocalization_method,string gene, bool finemap)
 {
   
   if(finemap==true)
    {
//      cout<<"Doing finemapping"<<endl;
      string ldFile1=gene+string("/ldFile");
      string ldFile2=gene+string("/ldFile2");

      ifstream check_ldFile(ldFile1.c_str());
      if(!check_ldFile)
       {
  //       cout<<"No LD infor provided"<<endl;
         ofstream outfile(ldFile1.c_str(), ios::out );
        outfile<<"Variant";
        vector<string> Test_variant;
        for(map<string, vector<double>>::iterator it=Res_fixed.begin();it!=Res_fixed.end();it++)
                Test_variant.push_back(it->first) ;
        
        for(int i=0;i<Test_variant.size();i++)
         {
           outfile<<" "<<Test_variant[i];
         }
        outfile<<endl;

        for(int i=0;i<Test_variant.size();i++)
         {
           outfile<<Test_variant[i];
           for(int j=0;j<Test_variant.size();j++)
            {
              if(i==j)
               {
                  outfile<<" "<<1;
               } else
               {
                  outfile<<" "<<0;
               }
            }
           outfile<<endl;
         }        

        outfile << endl;
        outfile.close();
       }



     string gwas_summary2=gene+string("/eQTL_summary1");
     string gwas_summary1=gene+string("/gwas_summary2");

     string GWAS_file=gene+string("/GWAS_test");

     fstream  FHOU(GWAS_file.c_str(),ios::out);
     if(!FHOU)
      {
        cout<<"Error for input or output"<<endl;
        return (1);
      }
     for(map<string,vector<double>>::iterator it=Res_fixed.begin();it!=Res_fixed.end();it++)
      {
        FHOU<<it->first<<"\t"<<Res_fixed[it->first][1]<<"\t1"<<endl;
      }
     FHOU.clear();
     FHOU.close();

     map<string, bool> index;
     vector<string> target;



     map<string, bool> MAP;
     importDataFirstColumn(ldFile1,MAP,1);

     prepare_gwas_summary(GWAS_file, gwas_summary1, Res_fixed, index,target,MAP);
     if(target.size()==0)
      {
//         cout<<"Unfortunately, there is no variant to explore"<<endl;
               char c[gwas_summary1.size() + 1];
               strcpy(c, gwas_summary1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
        return (1);
      }
     prepare_gwas_summary(Res_fixed,  gwas_summary2, index,target);

     prepare_ld(ldFile1,ldFile2,index);



     string outputFileName=output_prefix+string("_caviar_finemapping_with_");
     int totalCausalSNP =1;
     double NCP=6.0;
     double rho=0.95;
     bool histFlag =false;


     CaviarModel gwasModel_1(ldFile2, gwas_summary1, outputFileName + "_1", totalCausalSNP, NCP, rho, histFlag);
     gwasModel_1.run();
     gwasModel_1.finishUp();
     return 0;
    }

   for(map<string,string>::iterator it=GWAS_file.begin();it!=GWAS_file.end();it++ )
    {

//     cout<<"Come into colocalization function"<<endl;
//     string ldFile=gene+string("/ldFile");
     string ldFile1=gene+string("/ldFile");
     string ldFile2=gene+string("/ldFile2");
 

     ifstream check_ldFile(ldFile1.c_str());
      if(!check_ldFile)
       {
  //       cout<<"No LD infor provided"<<endl;
         ofstream outfile(ldFile1.c_str(), ios::out );
        outfile<<"Variant";
        vector<string> Test_variant;
        for(map<string, vector<double>>::iterator it=Res_fixed.begin();it!=Res_fixed.end();it++)
                Test_variant.push_back(it->first) ;

        for(int i=0;i<Test_variant.size();i++)
         {
           outfile<<" "<<Test_variant[i];
         }
        outfile<<endl;

        for(int i=0;i<Test_variant.size();i++)
         {
           outfile<<Test_variant[i];
           for(int j=0;j<Test_variant.size();j++)
            {
              if(i==j)
               {
                  outfile<<" "<<1;
               } else
               {
                  outfile<<" "<<0;
               }
            }
           outfile<<endl;
         }

        outfile << endl;
        outfile.close();
       }
     


 
//     string ldFile2=gene+string("//ldFile2");
     string gwas_summary2=gene+string("/eQTL_summary1");
     string gwas_summary1=gene+string("/gwas_summary2");
     

     map<string, bool> index;
     vector<string> target;

     
//     prepare_gwas_summary(Res_fixed, gwas_summary1, index, vector<string> &target);

//     prepare_gwas_summary(it->second,gwas_summary2,Res_fixed,index,target);
     map<string, bool> MAP;
//     cout<<"Extract variants with LD profile available"<<endl;
     importDataFirstColumn(ldFile1,MAP,1);

//     cout<<"Prepare gwas summary results"<<endl;
     prepare_gwas_summary(it->second, gwas_summary1, Res_fixed, index,target,MAP);     
     if(target.size()==0)
      {
//        cout<<"There is no explored variant"<<endl;
//        string file1=gene+string("/affected_allele_tissue_")+to_string(i);
               char c[gwas_summary1.size() + 1];
               strcpy(c, gwas_summary1.c_str());
               if(remove(c))
                {
                  cout<<"Successfully delete"<<endl;
                }
        continue;
      }
//     cout<<"Number of explored variants is: "<<target.size()<<endl;
//     cout<<"Prepare eQTL summary results"<<endl;
     prepare_gwas_summary(Res_fixed,  gwas_summary2, index,target);

//     cout<<"Prepare LD structure matrix"<<endl;
     prepare_ld(ldFile1,ldFile2,index);     



     string outputFileName=output_prefix+string("_ecaviar_colocalized_with_")+it->first;
     int totalCausalSNP =1;
     double NCP=6.0;
     double rho=0.95;
     bool histFlag =false;
//     cout<<"First, perform finemapping for gwas summary result"<<endl;


     CaviarModel gwasModel_1(ldFile2, gwas_summary1, outputFileName + "_1", totalCausalSNP, NCP, rho, histFlag);


//     cout<<"Secondly, perform finemapping for eQTL summary result"<<endl;


     CaviarModel gwasModel_2(ldFile2, gwas_summary2, outputFileName + "_2", totalCausalSNP, NCP, rho, histFlag);



     gwasModel_1.run();
     gwasModel_1.finishUp();

     gwasModel_2.run();
     gwasModel_2.finishUp();
     
     int snpCount  = gwasModel_1.snpCount;
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
     for(int i = 0; i < snpCount; i++) 
      {
        outfile << snpNames[i] << "\t" << (stat1[i] * stat2[i])/sumCLPP << "\t" << stat1[i] * stat2[i] << endl;
      }
     outfile.close();

     string file1=outputFileName+"_1"+"_post";
     string file2=outputFileName+"_2"+"_post";
     string file3=outputFileName+"_1"+"_set";
     string file4=outputFileName+"_2"+"_set";
     
               char c1[file1.size() + 1];
               strcpy(c1, file1.c_str());
               if(remove(c1))
                {
                  cout<<"Successfully delete"<<endl;
                }
    char c2[file2.size() + 1];
               strcpy(c2, file2.c_str());
               if(remove(c2))
                {
                  cout<<"Successfully delete"<<endl;
                }

    char c3[file3.size() + 1];
               strcpy(c3, file3.c_str());
               if(remove(c3))
                {
                  cout<<"Successfully delete"<<endl;
                }

    char c4[file4.size() + 1];
               strcpy(c4, file4.c_str());
               if(remove(c4))
                {
                  cout<<"Successfully delete"<<endl;
                }
     



/*

      map<string, vector<double>> summary;
      ReadSummary(it->second, summary);

      map<string, vector<double>> H4;
      for(int i=0;i<variant_[0].size();i++)
       {
         if(summary.find(variant_[0][i])!=summary.end())
          {
            double r1;
            double r2;
            double z1;
            double z2;
            vector<double> res;
            double prior1=0.0001;
            double prior2=0.0001;
            double prior12=0.000001;
            double beta1=summary[variant_[0][i]][0];
            z1=beta1/summary[variant_[0][i]][1];
            double beta2=Res_fixed[variant_[0][i]][0];
            z2=Res_fixed[variant_[0][i]][1];
            if(beta2<0 && z2<0)
             {
               continue;
             }
            double var1 = (beta1/z1)*(beta1/z1);
            double var2 = (beta2/z2)*(beta2/z2);

            r1 = 0.15*0.15/(0.15*0.15+var1);
            r2 = 0.2* 0.2/(0.2*0.2+var2);

            double abf1=0.5*(log(1-r1)+r1*z1*z1);
            double abf2=0.5*(log(1-r2)+r2*z2*z2);
            double h0=0.0;
            double h1=abf1+log(prior1);
            double h2=abf2+log(prior2);
            double h3=abf1+abf2+log(prior1)+log(prior2);
            double h4=abf1+abf2+log(prior12);
            double sum=(exp(h0)+exp(h1)+exp(h2)+exp(h3)+exp(h4));
            h0=exp(h0)/sum;
            h1=exp(h1)/sum;
            h2=exp(h2)/sum;
            h3=exp(h3)/sum;
            h4=exp(h4)/sum;
            res.push_back(h0);res.push_back(h1);res.push_back(h2);res.push_back(h3);res.push_back(h4);
            H4.insert(pair<string, vector<double>>(variant_[0][i],res));
          }
       }

      if(H4.size()==0)
       {
         continue;
       }
      string Output=output_prefix+string("_colocalized_with_")+it->first;
      fstream  FHOU(Output.c_str(),ios::out);
      if(!FHOU)
       {
         cout<<"Error for input or output"<<endl;
         return (1);
       }
     FHOU<<"variant     H0      H1      H2      H3      H4"<<endl;
     for(map<string,vector<double>>::iterator T=H4.begin();T!=H4.end();T++ )
      {
        FHOU<<T->first<<"       "<<H4[T->first][0]<<"   "<<H4[T->first][1]<<"   "<<H4[T->first][2]<<"   "<<H4[T->first][3]<<"   "<<H4[T->first][4]<<endl;
      }
     FHOU.clear();
     FHOU.close();
*/
   }
  return 0;
 }
















double MetaAnalysis(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > >&beta_, vector<vector<vector<double> > >&sd_,int index, mat &Res_fixed, mat &Res_random, string meta_mode, string output, vector<vector<string> > &variant_, map<string, bool> &pos, int colocalization, map<string, string> &GWAS_file, map<string, double > &MAF, vector<double > &VARIANCE, map<int,bool> &data_index,string gene, bool INDEX, bool finemap, map<string, string> & ALLELE, double p_cutoff,string cor_file, bool write_cor, bool Han)
 {
//   cout<<"Come into MetaAnalysis function to perform final meta-analysis"<<endl;
//   cout<<"Come into MetaAnalysis function to perform final meta-analysi"<<endl;
   string output1= output+"_statistical_signal";
   map <string, bool> pos1;
   map <string, vector<double>> Res;
   pos1=pos;
   vector<double> test1;
   vector<double> test2;
   vector<double> test3;
//   cout<<"Start.............."<<endl;
//   cout<<"Infor1 is "<<z_score_.size()<<" "<<beta_.size()<<" "<<sd_.size()<<" "<<variant_.size()<<" "<<data_index.size()<<endl;
//   cout<<"Infor2 is "<<z_score_[0][0].size()<<" "<<beta_[0][0].size()<<" "<<sd_[0][0].size()<<" "<<variant_[0].size()<<endl; 


   double P=ConstructCovRandom(z_score_, beta_, sd_ ,index , Res_fixed, Res_random, meta_mode, output1, variant_, pos, pos1, Res, MAF, VARIANCE, data_index, test1, test2,test3, ALLELE, cor_file, write_cor, Han);
//   cout<<"P is: "<<P<<endl;
   if(P<0)
    {
      cout<<"Abnormal, as resulted P value is "<<P<<" , skip"<<endl;
      return P;
    }
//   cout<<"Final meta analysis is over"<<endl;
   string output2= output+"_Mvalue";
   string localization_method="ecaviar";
//   cout<<"Meta analysis is over, and I will determine perform localization or not"<<endl;
//   cout.precision(17);
//   cout<<"Resulted p-value from meta-analysis is: "<<P<<endl;
//   cout<<"colocalization is: "<<colocalization<<", finemap is: "<<finemap<<", P is: "<<P<<", INDEX is: "<<INDEX<<endl; 
   if((colocalization==1 || finemap)  && P<=p_cutoff && INDEX)
    {
//      cout<<"Perform colocalization analysis"<<endl;
//      localization(variant_,Res, GWAS_file, output1, MAF);
//      localization(variant_,Res, GWAS_file, output1, MAF,localization_method,gene);   //Changed by Biao Zeng 01/24/2020
      localization(variant_,Res, GWAS_file, output1, MAF,localization_method,gene,finemap);
//      cout<<"colocalization analysis is over"<<endl;
//      localization(vector<vector<string> > &variant_, map <string, vector<double>>  &Res_fixed, map<string, string> &GWAS_file, string output_prefix, map<string, double > &MAF);
//
//
//
//      localization(variant_,Res, GWAS_file, output_prefix, map<string, double > &MAF, string colocalization_method,string gene)
    }
  //Changed by Biao Zeng, 11/02/2020, 5:29 pm
   INDEX=false;
   if(INDEX && z_score_.size()>1)
    {
//   cout<<"Come to estimate Mvalue"<<endl;
   calculate_Mvalue(z_score_, beta_, sd_ ,index , output2,variant_, pos, pos1, data_index);  
    }
   return P;
 }
















double MetaAnalysis(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > >&beta_, vector<vector<vector<double> > >&sd_,int index, mat &Res_fixed, mat &Res_random, string meta_mode, string output, vector<vector<string> > &variant_, map<string, bool> &pos, int colocalization, map<string, string> &GWAS_file, map<string, double > &MAF, vector<double > &VARIANCE, map<int,bool> &data_index,string gene,vector<double> &test1_1, vector<double> &test2_1, vector<double> &test3_1, map<string, string> &ALLELE, string cor_file, bool write_cor, bool Han)
 {
//   cout<<"Come into MetaAnalysis function to perform final meta-analysis"<<endl;
   string output1= output+"_statistical_signal";
   map <string, bool> pos1;
   map <string, vector<double>> Res;
   pos1=pos;
//   cout<<"Perform real meta-analysis"<<endl;
   
   double P=ConstructCovRandom(z_score_, beta_, sd_ ,index , Res_fixed, Res_random, meta_mode, output1, variant_, pos, pos1, Res, MAF, VARIANCE, data_index,test1_1,test2_1, test3_1, ALLELE, cor_file, write_cor, Han);
   
   
   return P;
 }


       int ReadInputFile(string input, vector<string> &output)
         {
//           cout<<"Come into function ReadInputFile "<<endl;
           int index = 0;
           ifstream fin(input.c_str(), std::ifstream::in);
           vector<double> words;
           string line;
           while(fin && getline(fin, line))
             {
               stringstream ss(line);
               string Str;
               ss >> Str;
  //             cout<<"Str is: "<<Str<<endl;
               output.push_back(Str);
             }
           return(0);
         }





       int ReadInputFile(string input, map<string,double> &output)
         {
           int index = 0;
           ifstream fin(input.c_str(), std::ifstream::in);
           vector<double> words;
           string line;
           while(fin && getline(fin, line))
             {
               stringstream ss(line);
               string Str;
               double Str1;
               ss >> Str;
               ss >> Str1;
               output.insert(pair<string, double>(Str, Str1));
             }
           return(0);
         }      
 
       int ReadInputFile(string input, map<string,string> &output)
         {
  //         cout<<"Come into function ReadInputFile"<<endl;
    //       cout<<"input file is: "<<input<<endl;
           int index = 0;
           ifstream fin(input.c_str(), std::ifstream::in);
           vector<double> words;
           string line;
           while(fin && getline(fin, line))
             {
               stringstream ss(line);
               string Str;
               string Str1;
               ss >> Str;
               ss >> Str1;
      //         cout<<"File infor: "<<Str<<" "<<Str1<<endl;
               output.insert(pair<string, string>(Str, Str1));
             }
           return(0);
         }


       int ReadInputFile(string input, vector<string> &output, vector<string> &allele_file)
         {
           int index = 0;
           ifstream fin(input.c_str(), std::ifstream::in);
           vector<double> words;
           string line;
           while(fin && getline(fin, line))
             {
               stringstream ss(line);
               string Str;
               ss >> Str;
               output.push_back(Str);
               ss >> Str;
               allele_file.push_back(Str);
             }
           return(0);
         }




void read_cis_summary(string file,vector<vector<double> > &z_score,vector<vector<double> > &beta,vector<vector<double> > &sd,vector<string> &variant, map<string, string> &allele)
 {

   

//   cout<<"Come into function read_cis_summary to extract eQTL signals"<<endl; 
   vector<double> beta_;
   vector<double> sd_;
   vector<double> z_;
    ifstream FHIN(file.c_str(),ios::in);
         if(!FHIN)
          {
            return ;
          }

         string A;
         int b=-1;
         while(FHIN && getline(FHIN, A))
          {
            b++;
            stringstream ss(A);
            vector<string> X;
            string B;
            string variant_;
            string allele_;
            ss>>variant_;
            ss>>allele_;
            double Beta_;
            ss>>Beta_;
            double Sd_;
            ss>>Sd_;
            double Z_;
            ss>>Z_;
            variant.push_back(variant_);
            allele.insert(pair<string,string>(variant_,allele_));
            z_.push_back(Z_);
            beta_.push_back(Beta_);
            sd_.push_back(Sd_);
         }
        z_score.push_back(z_);
        beta.push_back(beta_);
        sd.push_back(sd_);
//        cout<<"Tasting is over"<<endl;
        FHIN.clear();
        FHIN.close();
 }







int MM_conditional_analysis(vector<string> &file_, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,vector<string> &covariate_, vector <string> &grm_file_,string outputFile,vector<string> &geno_file_, string meta_mode, vector<string> &allele_file_, map<string, bool> &pos, int colocalization, map<string, string> &GWAS_file, map<string, double> &MAF,map<int,bool> &data_index, int nested, bool finemap, double p_cutoff, bool cis_summary, string cor_file, vector<string> &causal, int eQTL_number, bool Han)
 {
   int index=1;
   int Peak =0;
   vector<string> removal;
   map<string, int> Var;
   
   for(int i=0;i<file_.size();i++)
    {
      if(file_[i] =="NA")
       { 
         continue ;
       }
//      cout<<"Extract explored variants in dataset "<<i<<endl;
     if(!cis_summary)
      {
        import_first_row(file_[i], Var);
      } else
      {
        import_first_col(file_[i], Var);
//        for(map<string, int>::iterator it; it!=)
        for(map<string,int>::iterator it=Var.begin();it!=Var.end();it++)
         {
           pos.insert(pair<string,bool>(it->first,true));
         }
      }
    }
//   cout<<"I am here1"<<endl;
   map<string, string> ALLELE;
   bool write_cor=true;
   while(index==1)
    {
      if(cor_file!="")
       {
         index=0;
       } 
      vector<vector<vector<double> > > z_score_;
      vector<vector<vector<double> > > beta_;
      vector<vector<vector<double> > > sd_;
      vector<vector<vector<double> > > direction_;
      vector<vector<vector<string> > > allele_;
    
      vector<vector<string> > variant_;
      vector<double> VARIANCE;
//      cout<<"Start to perform eQTL or extract eQTL signal"<<endl;
      for(int I=0, indx=0;I<file_.size();)
       {
         vector<vector<double> > z_score;
         vector<vector<double> > beta;
         vector<vector<double> > sd;
         map<string,string>  allele;
         vector<string> variant;
         if(file_[I] =="NA")
          {
            data_index[I]=false;
            I++;
            continue;
          }
         string OutputFile=outputFile+string("_Tissue_")+to_string(I);
         vector<string> removal1;
         removal1=removal;
         double Variance=0.0;
         int conditional_index=0;
         if(!cis_summary)
          {


         conditional_index=conditional_analysis(file_[I], gene, totalCausalSNP, rho,histFlag,  gamma, weight, nthread, covariate_[I], grm_file_[I], OutputFile, geno_file_[I],z_score, beta, sd,1, variant, removal1, I, colocalization, Variance);
            VARIANCE.push_back(Variance);
         read_allele(allele_file_[I],allele);
          } else
          {
             read_cis_summary(file_[I],z_score,beta,sd,variant,allele);
             VARIANCE.push_back(1);
          }
//         VARIANCE.push_back(Variance);
//         read_allele(allele_file_[I],allele);
         map<string,int> Test;
//         cout<<"Value of conditional_index is "<<conditional_index<<endl;
         if(conditional_index!=-1)
          {
         for(int i=0;i<variant.size();i++)
          {
            Test.insert(pair<string,int>(variant[i],i));
          }
         for(map<string,int>::iterator it=Var.begin();it!=Var.end();it++)
          {
            if(ALLELE.find(it->first) != ALLELE.end())
             {
             } else
             {
               if(allele.find(it->first) != allele.end())
                {
                  ALLELE.insert(pair<string,string>(it->first,allele[it->first]));
                }
               
             }
          }
         }
         vector<vector<double> > z_score_1;
         vector<vector<double> > beta_1;
         vector<vector<double> > sd_1;
         vector<vector<double> > direction_1;
         vector<double> z_score_1_1;
         vector<double>  beta_1_1;
         vector<double>  direction_1_1;
         vector<double > sd_1_1;
         vector<string> variant_1;
//         cout<<"Start to combine the eQTL signals"<<endl;
         for(map<string,int>::iterator it=Var.begin();it!=Var.end();it++)
          {
//            cout<<"Test for variant "<<it->first<<endl;
            if(Test.find(it->first) != Test.end() && conditional_index!=-1)
             {
//               cout<<"OK, Get eQTL signal"<<endl;
//               cout<<"allele is: "<<allele[it->first]<<endl;
//               cout<<"In ALLELE, it is: "<<ALLELE[it->first]<<endl;
               if(allele[it->first]!=ALLELE[it->first])
                {
                  
                  double BETA=0.0-beta[0][Test[it->first]];
                  double z   =0.0-z_score[0][Test[it->first]];
                  z_score_1_1.push_back(z);
                  beta_1_1.push_back(BETA);
                  sd_1_1.push_back(sd[0][Test[it->first]]);
                  direction_1_1.push_back(-1.0);
                  variant_1.push_back(it->first);
                } else
                {
                  z_score_1_1.push_back(z_score[0][Test[it->first]]);
                  beta_1_1.push_back(beta[0][Test[it->first]]); 
                  sd_1_1.push_back(sd[0][Test[it->first]]);     
                  direction_1_1.push_back(1.0);
                  variant_1.push_back(it->first);
                }
             } else
             {
               z_score_1_1.push_back(0);
               beta_1_1.push_back(0);
               sd_1_1.push_back(0);
               direction_1_1.push_back(0);
               variant_1.push_back(it->first);
             }
         }
         z_score_1.push_back(z_score_1_1);
         beta_1.push_back(beta_1_1);
         sd_1.push_back(sd_1_1);
         direction_1.push_back(direction_1_1);
         
//         cout<<"Recording QTL signal is over"<<endl;
         z_score_.push_back(z_score_1);
         beta_.push_back(beta_1);
         sd_.push_back(sd_1);
         direction_.push_back(direction_1);
         variant_.push_back(variant_1);
         I++;
         indx++;
//         cout<<"Here1"<<endl;
       }
     Peak++;
     int N_peak=9999;
     for(int x=0;x<beta_.size();x++)
      {
        if(N_peak>beta_[x].size())
         {
           N_peak=beta_[x].size();
         }   
      }
//     cout<<"Here2"<<endl;
     for(int x=0;x<N_peak;x++)
      {
         double P;
         int nested_1=nested;
         if(nested>0)
          {
         for(int indx=0;indx<nested;indx++)
           {
             if(file_[indx]=="NA")
              {
                nested_1--;
              }
           } 
          }        
//         cout<<"Here3_1"<<endl; 
         if(nested_1<0)
          {
            nested_1=1;
          }
         mat Res_fixed_2  = mat(z_score_.size()-nested_1+1,3,fill::zeros);
         mat Res_random_2 = mat(z_score_.size()-nested_1+1,3,fill::zeros);
         if(nested>0 && nested_1 >1 && file_.size()>nested)
          {
            vector<vector<vector<double> > > z_score_nested;
            vector<vector<vector<double> > > beta_nested;
            vector<vector<vector<double> > > sd_nested;
            vector<vector<string> > variant_nested;
            mat Res_fixed  = mat(nested_1,3,fill::zeros);
            mat Res_random = mat(nested_1,3,fill::zeros);
            for(int indx1=0;indx1<nested_1;indx1++)
             {
                variant_nested.push_back(variant_[indx1]);
                z_score_nested.push_back(z_score_[indx1]);
                beta_nested.push_back(beta_[indx1]);
                sd_nested.push_back(sd_[indx1]);
             }
            vector<double> test1_1;
            vector<double> test2_1;
            vector<double> test3_1;
            vector<vector<double>> test1;
            vector<vector<double>> test2;
            vector<vector<double>> test3;
            string Output1 = outputFile + "_temp1_peak_" + to_string(Peak);  
            map<int, bool> data_index1;
            for(int indx=0,indx1=0;indx1<nested;indx1++)
             {
                data_index1.insert(pair<int,bool>(indx1, data_index[indx1]));
             }
            string test_file=""; 
            MetaAnalysis(z_score_nested,  beta_nested, sd_nested, x, Res_fixed, Res_random, "fixed", Output1, variant_nested, pos, 0,GWAS_file, MAF, VARIANCE, data_index1, gene, test1_1, test2_1, test3_1, ALLELE, test_file,false, Han);
           
            test1.push_back(test1_1);
            test2.push_back(test2_1);
            test3.push_back(test3_1); 
            vector<vector<vector<double> > > z_score_nested2;
            vector<vector<vector<double> > > beta_nested2;
            vector<vector<vector<double> > > sd_nested2;
            vector<vector<string> > variant_nested2;
            z_score_nested2.push_back(test1);
            beta_nested2.push_back(test2);
            sd_nested2.push_back(test3);
            variant_nested2.push_back(variant_[0]);
            for(int indx1=nested_1;indx1<z_score_.size();indx1++)
             {
               z_score_nested2.push_back(z_score_[indx1]);
               beta_nested2.push_back(beta_[indx1]);
               sd_nested2.push_back(sd_[indx1]);
               variant_nested2.push_back(variant_[indx1]);
             }
            string Output2 = outputFile + "_nested_peak_" + to_string(Peak);
            map<int, bool> data_index2;
            vector<double> VARIANCE2;
            data_index2.insert(pair<int,bool>(0,true));
            VARIANCE2.push_back(1);
            for(int indx1=nested,indx2=0;indx1<data_index.size();indx1++)
             {
               if(!data_index[indx1])
                {
                  continue;
                }
               indx2++;
               data_index2.insert(pair<int,bool>(indx2, data_index[indx1]));
               VARIANCE2.push_back(VARIANCE[nested_1+indx2-1]);
               
             }
            bool INDEX=true;
            P=MetaAnalysis(z_score_nested2,  beta_nested2, sd_nested2, x, Res_fixed_2, Res_random_2, "random", Output2, variant_nested2, pos, 0,GWAS_file, MAF, VARIANCE2, data_index2, gene, INDEX, finemap, ALLELE,p_cutoff,cor_file,false, Han);
          }
//          cout<<"Here3_2"<<endl; 
          mat Res_fixed = mat(z_score_.size(),3,fill::zeros);
          mat Res_random = mat(z_score_.size(),3,fill::zeros);
          string Output = outputFile + "_peak_" + to_string(Peak);
          double P1;
          bool INDEX=false;
          if(nested<0)
           {
             INDEX=true;
           }
//          INDEX=false;
//          cout<<"Do the fuc now"<<endl;
//          cout<<"Dimension for the summary result is: "<<endl;
//          cout<<z_score_.size()<<endl;
//          cout<<z_score_[0].size()<<endl;
//          cout<<z_score_[0][0].size()<<endl;
          P1=MetaAnalysis(z_score_,  beta_, sd_, x, Res_fixed, Res_random, meta_mode, Output, variant_, pos, colocalization,GWAS_file, MAF, VARIANCE, data_index, gene, INDEX, finemap, ALLELE,p_cutoff,cor_file, write_cor, Han);
          write_cor=false;
//          cout<<"Here3_3"<<endl;
          if(cis_summary)
           {
             return P1;
           }
          if(nested<0 || nested_1<=1)
           {
             P=P1;
           }
          if(P<0)
           {
             return P;
           }
//        cout<<"Resulted p value is: "<<P<<endl;
        if(P<=p_cutoff)
        {
         if(meta_mode =="random")
          {
            if(nested>0 && nested_1>1)
             {
               string PEAK=variant_[0][int(Res_random_2(0,0))];
               removal.push_back(PEAK);
               causal.push_back(PEAK); 
             } else
             {
               string PEAK=variant_[0][int(Res_random(0,0))];
               removal.push_back(PEAK);
               causal.push_back(PEAK);
             }
          }else
          {
            if(nested>0 && nested_1>1)
             {
               string PEAK=variant_[0][int(Res_fixed_2(0,0))];
               removal.push_back(PEAK);
               causal.push_back(PEAK);
             } else
             {
               string PEAK=variant_[0][int(Res_fixed(0,0))];
               removal.push_back(PEAK);
               causal.push_back(PEAK);
             }
          }
        
//	 cout<<"Let's update phenotype"<<endl;
         for(int i=0,indx=-1;i<file_.size();i++)
          {
            if(file_[i]=="NA")
             {
               continue;
             }
            indx++;
//            cout<<"Try to update phenotype in tissue"<<i<<endl;
            vector< vector<double> >  Data;
            vector<string> ind;
            vector<string> variant;
            map<string, int> VAR;
//            cout<<"Read input file"<<endl;
//            cout<<"file is: "<<file_[i]<<endl;
//            cout<<"explored gene is: "<<gene<<endl;
            read_table(file_[i],Data,ind,variant, VAR);
//            cout<<"Output explored variants"<<endl;


      //      for (map<string,int>::iterator it=VAR.begin(); it!=VAR.end(); ++it)
        //     {
          //     cout << it->first << " => " << it->second <<endl;
 	    // }

            mat data=mat(Data.size(),Data[0].size(),fill::zeros);
            for(int indx=0;indx<Data.size();indx++)
             {
               for(int j=0;j<Data[i].size();j++)
                {
                  data(indx,j)=Data[indx][j];
                }
             }
            vector<double>  pheno;
            vector < vector<string> > fam_infor;
            readfam(file_[i],pheno,fam_infor,1);

//          cout<<"i is: "<<i<<endl;
            

//            cout<<"Peak variant is: "<<" "<<Res_fixed(i,0)<<" "<<variant_[i][int(Res_fixed(i,0))]<<endl;  
            
//            cout<<"Results of Res_random is: "<<Res_random<<endl;
            if(meta_mode == "random")
             {
//               cout<<"In random model, beta value for tissue"<<i<<" is: "<<Res_random(i,1)<<endl;
               string PEAK;
               if(nested>0 && nested_1>1)
                {
                  PEAK=variant_[indx][int(Res_random_2(0,0))];
//                  cout<<"Run in nested manner, and peak variant is: "<<PEAK<<endl;
                }else
                {
                  PEAK=variant_[indx][int(Res_random(0,0))]; 
//                  cout<<"Run in standard manner, and peak variant is: "<<PEAK<<endl;
                }
              // removal.push_back(PEAK);

               map<string,int>::iterator Ite;

               Ite=VAR.find(PEAK);
               if(Ite==VAR.end())
                   {
//                   cout<<"Fuck, there is no genotype infor of causal variant in this data set"<<i<<endl;
                     char c[file_[i].size() + 1];
                     strcpy(c, file_[i].c_str());
                     if(remove(c))
                      {
  //                      cout<<"Successfully delete"<<endl;
                      }
                     file_[i]=string("NA");

                     
                      
                     continue;
                   }
               int COL=VAR[PEAK];
                
//             cout<<"col for causal variant is: "<<COL<<endl;
               if(nested>0 && nested_1>1)
                {
//                  data.col(0)=data.col(0)- direction_[i][0][int(Res_random_2(i,0))]*Res_random(i,1)*data.col(COL); 
                  data.col(0)=data.col(0)- direction_[indx][0][int(Res_random_2(0,0))]*beta_[indx][0][int(Res_random_2(0,0))]*data.col(COL); 
                   
                } else
                { 
//                cout<<"Estimated alleleic effect size is: "<<beta_[indx][0][int(Res_random(0,0))]<<", and direction of the effect is: "<<direction_[indx][0][int(Res_random(0,0))] <<endl;
//                cout<<"Genotype is: "<<endl;
//                cout<<data.col(COL)<<endl;
//                  data.col(0)=data.col(0)- direction_[i][0][int(Res_random(i,0))]*Res_random(i,1)*data.col(COL);
//                cout<<"The effect to remove is: "<<endl;
//                cout<<direction_[indx][0][int(Res_random(0,0))]*beta_[indx][0][int(Res_random(0,0))]*data.col(COL)<<endl;
                  data.col(0)=data.col(0)- direction_[indx][0][int(Res_random(0,0))]*beta_[indx][0][int(Res_random(0,0))]*data.col(COL);
                }
             } else
             {
//             cout<<"In fixed model, beta value for tissue"<<i<<" is: "<<Res_fixed(i,1)<<endl;
              //Changed by Biao Zeng, 04/15/2021, 5:40 pm 
//               string PEAK=variant_[indx][int(Res_fixed_2(0,0))];
//
               string PEAK;
               if(nested>0 && nested_1>1)
                {
                  PEAK=variant_[indx][int(Res_fixed_2(0,0))];
                  cout<<"Run in nested manner, and peak variant is: "<<PEAK<<endl;
                }else
                {
                  PEAK=variant_[indx][int(Res_fixed(0,0))];
                }
           //    removal.push_back(PEAK);
               map<string,int>::iterator Ite;
//             cout<<"Try to find the index of causal variant"<<PEAK<<endl;
               Ite=VAR.find(PEAK);
               if(Ite==VAR.end())
                   {
  //                 cout<<"Fuck, there is no genotype infor of causal variant in this data set"<<i<<endl;
                     file_[i]=string("NA");
                     continue;
                   }

               int COL=VAR[PEAK];
    //         cout<<"col for causal variant is: "<<COL<<endl;
               if(nested>0 && nested_1>1)
                {

      //          cout<<"Estimated alleleic effect size is: "<<beta_[indx][0][int(Res_fixed_2(0,0))]<<", and direction of the effect is: "<<direction_[indx][0][int(Res_fixed_2(0,0))] <<endl;
                  data.col(0)=data.col(0)- direction_[indx][0][int(Res_fixed_2(0,0))]*beta_[indx][0][int(Res_fixed_2(0,0))]*data.col(COL);
                } else
                {
                  data.col(0)=data.col(0)- direction_[indx][0][int(Res_fixed(0,0))]*beta_[indx][0][int(Res_fixed(0,0))]*data.col(COL);
        //        cout<<"Estimated alleleic effect size is: "<<beta_[indx][0][int(Res_fixed(0,0))]<<", and direction of the effect is: "<<direction_[indx][0][int(Res_fixed(0,0))] <<endl;
                }
             }

//          cout<<"Mean is: "<<mean(data.col(0))<<endl;
//          cout<<"Residual value is: "<<endl;
//          cout<<data.col(0)<<endl;
            data.col(0)=data.col(0)-mean(data.col(0)); 

//          cout<<"New residual value is: "<<endl;

//          cout<<data.col(0)<<endl;

//          mat Pheno =data.col(0);
//          cout<<"Come into updatefile function to update phenotype file"<<endl;
//          cout<<"Phenotype file is: "<<file_[i]<<endl;
            string updated =file_[i]+"_updated";


             
        

 
            updatefile(fam_infor,data,updated);
//          cout<<"Update is over"<<endl;
//            string file1=gene+string("/affected_allele_tissue_")+to_string(i);
            char c[file_[i].size() + 1];
            strcpy(c, file_[i].c_str());
            if(remove(c))  
             {
    //           cout<<"Successfully delete"<<endl;
             } 
            

            file_[i]=updated;
          }
         }
   
//FOR TEST
//    P=0.01;
       
         if(eQTL_number>=0   && causal.size()>=eQTL_number)
          {
            index=0;
          }

         if(P>p_cutoff)
          {
            index=0;
          }
         
      }
    }
   return 0;
 }


//int conditional_analysis(string file, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,string covariate, string grm_file,string outputFile,string geno_file,vector<vector<double> > &z_score, vector<vector<double> > &beta);



int create_bimbam(string file,string str)
 { 
   string str_phe=str+"_phe";    
   vector< vector<double> >  Data;
   vector<string> ind;
   vector<string> variant;
//   cout<<"Read input file"<<endl;
//   cout<<"file is: "<<file<<endl;
 //  cout<<"explored gene is: "<<gene<<endl;
   read_table(file,Data,ind,variant);
   fstream  FHOU(str.c_str(),ios::out);
   if(!FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }

   for(int i=1;i<variant.size();i++)
    {
      FHOU<<variant[i];
      FHOU<<",A,B";
      for(int j=0;j<Data.size();j++)
       {
//         Changed by Biao Zeng, 10/26/2020, 2:51 pm
//          Data[j][i]+=0.001*randn();
          FHOU<<","<<Data[j][i];
       }
      FHOU<<endl;
    }
   FHOU.clear();
   FHOU.close();
   

   fstream  FHOU1(str_phe.c_str(),ios::out);
   if(!FHOU1)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }

//   for(int i=0;i<variant.size();i++)
//    {
//      FHOU<<variant[i];
//      FHOU<<",A,B";
      for(int j=0;j<Data.size();j++)
       {

         FHOU1<<Data[j][0]<<endl;
       }

     FHOU1.clear();
     FHOU1.close();
  }


int conditional_analysis(string file, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,string covariate, string grm_file,string outputFile,string geno_file,vector<vector<double> > &z_score, vector<vector<double> > &beta, vector<vector<double> > &sd,int mode, vector<string> &Variant, vector<string> &removal1, int I, int colocalization, double &Variance)
  {
//   cout<<"Come into conditional_analysis function to perform real QTL detection"<<endl;
   vector< vector<double> >  Data;
   vector<string> ind;
   vector<string> variant;
   mode=1;
   vector<string> removal=removal1;

   cout<<"file is: "<<file<<", and the corresponding grm file is: "<<grm_file<<endl;

   if(removal.size()>0)
    {
//   cout<<"Variants detected in previous step"<<endl;
   for(int indx=0;indx<removal.size();indx++)
    {
      cout<<removal[indx]<<endl;
    }
    }
//   cout<<"Read input file"<<endl;
//   cout<<"file is: "<<file<<endl;
//   cout<<"explored gene is: "<<gene<<endl;
   read_table(file,Data,ind,variant);
//   cout<<"1: Read genotype and phenotype is over"<<endl;
   Variant.assign(variant.begin(), variant.end());
//   cout<<"Store variant infor"<<endl;
//   cout<<"Reading input file is finished"<<endl;
   mat data=mat(Data.size(),Data[0].size(),fill::zeros);
   string out_test="";
   string fam_file_test="";
//   cout<<"size of Data is: "<<Data.size()<<" "<<Data[0].size()<<endl;
   for(int i=0;i<Data.size();i++)
    {
      for(int j=0;j<Data[i].size();j++)
       {
         data(i,j)=Data[i][j];
       }
    }
 
   for(int i=0;i<data.n_cols;i++)
     {
       double Mean=0;
       int num=0;
       for(int j=0;j<data.n_rows;j++)
        {
          if(data(j,i)!=-9)
           {
             Mean+=data(j,i);
             num++;
           }
        }
       Mean/=num;
       for(int j=0;j<data.n_rows;j++)
        {
          if(data(j,i)==-9)
           {
             data(j,i)=Mean;
           }
        }
     }


   Variance = var(data.col(0));
   
//   cout<<"Convert input vector into mat finished"<<endl;
   int len=data.n_cols;
   int dim=data.n_rows;
   string out_text="";
   vector<string> col;
   for(int i=0;i<variant.size();i++)
    {
      col.push_back(variant[i]);
    }
   int j=0;
   for(int fea=0;fea<col.size();fea++)
    {
      if (regex_match (col[fea],regex("(IL)(.*)") ) || regex_match (col[fea],regex("(phe)(.*)") ) )
       {
         j=j+1;
       }
    }
   double cutoff=1.0;
//   cout<<"It is safe"<<endl;
//   cout<<"Output individuals"<<endl;

//   cout<<"j is: "<<j<<endl;
   for(int t=0;t<j;t++)
    {
//      cout<<"Step1"<<endl;
//      cout<<"I am a"<<endl;
      int x=0;
//      int bi=removal.size();
      int bi=1;
      vector<string> left;
      left=removal;
      int num_rev=0;
//      vector<string> removal;
      int target=0;
      double p_target;
      mat data1=data;
      vector<string> variant1;
      for(int i=0;i<variant.size();i++)
     {
         variant1.push_back(variant[i]);
       }
      string bed_file=string("probe")+to_string(t);
     if(geno_file!="NA")
      {
//      cout<<"Copy genotype to each probe"<<endl;
//      cout<<"Copy bed file"<<endl;
      copybed(geno_file+".bed",gene+"/"+bed_file+".bed");
//      cout<<"Copy fam file"<<endl;
      copyfam(geno_file+".fam",gene+"/"+bed_file+".fam");
//      cout<<"Copy bim file"<<endl;
      copybim(geno_file+".bim",gene+"/"+bed_file+".bim");

      string fam=gene+"/"+bed_file+".fam";
      vector<double>  pheno;
      vector < vector<string> > fam_infor;
      readfam(fam,pheno,fam_infor);
      mat fam_file=mat(pheno.size(),1,fill::zeros);
      for(int i=0;i<pheno.size();i++)
       {
         fam_file(i,0)=pheno[i];
       }
      int BIAO=fam_file.n_rows;
      BIAO=data.n_rows;
      fam_file.col(0)=data.col(t);
      string fam_file_out=fam+"_adjust";
//      cout<<"updatefam"<<endl;
      updatefam(fam_infor,fam_file,fam);
     }
//      cout<<"I am Ok here"<<endl;
      for(int k=j;k<len;k++)
       {
         vector<string> col1;
         for(int i=0;i<variant.size();i++)
          {
            col1.push_back(variant1[i]);
          }
         int len1=col1.size();
         data1.col(t)=data.col(t);
         target=0;
         p_target=1;
         int alle_target=0;
         string file_gemma_out=string("output_of_GEMMA_")+to_string(I)+"_"+gene;
         PARAM cPar;
         GEMMA cGemma;
//         cout<<"Inititation for GEMMA parameter and GEMMA core Over"<<endl;
         string str=gene+"/"+bed_file;
         if(geno_file!="NA")
          {
            cPar.file_bfile=str;
          } else
          {
            create_bimbam(file,str);
            cPar.file_geno=str;
            string str_phe=str+"_phe";
            cPar.file_pheno=str_phe;
          }
         int index_mixed=1;
         if(grm_file=="NA")
          {
            index_mixed=0;
          }
       //  string str;
         str.assign(grm_file);
         if(index_mixed==1)
          {
            cPar.file_kin=str;
          }
         str.assign(file_gemma_out);
         cPar.file_out=str;
//         str.clear();
         if(covariate!="NA")
          {
            str.clear();
            str.assign(covariate);
            cPar.file_cvt=str;          
            cPar.maf_level=-1;
          }
         if(grm_file!="NA")
          {
            cPar.a_mode=1;
          }else
          {  
            cPar.a_mode=51;
          }
//         cout<<"Perform parameter checking"<<endl;
         cPar.CheckParam();
//         cout<<"Parameter check is over"<<endl;
         cGemma.BatchRun(cPar);
//         cout<<"GEMMA running is over"<<endl;



         string file_gemma=string("./output")+"/output_of_GEMMA_"+to_string(I)+"_"+gene+".assoc.txt";
//         cout<<"Output of gemma is: "<<file_gemma<<endl;
         vector< vector<double> > Out_gemma;
         vector<string> Out_gemma_variant;
         readgemma(file_gemma,Out_gemma,Out_gemma_variant,z_score, beta,sd, index_mixed);
         if(Out_gemma_variant.size()<1)
          {
            return -1;
          }
         Variant=Out_gemma_variant;
         vector<string> Out_gemma_variant1=Out_gemma_variant;
//         cout<<"I am here11"<<endl;
         mat out_gemma=mat(Out_gemma.size(),Out_gemma[0].size(),fill::zeros);
//         cout<<"I am here22"<<endl;
         for(int i=0;i<Out_gemma.size();i++)
          {
            for(int j=0;j<Out_gemma[i].size();j++)
             {
               out_gemma(i,j)=Out_gemma[i][j];
             }
          }
//         cout<<"I am here33"<<endl;
         int target_gemma  =index_max(out_gemma.col(2));
//         cout<<"position for detected SNP is: "<<target_gemma<<endl;
//         cout<<"detected SNP in GEMMA is: "<<Out_gemma_variant[target_gemma]<<endl;
//         cout<<"Infor is: "<<out_gemma.row(target_gemma)<<endl;
         target=match(col1,Out_gemma_variant[target_gemma]);


         if(target<0)
          {
            return -1;
          }

         p_target=out_gemma(target_gemma,1);
         double beta_target =out_gemma(target_gemma,0);
//         cout<<"I am here4"<<endl;
//         cout<<"beta value is: "<<beta_target<<endl;
 //        for(int X=0;X<removal.size();X++)
 //         {
 //           cout<<"removal is: "<<removal[X]<<endl;
 //         }
  //       cout<<"p_target is: "<<p_target<<", target is: "<<target<<" and SNP is: "<<Out_gemma_variant[target_gemma]<<endl;
    //     cout<<"If there are some peak signals detected, remove"<<endl;
           
          if(removal.size()>0)
           { 
           }
 
         while(removal.size()>0 && match(removal,Out_gemma_variant[target_gemma])!=-1) 
          {
             out_gemma.shed_row(target_gemma);
             Out_gemma_variant.erase(Out_gemma_variant.begin()+target_gemma);
             col1.erase(col1.begin()+target);
             target_gemma  =index_max(out_gemma.col(2));
             target=match(col1,Out_gemma_variant[target_gemma]);
             p_target=out_gemma(target_gemma,1);
             if(p_target>cutoff)
              {
                break;
              }
             beta_target =out_gemma(target_gemma,0);
          }
      //   cout<<"p_target is:"<<p_target<<", target is: "<<target<<" and SNP is: "<<Out_gemma_variant[target_gemma]<<endl;
         if(p_target<=cutoff)
          {
        //     cout<<"bi is: "<<bi<<endl;
             if(bi>0)
              {
                int y=0;
                vector <string> left_test;
                for(int test=0;test<left.size();test++)
                 {
                   left_test.push_back(left[test]);
                 }
                vector<string> left_test_temp1=left_test;
                if(match(left_test,col1[target])<0)
                 {
                left_test_temp1.push_back(col1[target]);
                vector<string>::iterator z;
                for(int Z=0;Z<left_test_temp1.size()-1;Z++)
                 {
                   vector<string> col_test;
                   double vif=calculate_VIF(data,variant,*(left_test_temp1.begin()+Z),col_test);
                   if(vif>=2.5)
                    {
                      y=1;
                      break;
                    }
                  }
                 }
                
                if(y==1)
                 {
                   num_rev++;
//                   removal.push_back(col1[target]);
                   data1.shed_col(target);
                   removal.push_back(col1[target]);
                   col1.erase(col1.begin()+target);
                   out_gemma.shed_row(target_gemma);
                   Out_gemma_variant.erase(Out_gemma_variant.begin()+target_gemma);
                   int z=0;
                   while(y==1 && z==0  && col1.size()>10)
                    {
                      y=0;
//                      cout<<"sat1"<<endl;
                      target_gemma=index_max(out_gemma.col(2));
  //                    cout<<"sat2"<<endl;
                      int target=match(col1,Out_gemma_variant[target_gemma]);
                      while(target<0)
                       {
                         out_gemma.shed_row(target_gemma);
                         Out_gemma_variant.erase(Out_gemma_variant.begin()+target);
    //                     cout<<"sat3"<<endl;
                         target_gemma=index_max(out_gemma.col(2));
      //                   cout<<"sat4"<<endl;
                         target=match(col1,Out_gemma_variant[target_gemma]);
                       }
                      p_target=out_gemma(target_gemma,1);
                      beta_target=out_gemma(target_gemma,0);
//                      out_gemma.shed_row(target_gemma);
//                      Out_gemma_variant.erase(Out_gemma_variant.begin()+target_gemma);
//                      if(p_target>cutoff)
//                       {
//                         z=1;
//                         break;
//                       }
//
/*                       
                      if(p_target<=cutoff)
                       {
			    vector<string> left_test_temp=left_test;                            
			    left_test.push_back(col1[target]);
                            vector<string>::iterator z;
                            for(int Z=0;Z<left_test.size()-1;Z++)
                             {
                               vector<string> col_test;
                               vector<string> left_test_1;
 			       left_test_1=left_test;
                               left_test_1.erase(left_test_1.begin()+Z);
                               for(int i=0;i<left_test_1.size();i++)
                                {
                                  col_test.push_back(left_test_1[i]);
                                }
                               double vif=calculate_VIF(data,variant,(*(left_test.begin()+Z)),col_test);
                               if(vif>=5)
                                {
                                  y=1;
                                }
                             }

                            if(y==1)
                             {
                               num_rev=num_rev+1;
                               removal.push_back(col1[target]);
			       data1.shed_col(target);
                               col1.erase(col1.begin()+target);
                             }
                       }
*/
        //             cout<<"Why am I here1"<<endl;
 		      if(p_target<=cutoff)
                       {
                            vector<string> left_test_temp=left_test;
                            left_test_temp.push_back(col1[target]);
                            vector<string>::iterator z;
                          if(match(left_test,col1[target])<0)
                           {
                            for(int Z=0;Z<left_test_temp.size()-1;Z++)
                             {
                               vector<string> col_test;
                               col_test.push_back(col1[target]);
                               double vif=calculate_VIF(data,variant,*(left_test_temp.begin()+Z),col_test);
                               if(vif>=2.5)
                                {
                                  y=1;
                                  break;
                                }
                             }
                            } else
                            {
                              y=1;
                            }
                            if(y==1)
                             {
                               num_rev=num_rev+1;
                               removal.push_back(col1[target]);
                               data1.shed_col(target);
                               col1.erase(col1.begin()+target);
                               out_gemma.shed_row(target_gemma);
                               Out_gemma_variant.erase(Out_gemma_variant.begin()+target_gemma);
                             }
                       }
          //            cout<<"Why am I here2"<<endl;
                    }
  

                    
               
                  if(z==1)
                   {
                     continue;
                   }
                  }
                int x=1;              
                if(x==1)
                 {
                  num_rev=num_rev+1;
            //      cout<<"Why am I here3"<<endl;

                  map<string,int> VAR;
                  for(int BIAO=0;BIAO<removal.size();BIAO++)
                   {
                      VAR.insert(pair<string, int>(removal[BIAO],1));
                   }
                  for(int BIAO=0;BIAO<Out_gemma_variant1.size();BIAO++)
                   {
//                      VAR.find(Out_gemma_variant[BIAO]);
                      if(VAR.find(Out_gemma_variant1[BIAO]) != VAR.end())
                       {
                         beta[0][BIAO]=0.0000001;
                         if(sd[0][BIAO]<=0.00001)
                          {
                            z_score[0][BIAO]   = 0.001;
                          } else
                          {
                             z_score[0][BIAO]   =0.0000001/sd[0][BIAO];
                          }
                       }
                   }
//                  cout<<"I am here now"<<endl;
  //                removal.push_back(col1[target]);
                  if(mode==1)

                   {
                      cout<<"eQTL detection is over, and return to parent function"<<endl;
                      return 0;
                   }
                  ostringstream streamObj;
                  streamObj << p_target;
                  string str_p_target = streamObj.str();
                  string res=col1[t]+" "+col1[target]+" "+str_p_target+" "+to_string(beta_target);
                  out_test=out_test+"\n"+res;
                  len--;
                  if(len<=0)
                   {
                     
                     break;
                   }                  
                }
              } else
              {
  //                cout<<"Why am I here4"<<endl;
                  num_rev=num_rev+1;
                  ostringstream streamObj;
                  streamObj << p_target;
                  string str_p_target = streamObj.str();
                  
                  string res=col1[t]+" "+col1[target]+" "+str_p_target+" "+to_string(beta_target);
                  out_test=res;
                  bi++;
                  len--;
                  if(len<=0)
                   {
                     break;
                   }
                  if(mode==1)
                   {
                     return 0;
                   }
              }                     
           }
        }           
//      cout<<"conditional analysis result is: "<<out_test<<endl;
//      fstream  FHOU(outputFile.c_str(),ios::out);
//      if(!FHOU)
//       {
//         cout<<"Error for input or output"<<endl;
//         return (1);
//       }
//      FHOU<<out_test<<endl;
//      FHOU.clear();
//      FHOU.close();
      
   }
//   cout<<"Single eQTL detection is over"<<endl;
   return (0); 
  }


//Perform conditional analysis to detect peak signal
int conditional_analysis(string file, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,string covariate, vector <string> &grm_file,string outputFile,string geno_file)
 {

   vector< vector<double> >  Data;
   vector<string> ind;
   vector<string> variant;
//   cout<<"Read input file"<<endl;
//   cout<<"file is: "<<file<<endl;
//   cout<<"explored gene is: "<<gene<<endl;
   read_table(file,Data,ind,variant);
//   cout<<"Reading input file is finished"<<endl;
   mat data=mat(Data.size(),Data[0].size(),fill::zeros);
   string out_test="";
   string fam_file_test="";
   for(int i=0;i<Data.size();i++)
    {
      for(int j=0;j<Data[i].size();j++)
       {
         data(i,j)=Data[i][j];
       }
    }
//   cout<<"Convert input vector into mat finished"<<endl;
   int len=data.n_cols;
   int dim=data.n_rows;
   string out_text="";
   vector<string> col;
   for(int i=0;i<variant.size();i++)
    {
      col.push_back(variant[i]);
    }
   int j=0;
   for(int fea=0;fea<col.size();fea++)
    {
      if (regex_match (col[fea],regex("(IL)(.*)") ) || regex_match (col[fea],regex("(phe)(.*)") ) )
       {
         j=j+1;
       }
    }
   double cutoff=0.00001;
//   cout<<"It is safe"<<endl;
//   cout<<"Output individuals"<<endl;

//   cout<<"j is: "<<j<<endl;
   for(int t=0;t<j;t++)
    {
  //    cout<<"Step1"<<endl;
    //  cout<<"I am a"<<endl;
      int x=0;
      int bi=0;
      vector<string> left;
      int num_rev=0;
      vector<string> removal;
      int target=0;
      double p_target;
      mat data1=data;
      vector<string> variant1;
      for(int i=0;i<variant.size();i++)
       {
         variant1.push_back(variant[i]);
       }
      string bed_file=string("probe")+to_string(t);

//      cout<<"Copy genotype to each probe"<<endl;
//      cout<<"Copy bed file"<<endl;
      copybed(geno_file+".bed",gene+"/"+bed_file+".bed");
//      cout<<"Copy fam file"<<endl;
      copyfam(geno_file+".fam",gene+"/"+bed_file+".fam");
//      cout<<"Copy bim file"<<endl;
      copybim(geno_file+".bim",gene+"/"+bed_file+".bim");

      string fam=gene+"/"+bed_file+".fam";
      vector<double>  pheno;
      vector < vector<string> > fam_infor;
      readfam(fam,pheno,fam_infor);
      mat fam_file=mat(pheno.size(),1,fill::zeros);
      for(int i=0;i<pheno.size();i++)
       {
         fam_file(i,0)=pheno[i];
       }
      int BIAO=fam_file.n_rows;
      BIAO=data.n_rows;
      fam_file.col(0)=data.col(t);
      string fam_file_out=fam+"_adjust";
//      cout<<"updatefam"<<endl;
      updatefam(fam_infor,fam_file,fam);
//      cout<<"I am Ok here"<<endl;
      for(int k=j;k<len;k++)
       {         
         vector<string> col1;
         for(int i=0;i<variant.size();i++)
          {
            col1.push_back(variant1[i]);
          }
         int len1=col1.size();
         data1.col(t)=data.col(t);
         target=0;
         p_target=1;
         int alle_target=0;
         string file_gemma_out=string("output_of_GEMMA_")+gene;
         PARAM cPar;
         GEMMA cGemma;
//         cout<<"Inititation for GEMMA parameter and GEMMA core Over"<<endl;
         string str=gene+"/"+bed_file;
         cPar.file_bfile=str;

         str.assign(grm_file[0]);
         cPar.file_kin=str;
         str.assign(file_gemma_out);
         cPar.file_out=str;
         cPar.a_mode=1;
//         cout<<"Perform parameter checking"<<endl;
         cPar.CheckParam();
//         cout<<"Parameter check is over"<<endl;
         cGemma.BatchRun(cPar);
//         cout<<"GEMMA running is over"<<endl;



         string file_gemma=string("./output")+"/output_of_GEMMA_"+gene+".assoc.txt";
//         cout<<"Output of gemma is: "<<file_gemma<<endl;
         vector< vector<double> > Out_gemma;
         vector<string> Out_gemma_variant;
         readgemma(file_gemma,Out_gemma,Out_gemma_variant);
         mat out_gemma=mat(Out_gemma.size(),Out_gemma[0].size(),fill::zeros);
         for(int i=0;i<Out_gemma.size();i++)
          {
            for(int j=0;j<Out_gemma[i].size();j++)
             {
               out_gemma(i,j)=Out_gemma[i][j];
             }
          }

         int target_gemma  =index_max(out_gemma.col(2));
//         cout<<"position for detected SNP is: "<<target_gemma<<endl;
//         cout<<"detected SNP in GEMMA is: "<<Out_gemma_variant[target_gemma]<<endl;
         target=match(col1,Out_gemma_variant[target_gemma]);
         p_target=out_gemma(target_gemma,1);
         double beta_target =out_gemma(target_gemma,0);
//         cout<<"beta value is: "<<beta_target<<endl;
         for(int X=0;X<removal.size();X++)
          {
//            cout<<"removal is: "<<removal[X]<<endl;
          }
//         cout<<"p_target is: "<<p_target<<", target is: "<<target<<" and SNP is: "<<Out_gemma_variant[target_gemma]<<endl;
//         cout<<"If there are some peak signals detected, remove"<<endl;
         while(removal.size()>1&& match(removal,Out_gemma_variant[target_gemma])!=-1)
          {
//             cout<<"Please skip this SNP "<<Out_gemma_variant[target_gemma]<<endl;
             out_gemma.shed_row(target_gemma);
             Out_gemma_variant.erase(Out_gemma_variant.begin()+target_gemma);
//             cout<<"Try to find the index of minimum p value"<<endl;
             target_gemma  =index_max(out_gemma.col(2));
             target=match(col1,Out_gemma_variant[target_gemma]);
             p_target=out_gemma(target_gemma,1);
             if(p_target>cutoff)
              {
                break;
              }
             beta_target =out_gemma(target_gemma,0);
          }
//         cout<<"p_target is:"<<p_target<<", target is: "<<target<<" and SNP is: "<<Out_gemma_variant[target_gemma]<<endl;
         if(p_target>cutoff)
          {
//            cout<<"No significant variant found"<<endl;
            break ;
          }
         if(p_target<=cutoff)
          {
  //           cout<<"bi is: "<<bi<<endl;
             if(bi>0)
              {
    //            cout<<"Some peak signals are already detected"<<endl;
                int y=0;
                vector <string> left_test;
                for(int test=0;test<left.size();test++)
                 {
                   left_test.push_back(left[test]);
                 }
                left_test.push_back(col1[target]);
                vector<string>::iterator z;
//                for(z=left_test.begin();z!=left_test.end();z++)
                for(int Z=left_test.size();Z<left_test.size();Z++)
                 {
                   vector<string> col_test;
                   vector<string> left_test_1;
                   left_test_1=left_test;
                   left_test_1.erase(left_test_1.begin()+Z);
                   for(int i=0;i<left_test_1.size();i++)
                     {
                       col_test.push_back(left_test_1[i]);
                     }
//                   cout<<"Come to calculate VIF"<<endl;
                   double vif=calculate_VIF(data,variant,(*(left_test.begin()+Z)),col_test);
//                   cout<<"VIF value is: "<<vif<<endl;
                   if(vif>=10)
                    {
                      y=1;
                    }
                 }
              
          
                if(y==1)
                 {
//                   cout<<"Sorry, skip this SNP, as it locates in high LD with detected variants"<<endl;
                   num_rev++;
                   removal.push_back(col1[target]);
                   data1.shed_row(target);
                   col1.erase(col1.begin()+target);
                   int z=0;
                   while(y==1)
                    {
                      y=0;
                      int test=index_min(out_gemma.col(1));
                      out_gemma.shed_row(test);
                      int target=match(col1,Out_gemma_variant[target_gemma]);
                      while(target<0)
                       {
                         out_gemma.shed_row(target_gemma);
                         target_gemma=index_min(out_gemma.col(1));
                         target=match(col1,Out_gemma_variant[target_gemma]);
//                         cout<<"Please skip this SNP, as it is in the list of removal"<<endl;
                       }
                      p_target=out_gemma(target_gemma,1);
                      beta_target=out_gemma(target_gemma,0);
                      if(p_target>cutoff)
                       {
                         z=1;
                         break;
                       }
                      if(p_target<=cutoff)
                       {
                            
                            y=0;
                            vector<string>::iterator z;
//                            for(z=left_test.begin();z!=left_test.end();z++)
 			    for(int Z=0;Z<left_test.size();Z++)
                             {
                               vector<string> col_test;
                               vector<string> left_test_1;
                               left_test_1=left_test;
                               left_test_1.erase(left_test.begin()+Z);
                               for(int i=0;i<left_test_1.size();i++)
                                {
                                  col_test.push_back(left_test_1[i]);
                                }
                               double vif=calculate_VIF(data,variant,(*(left_test.begin()+Z)),col_test);
//                               cout<<"VIF value is: "<<vif<<endl;
                               if(vif>=10)
                                {
                                  y=1;
                                }
                             }

                            if(y==1)
                             {
//                               cout<<"Sorry, skip this SNP, as it locates in high LD with detected variants"<<endl;
                               num_rev=num_rev+1;
                               removal.push_back(col1[target]);
                               data1.shed_row(target);
                               col1.erase(col1.begin()+target);
                             }
                       }
                    }
                  if(z==1)
                   {
                     continue;
                   }
                  }
                int x=1;              
                if(x==1)
                 {
//                  cout<<"Conglatulations, another e-signal detected"<<endl;
                  num_rev=num_rev+1;
                  removal.push_back(col1[target]);
                 
                  ostringstream streamObj;
                  streamObj << p_target;
                  string str_p_target = streamObj.str();

                  string res=col1[t]+" "+col1[target]+" "+str_p_target+" "+to_string(beta_target);
                  extractvariant(gene+"/"+bed_file+".bed", gene+"/"+bed_file+"_high_LD.bed", 0.3, col1[target]);
                  string yFile=gene+"/"+bed_file+"_high_LD.phe";
                  string outputFileName=gene+"/"+col1[t]+"_"+to_string(bi)+"_high_LD.PolyQTL.output";
                  string X_file=gene+"/"+bed_file+"_high_LD.geno";
                  prepare_geno_phe(gene+"/"+bed_file+"_high_LD",X_file,yFile);
                  MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
                  MeQTLPoly.run();
                  MeQTLPoly.finishUp();


                  out_test=out_test+"\n"+res;
//                  cout<<"Get residuals and save it to fam file"<<endl;
//                  cout<<"beta value is: "<<beta_target<<endl;
                  fam_file.col(0)=fam_file.col(0)- beta_target*data1.col(target);
                  updatefam(fam_infor,fam_file,fam);
                  string fam_file_out=fam+"_adjust";
                  fam_file_test=fam+"_"+to_string(num_rev);
                  data1.shed_col(target);
                  bi++;
                  left.push_back(col1[target]);
                  col1.erase(col1.begin()+target);

                  len--;
                  if(len<=0)
                   {
                     break;
                   }                  
                }
//              cout<<"Next running"<<endl;
              } else
              {
                  num_rev=num_rev+1;
                  removal.push_back(col1[target]);
                  ostringstream streamObj;
                  streamObj << p_target;
                  string str_p_target = streamObj.str();
                  
                  string res=col1[t]+" "+col1[target]+" "+str_p_target+" "+to_string(beta_target);
                  extractvariant(gene+"/"+bed_file+".bed", gene+"/"+bed_file+"_high_LD.bed", 0.3, col1[target]);
                  string yFile=gene+"/"+bed_file+"_high_LD.phe";
                  string outputFileName=gene+"/"+col1[t]+"_"+to_string(bi)+"_high_LD.PolyQTL.output";
                  string X_file=gene+"/"+bed_file+"_high_LD.geno";
                  prepare_geno_phe(gene+"/"+bed_file+"_high_LD",X_file,yFile);
                  MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
                  MeQTLPoly.run();
                  MeQTLPoly.finishUp(); 
                  out_test=res;
//                  cout<<"Get residuals and save it to fam file"<<endl;
//                  cout<<"beta value is: "<<beta_target<<endl;

                  fam_file.col(0)=fam_file.col(0)- beta_target*data1.col(target);
                  updatefam(fam_infor,fam_file,fam);
                  string fam_file_out=fam+"_adjust";
                  fam_file_test=fam+"_"+to_string(num_rev);

                  data1.shed_col(target);
                  bi++;
                  left.push_back(col1[target]);
                  col1.erase(col1.begin()+target);
                  len--;
                  if(len<=0)
                   {
                     break;
                   }
              }                     
           }
        }           
//      cout<<"conditional analysis result is: "<<out_test<<endl;
      fstream  FHOU(outputFile.c_str(),ios::out);
      if(!FHOU)
       {
         cout<<"Error for input or output"<<endl;
         return (1);
       }
      FHOU<<out_test<<endl;
      FHOU.clear();
      FHOU.close();
      
   }
 }
double Han_adjust(int nstudy, mat &Cov1, double Z_score)
 {
  vec Han_data1={
0.9999999973,0.9999999944,0.9999999937,0.9999999922,0.9999999910,0.9999999911,0.9999999887,0.9999999898,0.9999999873,0.9999999878,0.9999999891,0.9999999881,0.9999999883,0.9999999886,0.9999999871,0.9999999859,0.9999999857,0.9999999840,0.9999999831,0.9999999836,0.9999999831,0.9999999849,0.9999999840,0.9999999818,0.9999999814,0.9999999822,0.9999999820,0.9999999799,0.9999999806,0.9999999825,0.9999999825,0.9999999817,0.9999999810,0.9999999770,0.9999999764,0.9999999771,0.9999999772,0.9999999761,0.9999999759,0.9999999738,0.9999999719,0.9999999765,0.9999999758,0.9999999756,0.9999999752,0.9999999750,0.9999999729,0.9999999716,0.9999999712,0.7794611342,0.7919015946,0.7994158204,0.8046048700,0.8084757110,0.8114953534,0.8139457442,0.8159911789,0.8177330453,0.8192375068,0.8205530101,0.8217130004,0.8227544457,0.8236859931,0.8245331997,0.8253093756,0.8260166552,0.8266725797,0.8272880358,0.8278437829,0.8283780744,0.8288622848,0.8293222676,0.8297499896,0.8301627204,0.8305528443,0.8309184155,0.8312673376,0.8315918620,0.8319084978,0.8322135471,0.8324946160,0.8327742759,0.8330307618,0.8332880068,0.8335259710,0.8337539369,0.8339753933,0.8341926987,0.8343988438,0.8346015373,0.8347984036,0.8349755665,0.8351529936,0.8353354597,0.8355066567,0.8356739711,0.8358266487,0.8359783251,0.6876640497,0.7028973271,0.7122335635,0.7187245494,0.7236039761,0.7274332778,0.7305612090,0.7331662804,0.7353980415,0.7373219423,0.7390173309,0.7405123302,0.7418533836,0.7430577297,0.7441521408,0.7451561296,0.7460753906,0.7469285156,0.7477217490,0.7484499339,0.7491463421,0.7497716683,0.7503733659,0.7509293191,0.7514701010,0.7519800017,0.7524606401,0.7529111682,0.7533368164,0.7537462473,0.7541475559,0.7545148359,0.7548811021,0.7552183582,0.7555522011,0.7558646513,0.7561709866,0.7564569970,0.7567421375,0.7570215454,0.7572787527,0.7575348074,0.7577781247,0.7580013332,0.7582446281,0.7584680400,0.7586859840,0.7588897512,0.7590871438,0.6189762512,0.6355195263,0.6457593210,0.6529263745,0.6583393158,0.6626078934,0.6661043542,0.6690167837,0.6715187430,0.6736854005,0.6755986863,0.6772783638,0.6787873316,0.6801533540,0.6813886350,0.6825270422,0.6835727278,0.6845325700,0.6854307189,0.6862591911,0.6870566013,0.6877618844,0.6884474876,0.6890804224,0.6896969405,0.6902706680,0.6908200887,0.6913349832,0.6918202632,0.6922842904,0.6927456494,0.6931622577,0.6935804807,0.6939634278,0.6943483958,0.6947034756,0.6950500159,0.6953835850,0.6957108105,0.6960222585,0.6963204241,0.6966161118,0.6968914574,0.6971518640,0.6974258801,0.6976816309,0.6979281887,0.6981637181,0.6983918877,0.5629047996,0.5800548683,0.5907510628,0.5982828152,0.6039958940,0.6085115296,0.6122155175,0.6153081092,0.6179765316,0.6202858680,0.6223246585,0.6241201352,0.6257362237,0.6271966610,0.6285218734,0.6297438947,0.6308589692,0.6318934371,0.6328626164,0.6337514509,0.6346028980,0.6353732717,0.6361044222,0.6367851730,0.6374457295,0.6380669645,0.6386600927,0.6392111022,0.6397386448,0.6402382804,0.6407390998,0.6411874461,0.6416385247,0.6420518566,0.6424673707,0.6428558814,0.6432266434,0.6435853824,0.6439370904,0.6442761083,0.6445966978,0.6449157587,0.6452130022,0.6454977554,0.6457936541,0.6460693492,0.6463420167,0.6465972415,0.6468403061,0.5152576063,0.5326012537,0.5434962001,0.5512059074,0.5570693124,0.5617162879,0.5655370332,0.5687338164,0.5714946590,0.5738874082,0.5760000168,0.5778640770,0.5795452095,0.5810628708,0.5824431418,0.5837178806,0.5848780049,0.5859567450,0.5869666871,0.5878931794,0.5887826313,0.5895861076,0.5903493143,0.5910638091,0.5917510439,0.5923998771,0.5930230117,0.5936003481,0.5941534274,0.5946757359,0.5951983521,0.5956690000,0.5961385162,0.5965741897,0.5970049425,0.5974157858,0.5978061127,0.5981822668,0.5985506533,0.5989044779,0.5992424470,0.5995734335,0.5998888513,0.6001872696,0.6005021764,0.6007879500,0.6010716974,0.6013408367,0.6015965468,0.4738079687,0.4910920011,0.5020120785,0.5097744484,0.5156953286,0.5203980059,0.5242739036,0.5275181217,0.5303255370,0.5327594909,0.5349148294,0.5368127587,0.5385322146,0.5400812018,0.5414903671,0.5427956516,0.5439795824,0.5450840353,0.5461177515,0.5470643928,0.5479759454,0.5488019982,0.5495850580,0.5503173734,0.5510229744,0.5516886242,0.5523305225,0.5529196582,0.5534875600,0.5540273359,0.5545596693,0.5550477009,0.5555276097,0.5559744210,0.5564224192,0.5568387273,0.5572426018,0.5576307306,0.5580115087,0.5583720809,0.5587192055,0.5590625508,0.5593878313,0.5596971997,0.5600167165,0.5603149479,0.5606040640,0.5608836024,0.5611454328,0.4371929299,0.4542472865,0.4650817286,0.4728124601,0.4787226105,0.4834244598,0.4873102349,0.4905659858,0.4933877692,0.4958331921,0.4980021564,0.4999170029,0.5016495982,0.5032096254,0.5046380017,0.5059524675,0.5071483509,0.5082596134,0.5093093195,0.5102654488,0.5111890866,0.5120194479,0.5128173446,0.5135598375,0.5142751861,0.5149508323,0.5155980184,0.5161965115,0.5167781269,0.5173151819,0.5178579499,0.5183549732,0.5188392856,0.5192967122,0.5197508392,0.5201735724,0.5205830217,0.5209742550,0.5213657298,0.5217296959,0.5220838560,0.5224326317,0.5227646208,0.5230778781,0.5234023691,0.5237031765,0.5239949264,0.5242818137,0.5245496558,0.4045012795,0.4212094099,0.4318791765,0.4395167064,0.4453659553,0.4500346121,0.4538928471,0.4571350476,0.4599434346,0.4623860442,0.4645438497,0.4664600101,0.4681858582,0.4697514317,0.4711772256,0.4724941748,0.4736944131,0.4748092282,0.4758565114,0.4768168446,0.4777403139,0.4785741308,0.4793735158,0.4801202374,0.4808367440,0.4815148030,0.4821650190,0.4827630615,0.4833515662,0.4838903697,0.4844333793,0.4849341353,0.4854189090,0.4858823035,0.4863358435,0.4867622121,0.4871749094,0.4875671013,0.4879616157,0.4883289769,0.4886838781,0.4890373862,0.4893663245,0.4896865806,0.4900098942,0.4903142116,0.4906063293,0.4908981064,0.4911689938,0.3750782719,0.3913603842,0.4018076248,0.4093087507,0.4150674687,0.4196670809,0.4234789476,0.4266877257,0.4294681049,0.4318847474,0.4340220868,0.4359220856,0.4376367744,0.4391924391,0.4406043333,0.4419141863,0.4431082241,0.4442226246,0.4452622314,0.4462162267,0.4471355343,0.4479644885,0.4487637205,0.4495072649,0.4502216888,0.4508950796,0.4515435373,0.4521376805,0.4527269254,0.4532674760,0.4538083105,0.4543070608,0.4547931301,0.4552533256,0.4557050199,0.4561310678,0.4565445880,0.4569375476,0.4573286807,0.4576992235,0.4580526759,0.4584026731,0.4587318807,0.4590535174,0.4593764190,0.4596791319,0.4599741386,0.4602632155,0.4605344380,0.3484267609,0.3642304014,0.3744191200,0.3817524232,0.3873942492,0.3919128834,0.3956537491,0.3988109471,0.4015497673,0.4039314690,0.4060351661,0.4079108261,0.4096042304,0.4111414077,0.4125353768,0.4138301922,0.4150138001,0.4161160882,0.4171430120,0.4180881638,0.4189965429,0.4198199900,0.4206122876,0.4213464040,0.4220565069,0.4227244101,0.4233666809,0.4239521542,0.4245398867,0.4250776637,0.4256124413,0.4261055908,0.4265888381,0.4270470365,0.4274987564,0.4279212508,0.4283295836,0.4287208646,0.4291063773,0.4294772922,0.4298270920,0.4301743424,0.4305016108,0.4308230038,0.4311415849,0.4314452535,0.4317389607,0.4320218619,0.4322927400,0.3241680507,0.3394631164,0.3493626705,0.3565087462,0.3620159306,0.3664315253,0.3700914747,0.3731828569,0.3758710657,0.3782079777,0.3802762660,0.3821188347,0.3837849832,0.3852953460,0.3866673061,0.3879389711,0.3891074873,0.3901921491,0.3912042291,0.3921321300,0.3930321771,0.3938419115,0.3946233353,0.3953520772,0.3960499662,0.3967070597,0.3973420437,0.3979203058,0.3984988002,0.3990342701,0.3995576251,0.4000477241,0.4005255583,0.4009798532,0.4014255044,0.4018407615,0.4022475606,0.4026303568,0.4030143343,0.4033786251,0.4037268609,0.4040703851,0.4043942218,0.4047137708,0.4050257802,0.4053323337,0.4056198459,0.4058972905,0.4061664330,0.3020010486,0.3167647736,0.3263583957,0.3332999842,0.3386581734,0.3429606551,0.3465307269,0.3495523007,0.3521795315,0.3544654334,0.3564885835,0.3582907267,0.3599247944,0.3614053127,0.3627525510,0.3639959948,0.3651441784,0.3662106774,0.3672034887,0.3681134875,0.3689988768,0.3697944060,0.3705632462,0.3712771904,0.3719677628,0.3726113233,0.3732336514,0.3738044416,0.3743753893,0.3749011926,0.3754136067,0.3758971420,0.3763688851,0.3768142879,0.3772558471,0.3776610787,0.3780633634,0.3784412182,0.3788195568,0.3791776858,0.3795210236,0.3798620724,0.3801783328,0.3804928141,0.3808004081,0.3811042955,0.3813865790,0.3816582434,0.3819254984,0.2816708635,0.2958974129,0.3051726990,0.3119011376,0.3170988249,0.3212808323,0.3247546101,0.3276994071,0.3302599130,0.3324893789,0.3344654764,0.3362245767,0.3378205949,0.3392697938,0.3405837754,0.3418019915,0.3429239086,0.3439680953,0.3449394837,0.3458311470,0.3466965958,0.3474751570,0.3482306485,0.3489302073,0.3496073272,0.3502363811,0.3508473961,0.3514116224,0.3519674418,0.3524852499,0.3529850367,0.3534617042,0.3539242278,0.3543609059,0.3547946008,0.3551911460,0.3555843938,0.3559604964,0.3563283620,0.3566823352,0.3570194637,0.3573513382,0.3576655401,0.3579761978,0.3582767557,0.3585747592,0.3588501768,0.3591198715,0.3593829551,0.2629797330,0.2766604550,0.2856123603,0.2921155372,0.2971536606,0.3012102997,0.3045809983,0.3074444135,0.3099333000,0.3121066098,0.3140264729,0.3157387153,0.3172954527,0.3187116726,0.3199922403,0.3211800978,0.3222772431,0.3232950295,0.3242426266,0.3251174517,0.3259601155,0.3267237819,0.3274575749,0.3281417896,0.3288038178,0.3294207955,0.3300169127,0.3305682996,0.3311123202,0.3316207367,0.3321086664,0.3325764659,0.3330268864,0.3334560519,0.3338799319,0.3342696299,0.3346554229,0.3350207736,0.3353832348,0.3357296956,0.3360598294,0.3363840869,0.3366946756,0.3369961681,0.3372906736,0.3375838467,0.3378534654,0.3381195688,0.3383752388,0.2457523505,0.2588916592,0.2675159500,0.2737926110,0.2786628396,0.2825902945,0.2858598502,0.2886341223,0.2910489196,0.2931619601,0.2950258452,0.2966938919,0.2982057849,0.2995820602,0.3008313132,0.3019882260,0.3030573414,0.3040453938,0.3049688549,0.3058239538,0.3066457782,0.3073873416,0.3081053496,0.3087721654,0.3094177409,0.3100215443,0.3106013976,0.3111394944,0.3116684152,0.3121678021,0.3126441312,0.3131011858,0.3135424852,0.3139610235,0.3143731511,0.3147554818,0.3151312367,0.3154885356,0.3158425282,0.3161800536,0.3165030185,0.3168206193,0.3171251714,0.3174198341,0.3177078360,0.3179938001,0.3182569238,0.3185155935,0.3187685339,0.2298404321,0.2424420872,0.2507393463,0.2567895999,0.2614912703,0.2652902186,0.2684511705,0.2711386687,0.2734788430,0.2755273122,0.2773343549,0.2789549871,0.2804212383,0.2817561097,0.2829713790,0.2840955368,0.2851350682,0.2860962815,0.2869946302,0.2878277762,0.2886267758,0.2893498147,0.2900469726,0.2906968457,0.2913251640,0.2919111497,0.2924786548,0.2930023572,0.2935174972,0.2940041956,0.2944669288,0.2949129354,0.2953417462,0.2957528556,0.2961528637,0.2965268128,0.2968911502,0.2972399473,0.2975868581,0.2979153051,0.2982304164,0.2985401818,0.2988357716,0.2991259362,0.2994040855,0.2996855830,0.2999398282,0.3001929053,0.3004413283,0.2151140052,0.2271883696,0.2351609426,0.2409906576,0.2455233097,0.2491873448,0.2522417497,0.2548431012,0.2571061163,0.2590896002,0.2608397235,0.2624081213,0.2638293137,0.2651236688,0.2663022263,0.2673949710,0.2684055252,0.2693371259,0.2702091208,0.2710179915,0.2717963649,0.2724948240,0.2731720960,0.2738059744,0.2744162124,0.2749883753,0.2755370971,0.2760488513,0.2765495517,0.2770220444,0.2774722817,0.2779055252,0.2783230309,0.2787243195,0.2791125863,0.2794768105,0.2798304076,0.2801735085,0.2805068813,0.2808305522,0.2811351222,0.2814356692,0.2817262954,0.2820090114,0.2822786742,0.2825543514,0.2827995752,0.2830463512,0.2832874714,0.2014660161,0.2130247893,0.2206784479,0.2262847705,0.2306496623,0.2341819291,0.2371313533,0.2396422593,0.2418290108,0.2437456223,0.2454410967,0.2469555080,0.2483323220,0.2495862286,0.2507298945,0.2517871805,0.2527661252,0.2536695027,0.2545145656,0.2553009462,0.2560545576,0.2567300260,0.2573895763,0.2580060292,0.2585952379,0.2591511765,0.2596855913,0.2601837470,0.2606700495,0.2611264171,0.2615634695,0.2619875385,0.2623913850,0.2627794108,0.2631571490,0.2635139074,0.2638547512,0.2641905559,0.2645145631,0.2648269915,0.2651229267,0.2654168394,0.2656979633,0.2659743332,0.2662376089,0.2665046264,0.2667454151,0.2669826887,0.2672168078,0.1887961087,0.1998537805,0.2071963508,0.2125839053,0.2167847210,0.2201839478,0.2230272184,0.2254478302,0.2275594316,0.2294088452,0.2310491291,0.2325126898,0.2338441470,0.2350564646,0.2361605659,0.2371868368,0.2381352820,0.2390063472,0.2398288758,0.2405865870,0.2413163881,0.2419735937,0.2426122943,0.2432090922,0.2437804110,0.2443204846,0.2448353518,0.2453214856,0.2457937259,0.2462360551,0.2466604245,0.2470723345,0.2474626588,0.2478392994,0.2482048079,0.2485530103,0.2488828388,0.2492090138,0.2495244778,0.2498271268,0.2501154750,0.2504005880,0.2506743349,0.2509417878,0.2511960502,0.2514532650,0.2516891489,0.2519193936,0.2521463535,0.1770219879,0.1875936684,0.1946322239,0.1998037221,0.2038413144,0.2071115890,0.2098472500,0.2121797278,0.2142160992,0.2160025785,0.2175842405,0.2189976626,0.2202830763,0.2214531950,0.2225240101,0.2235159177,0.2244306514,0.2252749162,0.2260710592,0.2268025095,0.2275071557,0.2281445013,0.2287652965,0.2293394120,0.2298936549,0.2304157795,0.2309141910,0.2313852235,0.2318439933,0.2322732601,0.2326847371,0.2330838747,0.2334605521,0.2338256881,0.2341804631,0.2345174683,0.2348395383,0.2351563321,0.2354619404,0.2357530234,0.2360348705,0.2363093304,0.2365760590,0.2368356082,0.2370819894,0.2373328957,0.2375599703,0.2377848892,0.2380033135,0.1660655820,0.1761667016,0.1829091340,0.1878688308,0.1917461217,0.1948897306,0.1975231727,0.1997689787,0.2017303217,0.2034526361,0.2049766043,0.2063429891,0.2075820775,0.2087119746,0.2097447836,0.2107032485,0.2115855967,0.2124027729,0.2131743019,0.2138790370,0.2145584599,0.2151741697,0.2157745423,0.2163331133,0.2168682202,0.2173721136,0.2178537354,0.2183111700,0.2187561449,0.2191687375,0.2195668716,0.2199549689,0.2203185577,0.2206727845,0.2210166802,0.2213422551,0.2216546933,0.2219601498,0.2222568435,0.2225410680,0.2228126809,0.2230781941,0.2233361000,0.2235893378,0.2238244510,0.2240689827,0.2242915264,0.2245078655,0.2247195200,0.1558615879,0.1655084018,0.1719613320,0.1767167451,0.1804384380,0.1834595700,0.1859885061,0.1881494984,0.1900385230,0.1916951716,0.1931630414,0.1944796135,0.1956758549,0.1967659055,0.1977630919,0.1986858967,0.1995374130,0.2003248578,0.2010705896,0.2017513544,0.2024086045,0.2030034847,0.2035819822,0.2041229395,0.2046394703,0.2051267800,0.2055918605,0.2060340446,0.2064630049,0.2068633342,0.2072475014,0.2076230503,0.2079757860,0.2083185553,0.2086518552,0.2089672530,0.2092685799,0.2095641835,0.2098513001,0.2101266564,0.2103905396,0.2106471188,0.2108968457,0.2111405774,0.2113689326,0.2116065083,0.2118206888,0.2120298114,0.2122349296,0.1463465633,0.1555567263,0.1617285663,0.1662854676,0.1698556549,0.1727553859,0.1751854906,0.1772619386,0.1790808936,0.1806733178,0.1820885148,0.1833548484,0.1845061523,0.1855586171,0.1865175962,0.1874078430,0.1882310407,0.1889878516,0.1897083229,0.1903657236,0.1909979913,0.1915730210,0.1921313369,0.1926538544,0.1931525773,0.1936246492,0.1940716668,0.1944985042,0.1949139397,0.1952999064,0.1956717068,0.1960351663,0.1963757990,0.1967060856,0.1970283476,0.1973328740,0.1976277571,0.1979115263,0.1981893075,0.1984553429,0.1987100152,0.1989585007,0.1991990707,0.1994358072,0.1996566527,0.1998878060,0.2000954904,0.2002961434,0.2004934516,0.1374687691,0.1462578407,0.1521587144,0.1565217119,0.1599450600,0.1627278992,0.1650604264,0.1670558472,0.1688036519,0.1703348291,0.1716964009,0.1729143628,0.1740226985,0.1750354641,0.1759594301,0.1768189898,0.1776096477,0.1783402798,0.1790329034,0.1796682928,0.1802778358,0.1808326165,0.1813715082,0.1818753337,0.1823567011,0.1828125538,0.1832440062,0.1836562335,0.1840568116,0.1844307960,0.1847903178,0.1851378974,0.1854674955,0.1857877499,0.1860964904,0.1863915386,0.1866763667,0.1869518685,0.1872192430,0.1874762527,0.1877224187,0.1879632460,0.1881949459,0.1884231760,0.1886358930,0.1888608932,0.1890615160,0.1892557974,0.1894469410,0.1291800634,0.1375623137,0.1432021170,0.1473805905,0.1506558789,0.1533272859,0.1555649728,0.1574796682,0.1591591315,0.1606297417,0.1619381157,0.1631097138,0.1641773688,0.1651515950,0.1660397920,0.1668677862,0.1676303820,0.1683338363,0.1690027336,0.1696122884,0.1701992515,0.1707347771,0.1712558491,0.1717422556,0.1722044562,0.1726410611,0.1730608384,0.1734559316,0.1738443089,0.1742040318,0.1745525184,0.1748855436,0.1752052768,0.1755138504,0.1758129167,0.1760953147,0.1763719717,0.1766379590,0.1768953705,0.1771421895,0.1773819509,0.1776149158,0.1778374670,0.1780579345,0.1782629897,0.1784795075,0.1786732895,0.1788626208,0.1790459897,0.1214348936,0.1294256631,0.1348145520,0.1388097441,0.1419479965,0.1445089886,0.1466525742,0.1484909917,0.1501020650,0.1515162655,0.1527726262,0.1538991246,0.1549258533,0.1558625940,0.1567187462,0.1575137902,0.1582466658,0.1589252788,0.1595678557,0.1601558778,0.1607221633,0.1612370953,0.1617393186,0.1622094783,0.1626535127,0.1630734551,0.1634770370,0.1638604289,0.1642329881,0.1645797664,0.1649142384,0.1652375333,0.1655463769,0.1658425829,0.1661327837,0.1664034228,0.1666702301,0.1669272117,0.1671746859,0.1674129871,0.1676445196,0.1678700798,0.1680847063,0.1682983508,0.1684961622,0.1687037269,0.1688917598,0.1690739677,0.1692527026,0.1141899901,0.1218071349,0.1269536087,0.1307739596,0.1337776403,0.1362322742,0.1382854651,0.1400496871,0.1415949709,0.1429512094,0.1441594356,0.1452412601,0.1462268598,0.1471283362,0.1479525551,0.1487167183,0.1494205403,0.1500753204,0.1506925297,0.1512582901,0.1518047030,0.1522987075,0.1527816931,0.1532337350,0.1536643337,0.1540682430,0.1544555228,0.1548255632,0.1551850569,0.1555189663,0.1558410009,0.1561533228,0.1564494393,0.1567357885,0.1570144830,0.1572753341,0.1575314290,0.1577808207,0.1580184310,0.1582494020,0.1584728611,0.1586893546,0.1588972683,0.1591031989,0.1592936633,0.1594929794,0.1596755223,0.1598504667,0.1600238494,0.1074122596,0.1146705023,0.1195827646,0.1232339242,0.1261091106,0.1284600223,0.1304279428,0.1321188394,0.1336016644,0.1349039524,0.1360627659,0.1371005317,0.1380478097,0.1389122866,0.1397051215,0.1404404743,0.1411161490,0.1417450038,0.1423397894,0.1428839963,0.1434086740,0.1438837895,0.1443491486,0.1447847390,0.1451987385,0.1455855603,0.1459585286,0.1463176884,0.1466625607,0.1469815337,0.1472950564,0.1475954816,0.1478800594,0.1481565903,0.1484247007,0.1486739541,0.1489223441,0.1491622761,0.1493911254,0.1496127665,0.1498289440,0.1500389244,0.1502381556,0.1504364711,0.1506202537,0.1508132356,0.1509877231,0.1511568799,0.1513248033,0.1010667666,0.1079797476,0.1126686730,0.1161583984,0.1189071083,0.1211574907,0.1230432462,0.1246640281,0.1260832685,0.1273340585,0.1284462398,0.1294426645,0.1303532343,0.1311803611,0.1319455313,0.1326479527,0.1332968414,0.1339031341,0.1344748488,0.1349979712,0.1355009259,0.1359588168,0.1364067991,0.1368247047,0.1372253617,0.1375963985,0.1379552461,0.1383010864,0.1386330333,0.1389396658,0.1392405898,0.1395294469,0.1398027913,0.1400710951,0.1403280692,0.1405677916,0.1408077267,0.1410387098,0.1412607192,0.1414736394,0.1416801868,0.1418838268,0.1420733224,0.1422666757,0.1424440476,0.1426282884,0.1427983730,0.1429619463,0.1431225976,0.0951231654,0.1017062551,0.1061802741,0.1095130816,0.1121409564,0.1142961706,0.1161008733,0.1176546597,0.1190128645,0.1202123557,0.1212789179,0.1222349032,0.1231086143,0.1239027740,0.1246378187,0.1253123792,0.1259360280,0.1265164799,0.1270667942,0.1275701638,0.1280521861,0.1284925356,0.1289216710,0.1293243974,0.1297102801,0.1300675267,0.1304130092,0.1307453998,0.1310649609,0.1313584899,0.1316493434,0.1319270947,0.1321906268,0.1324473284,0.1326947838,0.1329251278,0.1331577054,0.1333793509,0.1335935669,0.1337972948,0.1339983177,0.1341906789,0.1343736673,0.1345611057,0.1347315735,0.1349095640,0.1350721867,0.1352302119,0.1353848020,0.0895517420,0.0958220112,0.1000876059,0.1032695791,0.1057831642,0.1078439623,0.1095717784,0.1110606735,0.1123601226,0.1135099877,0.1145335236,0.1154484319,0.1162881147,0.1170499693,0.1177546664,0.1184026419,0.1190015967,0.1195580970,0.1200861827,0.1205681540,0.1210322721,0.1214567628,0.1218686071,0.1222555690,0.1226257745,0.1229702519,0.1233012758,0.1236219321,0.1239275871,0.1242109677,0.1244912113,0.1247588518,0.1250116343,0.1252584777,0.1254959671,0.1257181985,0.1259410132,0.1261532198,0.1263605410,0.1265585214,0.1267508554,0.1269348954,0.1271113832,0.1272919405,0.1274563965,0.1276268454,0.1277838824,0.1279347798,0.1280855593,0.0843271583,0.0902978744,0.0943646699,0.0974037646,0.0998051097,0.1017750206,0.1034292700,0.1048528959,0.1060996968,0.1071997606,0.1081814997,0.1090573790,0.1098646317,0.1105939285,0.1112695582,0.1118922281,0.1124662439,0.1129990330,0.1135069945,0.1139700495,0.1144149142,0.1148230976,0.1152194036,0.1155906720,0.1159470968,0.1162770810,0.1165940926,0.1169051797,0.1171966155,0.1174696327,0.1177406324,0.1179962523,0.1182381650,0.1184755476,0.1187052848,0.1189184273,0.1191316212,0.1193363390,0.1195370750,0.1197258283,0.1199115571,0.1200899596,0.1202586189,0.1204325921,0.1205906761,0.1207554797,0.1209058740,0.1210516824,0.1211948813,0.0794262471,0.0851108997,0.0889877537,0.0918876221,0.0941806779,0.0960654927,0.0976480004,0.0990095242,0.1002024135,0.1012583121,0.1021966404,0.1030375422,0.1038115361,0.1045093805,0.1051590324,0.1057562334,0.1063055015,0.1068170741,0.1073039327,0.1077489961,0.1081763044,0.1085675211,0.1089488502,0.1093054493,0.1096480404,0.1099643343,0.1102678038,0.1105655283,0.1108470578,0.1111091905,0.1113695965,0.1116167651,0.1118469888,0.1120763997,0.1122971835,0.1125017022,0.1127068600,0.1129020719,0.1130974975,0.1132786322,0.1134575328,0.1136275651,0.1137908034,0.1139565422,0.1141105343,0.1142670516,0.1144123920,0.1145524699,0.1146903260,0.0748282976,0.0802364492,0.0839329075,0.0867008332,0.0888912950,0.0906930365,0.0922055971,0.0935082511,0.0946493257,0.0956606099,0.0965586283,0.0973650095,0.0981054511,0.0987743280,0.0993978456,0.0999704736,0.1004959764,0.1009893168,0.1014550298,0.1018831270,0.1022913650,0.1026667052,0.1030338471,0.1033755513,0.1037043319,0.1040077494,0.1042991342,0.1045852505,0.1048561801,0.1051083033,0.1053572913,0.1055945824,0.1058157449,0.1060359636,0.1062484031,0.1064448206,0.1066421375,0.1068300907,0.1070171957,0.1071907578,0.1073635874,0.1075265944,0.1076843410,0.1078430328,0.1079910756,0.1081418380,0.1082807400,0.1084161036,0.1085481263,0.0705099876,0.0756574296,0.0791795987,0.0818216599,0.0839115582,0.0856333550,0.0870786925,0.0883268158,0.0894180970,0.0903855424,0.0912470468,0.0920175050,0.0927263230,0.0933676870,0.0939655585,0.0945141128,0.0950188806,0.0954902813,0.0959372100,0.0963493576,0.0967405402,0.0971007007,0.0974520092,0.0977793416,0.0980946382,0.0983860033,0.0986672051,0.0989406187,0.0992018411,0.0994432466,0.0996813689,0.0999111768,0.1001234831,0.1003342624,0.1005381460,0.1007265061,0.1009154665,0.1010980781,0.1012763140,0.1014426260,0.1016089157,0.1017664765,0.1019170339,0.1020698097,0.1022123979,0.1023582101,0.1024899036,0.1026219884,0.1027472322,0.0664565363,0.0713533212,0.0747091878,0.0772299952,0.0792249127,0.0808689276,0.0822501614,0.0834440316,0.0844875319,0.0854129078,0.0862378452,0.0869759060,0.0876552379,0.0882697077,0.0888423723,0.0893678337,0.0898523828,0.0903031628,0.0907310431,0.0911267299,0.0915020306,0.0918479647,0.0921828236,0.0924983333,0.0927993064,0.0930811934,0.0933509703,0.0936116977,0.0938624331,0.0940942744,0.0943236406,0.0945449462,0.0947468608,0.0949505051,0.0951454581,0.0953266261,0.0955078119,0.0956834079,0.0958556175,0.0960141117,0.0961737323,0.0963241251,0.0964696941,0.0966166396,0.0967541517,0.0968933476,0.0970190320,0.0971468594,0.0972670662,0.0626464832,0.0673050400,0.0705024175,0.0729064817,0.0748101720,0.0763812589,0.0777016519,0.0788421708,0.0798408092,0.0807257912,0.0815139652,0.0822202748,0.0828729318,0.0834606787,0.0840082628,0.0845112908,0.0849761363,0.0854078697,0.0858183531,0.0861980365,0.0865580404,0.0868884945,0.0872101182,0.0875125402,0.0878009586,0.0880724993,0.0883300509,0.0885810561,0.0888203565,0.0890419650,0.0892635767,0.0894751593,0.0896694138,0.0898652520,0.0900516222,0.0902258675,0.0903993780,0.0905685777,0.0907321514,0.0908857662,0.0910383469,0.0911826133,0.0913231497,0.0914635891,0.0915968807,0.0917310629,0.0918513860,0.0919747273,0.0920889687,0.0590654952,0.0634967662,0.0665432994,0.0688368287,0.0706529321,0.0721541873,0.0734142459,0.0745047830,0.0754586212,0.0763052272,0.0770598187,0.0777350207,0.0783604042,0.0789222465,0.0794471483,0.0799298759,0.0803740518,0.0807877704,0.0811812097,0.0815443325,0.0818909886,0.0822058901,0.0825144009,0.0828045879,0.0830809010,0.0833420298,0.0835881956,0.0838293869,0.0840588093,0.0842707229,0.0844835419,0.0846865891,0.0848723640,0.0850610321,0.0852397628,0.0854065619,0.0855753349,0.0857361571,0.0858932937,0.0860408179,0.0861860797,0.0863246415,0.0864604295,0.0865959099,0.0867224315,0.0868511309,0.0869670358,0.0870857626,0.0871961416,0.0557005062,0.0599136650,0.0628157868,0.0650028245,0.0667356040,0.0681686422,0.0693734852,0.0704139539,0.0713273651,0.0721357305,0.0728573928,0.0735035800,0.0741024315,0.0746397188,0.0751427557,0.0756056511,0.0760307405,0.0764268455,0.0768029956,0.0771507093,0.0774837930,0.0777852437,0.0780795567,0.0783572836,0.0786243808,0.0788734929,0.0791109652,0.0793414040,0.0795602691,0.0797651752,0.0799682195,0.0801633428,0.0803404522,0.0805207196,0.0806931020,0.0808528454,0.0810155025,0.0811688937,0.0813203680,0.0814608795,0.0816016452,0.0817347201,0.0818636928,0.0819948702,0.0821167416,0.0822395367,0.0823508728,0.0824641426,0.0825699377,0.0525356066,0.0565436057,0.0593063408,0.0613916665,0.0630445435,0.0644120370,0.0655627226,0.0665567803,0.0674283223,0.0682019818,0.0688939005,0.0695115023,0.0700830773,0.0705973895,0.0710791292,0.0715209752,0.0719294129,0.0723094488,0.0726685824,0.0730015535,0.0733202952,0.0736095565,0.0738909666,0.0741575088,0.0744142958,0.0746531104,0.0748800444,0.0751007761,0.0753102406,0.0755064170,0.0757020468,0.0758879990,0.0760563376,0.0762316408,0.0763954309,0.0765505317,0.0767058059,0.0768516180,0.0769980110,0.0771332018,0.0772689117,0.0773959660,0.0775192392,0.0776449119,0.0777623331,0.0778797572,0.0779855904,0.0780944788,0.0781964517,0.0495592721,0.0533703410,0.0560010254,0.0579881008,0.0595649208,0.0608697770,0.0619679220,0.0629188514,0.0637513935,0.0644902980,0.0651516087,0.0657420413,0.0662896938,0.0667815960,0.0672420352,0.0676647791,0.0680571136,0.0684201260,0.0687640439,0.0690840244,0.0693872878,0.0696659239,0.0699338289,0.0701895055,0.0704348227,0.0706642967,0.0708823528,0.0710942197,0.0712949243,0.0714820294,0.0716704772,0.0718482910,0.0720088467,0.0721780760,0.0723334487,0.0724831787,0.0726330176,0.0727737618,0.0729122362,0.0730415741,0.0731713651,0.0732936743,0.0734115919,0.0735322459,0.0736452737,0.0737584588,0.0738601597,0.0739637395,0.0740621426,0.0467588969,0.0503822986,0.0528865937,0.0547808073,0.0562838245,0.0575289260,0.0585770037,0.0594867015,0.0602807487,0.0609871486,0.0616179144,0.0621842415,0.0627071146,0.0631784024,0.0636186109,0.0640228738,0.0643985772,0.0647458533,0.0650744018,0.0653828643,0.0656721950,0.0659378917,0.0661949102,0.0664396382,0.0666746128,0.0668951622,0.0671030532,0.0673062454,0.0674989786,0.0676783399,0.0678591163,0.0680295919,0.0681848780,0.0683452929,0.0684956491,0.0686373454,0.0687813240,0.0689168995,0.0690501129,0.0691738446,0.0692972701,0.0694147027,0.0695291745,0.0696440715,0.0697524021,0.0698611654,0.0699585143,0.0700579381,0.0701513608,0.0441224966,0.0475672611,0.0499522896,0.0517569243,0.0531899166,0.0543785260,0.0553781582,0.0562475223,0.0570057930,0.0576809426,0.0582838203,0.0588252385,0.0593254684,0.0597752135,0.0601966261,0.0605827884,0.0609424076,0.0612741835,0.0615881722,0.0618841602,0.0621619136,0.0624154827,0.0626622464,0.0628966351,0.0631211051,0.0633317278,0.0635311133,0.0637271751,0.0639105660,0.0640823256,0.0642560525,0.0644198795,0.0645674863,0.0647218499,0.0648646446,0.0650020875,0.0651395591,0.0652688779,0.0653967157,0.0655154489,0.0656340156,0.0657471207,0.0658563400,0.0659656780,0.0660699001,0.0661735121,0.0662680597,0.0663632422,0.0664520331,0.0416401758,0.0449164216,0.0471863229,0.0489064867,0.0502729149,0.0514053830,0.0523598104,0.0531894812,0.0539145039,0.0545587988,0.0551346947,0.0556525900,0.0561310292,0.0565608339,0.0569639881,0.0573331480,0.0576766617,0.0579942961,0.0582949796,0.0585775813,0.0588436180,0.0590865041,0.0593234170,0.0595478349,0.0597632352,0.0599647060,0.0601548832,0.0603425352,0.0605187024,0.0606827000,0.0608490337,0.0610058919,0.0611486128,0.0612944587,0.0614328643,0.0615646818,0.0616951157,0.0618195551,0.0619416395,0.0620551120,0.0621697662,0.0622775535,0.0623810374,0.0624875871,0.0625880704,0.0626871269,0.0627778301,0.0628684748,0.0629535074,0.0393036187,0.0424185709,0.0445788370,0.0462169339,0.0475196320,0.0485994411,0.0495114082,0.0503031838,0.0509955569,0.0516117411,0.0521607786,0.0526557415,0.0531129383,0.0535232801,0.0539094910,0.0542619286,0.0545914564,0.0548949728,0.0551830135,0.0554528180,0.0557072742,0.0559405881,0.0561665428,0.0563813024,0.0565877838,0.0567798380,0.0569633484,0.0571424159,0.0573110266,0.0574681961,0.0576268604,0.0577775642,0.0579145679,0.0580531822,0.0581868213,0.0583127155,0.0584385968,0.0585567589,0.0586746802,0.0587839609,0.0588921954,0.0589960083,0.0590955594,0.0591967041,0.0592932108,0.0593886328,0.0594758492,0.0595632465,0.0596438308,0.0371031088,0.0400638203,0.0421201267,0.0436811575,0.0449228207,0.0459526880,0.0468228120,0.0475791176,0.0482388812,0.0488271842,0.0493513867,0.0498241516,0.0502620339,0.0506543880,0.0510230126,0.0513609127,0.0516751153,0.0519644617,0.0522417374,0.0524985570,0.0527421648,0.0529663775,0.0531831793,0.0533876089,0.0535846178,0.0537692231,0.0539455556,0.0541161110,0.0542778655,0.0544285654,0.0545811775,0.0547247761,0.0548558830,0.0549873325,0.0551170410,0.0552364248,0.0553565621,0.0554708850,0.0555837211,0.0556886549,0.0557917248,0.0558906635,0.0559870892,0.0560827241,0.0561757285,0.0562670760,0.0563511259,0.0564342993,0.0565117370,0.0350306850,0.0378438680,0.0398021829,0.0412886860,0.0424709063,0.0434536596,0.0442840256,0.0450054060,0.0456355169,0.0461974285,0.0466988580,0.0471495033,0.0475684058,0.0479420874,0.0482956112,0.0486186294,0.0489185226,0.0491962398,0.0494600904,0.0497057962,0.0499396919,0.0501527807,0.0503608070,0.0505569104,0.0507454445,0.0509216714,0.0510906972,0.0512529742,0.0514094611,0.0515532519,0.0516988048,0.0518366730,0.0519625418,0.0520882769,0.0522111294,0.0523267922,0.0524410042,0.0525506424,0.0526585095,0.0527598773,0.0528593919,0.0529535244,0.0530454935,0.0531367593,0.0532267215,0.0533136398,0.0533946648,0.0534730493,0.0535474995,0.0330786114,0.0357516122,0.0376156325,0.0390310623,0.0401567397,0.0410953465,0.0418865847,0.0425756095,0.0431758892,0.0437135307,0.0441914007,0.0446223181,0.0450216687,0.0453798497,0.0457177393,0.0460254313,0.0463132011,0.0465774139,0.0468301622,0.0470657448,0.0472882178,0.0474934966,0.0476923209,0.0478799072,0.0480607111,0.0482296038,0.0483902659,0.0485461448,0.0486945595,0.0488329949,0.0489723207,0.0491047204,0.0492252965,0.0493450757,0.0494619285,0.0495746830,0.0496822542,0.0497880215,0.0498908880,0.0499882095,0.0500840693,0.0501742632,0.0502616538,0.0503500604,0.0504353195,0.0505182076,0.0505967778,0.0506717103,0.0507428357,0.0312387547,0.0337793417,0.0355527820,0.0368991672,0.0379728992,0.0388677437,0.0396226353,0.0402795917,0.0408517136,0.0413665287,0.0418227876,0.0422343069,0.0426148958,0.0429582039,0.0432807253,0.0435743249,0.0438496245,0.0441020444,0.0443437655,0.0445694192,0.0447816452,0.0449780328,0.0451684026,0.0453479212,0.0455205219,0.0456819157,0.0458354130,0.0459856478,0.0461273654,0.0462595637,0.0463923000,0.0465196845,0.0466345055,0.0467492120,0.0468626917,0.0469690251,0.0470726213,0.0471749100,0.0472716129,0.0473659583,0.0474572597,0.0475440543,0.0476279372,0.0477111354,0.0477930853,0.0478726336,0.0479483825,0.0480191070,0.0480870533,0.0295047679,0.0319183576,0.0336066452,0.0348880545,0.0359110283,0.0367647914,0.0374838803,0.0381105947,0.0386566370,0.0391487876,0.0395833612,0.0399770599,0.0403402549,0.0406683246,0.0409760203,0.0412578315,0.0415207371,0.0417620343,0.0419924152,0.0422077540,0.0424115401,0.0425986525,0.0427807658,0.0429521660,0.0431173817,0.0432713490,0.0434188078,0.0435620711,0.0436982956,0.0438248882,0.0439510232,0.0440730952,0.0441835268,0.0442934203,0.0444022452,0.0445036915,0.0446023503,0.0447000995,0.0447924692,0.0448830573,0.0449714510,0.0450547774,0.0451341533,0.0452138506,0.0452927805,0.0453685649,0.0454414281,0.0455087769,0.0455751417,0.0278701526,0.0301630650,0.0317693327,0.0329892654,0.0339643455,0.0347784816,0.0354642124,0.0360614286,0.0365820029,0.0370526695,0.0374672654,0.0378433969,0.0381899939,0.0385030566,0.0387976502,0.0390663386,0.0393178351,0.0395485003,0.0397686798,0.0399744564,0.0401683313,0.0403482382,0.0405222057,0.0406861722,0.0408434313,0.0409907361,0.0411326367,0.0412688774,0.0413998499,0.0415210726,0.0416414676,0.0417585342,0.0418644467,0.0419695723,0.0420730116,0.0421697067,0.0422643223,0.0423587737,0.0424459800,0.0425337871,0.0426172997,0.0426977726,0.0427741614,0.0428497496,0.0429256507,0.0429978825,0.0430684747,0.0431324822,0.0431965818,0.0263281200,0.0285075182,0.0300356292,0.0311970833,0.0321259013,0.0329022611,0.0335557472,0.0341257577,0.0346217150,0.0350713523,0.0354670389,0.0358254630,0.0361574783,0.0364561960,0.0367373738,0.0369925590,0.0372349952,0.0374541781,0.0376649614,0.0378628009,0.0380468588,0.0382187781,0.0383862111,0.0385415476,0.0386933922,0.0388332026,0.0389687418,0.0390993782,0.0392249060,0.0393409257,0.0394563623,0.0395677297,0.0396696157,0.0397693478,0.0398688077,0.0399618949,0.0400520269,0.0401419251,0.0402259025,0.0403090138,0.0403896674,0.0404670463,0.0405398956,0.0406121350,0.0406846420,0.0407534473,0.0408213475,0.0408818648,0.0409434816,0.0248743883,0.0269455878,0.0283992049,0.0295053111,0.0303893744,0.0311294405,0.0317518535,0.0322959099,0.0327693523,0.0331985787,0.0335752527,0.0339179343,0.0342340845,0.0345209661,0.0347888885,0.0350323278,0.0352639227,0.0354730010,0.0356748490,0.0358628793,0.0360393975,0.0362041725,0.0363642034,0.0365128782,0.0366571825,0.0367914762,0.0369208178,0.0370464339,0.0371657679,0.0372775013,0.0373869435,0.0374938209,0.0375919634,0.0376868336,0.0377821903,0.0378717034,0.0379565058,0.0380431109,0.0381233564,0.0382022582,0.0382798237,0.0383541500,0.0384231781,0.0384931325,0.0385622742,0.0386290143,0.0386925285,0.0387509458,0.0388105709,0.0235029591,0.0254722567,0.0268542663,0.0279071383,0.0287495698,0.0294545456,0.0300478335,0.0305662615,0.0310176807,0.0314269938,0.0317868691,0.0321146723,0.0324159384,0.0326901974,0.0329460838,0.0331790209,0.0333988923,0.0335992210,0.0337921469,0.0339718876,0.0341407180,0.0342978762,0.0344505503,0.0345928058,0.0347314168,0.0348593154,0.0349827318,0.0351025418,0.0352169698,0.0353234976,0.0354287910,0.0355308986,0.0356249012,0.0357153904,0.0358061559,0.0358923457,0.0359724954,0.0360563106,0.0361329561,0.0362079886,0.0362825676,0.0363534196,0.0364202321,0.0364869883,0.0365531716,0.0366162933,0.0366769315,0.0367337285,0.0367900987,0.0222105785,0.0240805231,0.0253955672,0.0263977950,0.0272005834,0.0278716342,0.0284371666,0.0289320300,0.0293625133,0.0297528523,0.0300961685,0.0304094416,0.0306963030,0.0309590152,0.0312023487,0.0314247589,0.0316350592,0.0318261117,0.0320104291,0.0321821480,0.0323432075,0.0324942175,0.0326392181,0.0327758376,0.0329086533,0.0330304163,0.0331483999,0.0332627217,0.0333727303,0.0334745305,0.0335749013,0.0336725150,0.0337622437,0.0338492028,0.0339354996,0.0340186403,0.0340954390,0.0341745621,0.0342472843,0.0343196965,0.0343914279,0.0344592750,0.0345239426,0.0345870234,0.0346507677,0.0347104030,0.0347684700,0.0348234037,0.0348767322,0.0209905183,0.0227669408,0.0240181764,0.0249720465,0.0257360451,0.0263758244,0.0269147985,0.0273866762,0.0277972036,0.0281696867,0.0284977766,0.0287953934,0.0290695588,0.0293202990,0.0295538410,0.0297656418,0.0299661889,0.0301482760,0.0303237553,0.0304884871,0.0306425909,0.0307868426,0.0309250430,0.0310556755,0.0311829542,0.0312999147,0.0314121470,0.0315209131,0.0316265442,0.0317233321,0.0318191182,0.0319132352,0.0319989150,0.0320813220,0.0321648086,0.0322445515,0.0323180044,0.0323927917,0.0324624117,0.0325320976,0.0326004250,0.0326651798,0.0327267146,0.0327879258,0.0328482760,0.0329054743,0.0329611169,0.0330142068,0.0330639620,0.0198397176,0.0215269969,0.0227171383,0.0236250341,0.0243538064,0.0249615477,0.0254752050,0.0259255542,0.0263168106,0.0266723140,0.0269851661,0.0272689783,0.0275311448,0.0277700525,0.0279932037,0.0281952394,0.0283867371,0.0285605713,0.0287279202,0.0288868802,0.0290328487,0.0291703798,0.0293027608,0.0294281581,0.0295494472,0.0296607752,0.0297675775,0.0298717330,0.0299730511,0.0300653791,0.0301575818,0.0302471852,0.0303288232,0.0304083854,0.0304872905,0.0305644948,0.0306339863,0.0307061117,0.0307722311,0.0308385114,0.0309048023,0.0309661002,0.0310251515,0.0310835708,0.0311414403,0.0311957134,0.0312490212,0.0313001717,0.0313484192,0.0187537946,0.0203565235,0.0214881658,0.0223530880,0.0230463386,0.0236259481,0.0241147040,0.0245444859,0.0249170274,0.0252563652,0.0255539094,0.0258253626,0.0260761826,0.0263028166,0.0265164407,0.0267089449,0.0268923285,0.0270581285,0.0272179042,0.0273700495,0.0275090108,0.0276409372,0.0277670079,0.0278872591,0.0280025026,0.0281082395,0.0282121535,0.0283106931,0.0284071261,0.0284957148,0.0285841851,0.0286699767,0.0287474949,0.0288238443,0.0288987182,0.0289730216,0.0290390994,0.0291079280,0.0291722587,0.0292355149,0.0292983862,0.0293568935,0.0294140194,0.0294691664,0.0295248645,0.0295763388,0.0296275096,0.0296761695,0.0297228905,0.0177284393,0.0192514564,0.0203270614,0.0211502925,0.0218110359,0.0223628826,0.0228289428,0.0232376306,0.0235928563,0.0239170182,0.0242013807,0.0244594890,0.0246986725,0.0249148564,0.0251189272,0.0253032335,0.0254779512,0.0256362920,0.0257884542,0.0259340839,0.0260670992,0.0261922033,0.0263134902,0.0264287369,0.0265381689,0.0266388621,0.0267391859,0.0268324248,0.0269246770,0.0270088662,0.0270939711,0.0271759971,0.0272502821,0.0273232910,0.0273952595,0.0274649640,0.0275288990,0.0275950615,0.0276563539,0.0277165842,0.0277754764,0.0278329604,0.0278874280,0.0279395250,0.0279935033,0.0280430551,0.0280910069,0.0281375974,0.0281824256,0.0167608696,0.0182070012,0.0192302665,0.0200132474,0.0206426537,0.0211685691,0.0216123267,0.0220022601,0.0223411262,0.0226499335,0.0229208342,0.0231669630,0.0233955583,0.0236010481,0.0237969263,0.0239725148,0.0241389322,0.0242900418,0.0244354921,0.0245754248,0.0247015603,0.0248209611,0.0249368787,0.0250474387,0.0251522153,0.0252479147,0.0253437981,0.0254327022,0.0255208150,0.0256014109,0.0256830691,0.0257609597,0.0258320234,0.0259018055,0.0259709948,0.0260368315,0.0260983692,0.0261613028,0.0262200202,0.0262781944,0.0263344536,0.0263889913,0.0264409392,0.0264912053,0.0265428094,0.0265892404,0.0266359470,0.0266795373,0.0267228057,0.0158471070,0.0172210565,0.0181942449,0.0189397430,0.0195388013,0.0200389407,0.0204622908,0.0208339649,0.0211566546,0.0214514209,0.0217094529,0.0219443249,0.0221619359,0.0223579583,0.0225447979,0.0227127700,0.0228712949,0.0230166668,0.0231546384,0.0232884679,0.0234086105,0.0235235555,0.0236335906,0.0237388142,0.0238391394,0.0239312723,0.0240219228,0.0241066577,0.0241914178,0.0242679932,0.0243464481,0.0244215761,0.0244888912,0.0245553712,0.0246216470,0.0246844859,0.0247435959,0.0248038863,0.0248595062,0.0249150618,0.0249693654,0.0250208749,0.0250708656,0.0251191661,0.0251682489,0.0252128749,0.0252570524,0.0252982632,0.0253400548,0.0149847775,0.0162900867,0.0172146478,0.0179240483,0.0184945862,0.0189707707,0.0193737826,0.0197286893,0.0200360600,0.0203170304,0.0205630491,0.0207867819,0.0209940584,0.0211821364,0.0213602159,0.0215209914,0.0216714973,0.0218101088,0.0219420454,0.0220700371,0.0221849584,0.0222943338,0.0223994923,0.0225002962,0.0225963985,0.0226840160,0.0227704928,0.0228513702,0.0229322183,0.0230054859,0.0230802003,0.0231516762,0.0232164874,0.0232794542,0.0233431569,0.0234033323,0.0234599015,0.0235175238,0.0235708710,0.0236241126,0.0236758815,0.0237250472,0.0237729654,0.0238187504,0.0238659382,0.0239080794,0.0239502617,0.0239900672,0.0240298778,0.0141702989,0.0154102272,0.0162894486,0.0169637275,0.0175070281,0.0179603433,0.0183447042,0.0186829623,0.0189754898,0.0192438093,0.0194783063,0.0196917098,0.0198891184,0.0200691626,0.0202394522,0.0203918438,0.0205359721,0.0206684316,0.0207944189,0.0209158108,0.0210264064,0.0211301712,0.0212309390,0.0213277623,0.0214192979,0.0215024083,0.0215849726,0.0216621281,0.0217394483,0.0218096077,0.0218807206,0.0219489324,0.0220112083,0.0220709024,0.0221322842,0.0221897758,0.0222435496,0.0222988103,0.0223497548,0.0224003745,0.0224501297,0.0224970341,0.0225425088,0.0225868346,0.0226317878,0.0226721283,0.0227131249,0.0227496135,0.0227883058,0.0134012149,0.0145784688,0.0154146450,0.0160564558,0.0165736186,0.0170046832,0.0173715678,0.0176932632,0.0179723423,0.0182276876,0.0184517942,0.0186553631,0.0188435066,0.0190151532,0.0191778261,0.0193233669,0.0194601814,0.0195872233,0.0197076420,0.0198236045,0.0199288396,0.0200279863,0.0201244517,0.0202170700,0.0203042065,0.0203832846,0.0204625405,0.0205360023,0.0206095481,0.0206766542,0.0207446206,0.0208100471,0.0208693135,0.0209261004,0.0209848090,0.0210404845,0.0210912640,0.0211441185,0.0211925326,0.0212414047,0.0212884953,0.0213335182,0.0213773200,0.0214191685,0.0214627443,0.0215011993,0.0215403441,0.0215749700,0.0216117388,0.0126743196,0.0137928820,0.0145882418,0.0151987694,0.0156909400,0.0161011756,0.0164508979,0.0167573526,0.0170227228,0.0172667801,0.0174801436,0.0176741565,0.0178536209,0.0180170738,0.0181725637,0.0183113676,0.0184421246,0.0185631277,0.0186780663,0.0187891180,0.0188896801,0.0189839985,0.0190761404,0.0191642407,0.0192482215,0.0193231677,0.0193990209,0.0194690762,0.0195388286,0.0196029527,0.0196680246,0.0197300913,0.0197872220,0.0198418289,0.0198977920,0.0199504183,0.0199997020,0.0200500497,0.0200960205,0.0201430249,0.0201880873,0.0202310048,0.0202735675,0.0203124434,0.0203542728,0.0203919337,0.0204291213,0.0204619405,0.0204965334,0.0119879771,0.0130501642,0.0138065933,0.0143877462,0.0148555337,0.0152464040,0.0155797407,0.0158714503,0.0161245901,0.0163568609,0.0165605758,0.0167457489,0.0169164833,0.0170721743,0.0172208677,0.0173534122,0.0174784901,0.0175939142,0.0177033934,0.0178090510,0.0179048586,0.0179954958,0.0180826664,0.0181673798,0.0182475821,0.0183186357,0.0183914334,0.0184585986,0.0185250428,0.0185865195,0.0186488722,0.0187074785,0.0187617336,0.0188143032,0.0188679012,0.0189182780,0.0189649024,0.0190129478,0.0190569760,0.0191021340,0.0191449607,0.0191862216,0.0192269702,0.0192639053,0.0193043244,0.0193399278,0.0193760338,0.0194068856,0.0194398593,0.0113392996,0.0123487445,0.0130674680,0.0136208888,0.0140655912,0.0144375600,0.0147553195,0.0150334075,0.0152741185,0.0154961674,0.0156899310,0.0158665101,0.0160292703,0.0161776400,0.0163199573,0.0164459466,0.0165655119,0.0166760828,0.0167801746,0.0168807606,0.0169727070,0.0170588317,0.0171426922,0.0172231939,0.0172989050,0.0173674689,0.0174362670,0.0175009403,0.0175639976,0.0176228522,0.0176824820,0.0177385599,0.0177908645,0.0178411829,0.0178924322,0.0179397800,0.0179844652,0.0180310933,0.0180722397,0.0181155234,0.0181568615,0.0181957195,0.0182352793,0.0182703141,0.0183093843,0.0183430228,0.0183769348,0.0184064494,0.0184380392,0.0107266118,0.0116852188,0.0123685796,0.0128949726,0.0133185967,0.0136724357,0.0139749684,0.0142407733,0.0144698725,0.0146809340,0.0148658702,0.0150338748,0.0151891542,0.0153308500,0.0154663586,0.0155871620,0.0157007506,0.0158067830,0.0159056029,0.0160019377,0.0160894209,0.0161719191,0.0162519247,0.0163289512,0.0164007363,0.0164670592,0.0165320772,0.0165937559,0.0166542131,0.0167102227,0.0167666501,0.0168206122,0.0168708224,0.0169190197,0.0169675515,0.0170126655,0.0170554330,0.0170998248,0.0171394310,0.0171804136,0.0172204566,0.0172572904,0.0172954771,0.0173285181,0.0173658033,0.0173977665,0.0174298925,0.0174585137,0.0174885157,0.0101477164,0.0110582640,0.0117079317,0.0122087966,0.0126121043,0.0129485347,0.0132369502,0.0134898710,0.0137082034,0.0139091897,0.0140856709,0.0142459570,0.0143938805,0.0145288309,0.0146582765,0.0147735866,0.0148818962,0.0149828028,0.0150774014,0.0151685923,0.0152529717,0.0153313187,0.0154070102,0.0154810981,0.0155499701,0.0156131943,0.0156753554,0.0157339656,0.0157917341,0.0158453863,0.0158994367,0.0159503803,0.0159979168,0.0160447422,0.0160914483,0.0161341397,0.0161745307,0.0162176463,0.0162550562,0.0162944746,0.0163323379,0.0163673855,0.0164047580,0.0164357142,0.0164711011,0.0165020750,0.0165330202,0.0165597196,0.0165884678,0.0096006165,0.0104655725,0.0110837623,0.0115594860,0.0119434614,0.0122634351,0.0125385603,0.0127792912,0.0129878657,0.0131788484,0.0133476955,0.0134996049,0.0136401880,0.0137692810,0.0138933721,0.0140027771,0.0141061805,0.0142027181,0.0142927839,0.0143797943,0.0144603130,0.0145353256,0.0146076853,0.0146776762,0.0147437031,0.0148042494,0.0148634956,0.0149196734,0.0149742906,0.0150258678,0.0150777700,0.0151255742,0.0151713511,0.0152159002,0.0152598164,0.0153017768,0.0153401220,0.0153815588,0.0154173082,0.0154543734,0.0154910676,0.0155240820,0.0155599306,0.0155889511,0.0156227413,0.0156531685,0.0156822708,0.0157081684,0.0157359720,0.0090839825,0.0099049703,0.0104930806,0.0109453031,0.0113109949,0.0116154735,0.0118779398,0.0121067142,0.0123052847,0.0124874325,0.0126483676,0.0127932004,0.0129270963,0.0130501953,0.0131684896,0.0132727653,0.0133713292,0.0134638282,0.0135495388,0.0136326275,0.0137094239,0.0137813694,0.0138498664,0.0139172712,0.0139797417,0.0140375803,0.0140940444,0.0141477995,0.0141998276,0.0142488273,0.0142987051,0.0143441076,0.0143880557,0.0144307986,0.0144722760,0.0145132099,0.0145494423,0.0145888480,0.0146230566,0.0146584422,0.0146935876,0.0147244474,0.0147588776,0.0147865424,0.0148191317,0.0148480186,0.0148762803,0.0149009782,0.0149271808,0.0085952993,0.0093749029,0.0099344001,0.0103643763,0.0107119867,0.0110024415,0.0112520140,0.0114697129,0.0116588233,0.0118326034,0.0119858546,0.0121239005,0.0122519040,0.0123690337,0.0124825204,0.0125812676,0.0126752439,0.0127628698,0.0128454119,0.0129246816,0.0129984130,0.0130666657,0.0131320504,0.0131961552,0.0132561233,0.0133115571,0.0133648101,0.0134160218,0.0134658986,0.0135124368,0.0135604879,0.0136035729,0.0136452592,0.0136866517,0.0137260955,0.0137651690,0.0137993500,0.0138375548,0.0138701287,0.0139039760,0.0139372342,0.0139670113,0.0139999673,0.0140264225,0.0140572290,0.0140845326,0.0141117682,0.0141356239,0.0141600040,0.0081335420,0.0088740846,0.0094052958,0.0098145030,0.0101453122,0.0104221041,0.0106594206,0.0108666942,0.0110471246,0.0112126476,0.0113594271,0.0114909786,0.0116119836,0.0117240077,0.0118319311,0.0119264079,0.0120158909,0.0120995569,0.0121786424,0.0122532439,0.0123246619,0.0123894941,0.0124518377,0.0125129000,0.0125704646,0.0126232072,0.0126742196,0.0127225837,0.0127703598,0.0128148989,0.0128612324,0.0129015593,0.0129414513,0.0129811292,0.0130187781,0.0130557884,0.0130883804,0.0131248691,0.0131561453,0.0131887339,0.0132205277,0.0132488623,0.0132800117,0.0133053457,0.0133344268,0.0133617233,0.0133873221,0.0134097389,0.0134331156,0.0076968031,0.0084004038,0.0089050272,0.0092946456,0.0096094485,0.0098727392,0.0100989175,0.0102959846,0.0104678858,0.0106256191,0.0107657853,0.0108910858,0.0110058460,0.0111131526,0.0112158712,0.0113064308,0.0113908308,0.0114713811,0.0115465523,0.0116174615,0.0116860464,0.0117476806,0.0118075393,0.0118654246,0.0119206204,0.0119704956,0.0120193789,0.0120654092,0.0121106291,0.0121534498,0.0121982312,0.0122365413,0.0122745277,0.0123122770,0.0123482202,0.0123832041,0.0124147053,0.0124495933,0.0124790677,0.0125103026,0.0125415362,0.0125680721,0.0125976234,0.0126218450,0.0126496260,0.0126759313,0.0127001106,0.0127219994,0.0127443324,0.0072844413,0.0079523377,0.0084322266,0.0088027899,0.0091020558,0.0093529242,0.0095681123,0.0097558329,0.0099196558,0.0100698219,0.0102033962,0.0103229658,0.0104320866,0.0105342005,0.0106326597,0.0107188548,0.0107987825,0.0108757300,0.0109481621,0.0110152587,0.0110806213,0.0111396265,0.0111966940,0.0112520923,0.0113046372,0.0113516934,0.0113994683,0.0114428371,0.0114860009,0.0115267574,0.0115691350,0.0116059933,0.0116422189,0.0116779806,0.0117127029,0.0117463696,0.0117762555,0.0118094485,0.0118374011,0.0118677647,0.0118970585,0.0119223370,0.0119507258,0.0119742378,0.0120001065,0.0120258821,0.0120489057,0.0120691849,0.0120907354,0.0068941411,0.0075289813,0.0079847116,0.0083374482,0.0086220640,0.0088604839,0.0090657418,0.0092443335,0.0094004035,0.0095435742,0.0096706076,0.0097847058,0.0098887811,0.0099861824,0.0100797508,0.0101620860,0.0102382843,0.0103112972,0.0103808991,0.0104445881,0.0105071418,0.0105638937,0.0106173738,0.0106706728,0.0107206510,0.0107657758,0.0108111743,0.0108525407,0.0108935147,0.0109327108,0.0109730471,0.0110081162,0.0110426295,0.0110769011,0.0111102105,0.0111422402,0.0111713787,0.0112022988,0.0112294649,0.0112582651,0.0112862298,0.0113104702,0.0113377189,0.0113598207,0.0113847676,0.0114091601,0.0114307544,0.0114502483,0.0114713290,0.0065248743,0.0071280384,0.0075612020,0.0078967497,0.0081676627,0.0083944849,0.0085899785,0.0087593505,0.0089084790,0.0090454349,0.0091661267,0.0092749581,0.0093739953,0.0094670460,0.0095559879,0.0096344640,0.0097073886,0.0097766875,0.0098433200,0.0099035057,0.0099635662,0.0100179236,0.0100687969,0.0101198763,0.0101667527,0.0102099473,0.0102535327,0.0102928875,0.0103324757,0.0103697718,0.0104084579,0.0104416457,0.0104743907,0.0105076504,0.0105389910,0.0105694012,0.0105974294,0.0106266884,0.0106524499,0.0106800171,0.0107070007,0.0107299836,0.0107561089,0.0107772481,0.0108009097,0.0108245209,0.0108449866,0.0108633325,0.0108837619,0.0061760347,0.0067491076,0.0071607061,0.0074796858,0.0077370862,0.0079535267,0.0081393290,0.0083006852,0.0084427926,0.0085732741,0.0086883344,0.0087917413,0.0088862910,0.0089746976,0.0090595542,0.0091344365,0.0092038454,0.0092698546,0.0093337281,0.0093909729,0.0094481669,0.0095001399,0.0095486869,0.0095977708,0.0096425671,0.0096833248,0.0097250060,0.0097625927,0.0098003413,0.0098362503,0.0098728562,0.0099044628,0.0099357210,0.0099678099,0.0099975202,0.0100263850,0.0100529935,0.0100807035,0.0101060491,0.0101322355,0.0101575302,0.0101798354,0.0102047630,0.0102248851,0.0102472173,0.0102697447,0.0102896947,0.0103065493,0.0103267069,0.0058459537,0.0063903887,0.0067816846,0.0070849304,0.0073299051,0.0075359736,0.0077127812,0.0078660497,0.0080017782,0.0081259390,0.0082356504,0.0083344104,0.0084243022,0.0085081333,0.0085892581,0.0086609728,0.0087269078,0.0087899520,0.0088507923,0.0089056340,0.0089598942,0.0090095672,0.0090560415,0.0091025235,0.0091452820,0.0091841282,0.0092238008,0.0092595301,0.0092962871,0.0093303161,0.0093649209,0.0093951343,0.0094250816,0.0094555850,0.0094837351,0.0095115391,0.0095370474,0.0095633494,0.0095877390,0.0096122711,0.0096370463,0.0096582170,0.0096819541,0.0097009566,0.0097218175,0.0097436100,0.0097630179,0.0097794215,0.0097980805,0.0055339955,0.0060510225,0.0064233572,0.0067112989,0.0069443370,0.0071402229,0.0073086921,0.0074547894,0.0075840490,0.0077025355,0.0078065211,0.0079013823,0.0079867257,0.0080666680,0.0081436564,0.0082119622,0.0082752479,0.0083352787,0.0083931930,0.0084456689,0.0084972049,0.0085442178,0.0085888629,0.0086334321,0.0086738143,0.0087108770,0.0087490110,0.0087828100,0.0088179442,0.0088508541,0.0088836362,0.0089122676,0.0089409907,0.0089702891,0.0089970603,0.0090232214,0.0090478048,0.0090730529,0.0090964695,0.0091198893,0.0091431678,0.0091637080,0.0091862834,0.0092045442,0.0092243093,0.0092453576,0.0092637179,0.0092787152,0.0092968532,0.0052387997,0.0057299791,0.0060839482,0.0063572526,0.0065794792,0.0067659477,0.0069263092,0.0070648539,0.0071884705,0.0073012985,0.0074002034,0.0074906278,0.0075719730,0.0076478240,0.0077214111,0.0077862719,0.0078469211,0.0079038303,0.0079593031,0.0080094676,0.0080584444,0.0081031051,0.0081457147,0.0081881629,0.0082271728,0.0082619542,0.0082986510,0.0083308284,0.0083645009,0.0083961595,0.0084273952,0.0084548551,0.0084817742,0.0085098448,0.0085352972,0.0085601836,0.0085837842,0.0086079946,0.0086303341,0.0086526704,0.0086748922,0.0086947564,0.0087162078,0.0087334758,0.0087522832,0.0087725779,0.0087897087,0.0088037942,0.0088215639,0.0049599825,0.0054259869,0.0057624477,0.0060223760,0.0062338891,0.0064113907,0.0065641157,0.0066957248,0.0068139852,0.0069210140,0.0070155964,0.0071016048,0.0071789434,0.0072511515,0.0073211970,0.0073831422,0.0074410377,0.0074951685,0.0075481825,0.0075961387,0.0076427276,0.0076854627,0.0077260259,0.0077662636,0.0078032997,0.0078368511,0.0078719536,0.0079025581,0.0079346171,0.0079646319,0.0079950102,0.0080209639,0.0080463574,0.0080732187,0.0080974979,0.0081214902,0.0081438910,0.0081671146,0.0081883186,0.0082095692,0.0082309673,0.0082493541,0.0082701962,0.0082867461,0.0083044949,0.0083241329,0.0083407642,0.0083533745,0.0083708477,0.0046958404,0.0051384436,0.0054582950,0.0057055090,0.0059070053,0.0060757790,0.0062211744,0.0063465784,0.0064591028,0.0065607766,0.0066507313,0.0067329311,0.0068065275,0.0068753643,0.0069419611,0.0070014272,0.0070560773,0.0071079311,0.0071586467,0.0072042012,0.0072489652,0.0072891433,0.0073279752,0.0073665137,0.0074016789,0.0074335358,0.0074673710,0.0074959850,0.0075270919,0.0075558093,0.0075846489,0.0076093576,0.0076336971,0.0076593149,0.0076825908,0.0077054019,0.0077268987,0.0077489048,0.0077694288,0.0077891011,0.0078100912,0.0078270794,0.0078470880,0.0078627992,0.0078798573,0.0078988512,0.0079142344,0.0079268856,0.0079431235,0.0044461245,0.0048663924,0.0051707114,0.0054057802,0.0055971564,0.0057577085,0.0058964212,0.0060155144,0.0061229885,0.0062197745,0.0063054957,0.0063834200,0.0064538152,0.0065192701,0.0065829503,0.0066392697,0.0066917854,0.0067411187,0.0067892341,0.0068328303,0.0068756511,0.0069138168,0.0069506032,0.0069874419,0.0070213755,0.0070515129,0.0070835937,0.0071107713,0.0071406210,0.0071682536,0.0071953156,0.0072191114,0.0072424029,0.0072667865,0.0072888730,0.0073105618,0.0073310943,0.0073523136,0.0073719836,0.0073908088,0.0074108659,0.0074267199,0.0074458942,0.0074608808,0.0074771413,0.0074953631,0.0075100969,0.0075223077,0.0075376686,0.0042097577,0.0046088423,0.0048983479,0.0051217387,0.0053038217,0.0054563997,0.0055888176,0.0057019283,0.0058042333,0.0058963595,0.0059781114,0.0060522578,0.0061193096,0.0061818209,0.0062425236,0.0062958923,0.0063462195,0.0063934134,0.0064390262,0.0064807098,0.0065216336,0.0065580173,0.0065929540,0.0066280508,0.0066607135,0.0066891447,0.0067199276,0.0067457246,0.0067744346,0.0068001126,0.0068262672,0.0068491857,0.0068713059,0.0068944714,0.0069155616,0.0069364235,0.0069558528,0.0069759694,0.0069948980,0.0070132262,0.0070318350,0.0070470690,0.0070654476,0.0070794284,0.0070950150,0.0071127628,0.0071269722,0.0071382989,0.0071528008,0.0039862782,0.0043652973,0.0046402582,0.0048529348,0.0050258356,0.0051710883,0.0052973675,0.0054047755,0.0055025025,0.0055900408,0.0056677444,0.0057385792,0.0058026753,0.0058617757,0.0059200226,0.0059706650,0.0060186957,0.0060635229,0.0061073035,0.0061468106,0.0061857060,0.0062207020,0.0062543182,0.0062869705,0.0063186777,0.0063454579,0.0063747556,0.0063991632,0.0064267859,0.0064515542,0.0064764629,0.0064981344,0.0065192596,0.0065416095,0.0065613598,0.0065810604,0.0066000371,0.0066194528,0.0066370436,0.0066547436,0.0066726233,0.0066870467,0.0067045402,0.0067177795,0.0067326968,0.0067499189,0.0067632775,0.0067741924,0.0067877913,0.0037747841,0.0041348996,0.0043958686,0.0045982002,0.0047624908,0.0049011751,0.0050210937,0.0051236557,0.0052163505,0.0052999979,0.0053738630,0.0054412546,0.0055021873,0.0055585568,0.0056141443,0.0056621858,0.0057082540,0.0057508933,0.0057927555,0.0058303413,0.0058673018,0.0059008843,0.0059329841,0.0059643280,0.0059940777,0.0060198118,0.0060477180,0.0060708814,0.0060972382,0.0061211657,0.0061447376,0.0061652705,0.0061857516,0.0062069214,0.0062256215,0.0062442300,0.0062623882,0.0062814304,0.0062977759,0.0063149967,0.0063319456,0.0063456236,0.0063623453,0.0063748374,0.0063890721,0.0064052930,0.0064186360,0.0064287016,0.0064415846,0.0035745979,0.0039163325,0.0041645383,0.0043571579,0.0045134179,0.0046453308,0.0047595191,0.0048571358,0.0049452561,0.0050252987,0.0050954458,0.0051595829,0.0052176034,0.0052712557,0.0053240229,0.0053700281,0.0054140375,0.0054545354,0.0054945835,0.0055300774,0.0055656519,0.0055976049,0.0056280244,0.0056582405,0.0056864765,0.0057111287,0.0057375471,0.0057596735,0.0057846749,0.0058075620,0.0058299106,0.0058497827,0.0058695701,0.0058890592,0.0059072275,0.0059249561,0.0059422071,0.0059606628,0.0059761339,0.0059922682,0.0060085413,0.0060214996,0.0060375349,0.0060495596,0.0060627652,0.0060786401,0.0060913907,0.0061008511,0.0061132297,0.0033853421,0.0037098154,0.0039458446,0.0041286446,0.0042774990,0.0044030874,0.0045116649,0.0046048179,0.0046883596,0.0047644911,0.0048316779,0.0048928599,0.0049476718,0.0049989955,0.0050491690,0.0050931790,0.0051350353,0.0051737638,0.0052119814,0.0052456891,0.0052793558,0.0053100625,0.0053388500,0.0053677444,0.0053944941,0.0054181200,0.0054435791,0.0054642659,0.0054884113,0.0055102957,0.0055315914,0.0055505828,0.0055693602,0.0055877197,0.0056053064,0.0056219719,0.0056386471,0.0056563643,0.0056708547,0.0056862134,0.0057021950,0.0057140632,0.0057294843,0.0057411763,0.0057534226,0.0057686972,0.0057807273,0.0057895858,0.0058017545,0.0032062153,0.0035144404,0.0037386266,0.0039125855,0.0040541845,0.0041735948,0.0042767370,0.0043654105,0.0044451245,0.0045176209,0.0045815497,0.0046397719,0.0046918055,0.0047411351,0.0047885722,0.0048307203,0.0048704361,0.0049073965,0.0049439574,0.0049759299,0.0050080512,0.0050371560,0.0050648575,0.0050922283,0.0051179454,0.0051405566,0.0051643252,0.0051844695,0.0052073201,0.0052280580,0.0052486163,0.0052665473,0.0052843507,0.0053023038,0.0053187772,0.0053346355,0.0053506667,0.0053676931,0.0053813721,0.0053957480,0.0054113291,0.0054227059,0.0054371516,0.0054485199,0.0054598930,0.0054745479,0.0054863861,0.0054948099,0.0055063133,0.0030366577,0.0033296291,0.0035422256,0.0037079398,0.0038428359,0.0039562611,0.0040542909,0.0041385427,0.0042145475,0.0042837473,0.0043443896,0.0044000676,0.0044496353,0.0044963268,0.0045417121,0.0045818523,0.0046198115,0.0046547328,0.0046899969,0.0047200576,0.0047508446,0.0047784039,0.0048050758,0.0048312551,0.0048554278,0.0048772601,0.0048997408,0.0049190302,0.0049407063,0.0049605425,0.0049802025,0.0049974322,0.0050144096,0.0050313489,0.0050470215,0.0050625182,0.0050772692,0.0050935610,0.0051071315,0.0051206204,0.0051353177,0.0051461254,0.0051598612,0.0051706590,0.0051817633,0.0051958071,0.0052069154,0.0052150216,0.0052258720,0.0028762062,0.0031545590,0.0033564314,0.0035141559,0.0036422285,0.0037502369,0.0038432325,0.0039237596,0.0039963150,0.0040618262,0.0041199835,0.0041727846,0.0042199195,0.0042642707,0.0043078664,0.0043459793,0.0043819268,0.0044154427,0.0044488751,0.0044775599,0.0045069227,0.0045331297,0.0045585330,0.0045835582,0.0046067072,0.0046275103,0.0046490096,0.0046674637,0.0046880119,0.0047067450,0.0047256044,0.0047421746,0.0047580595,0.0047745151,0.0047892332,0.0048039871,0.0048186011,0.0048335328,0.0048462762,0.0048594951,0.0048735747,0.0048839357,0.0048967977,0.0049072523,0.0049179195,0.0049314467,0.0049417421,0.0049497347,0.0049597447,0.0027242540,0.0029886717,0.0031805198,0.0033303972,0.0034521691,0.0035549867,0.0036433721,0.0037199352,0.0037893755,0.0038514752,0.0039070757,0.0039571792,0.0040021312,0.0040442304,0.0040858989,0.0041220294,0.0041563410,0.0041887065,0.0042202101,0.0042478192,0.0042756642,0.0043005927,0.0043248404,0.0043489447,0.0043708424,0.0043906013,0.0044112171,0.0044286040,0.0044482026,0.0044658349,0.0044843723,0.0044997801,0.0045151336,0.0045310020,0.0045447686,0.0045591528,0.0045729548,0.0045867051,0.0045990702,0.0046118530,0.0046249864,0.0046352930,0.0046471987,0.0046573284,0.0046676495,0.0046806390,0.0046903029,0.0046979489,0.0047075630,0.0025804090,0.0028315539,0.0030140584,0.0031562677,0.0032721636,0.0033698433,0.0034540396,0.0035267983,0.0035929983,0.0036521669,0.0037051302,0.0037528256,0.0037955452,0.0038359281,0.0038754212,0.0039099695,0.0039426651,0.0039734125,0.0040033396,0.0040297820,0.0040563846,0.0040801348,0.0041033857,0.0041261570,0.0041469504,0.0041658933,0.0041855884,0.0042020861,0.0042207457,0.0042376721,0.0042554005,0.0042699127,0.0042846769,0.0042998635,0.0043126796,0.0043265036,0.0043400485,0.0043527283,0.0043647991,0.0043767420,0.0043893262,0.0043990174,0.0044106128,0.0044201779,0.0044299950,0.0044426520,0.0044516870,0.0044590792,0.0044680044,0.0024443780,0.0026830043,0.0028560998,0.0029913929,0.0031018818,0.0031946349,0.0032746285,0.0033439385,0.0034068019,0.0034633810,0.0035135644,0.0035591584,0.0035999580,0.0036383288,0.0036757870,0.0037088934,0.0037399569,0.0037691164,0.0037978523,0.0038231284,0.0038481852,0.0038711960,0.0038933608,0.0039149041,0.0039345307,0.0039525427,0.0039715098,0.0039871439,0.0040047823,0.0040210625,0.0040377439,0.0040521124,0.0040656606,0.0040802673,0.0040928363,0.0041058941,0.0041186568,0.0041309457,0.0041420552,0.0041538265,0.0041659554,0.0041751232,0.0041861585,0.0041950219,0.0042045789,0.0042166524,0.0042253800,0.0042323916,0.0042407487,0.0023154280,0.0025420554,0.0027066421,0.0028352343,0.0029403553,0.0030286924,0.0031047433,0.0031707480,0.0032305584,0.0032844881,0.0033321472,0.0033755213,0.0034144777,0.0034510558,0.0034866683,0.0035181826,0.0035478853,0.0035757473,0.0036028096,0.0036274961,0.0036508394,0.0036728530,0.0036938924,0.0037144290,0.0037333013,0.0037502403,0.0037683766,0.0037836167,0.0038000372,0.0038157542,0.0038316932,0.0038452234,0.0038580722,0.0038721802,0.0038841158,0.0038965514,0.0039086304,0.0039204270,0.0039313291,0.0039423174,0.0039537307,0.0039628553,0.0039729757,0.0039818370,0.0039907947,0.0040021954,0.0040103767,0.0040171433,0.0040255339,0.0021934503,0.0024087153,0.0025649686,0.0026874668,0.0027874105,0.0028714195,0.0029440039,0.0030065584,0.0030635639,0.0031149497,0.0031600741,0.0032014530,0.0032385532,0.0032735044,0.0033074279,0.0033372924,0.0033657827,0.0033924683,0.0034179908,0.0034416537,0.0034636334,0.0034847870,0.0035048473,0.0035243458,0.0035425744,0.0035586683,0.0035755993,0.0035904010,0.0036058339,0.0036210989,0.0036360999,0.0036489705,0.0036614994,0.0036748197,0.0036859264,0.0036978861,0.0037095106,0.0037209338,0.0037314117,0.0037416317,0.0037525334,0.0037612208,0.0037711136,0.0037791724,0.0037881759,0.0037987634,0.0038065861,0.0038129305,0.0038210816,0.0020778679,0.0022823477,0.0024310862,0.0025474381,0.0026424074,0.0027222391,0.0027914980,0.0028509678,0.0029052443,0.0029545066,0.0029970604,0.0030364877,0.0030717833,0.0031049851,0.0031375093,0.0031659545,0.0031930075,0.0032184125,0.0032429035,0.0032653647,0.0032864051,0.0033063773,0.0033256294,0.0033440466,0.0033613491,0.0033767977,0.0033928394,0.0034068186,0.0034216984,0.0034362880,0.0034506277,0.0034629115,0.0034748652,0.0034874453,0.0034980876,0.0035094028,0.0035206379,0.0035315364,0.0035414536,0.0035510936,0.0035617416,0.0035700333,0.0035793081,0.0035871859,0.0035954971,0.0036057152,0.0036132191,0.0036194961,0.0036268374,0.0019685916,0.0021628191,0.0023039775,0.0024149165,0.0025051544,0.0025807553,0.0026469689,0.0027035366,0.0027551980,0.0028018293,0.0028425133,0.0028799960,0.0029134216,0.0029452541,0.0029762791,0.0030033299,0.0030292012,0.0030533643,0.0030766909,0.0030981335,0.0031180118,0.0031371519,0.0031558450,0.0031730220,0.0031896076,0.0032043458,0.0032194572,0.0032330340,0.0032472113,0.0032610314,0.0032747215,0.0032864106,0.0032978019,0.0033097394,0.0033199981,0.0033305798,0.0033411912,0.0033517625,0.0033614198,0.0033705532,0.0033806158,0.0033885364,0.0033972501,0.0034049024,0.0034127595,0.0034225037,0.0034298696,0.0034358403,0.0034422601,0.0018650155,0.0020496484,0.0021836314,0.0022893448,0.0023749218,0.0024468348,0.0025096524,0.0025635970,0.0026128488,0.0026572952,0.0026959551,0.0027317646,0.0027636317,0.0027938404,0.0028236072,0.0028491091,0.0028736435,0.0028968617,0.0029189953,0.0029395343,0.0029583083,0.0029768142,0.0029944928,0.0030108676,0.0030269670,0.0030404435,0.0030549497,0.0030681035,0.0030815709,0.0030947942,0.0031076927,0.0031188289,0.0031299520,0.0031410700,0.0031508950,0.0031610791,0.0031712534,0.0031811373,0.0031905300,0.0031992064,0.0032088830,0.0032162576,0.0032248405,0.0032320608,0.0032396741,0.0032487270,0.0032554990,0.0032615167,0.0032676572,0.0017671605,0.0019422771,0.0020697268,0.0021702557,0.0022514718,0.0023197214,0.0023798501,0.0024310907,0.0024779349,0.0025203309,0.0025571458,0.0025910811,0.0026212306,0.0026501595,0.0026783770,0.0027028090,0.0027265078,0.0027483698,0.0027696024,0.0027887432,0.0028068277,0.0028244900,0.0028413149,0.0028570802,0.0028725992,0.0028853673,0.0028989286,0.0029118326,0.0029242457,0.0029369469,0.0029493653,0.0029599124,0.0029706860,0.0029810312,0.0029904232,0.0030001475,0.0030100483,0.0030193431,0.0030283073,0.0030366716,0.0030457319,0.0030527990,0.0030611484,0.0030679653,0.0030751739,0.0030839120,0.0030903793,0.0030959158,0.0031019335,0.0016743073,0.0018406832,0.0019616820,0.0020574090,0.0021346773,0.0021997904,0.0022564476,0.0023056114,0.0023500513,0.0023905139,0.0024256327,0.0024576700,0.0024864223,0.0025140178,0.0025408440,0.0025644102,0.0025867465,0.0026076153,0.0026278277,0.0026461326,0.0026632636,0.0026800769,0.0026960949,0.0027109134,0.0027258062,0.0027382471,0.0027511004,0.0027634639,0.0027750395,0.0027872395,0.0027993294,0.0028091701,0.0028194592,0.0028294750,0.0028383036,0.0028472891,0.0028570168,0.0028659211,0.0028745052,0.0028823326,0.0028909001,0.0028976791,0.0029056599,0.0029122646,0.0029192596,0.0029275014,0.0029336823,0.0029388331,0.0029446355,0.0015863610,0.0017444662,0.0018594712,0.0019503192,0.0020241036,0.0020858614,0.0021397879,0.0021867156,0.0022289895,0.0022673583,0.0023009730,0.0023315878,0.0023585398,0.0023848899,0.0024105448,0.0024331702,0.0024543573,0.0024740951,0.0024931603,0.0025108084,0.0025270391,0.0025432404,0.0025583644,0.0025726294,0.0025865946,0.0025985674,0.0026107662,0.0026226298,0.0026336972,0.0026450623,0.0026566275,0.0026658982,0.0026759164,0.0026853854,0.0026938085,0.0027024911,0.0027117595,0.0027202924,0.0027284341,0.0027359655,0.0027439473,0.0027504628,0.0027580300,0.0027644434,0.0027709449,0.0027789110,0.0027849160,0.0027896011,0.0027952174,0.0015031537,0.0016533484,0.0017625141,0.0018489120,0.0019189299,0.0019778693,0.0020292507,0.0020739214,0.0021141328,0.0021505977,0.0021826048,0.0022118220,0.0022374480,0.0022625272,0.0022869680,0.0023088030,0.0023287710,0.0023473572,0.0023657108,0.0023824719,0.0023979782,0.0024134063,0.0024275692,0.0024414108,0.0024546234,0.0024661664,0.0024776605,0.0024890582,0.0024995443,0.0025102833,0.0025212914,0.0025303934,0.0025399610,0.0025489648,0.0025568334,0.0025651293,0.0025738112,0.0025821787,0.0025899486,0.0025971835,0.0026046393,0.0026104946,0.0026182732,0.0026243222,0.0026303899,0.0026379979,0.0026436223,0.0026483049,0.0026536630,0.0014244348,0.0015669563,0.0016708312,0.0017527878,0.0018195190,0.0018756612,0.0019245052,0.0019668451,0.0020052937,0.0020398149,0.0020703475,0.0020983365,0.0021226110,0.0021464373,0.0021697561,0.0021903916,0.0022095573,0.0022271497,0.0022447964,0.0022606060,0.0022754224,0.0022900228,0.0023038182,0.0023169585,0.0023294494,0.0023403879,0.0023514093,0.0023623305,0.0023720300,0.0023823689,0.0023930133,0.0024016947,0.0024107945,0.0024193901,0.0024269577,0.0024347399,0.0024429801,0.0024511495,0.0024583866,0.0024652866,0.0024723928,0.0024779103,0.0024853424,0.0024909837,0.0024971239,0.0025041955,0.0025096100,0.0025140266,0.0025192977,0.0013499516,0.0014852128,0.0015839246,0.0016617280,0.0017252286,0.0017787812,0.0018251235,0.0018652427,0.0019019779,0.0019345810,0.0019640141,0.0019905219,0.0020136529,0.0020363836,0.0020586942,0.0020784746,0.0020964022,0.0021132313,0.0021298144,0.0021451400,0.0021592148,0.0021730078,0.0021862159,0.0021987554,0.0022107582,0.0022212708,0.0022316664,0.0022420561,0.0022513572,0.0022610972,0.0022712955,0.0022796839,0.0022882732,0.0022964643,0.0023035208,0.0023110538,0.0023188751,0.0023265431,0.0023334861,0.0023402459,0.0023469380,0.0023521523,0.0023592509,0.0023647456,0.0023705113,0.0023771934,0.0023825735,0.0023867115,0.0023914467,0.0012792701,0.0014077644,0.0015014436,0.0015755093,0.0016358546,0.0016867676,0.0017308773,0.0017693036,0.0018040955,0.0018351406,0.0018630482,0.0018883215,0.0019101863,0.0019319576,0.0019531798,0.0019719724,0.0019892409,0.0020051079,0.0020207841,0.0020355895,0.0020489935,0.0020623382,0.0020745619,0.0020865012,0.0020980733,0.0021080321,0.0021179654,0.0021279192,0.0021369634,0.0021459139,0.0021555464,0.0021637165,0.0021719467,0.0021797110,0.0021863807,0.0021937156,0.0022010656,0.0022083428,0.0022150544,0.0022215829,0.0022279375,0.0022329105,0.0022395448,0.0022447785,0.0022503101,0.0022567484,0.0022618079,0.0022658596,0.0022703157,0.0012125171,0.0013345629,0.0014234105,0.0014938330,0.0015511077,0.0015995984,0.0016414893,0.0016779280,0.0017112364,0.0017408520,0.0017673174,0.0017914496,0.0018122776,0.0018329743,0.0018533292,0.0018710819,0.0018876350,0.0019026840,0.0019175441,0.0019317336,0.0019443535,0.0019571620,0.0019688212,0.0019802825,0.0019911381,0.0020005744,0.0020101055,0.0020197985,0.0020283593,0.0020366257,0.0020458762,0.0020538575,0.0020614795,0.0020690890,0.0020753504,0.0020822315,0.0020892578,0.0020960985,0.0021027838,0.0021090672,0.0021149447,0.0021196523,0.0021259834,0.0021310593,0.0021363165,0.0021425680,0.0021471789,0.0021512246,0.0021555426,0.0011489956,0.0012650382,0.0013494867,0.0014163520,0.0014708967,0.0015170510,0.0015567937,0.0015914174,0.0016230441,0.0016514262,0.0016765359,0.0016995341,0.0017195077,0.0017389742,0.0017584679,0.0017751945,0.0017913275,0.0018054061,0.0018195810,0.0018329059,0.0018451139,0.0018571965,0.0018684933,0.0018794285,0.0018896432,0.0018987255,0.0019076900,0.0019170671,0.0019252486,0.0019331697,0.0019419895,0.0019492255,0.0019566788,0.0019638765,0.0019700090,0.0019764612,0.0019832070,0.0019899763,0.0019960922,0.0020021647,0.0020078197,0.0020123306,0.0020184362,0.0020231264,0.0020280971,0.0020341216,0.0020385127,0.0020423437,0.0020464748,0.0010890816,0.0011991306,0.0012793356,0.0013429986,0.0013948241,0.0014388091,0.0014765330,0.0015095303,0.0015395472,0.0015665251,0.0015905907,0.0016124595,0.0016315785,0.0016500570,0.0016683677,0.0016845688,0.0016998003,0.0017132309,0.0017266711,0.0017393103,0.0017509880,0.0017625199,0.0017733905,0.0017834814,0.0017934007,0.0018022689,0.0018105816,0.0018194992,0.0018274217,0.0018350290,0.0018433793,0.0018500978,0.0018571434,0.0018641444,0.0018699222,0.0018762316,0.0018827091,0.0018890655,0.0018948593,0.0019007916,0.0019060158,0.0019102294,0.0019162144,0.0019208398,0.0019254201,0.0019309612,0.0019353242,0.0019390128,0.0019429545,0.0010320568,0.0011367692,0.0012128101,0.0012733826,0.0013226407,0.0013647066,0.0014004330,0.0014319117,0.0014602834,0.0014861067,0.0015089499,0.0015298407,0.0015478874,0.0015656203,0.0015829805,0.0015983611,0.0016128827,0.0016258187,0.0016386399,0.0016505072,0.0016617589,0.0016727385,0.0016829559,0.0016926449,0.0017021143,0.0017106371,0.0017184227,0.0017270634,0.0017345970,0.0017418174,0.0017496673,0.0017560754,0.0017630727,0.0017694633,0.0017749908,0.0017809143,0.0017870784,0.0017931860,0.0017989277,0.0018043988,0.0018095174,0.0018133631,0.0018191734,0.0018233284,0.0018277831,0.0018331395,0.0018373901,0.0018408777,0.0018446312,0.0009782315,0.0010777320,0.0011499328,0.0012073517,0.0012543121,0.0012943433,0.0013282042,0.0013581694,0.0013850673,0.0014098176,0.0014314747,0.0014513190,0.0014687731,0.0014854394,0.0015019088,0.0015166339,0.0015305170,0.0015428362,0.0015552245,0.0015663109,0.0015770819,0.0015874615,0.0015972155,0.0016064201,0.0016154528,0.0016235446,0.0016311370,0.0016393851,0.0016465444,0.0016533318,0.0016607086,0.0016670377,0.0016736374,0.0016797415,0.0016848611,0.0016905446,0.0016964933,0.0017023504,0.0017078006,0.0017129399,0.0017178474,0.0017215138,0.0017272387,0.0017311160,0.0017352345,0.0017402739,0.0017444208,0.0017477003,0.0017514608,0.0009271528,0.0010216441,0.0010902209,0.0011448839,0.0011894919,0.0012274860,0.0012598542,0.0012883072,0.0013139085,0.0013374737,0.0013580043,0.0013771047,0.0013936828,0.0014095788,0.0014250535,0.0014390854,0.0014523743,0.0014640928,0.0014759224,0.0014864956,0.0014967188,0.0015066819,0.0015159443,0.0015246878,0.0015331001,0.0015408677,0.0015482195,0.0015561761,0.0015628801,0.0015694921,0.0015763536,0.0015824317,0.0015885276,0.0015946689,0.0015993405,0.0016048292,0.0016104159,0.0016159697,0.0016212590,0.0016261043,0.0016308456,0.0016343701,0.0016397880,0.0016435761,0.0016475683,0.0016521397,0.0016562854,0.0016592107,0.0016629010,0.0008787943,0.0009685757,0.0010337295,0.0010856816,0.0011281153,0.0011642648,0.0011949681,0.0012221980,0.0012464015,0.0012689350,0.0012882947,0.0013064580,0.0013222515,0.0013375417,0.0013523023,0.0013656596,0.0013781257,0.0013892935,0.0014007934,0.0014106748,0.0014204331,0.0014298424,0.0014386747,0.0014471135,0.0014550066,0.0014625336,0.0014694875,0.0014772108,0.0014836208,0.0014898606,0.0014963442,0.0015022248,0.0015080528,0.0015138114,0.0015183325,0.0015235175,0.0015289289,0.0015342184,0.0015391650,0.0015437928,0.0015482184,0.0015514350,0.0015567544,0.0015602920,0.0015641730,0.0015684670,0.0015725475,0.0015752835,0.0015789524,0.0008329873,0.0009183142,0.0009803429,0.0010294266,0.0010698820,0.0011042495,0.0011334173,0.0011593585,0.0011824528,0.0012038200,0.0012223924,0.0012396397,0.0012545784,0.0012690419,0.0012830814,0.0012958931,0.0013076199,0.0013183832,0.0013294608,0.0013388876,0.0013479860,0.0013570880,0.0013654922,0.0013734657,0.0013810370,0.0013881987,0.0013948170,0.0014021179,0.0014082882,0.0014142578,0.0014202970,0.0014261789,0.0014315134,0.0014369941,0.0014412766,0.0014463292,0.0014514223,0.0014565002,0.0014611896,0.0014655696,0.0014699504,0.0014728993,0.0014780123,0.0014815487,0.0014849042,0.0014892468,0.0014930090,0.0014956616,0.0014991248,0.0007895405,0.0008706512,0.0009295588,0.0009763361,0.0010147015,0.0010473502,0.0010751462,0.0010996997,0.0011217308,0.0011421234,0.0011597981,0.0011761354,0.0011904254,0.0012040906,0.0012175639,0.0012296696,0.0012409069,0.0012512176,0.0012617249,0.0012706056,0.0012793790,0.0012879285,0.0012963394,0.0013037399,0.0013109840,0.0013175832,0.0013238189,0.0013308689,0.0013369220,0.0013423857,0.0013482212,0.0013538887,0.0013590850,0.0013642243,0.0013681053,0.0013730705,0.0013779394,0.0013827397,0.0013873834,0.0013913746,0.0013956044,0.0013983278,0.0014033354,0.0014066824,0.0014098775,0.0014139255,0.0014175696,0.0014201533,0.0014235514,0.0007483274,0.0008254029,0.0008814229,0.0009259633,0.0009624474,0.0009934606,0.0010198235,0.0010431811,0.0010639780,0.0010834615,0.0011002621,0.0011158113,0.0011294658,0.0011425488,0.0011554521,0.0011669116,0.0011776236,0.0011874343,0.0011975012,0.0012058356,0.0012143913,0.0012225316,0.0012303474,0.0012375245,0.0012444049,0.0012506556,0.0012566742,0.0012634512,0.0012691616,0.0012743044,0.0012798784,0.0012852871,0.0012900905,0.0012949479,0.0012988799,0.0013034912,0.0013080633,0.0013127734,0.0013171572,0.0013209305,0.0013250778,0.0013275057,0.0013323637,0.0013354147,0.0013386333,0.0013423650,0.0013458484,0.0013485643,0.0013515592,0.0007094651,0.0007825686,0.0008358258,0.0008780495,0.0009129154,0.0009423696,0.0009673249,0.0009896545,0.0010092772,0.0010281556,0.0010440293,0.0010587086,0.0010716954,0.0010842364,0.0010964239,0.0011074207,0.0011175851,0.0011268317,0.0011364167,0.0011444504,0.0011526225,0.0011604475,0.0011679316,0.0011747406,0.0011812610,0.0011872576,0.0011928212,0.0011994029,0.0012047040,0.0012097335,0.0012150199,0.0012201311,0.0012247328,0.0012293101,0.0012331400,0.0012375287,0.0012417526,0.0012462901,0.0012504985,0.0012540319,0.0012580564,0.0012603059,0.0012650597,0.0012678597,0.0012710310,0.0012747592,0.0012778880,0.0012803229,0.0012833072,0.0006725499,0.0007420300,0.0007926248,0.0008327387,0.0008659483,0.0008939110,0.0009175482,0.0009388878,0.0009575772,0.0009754324,0.0009905582,0.0010046063,0.0010168776,0.0010289160,0.0010404814,0.0010508179,0.0010606859,0.0010695449,0.0010786470,0.0010862493,0.0010938584,0.0011014860,0.0011084375,0.0011152134,0.0011211941,0.0011268012,0.0011322027,0.0011385529,0.0011435562,0.0011483012,0.0011532837,0.0011583711,0.0011626537,0.0011670556,0.0011708126,0.0011748330,0.0011788932,0.0011832188,0.0011872749,0.0011907089,0.0011945929,0.0011966700,0.0012010912,0.0012036988,0.0012067231,0.0012104881,0.0012134469,0.0012156152,0.0012185935,0.0006375201,0.0007036311,0.0007515773,0.0007896721,0.0008213566,0.0008479645,0.0008704140,0.0008906458,0.0009085437,0.0009255287,0.0009397789,0.0009531860,0.0009649140,0.0009764479,0.0009872916,0.0009972501,0.0010066040,0.0010150935,0.0010237964,0.0010310163,0.0010382347,0.0010455408,0.0010520432,0.0010585665,0.0010642582,0.0010695594,0.0010748397,0.0010807911,0.0010855404,0.0010900164,0.0010948924,0.0010994281,0.0011038477,0.0011079958,0.0011114358,0.0011153529,0.0011192209,0.0011233506,0.0011272497,0.0011305245,0.0011342015,0.0011361439,0.0011404467,0.0011428995,0.0011458883,0.0011493932,0.0011522034,0.0011543283,0.0011570306,0.0006043583,0.0006672437,0.0007126862,0.0007490056,0.0007791304,0.0008043136,0.0008256817,0.0008450708,0.0008618986,0.0008780984,0.0008917574,0.0009044353,0.0009155622,0.0009267285,0.0009368967,0.0009464873,0.0009551742,0.0009633191,0.0009717475,0.0009786306,0.0009855033,0.0009923582,0.0009985706,0.0010048213,0.0010101999,0.0010152678,0.0010203528,0.0010259819,0.0010305494,0.0010348487,0.0010393482,0.0010437707,0.0010478999,0.0010519546,0.0010553068,0.0010589208,0.0010624730,0.0010663662,0.0010702360,0.0010733838,0.0010768482,0.0010787892,0.0010828491,0.0010849921,0.0010880339,0.0010914017,0.0010940134,0.0010962839,0.0010984941,0.0005729697,0.0006325932,0.0006757588,0.0007104285,0.0007390082,0.0007628779,0.0007833492,0.0008016446,0.0008178054,0.0008330478,0.0008462107,0.0008581529,0.0008689346,0.0008793676,0.0008891583,0.0008981769,0.0009066277,0.0009142853,0.0009222669,0.0009289159,0.0009354767,0.0009419767,0.0009478155,0.0009538036,0.0009589328,0.0009637106,0.0009686238,0.0009738540,0.0009783829,0.0009824627,0.0009867156,0.0009909041,0.0009947179,0.0009986455,0.0010020386,0.0010053158,0.0010087451,0.0010123991,0.0010161565,0.0010193786,0.0010224474,0.0010243003,0.0010281951,0.0010301948,0.0010330942,0.0010362993,0.0010388814,0.0010408987,0.0010430015,0.0005431409,0.0005998314,0.0006408317,0.0006738062,0.0007010072,0.0007236881,0.0007430656,0.0007605174,0.0007760374,0.0007903811,0.0008030235,0.0008143226,0.0008244190,0.0008346535,0.0008438108,0.0008525065,0.0008604819,0.0008677872,0.0008753197,0.0008818033,0.0008879339,0.0008940700,0.0008996900,0.0009054633,0.0009102875,0.0009149808,0.0009194980,0.0009244902,0.0009288999,0.0009326000,0.0009367601,0.0009408287,0.0009443919,0.0009481391,0.0009513178,0.0009544989,0.0009577495,0.0009612678,0.0009647145,0.0009678129,0.0009708461,0.0009724834,0.0009763008,0.0009781490,0.0009809451,0.0009839662,0.0009864606,0.0009884575,0.0009905084,0.0005149864,0.0005687383,0.0006077606,0.0006390311,0.0006649356,0.0006866019,0.0007050595,0.0007215465,0.0007363455,0.0007499561,0.0007620239,0.0007728112,0.0007823901,0.0007921513,0.0008007538,0.0008090687,0.0008165741,0.0008235543,0.0008307707,0.0008371715,0.0008428974,0.0008485677,0.0008540169,0.0008595723,0.0008641650,0.0008684846,0.0008729391,0.0008776026,0.0008819097,0.0008854069,0.0008892893,0.0008932832,0.0008966540,0.0009002265,0.0009032041,0.0009062995,0.0009093528,0.0009127189,0.0009160752,0.0009189639,0.0009218007,0.0009235128,0.0009270296,0.0009287248,0.0009314437,0.0009343903,0.0009367513,0.0009385295,0.0009404843,0.0004882234,0.0005393377,0.0005763640,0.0006061095,0.0006306881,0.0006512844,0.0006687806,0.0006845633,0.0006988058,0.0007116140,0.0007231848,0.0007333458,0.0007424361,0.0007517260,0.0007600154,0.0007679234,0.0007749626,0.0007815724,0.0007885275,0.0007945196,0.0008000464,0.0008055558,0.0008105910,0.0008160293,0.0008202051,0.0008244252,0.0008287617,0.0008331350,0.0008371809,0.0008405811,0.0008443172,0.0008479802,0.0008513064,0.0008546298,0.0008576193,0.0008605908,0.0008634144,0.0008666444,0.0008698337,0.0008725555,0.0008753168,0.0008768942,0.0008801647,0.0008818871,0.0008843769,0.0008872668,0.0008896555,0.0008912178,0.0008930608,0.0004628785,0.0005113882,0.0005466153,0.0005749166,0.0005982376,0.0006178463,0.0006344885,0.0006494571,0.0006631031,0.0006752511,0.0006862287,0.0006958902,0.0007045411,0.0007133912,0.0007214252,0.0007288344,0.0007355035,0.0007417289,0.0007484640,0.0007541065,0.0007594957,0.0007647184,0.0007694511,0.0007746690,0.0007785847,0.0007826728,0.0007867768,0.0007909244,0.0007947869,0.0007980962,0.0008015503,0.0008050909,0.0008082578,0.0008114102,0.0008143198,0.0008170398,0.0008197705,0.0008227961,0.0008258212,0.0008285241,0.0008311118,0.0008326236,0.0008357614,0.0008374601,0.0008398370,0.0008425564,0.0008447132,0.0008463112,0.0008480463,0.0004388764,0.0004849174,0.0005183475,0.0005452617,0.0005675236,0.0005862285,0.0006019638,0.0006161856,0.0006291964,0.0006407177,0.0006512191,0.0006603954,0.0006686151,0.0006769768,0.0006845582,0.0006918216,0.0006981482,0.0007040126,0.0007104470,0.0007156950,0.0007209819,0.0007258622,0.0007303995,0.0007353361,0.0007392425,0.0007430123,0.0007470702,0.0007509057,0.0007545877,0.0007576844,0.0007609150,0.0007644933,0.0007674182,0.0007703584,0.0007731450,0.0007757980,0.0007784772,0.0007813058,0.0007840779,0.0007867074,0.0007892311,0.0007905067,0.0007935468,0.0007951502,0.0007975233,0.0007999761,0.0008021118,0.0008036770,0.0008054309,0.0004161080,0.0004598677,0.0004915962,0.0005171180,0.0005383387,0.0005560798,0.0005710629,0.0005846520,0.0005970683,0.0006081370,0.0006179074,0.0006267237,0.0006345647,0.0006424104,0.0006496945,0.0006566017,0.0006627063,0.0006682351,0.0006744137,0.0006794287,0.0006843953,0.0006890254,0.0006932891,0.0006980481,0.0007017997,0.0007054699,0.0007092630,0.0007128179,0.0007164905,0.0007193849,0.0007224582,0.0007258364,0.0007285775,0.0007314767,0.0007340764,0.0007366277,0.0007391541,0.0007418397,0.0007444530,0.0007470051,0.0007493490,0.0007505962,0.0007535597,0.0007551383,0.0007574452,0.0007596076,0.0007617477,0.0007631770,0.0007648572,0.0003945285,0.0004359956,0.0004662214,0.0004905186,0.0005106931,0.0005274322,0.0005418458,0.0005547313,0.0005665441,0.0005771202,0.0005864405,0.0005947801,0.0006022838,0.0006096595,0.0006166130,0.0006232082,0.0006290480,0.0006343312,0.0006401269,0.0006449292,0.0006497459,0.0006541277,0.0006580990,0.0006625791,0.0006662345,0.0006697930,0.0006733804,0.0006766876,0.0006802572,0.0006830027,0.0006858786,0.0006891488,0.0006917524,0.0006945980,0.0006970351,0.0006994324,0.0007017823,0.0007044598,0.0007067469,0.0007092928,0.0007115445,0.0007127012,0.0007155537,0.0007169992,0.0007192937,0.0007213385,0.0007233870,0.0007247077,0.0007262148,0.0003740573,0.0004135004,0.0004422415,0.0004653038,0.0004844167,0.0005004170,0.0005141335,0.0005263517,0.0005376298,0.0005476543,0.0005565002,0.0005644348,0.0005716866,0.0005786364,0.0005852980,0.0005914862,0.0005971567,0.0006021209,0.0006076470,0.0006122912,0.0006167394,0.0006210069,0.0006247378,0.0006289174,0.0006325630,0.0006358359,0.0006392380,0.0006424531,0.0006458797,0.0006484335,0.0006511214,0.0006544162,0.0006567702,0.0006594664,0.0006617691,0.0006641965,0.0006664085,0.0006688614,0.0006710492,0.0006735412,0.0006756313,0.0006767342,0.0006794221,0.0006808802,0.0006831210,0.0006849511,0.0006869656,0.0006882669,0.0006896758,0.0003546812,0.0003920611,0.0004194517,0.0004414166,0.0004595043,0.0004747969,0.0004877935,0.0004994634,0.0005101472,0.0005197929,0.0005280561,0.0005356441,0.0005425728,0.0005491952,0.0005555137,0.0005614760,0.0005668684,0.0005715585,0.0005767954,0.0005811662,0.0005854493,0.0005895203,0.0005932066,0.0005969848,0.0006005183,0.0006036774,0.0006069951,0.0006099481,0.0006132985,0.0006157548,0.0006181323,0.0006213374,0.0006235567,0.0006261549,0.0006283715,0.0006305979,0.0006327516,0.0006350039,0.0006370526,0.0006395573,0.0006415996,0.0006425450,0.0006452606,0.0006464877,0.0006486985,0.0006504538,0.0006522788,0.0006536321,0.0006549361,0.0003363293,0.0003718637,0.0003979146,0.0004186784,0.0004359237,0.0004504622,0.0004628591,0.0004738773,0.0004841370,0.0004932274,0.0005011577,0.0005083344,0.0005149067,0.0005212553,0.0005273388,0.0005329839,0.0005381053,0.0005425687,0.0005475994,0.0005516972,0.0005557281,0.0005596634,0.0005631309,0.0005668431,0.0005700464,0.0005731840,0.0005763657,0.0005791697,0.0005822750,0.0005846388,0.0005868989,0.0005898092,0.0005921284,0.0005944975,0.0005966488,0.0005988664,0.0006007989,0.0006030334,0.0006048829,0.0006073207,0.0006091892,0.0006101124,0.0006127978,0.0006139038,0.0006160622,0.0006177095,0.0006194754,0.0006207218,0.0006219287,0.0003189975,0.0003526135,0.0003774335,0.0003971470,0.0004135670,0.0004274315,0.0004391712,0.0004495384,0.0004593956,0.0004679917,0.0004756787,0.0004824878,0.0004886575,0.0004946895,0.0005005145,0.0005059165,0.0005107330,0.0005150467,0.0005198849,0.0005236722,0.0005275918,0.0005313205,0.0005345272,0.0005382284,0.0005411729,0.0005441523,0.0005471798,0.0005498822,0.0005529160,0.0005551272,0.0005572803,0.0005599500,0.0005621702,0.0005645293,0.0005665800,0.0005686083,0.0005704893,0.0005726033,0.0005744278,0.0005767406,0.0005784877,0.0005793087,0.0005818754,0.0005830423,0.0005851094,0.0005866126,0.0005882703,0.0005895173,0.0005906540,0.0003024613,0.0003343787,0.0003579752,0.0003767378,0.0003923087,0.0004055313,0.0004167707,0.0004265943,0.0004360211,0.0004440644,0.0004514241,0.0004579840,0.0004638207,0.0004695224,0.0004751389,0.0004801641,0.0004846847,0.0004888984,0.0004935605,0.0004971255,0.0005008600,0.0005043763,0.0005074978,0.0005110085,0.0005137972,0.0005165792,0.0005194872,0.0005221770,0.0005250146,0.0005270585,0.0005291526,0.0005316748,0.0005338041,0.0005360111,0.0005378563,0.0005399668,0.0005417875,0.0005438424,0.0005455216,0.0005476591,0.0005493672,0.0005501585,0.0005525446,0.0005536680,0.0005555738,0.0005571775,0.0005585599,0.0005598100,0.0005608602,0.0002868267,0.0003171544,0.0003395918,0.0003573861,0.0003721841,0.0003847637,0.0003954826,0.0004048732,0.0004138368,0.0004214051,0.0004285048,0.0004346287,0.0004401946,0.0004456235,0.0004509210,0.0004558344,0.0004600222,0.0004641071,0.0004686206,0.0004719676,0.0004754222,0.0004788347,0.0004817406,0.0004851382,0.0004877620,0.0004905089,0.0004931518,0.0004957914,0.0004985207,0.0005004655,0.0005024751,0.0005048090,0.0005068482,0.0005089958,0.0005106890,0.0005127773,0.0005144267,0.0005164014,0.0005180705,0.0005200246,0.0005216861,0.0005224660,0.0005248014,0.0005257983,0.0005275837,0.0005290574,0.0005304060,0.0005315883,0.0005325689,0.0002720905,0.0003007761,0.0003220830,0.0003390039,0.0003531055,0.0003650979,0.0003752734,0.0003841616,0.0003926333,0.0003998376,0.0004066997,0.0004126011,0.0004177811,0.0004229234,0.0004279613,0.0004326836,0.0004366753,0.0004405257,0.0004448682,0.0004480706,0.0004514231,0.0004546378,0.0004572713,0.0004605997,0.0004631063,0.0004657011,0.0004682455,0.0004707003,0.0004733403,0.0004751212,0.0004769925,0.0004792258,0.0004812099,0.0004832827,0.0004848922,0.0004868675,0.0004885317,0.0004904498,0.0004919811,0.0004938509,0.0004954639,0.0004961060,0.0004984593,0.0004993255,0.0005010673,0.0005024029,0.0005037390,0.0005048006,0.0005057723,0.0002580657,0.0002853334,0.0003055056,0.0003215223,0.0003350156,0.0003464599,0.0003561193,0.0003644902,0.0003726248,0.0003794107,0.0003860375,0.0003916587,0.0003965569,0.0004013975,0.0004061276,0.0004107473,0.0004145584,0.0004182539,0.0004224054,0.0004253161,0.0004285272,0.0004316214,0.0004341175,0.0004373214,0.0004396789,0.0004421001,0.0004446183,0.0004469362,0.0004493414,0.0004510467,0.0004529391,0.0004550287,0.0004569411,0.0004588909,0.0004604780,0.0004622717,0.0004639214,0.0004656015,0.0004671953,0.0004689674,0.0004705463,0.0004711571,0.0004733128,0.0004741897,0.0004757828,0.0004771232,0.0004784051,0.0004794094,0.0004803085,0.0002447323,0.0002706912,0.0002897497,0.0003050530,0.0003178847,0.0003287830,0.0003379673,0.0003459013,0.0003536454,0.0003600977,0.0003663728,0.0003717559,0.0003763975,0.0003809507,0.0003854737,0.0003898238,0.0003935062,0.0003970007,0.0004009540,0.0004039119,0.0004068700,0.0004098121,0.0004121707,0.0004152654,0.0004174903,0.0004196894,0.0004221692,0.0004244467,0.0004266058,0.0004282624,0.0004299958,0.0004321461,0.0004340033,0.0004356903,0.0004371899,0.0004390173,0.0004405657,0.0004421121,0.0004437831,0.0004453594,0.0004468593,0.0004474967,0.0004494231,0.0004502571,0.0004518305,0.0004531571,0.0004543309,0.0004553280,0.0004561205,0.0002320374,0.0002567102,0.0002748869,0.0002894073,0.0003015969,0.0003119938,0.0003207093,0.0003282365,0.0003355105,0.0003418066,0.0003477147,0.0003528379,0.0003572866,0.0003615781,0.0003659190,0.0003700898,0.0003735215,0.0003768214,0.0003805845,0.0003833947,0.0003862886,0.0003890368,0.0003912935,0.0003942701,0.0003963534,0.0003984101,0.0004008799,0.0004029506,0.0004051680,0.0004066366,0.0004083564,0.0004103901,0.0004120811,0.0004136853,0.0004152237,0.0004168253,0.0004183222,0.0004199236,0.0004214099,0.0004229480,0.0004243391,0.0004249901,0.0004267795,0.0004275334,0.0004291067,0.0004303063,0.0004314844,0.0004324546,0.0004331968,0.0002200850,0.0002434460,0.0002607403,0.0002745862,0.0002861203,0.0002959905,0.0003043163,0.0003114774,0.0003183126,0.0003244230,0.0003300624,0.0003348602,0.0003391428,0.0003432222,0.0003473954,0.0003513014,0.0003546685,0.0003577083,0.0003613450,0.0003640163,0.0003667466,0.0003693137,0.0003715246,0.0003743225,0.0003763421,0.0003783454,0.0003806285,0.0003826704,0.0003846967,0.0003861463,0.0003877756,0.0003897195,0.0003912899,0.0003928436,0.0003943103,0.0003957915,0.0003971620,0.0003988521,0.0004001826,0.0004017166,0.0004028745,0.0004035185,0.0004052107,0.0004059522,0.0004074694,0.0004086615,0.0004098024,0.0004107652,0.0004113832,0.0002087578,0.0002309532,0.0002473614,0.0002604811,0.0002714540,0.0002807910,0.0002887631,0.0002955859,0.0003021356,0.0003078851,0.0003133358,0.0003178940,0.0003218764,0.0003258729,0.0003297065,0.0003335193,0.0003367172,0.0003395902,0.0003430381,0.0003456660,0.0003482283,0.0003506273,0.0003527464,0.0003553947,0.0003573257,0.0003592757,0.0003614893,0.0003633535,0.0003652776,0.0003667374,0.0003681828,0.0003700409,0.0003715603,0.0003730197,0.0003744901,0.0003758243,0.0003770968,0.0003788163,0.0003800700,0.0003814705,0.0003826104,0.0003831236,0.0003848696,0.0003854342,0.0003869812,0.0003881799,0.0003893011,0.0003901282,0.0003907278,0.0001979829,0.0002190216,0.0002346906,0.0002471006,0.0002575211,0.0002664329,0.0002740115,0.0002805305,0.0002867448,0.0002921849,0.0002973490,0.0003017943,0.0003055213,0.0003093277,0.0003129488,0.0003166026,0.0003196562,0.0003223481,0.0003256808,0.0003281848,0.0003306386,0.0003329797,0.0003348948,0.0003374652,0.0003392749,0.0003411351,0.0003432075,0.0003450181,0.0003468676,0.0003481792,0.0003495793,0.0003514197,0.0003527776,0.0003541740,0.0003556305,0.0003568910,0.0003581836,0.0003597124,0.0003609917,0.0003622451,0.0003633346,0.0003638897,0.0003655615,0.0003660451,0.0003674641,0.0003686685,0.0003696619,0.0003705935,0.0003711061,0.0001877728,0.0002077078,0.0002226446,0.0002344802,0.0002443651,0.0002528296,0.0002599882,0.0002662578,0.0002721329,0.0002772754,0.0002822366,0.0002864641,0.0002899875,0.0002936146,0.0002970520,0.0003005574,0.0003034552,0.0003060419,0.0003091479,0.0003116326,0.0003138961,0.0003161351,0.0003179264,0.0003204155,0.0003221626,0.0003239251,0.0003258939,0.0003276014,0.0003293154,0.0003306375,0.0003320221,0.0003337010,0.0003350468,0.0003363292,0.0003377077,0.0003388898,0.0003400810,0.0003415977,0.0003428446,0.0003440544,0.0003450111,0.0003456822,0.0003472002,0.0003476930,0.0003489557,0.0003501719,0.0003511445,0.0003519930,0.0003524543,0.0001780962,0.0001970273,0.0002112867,0.0002224432,0.0002318148,0.0002399521,0.0002467587,0.0002526215,0.0002582632,0.0002632451,0.0002678753,0.0002719504,0.0002752497,0.0002787266,0.0002819767,0.0002854071,0.0002880619,0.0002905686,0.0002934935,0.0002957944,0.0002980004,0.0003001623,0.0003017720,0.0003041915,0.0003059799,0.0003075023,0.0003094626,0.0003111048,0.0003127288,0.0003139546,0.0003152087,0.0003169576,0.0003181482,0.0003193317,0.0003206158,0.0003219152,0.0003229915,0.0003243431,0.0003256243,0.0003267589,0.0003276260,0.0003283186,0.0003297726,0.0003301998,0.0003313867,0.0003325793,0.0003335177,0.0003343011,0.0003347438,0.0001689266,0.0001869081,0.0002004663,0.0002111225,0.0002199673,0.0002276937,0.0002341450,0.0002396819,0.0002450881,0.0002497944,0.0002542033,0.0002581657,0.0002612884,0.0002645583,0.0002676478,0.0002709341,0.0002735042,0.0002758644,0.0002786537,0.0002807278,0.0002829223,0.0002850096,0.0002865657,0.0002887669,0.0002906219,0.0002919718,0.0002939345,0.0002953794,0.0002968892,0.0002980981,0.0002992311,0.0003009683,0.0003021589,0.0003032219,0.0003044493,0.0003056332,0.0003066796,0.0003080731,0.0003092795,0.0003102975,0.0003111382,0.0003118274,0.0003132452,0.0003135891,0.0003147117,0.0003158182,0.0003167255,0.0003175286,0.0003179206,0.0001602549,0.0001772800,0.0001901974,0.0002003229,0.0002087869,0.0002160840,0.0002221747,0.0002275392,0.0002326204,0.0002371016,0.0002413075,0.0002450588,0.0002480654,0.0002510809,0.0002540694,0.0002572095,0.0002595627,0.0002619423,0.0002645958,0.0002665613,0.0002686521,0.0002706818,0.0002721735,0.0002741839,0.0002759191,0.0002772099,0.0002790837,0.0002805082,0.0002819231,0.0002830267,0.0002841274,0.0002857501,0.0002869335,0.0002879222,0.0002892042,0.0002902024,0.0002912904,0.0002925393,0.0002937459,0.0002947058,0.0002954533,0.0002961646,0.0002975507,0.0002978207,0.0002988846,0.0002999618,0.0003008460,0.0003015336,0.0003020072,0.0001519991,0.0001681744,0.0001804651,0.0001900918,0.0001981450,0.0002050458,0.0002108311,0.0002159638,0.0002207367,0.0002250099,0.0002290647,0.0002326013,0.0002354368,0.0002383898,0.0002411977,0.0002442338,0.0002464010,0.0002486280,0.0002511857,0.0002530710,0.0002550588,0.0002569851,0.0002583169,0.0002603892,0.0002619772,0.0002633226,0.0002650038,0.0002663063,0.0002677241,0.0002687862,0.0002697584,0.0002713419,0.0002724687,0.0002733984,0.0002746103,0.0002756218,0.0002766274,0.0002778425,0.0002789954,0.0002799442,0.0002805960,0.0002812879,0.0002825898,0.0002828375,0.0002839166,0.0002848517,0.0002856837,0.0002864027,0.0002868414,0.0001441564,0.0001595328,0.0001712263,0.0001803568,0.0001880180,0.0001945482,0.0002000628,0.0002049255,0.0002095272,0.0002136331,0.0002174414,0.0002208025,0.0002234562,0.0002262822,0.0002289939,0.0002318578,0.0002339135,0.0002360327,0.0002384822,0.0002402866,0.0002421682,0.0002439420,0.0002452461,0.0002472305,0.0002487599,0.0002500838,0.0002516831,0.0002529180,0.0002542656,0.0002552711,0.0002562133,0.0002576394,0.0002588054,0.0002596368,0.0002607937,0.0002618073,0.0002627658,0.0002638408,0.0002649765,0.0002659406,0.0002664606,0.0002671584,0.0002683915,0.0002685520,0.0002696444,0.0002706153,0.0002713344,0.0002720321,0.0002723807,0.0001367176,0.0001513433,0.0001624783,0.0001711499,0.0001783955,0.0001846751,0.0001898681,0.0001944627,0.0001988198,0.0002027637,0.0002063930,0.0002096096,0.0002121213,0.0002148666,0.0002174059,0.0002200778,0.0002220241,0.0002240929,0.0002264378,0.0002281334,0.0002299298,0.0002316090,0.0002329002,0.0002348559,0.0002361922,0.0002375001,0.0002390034,0.0002401591,0.0002414054,0.0002424241,0.0002432681,0.0002446471,0.0002458312,0.0002466191,0.0002477141,0.0002486537,0.0002495268,0.0002506173,0.0002516085,0.0002525924,0.0002530691,0.0002537689,0.0002549243,0.0002550295,0.0002560623,0.0002570636,0.0002576917,0.0002583893,0.0002587045,0.0001296741,0.0001435747,0.0001541891,0.0001624039,0.0001693020,0.0001752354,0.0001801623,0.0001846072,0.0001887032,0.0001924570,0.0001958940,0.0001990055,0.0002013534,0.0002040282,0.0002064013,0.0002089801,0.0002108165,0.0002127331,0.0002150459,0.0002166742,0.0002183324,0.0002198564,0.0002211850,0.0002229738,0.0002242436,0.0002254737,0.0002268910,0.0002280572,0.0002292446,0.0002302359,0.0002310402,0.0002322430,0.0002334918,0.0002342368,0.0002352369,0.0002361433,0.0002369876,0.0002380120,0.0002389460,0.0002398908,0.0002403774,0.0002410803,0.0002420615,0.0002422921,0.0002431763,0.0002441569,0.0002447407,0.0002453752,0.0002456905,0.0001229729,0.0001362164,0.0001462771,0.0001541107,0.0001606845,0.0001662914,0.0001709835,0.0001752287,0.0001791394,0.0001826568,0.0001859675,0.0001888974,0.0001911281,0.0001937641,0.0001959124,0.0001983561,0.0002001884,0.0002019761,0.0002041949,0.0002057511,0.0002072403,0.0002087642,0.0002099903,0.0002117049,0.0002128978,0.0002140768,0.0002154296,0.0002165943,0.0002176994,0.0002186274,0.0002194234,0.0002205234,0.0002217503,0.0002224712,0.0002233613,0.0002243189,0.0002251266,0.0002260865,0.0002269175,0.0002278298,0.0002283237,0.0002289823,0.0002298697,0.0002301221,0.0002309556,0.0002318793,0.0002324515,0.0002330889,0.0002333652,0.0001166300,0.0001292127,0.0001388026,0.0001462034,0.0001524922,0.0001578160,0.0001622352,0.0001663180,0.0001699851,0.0001733619,0.0001765168,0.0001793255,0.0001814457,0.0001839707,0.0001859859,0.0001883163,0.0001900654,0.0001917935,0.0001939088,0.0001953990,0.0001967179,0.0001982223,0.0001993856,0.0002009947,0.0002021425,0.0002033231,0.0002045891,0.0002057424,0.0002067016,0.0002076117,0.0002083262,0.0002093822,0.0002105769,0.0002112783,0.0002120766,0.0002130193,0.0002137792,0.0002147196,0.0002155067,0.0002163682,0.0002168055,0.0002175400,0.0002183393,0.0002185969,0.0002193855,0.0002202322,0.0002207536,0.0002213804,0.0002216405,0.0001106322,0.0001225831,0.0001316757,0.0001386987,0.0001447386,0.0001497805,0.0001539454,0.0001578421,0.0001613601,0.0001645289,0.0001675931,0.0001702001,0.0001722383,0.0001746203,0.0001764938,0.0001787650,0.0001804420,0.0001820501,0.0001840846,0.0001854759,0.0001867410,0.0001881739,0.0001892989,0.0001908546,0.0001920006,0.0001930620,0.0001942573,0.0001954007,0.0001963003,0.0001971726,0.0001978167,0.0001988539,0.0002000392,0.0002007214,0.0002014462,0.0002023113,0.0002030729,0.0002039692,0.0002047027,0.0002055256,0.0002059322,0.0002066955,0.0002073728,0.0002076360,0.0002083572,0.0002092046,0.0002096779,0.0002102488,0.0002105374,0.0001049568,0.0001163164,0.0001248965,0.0001316130,0.0001373476,0.0001421532,0.0001461164,0.0001498005,0.0001531941,0.0001561708,0.0001590898,0.0001615431,0.0001635221,0.0001657933,0.0001675616,0.0001696927,0.0001712883,0.0001727673,0.0001747577,0.0001761099,0.0001772932,0.0001786907,0.0001797217,0.0001812041,0.0001823607,0.0001833326,0.0001844490,0.0001855764,0.0001864779,0.0001872467,0.0001878499,0.0001888484,0.0001899745,0.0001906466,0.0001913117,0.0001921365,0.0001928402,0.0001937542,0.0001943959,0.0001951916,0.0001955638,0.0001963406,0.0001970052,0.0001972175,0.0001978816,0.0001987240,0.0001992063,0.0001997546,0.0002000067,0.0000995737,0.0001103031,0.0001184708,0.0001249273,0.0001303456,0.0001349300,0.0001386946,0.0001422061,0.0001454107,0.0001482610,0.0001510307,0.0001533303,0.0001552660,0.0001573868,0.0001591191,0.0001611121,0.0001625959,0.0001640304,0.0001659395,0.0001672289,0.0001683635,0.0001696871,0.0001706073,0.0001720989,0.0001731894,0.0001740470,0.0001751598,0.0001762828,0.0001771009,0.0001778145,0.0001784087,0.0001793695,0.0001804334,0.0001810427,0.0001817260,0.0001825130,0.0001831304,0.0001840329,0.0001846152,0.0001854122,0.0001857232,0.0001865114,0.0001871407,0.0001872756,0.0001879131,0.0001887626,0.0001892042,0.0001897159,0.0001900195,0.0000944527,0.0001046267,0.0001123926,0.0001185418,0.0001237348,0.0001280494,0.0001315755,0.0001349625,0.0001380409,0.0001407750,0.0001434218,0.0001455532,0.0001473697,0.0001494488,0.0001511031,0.0001530159,0.0001543886,0.0001557314,0.0001575971,0.0001587400,0.0001598351,0.0001611211,0.0001619450,0.0001634711,0.0001644646,0.0001652940,0.0001663151,0.0001673267,0.0001682432,0.0001688471,0.0001694502,0.0001703418,0.0001714291,0.0001719782,0.0001725847,0.0001733469,0.0001739601,0.0001747941,0.0001753276,0.0001760954,0.0001763695,0.0001771264,0.0001776802,0.0001779048,0.0001784673,0.0001792807,0.0001796662,0.0001802315,0.0001804621,0.0000896004,0.0000992859,0.0001066565,0.0001125248,0.0001174488,0.0001215445,0.0001249081,0.0001280927,0.0001310194,0.0001336172,0.0001361284,0.0001382058,0.0001399041,0.0001419475,0.0001434733,0.0001452650,0.0001465589,0.0001478722,0.0001496753,0.0001507498,0.0001517897,0.0001529642,0.0001537710,0.0001552802,0.0001561726,0.0001569578,0.0001579844,0.0001589343,0.0001597713,0.0001603671,0.0001608959,0.0001617690,0.0001628840,0.0001633124,0.0001639431,0.0001645654,0.0001652192,0.0001660781,0.0001665414,0.0001673180,0.0001675429,0.0001682068,0.0001687583,0.0001690050,0.0001695146,0.0001703315,0.0001706521,0.0001711996,0.0001713649,0.0000849998,0.0000941834,0.0001012210,0.0001067262,0.0001114737,0.0001153824,0.0001185485,0.0001215845,0.0001243631,0.0001268183,0.0001292651,0.0001312835,0.0001327772,0.0001347725,0.0001362512,0.0001378830,0.0001391602,0.0001404654,0.0001421544,0.0001431832,0.0001441431,0.0001452362,0.0001460196,0.0001474509,0.0001483526,0.0001490683,0.0001500463,0.0001509480,0.0001517479,0.0001522353,0.0001528073,0.0001535957,0.0001547176,0.0001551355,0.0001556880,0.0001563024,0.0001569101,0.0001576848,0.0001581769,0.0001588897,0.0001591187,0.0001597112,0.0001602954,0.0001605282,0.0001610383,0.0001618422,0.0001620856,0.0001626655,0.0001627998,0.0000806601,0.0000893911,0.0000960369,0.0001012471,0.0001057703,0.0001095049,0.0001125551,0.0001154272,0.0001180528,0.0001203503,0.0001226962,0.0001246230,0.0001260520,0.0001279500,0.0001293660,0.0001309020,0.0001321126,0.0001333175,0.0001350099,0.0001359914,0.0001368321,0.0001379364,0.0001386184,0.0001400702,0.0001408872,0.0001415566,0.0001425465,0.0001433344,0.0001441318,0.0001445674,0.0001450919,0.0001459024,0.0001469261,0.0001473564,0.0001478004,0.0001484462,0.0001490452,0.0001497713,0.0001502410,0.0001509282,0.0001510816,0.0001516528,0.0001522381,0.0001524811,0.0001530188,0.0001536667,0.0001539418,0.0001544930,0.0001546746,0.0000765263,0.0000847949,0.0000911425,0.0000960877,0.0001003527,0.0001039277,0.0001067941,0.0001095765,0.0001121054,0.0001142487,0.0001164664,0.0001182900,0.0001196480,0.0001214673,0.0001228237,0.0001242997,0.0001254365,0.0001265971,0.0001282173,0.0001291468,0.0001299319,0.0001310152,0.0001316474,0.0001330579,0.0001338036,0.0001344576,0.0001353845,0.0001361163,0.0001368824,0.0001372800,0.0001378344,0.0001385573,0.0001395613,0.0001399569,0.0001403627,0.0001409370,0.0001416111,0.0001422562,0.0001426747,0.0001433163,0.0001435086,0.0001440295,0.0001446459,0.0001447852,0.0001452884,0.0001459320,0.0001462047,0.0001467012,0.0001469368,0.0000726201,0.0000804959,0.0000864904,0.0000912080,0.0000952548,0.0000986276,0.0001013827,0.0001040068,0.0001064386,0.0001084884,0.0001105917,0.0001123180,0.0001135776,0.0001152857,0.0001166374,0.0001180267,0.0001190910,0.0001201881,0.0001217408,0.0001226210,0.0001234100,0.0001244455,0.0001250383,0.0001263752,0.0001270558,0.0001276949,0.0001285994,0.0001292558,0.0001299993,0.0001303729,0.0001309077,0.0001316028,0.0001325543,0.0001329324,0.0001333250,0.0001338520,0.0001344932,0.0001351196,0.0001355415,0.0001361550,0.0001363054,0.0001368436,0.0001373981,0.0001375412,0.0001380057,0.0001385967,0.0001388228,0.0001393462,0.0001396180,0.0000689114,0.0000763671,0.0000820807,0.0000865729,0.0000904089,0.0000935994,0.0000962539,0.0000987400,0.0001010327,0.0001030368,0.0001049790,0.0001066795,0.0001078497,0.0001094536,0.0001107404,0.0001120756,0.0001131006,0.0001140948,0.0001156017,0.0001164406,0.0001171775,0.0001181583,0.0001187487,0.0001199921,0.0001206236,0.0001212788,0.0001221375,0.0001227522,0.0001234624,0.0001238269,0.0001242844,0.0001249827,0.0001258969,0.0001262335,0.0001266135,0.0001271303,0.0001277264,0.0001283192,0.0001287547,0.0001293398,0.0001294936,0.0001299743,0.0001305292,0.0001306084,0.0001310818,0.0001316381,0.0001318744,0.0001323739,0.0001326338,0.0000653719,0.0000724980,0.0000779098,0.0000821739,0.0000858108,0.0000888339,0.0000913611,0.0000937420,0.0000959068,0.0000978022,0.0000996733,0.0001012633,0.0001023946,0.0001039387,0.0001051628,0.0001064314,0.0001074105,0.0001083477,0.0001097233,0.0001105482,0.0001112648,0.0001122076,0.0001127685,0.0001139757,0.0001145486,0.0001151868,0.0001160183,0.0001165741,0.0001172554,0.0001175665,0.0001180342,0.0001187091,0.0001195598,0.0001199005,0.0001202303,0.0001207724,0.0001212919,0.0001218544,0.0001222902,0.0001228514,0.0001229647,0.0001234818,0.0001239821,0.0001240483,0.0001245401,0.0001249939,0.0001252680,0.0001257187,0.0001259416,0.0000620026,0.0000687925,0.0000738968,0.0000780060,0.0000814521,0.0000843013,0.0000867588,0.0000889693,0.0000910803,0.0000928354,0.0000946251,0.0000961147,0.0000971943,0.0000986976,0.0000998736,0.0001010804,0.0001019787,0.0001029004,0.0001041734,0.0001049865,0.0001056745,0.0001065624,0.0001070885,0.0001082478,0.0001087554,0.0001094172,0.0001101670,0.0001107073,0.0001114125,0.0001116485,0.0001121020,0.0001127557,0.0001135586,0.0001138536,0.0001142232,0.0001147247,0.0001151983,0.0001157069,0.0001161383,0.0001166674,0.0001167980,0.0001172901,0.0001177693,0.0001178521,0.0001183067,0.0001187401,0.0001189655,0.0001194258,0.0001196101,0.0000588248,0.0000652810,0.0000701064,0.0000740120,0.0000773301,0.0000799936,0.0000823808,0.0000844746,0.0000864793,0.0000881394,0.0000898252,0.0000912549,0.0000922327,0.0000937665,0.0000948385,0.0000959755,0.0000968497,0.0000977174,0.0000988978,0.0000997414,0.0001003832,0.0001011903,0.0001016551,0.0001028072,0.0001032967,0.0001038938,0.0001046551,0.0001051179,0.0001058101,0.0001060709,0.0001064716,0.0001071066,0.0001078561,0.0001081315,0.0001084869,0.0001089561,0.0001093883,0.0001098848,0.0001103125,0.0001108039,0.0001108913,0.0001114240,0.0001118596,0.0001119553,0.0001123670,0.0001127424,0.0001130119,0.0001134635,0.0001136124,0.0000557925,0.0000619400,0.0000665047,0.0000702657,0.0000734100,0.0000759260,0.0000782419,0.0000801911,0.0000820800,0.0000836511,0.0000852830,0.0000866883,0.0000875514,0.0000890418,0.0000900422,0.0000911342,0.0000920005,0.0000927849,0.0000938949,0.0000947045,0.0000952994,0.0000960941,0.0000965220,0.0000976358,0.0000980715,0.0000986914,0.0000994035,0.0000998367,0.0001004899,0.0001007655,0.0001011222,0.0001017072,0.0001024612,0.0001026929,0.0001030379,0.0001034670,0.0001038905,0.0001043526,0.0001047802,0.0001052588,0.0001053027,0.0001058623,0.0001062753,0.0001063600,0.0001066929,0.0001070664,0.0001073225,0.0001078048,0.0001079138,0.0000529531,0.0000587811,0.0000631162,0.0000666838,0.0000696923,0.0000720988,0.0000742826,0.0000761491,0.0000779316,0.0000794240,0.0000809654,0.0000823464,0.0000831304,0.0000845332,0.0000854858,0.0000865328,0.0000873218,0.0000881531,0.0000891436,0.0000899025,0.0000905010,0.0000912388,0.0000916668,0.0000927260,0.0000931099,0.0000937182,0.0000943820,0.0000948306,0.0000954408,0.0000957194,0.0000960678,0.0000966269,0.0000973265,0.0000975502,0.0000978799,0.0000982719,0.0000986872,0.0000991212,0.0000994536,0.0000999497,0.0000999944,0.0001005526,0.0001009562,0.0001010879,0.0001013122,0.0001016992,0.0001019176,0.0001024209,0.0001024801,0.0000502409,0.0000558110,0.0000599254,0.0000632983,0.0000661555,0.0000684194,0.0000705107,0.0000722969,0.0000739570,0.0000754423,0.0000769224,0.0000782010,0.0000789477,0.0000802726,0.0000811990,0.0000821848,0.0000829452,0.0000837039,0.0000846397,0.0000853666,0.0000859333,0.0000866237,0.0000870528,0.0000880881,0.0000884394,0.0000889737,0.0000896113,0.0000900529,0.0000906578,0.0000908948,0.0000912320,0.0000918128,0.0000924416,0.0000926292,0.0000929755,0.0000933143,0.0000937430,0.0000941601,0.0000944301,0.0000949549,0.0000949820,0.0000954888,0.0000958823,0.0000960292,0.0000962394,0.0000965890,0.0000968207,0.0000973057,0.0000973555,0.0000476835,0.0000529600,0.0000568548,0.0000600883,0.0000627790,0.0000649434,0.0000669499,0.0000686313,0.0000702193,0.0000716386,0.0000729894,0.0000742395,0.0000749759,0.0000762389,0.0000771495,0.0000780161,0.0000787794,0.0000794740,0.0000804014,0.0000810756,0.0000815879,0.0000822244,0.0000826578,0.0000836340,0.0000839816,0.0000845155,0.0000851194,0.0000855395,0.0000861121,0.0000862901,0.0000866531,0.0000871778,0.0000878385,0.0000879731,0.0000883011,0.0000886267,0.0000890396,0.0000894439,0.0000897437,0.0000901935,0.0000901993,0.0000907145,0.0000910430,0.0000912214,0.0000913963,0.0000917359,0.0000919811,0.0000924338,0.0000924887,0.0000452297,0.0000502575,0.0000539889,0.0000570290,0.0000595946,0.0000616874,0.0000635359,0.0000651615,0.0000666515,0.0000679989,0.0000692528,0.0000705029,0.0000711841,0.0000723806,0.0000732785,0.0000740852,0.0000748100,0.0000754549,0.0000763389,0.0000770071,0.0000775038,0.0000780799,0.0000784969,0.0000794389,0.0000797802,0.0000802664,0.0000808195,0.0000812499,0.0000817933,0.0000819959,0.0000823037,0.0000827740,0.0000834258,0.0000835837,0.0000838839,0.0000841773,0.0000846192,0.0000849699,0.0000852656,0.0000856612,0.0000856755,0.0000861631,0.0000865059,0.0000866717,0.0000868016,0.0000871557,0.0000873221,0.0000878112,0.0000878007,0.0000429233,0.0000476943,0.0000512557,0.0000541171,0.0000565738,0.0000585831,0.0000603304,0.0000618612,0.0000632875,0.0000645597,0.0000657348,0.0000669643,0.0000675829,0.0000687587,0.0000695488,0.0000703544,0.0000710314,0.0000716470,0.0000725097,0.0000731116,0.0000736110,0.0000741351,0.0000745384,0.0000754300,0.0000757866,0.0000762515,0.0000767625,0.0000771568,0.0000776863,0.0000778880,0.0000781737,0.0000786286,0.0000792662,0.0000793974,0.0000796739,0.0000799472,0.0000803523,0.0000806839,0.0000809530,0.0000813685,0.0000813875,0.0000818587,0.0000822012,0.0000823276,0.0000824426,0.0000827953,0.0000829280,0.0000834360,0.0000834005,0.0000407264,0.0000452569,0.0000486493,0.0000513793,0.0000537178,0.0000556260,0.0000572819,0.0000587329,0.0000601168,0.0000612840,0.0000623980,0.0000635894,0.0000642042,0.0000652928,0.0000660293,0.0000668256,0.0000674896,0.0000680550,0.0000688915,0.0000694486,0.0000699248,0.0000704412,0.0000707903,0.0000716613,0.0000719913,0.0000723891,0.0000729117,0.0000732720,0.0000737913,0.0000739675,0.0000742524,0.0000746556,0.0000752944,0.0000753849,0.0000756815,0.0000759546,0.0000763183,0.0000766769,0.0000768790,0.0000773155,0.0000772945,0.0000777506,0.0000781049,0.0000782241,0.0000782878,0.0000786390,0.0000787692,0.0000792777,0.0000792510,0.0000386420,0.0000429495,0.0000461880,0.0000487919,0.0000509905,0.0000528352,0.0000543534,0.0000557472,0.0000570586,0.0000581827,0.0000592897,0.0000603806,0.0000609445,0.0000619914,0.0000627015,0.0000634322,0.0000640530,0.0000646146,0.0000654144,0.0000659356,0.0000664104,0.0000668946,0.0000672293,0.0000680705,0.0000683748,0.0000687842,0.0000692414,0.0000695769,0.0000700890,0.0000702438,0.0000705026,0.0000709283,0.0000715188,0.0000715929,0.0000718985,0.0000721408,0.0000725343,0.0000728299,0.0000729709,0.0000734238,0.0000734321,0.0000738672,0.0000741820,0.0000743196,0.0000743400,0.0000747384,0.0000748298,0.0000753037,0.0000752580,0.0000366909,0.0000407497,0.0000438433,0.0000463082,0.0000484092,0.0000501456,0.0000515964,0.0000529263,0.0000541965,0.0000552535,0.0000562741,0.0000573668,0.0000578568,0.0000588654,0.0000595568,0.0000602354,0.0000608304,0.0000613525,0.0000621594,0.0000626121,0.0000630466,0.0000635130,0.0000637977,0.0000646529,0.0000649426,0.0000653096,0.0000657627,0.0000660870,0.0000665701,0.0000667391,0.0000669629,0.0000673731,0.0000679211,0.0000680086,0.0000682764,0.0000684634,0.0000688792,0.0000691952,0.0000692917,0.0000697697,0.0000697349,0.0000701303,0.0000704696,0.0000706033,0.0000705692,0.0000709964,0.0000711015,0.0000715079,0.0000714674,0.0000348145,0.0000386766,0.0000416070,0.0000439406,0.0000459491,0.0000476028,0.0000489845,0.0000502425,0.0000514498,0.0000524586,0.0000534494,0.0000544700,0.0000549041,0.0000558987,0.0000565716,0.0000572342,0.0000577695,0.0000582663,0.0000590007,0.0000594508,0.0000598826,0.0000602911,0.0000605624,0.0000614342,0.0000616618,0.0000620218,0.0000624635,0.0000627711,0.0000632203,0.0000634113,0.0000636029,0.0000640130,0.0000644979,0.0000646252,0.0000648634,0.0000649689,0.0000654176,0.0000657289,0.0000658150,0.0000662763,0.0000662680,0.0000666414,0.0000669388,0.0000670626,0.0000670174,0.0000674481,0.0000675594,0.0000679207,0.0000679121,0.0000330207,0.0000367048,0.0000394964,0.0000416815,0.0000436078,0.0000451884,0.0000465247,0.0000477225,0.0000488825,0.0000498252,0.0000507526,0.0000517084,0.0000521273,0.0000530814,0.0000537121,0.0000543481,0.0000548535,0.0000552851,0.0000560317,0.0000564511,0.0000568560,0.0000572271,0.0000575303,0.0000583614,0.0000585322,0.0000589320,0.0000593343,0.0000596309,0.0000600465,0.0000602136,0.0000603955,0.0000608279,0.0000612921,0.0000613615,0.0000615943,0.0000617126,0.0000621404,0.0000624425,0.0000625103,0.0000629619,0.0000629508,0.0000633176,0.0000636270,0.0000637177,0.0000636532,0.0000640518,0.0000641555,0.0000645361,0.0000645064,0.0000313500,0.0000348522,0.0000374701,0.0000395596,0.0000413805,0.0000429021,0.0000441489,0.0000452993,0.0000464138,0.0000472872,0.0000482027,0.0000491139,0.0000494913,0.0000503774,0.0000509915,0.0000516088,0.0000521073,0.0000525170,0.0000531936,0.0000536156,0.0000540002,0.0000543603,0.0000546353,0.0000554143,0.0000555992,0.0000559568,0.0000563782,0.0000566219,0.0000570471,0.0000571990,0.0000573529,0.0000578087,0.0000581999,0.0000582732,0.0000584921,0.0000586465,0.0000589921,0.0000593241,0.0000593939,0.0000598020,0.0000598118,0.0000601446,0.0000604296,0.0000605240,0.0000604893,0.0000608284,0.0000609483,0.0000613054,0.0000612953,0.0000297534,0.0000330578,0.0000355456,0.0000375466,0.0000392681,0.0000407489,0.0000419361,0.0000430296,0.0000440374,0.0000448954,0.0000457717,0.0000466203,0.0000470128,0.0000478426,0.0000484156,0.0000490373,0.0000494994,0.0000498711,0.0000505252,0.0000509293,0.0000512814,0.0000516322,0.0000519074,0.0000526407,0.0000527889,0.0000531446,0.0000535502,0.0000537726,0.0000541497,0.0000543078,0.0000544985,0.0000548865,0.0000552795,0.0000553577,0.0000555753,0.0000557192,0.0000559747,0.0000563633,0.0000564085,0.0000568029,0.0000568124,0.0000571380,0.0000574207,0.0000575027,0.0000574977,0.0000577906,0.0000579069,0.0000582148,0.0000582270,0.0000282289,0.0000313723,0.0000337209,0.0000356479,0.0000372593,0.0000386611,0.0000398166,0.0000408546,0.0000418118,0.0000426175,0.0000434515,0.0000442837,0.0000446402,0.0000454351,0.0000459978,0.0000465492,0.0000470171,0.0000473935,0.0000480099,0.0000483890,0.0000487040,0.0000490287,0.0000492984,0.0000499908,0.0000501277,0.0000504651,0.0000508318,0.0000510859,0.0000514356,0.0000515582,0.0000517890,0.0000521418,0.0000525036,0.0000525846,0.0000527926,0.0000528976,0.0000531882,0.0000535325,0.0000535530,0.0000539803,0.0000539631,0.0000542939,0.0000545233,0.0000546163,0.0000546306,0.0000548934,0.0000549899,0.0000553187,0.0000553056,0.0000267820,0.0000297822,0.0000319932,0.0000338275,0.0000353664,0.0000366832,0.0000378320,0.0000388030,0.0000396647,0.0000404821,0.0000412458,0.0000420609,0.0000424215,0.0000431527,0.0000436796,0.0000442174,0.0000446351,0.0000450018,0.0000455818,0.0000459437,0.0000462665,0.0000465571,0.0000468081,0.0000474566,0.0000476133,0.0000479452,0.0000482920,0.0000485236,0.0000488846,0.0000489811,0.0000492150,0.0000495446,0.0000498257,0.0000499523,0.0000501475,0.0000502656,0.0000505388,0.0000508764,0.0000508884,0.0000512929,0.0000512782,0.0000515537,0.0000518171,0.0000518541,0.0000518435,0.0000521196,0.0000522258,0.0000525296,0.0000525684,0.0000254217,0.0000282468,0.0000303628,0.0000320980,0.0000335721,0.0000348285,0.0000359301,0.0000368349,0.0000376916,0.0000384534,0.0000391595,0.0000399338,0.0000402860,0.0000409697,0.0000414847,0.0000420071,0.0000423813,0.0000427711,0.0000432966,0.0000436306,0.0000439457,0.0000442458,0.0000444473,0.0000450552,0.0000452236,0.0000455292,0.0000458529,0.0000460855,0.0000463996,0.0000465101,0.0000467680,0.0000470817,0.0000473331,0.0000474439,0.0000476368,0.0000477510,0.0000480023,0.0000483381,0.0000483375,0.0000487027,0.0000487316,0.0000489471,0.0000492221,0.0000492795,0.0000492586,0.0000495091,0.0000495917,0.0000498917,0.0000499190,0.0000241174,0.0000268421,0.0000288154,0.0000304638,0.0000318585,0.0000330680,0.0000341250,0.0000349558,0.0000357789,0.0000365180,0.0000371870,0.0000379048,0.0000382355,0.0000389155,0.0000394075,0.0000398883,0.0000402658,0.0000406418,0.0000411126,0.0000414617,0.0000417184,0.0000420428,0.0000422291,0.0000427778,0.0000429535,0.0000432473,0.0000435557,0.0000437592,0.0000440608,0.0000441731,0.0000444287,0.0000447097,0.0000449651,0.0000450795,0.0000452331,0.0000453656,0.0000456050,0.0000459287,0.0000459197,0.0000462392,0.0000462821,0.0000465134,0.0000467407,0.0000468319,0.0000467905,0.0000470185,0.0000471160,0.0000473944,0.0000474106,0.0000228879,0.0000254879,0.0000273585,0.0000289336,0.0000302482,0.0000313646,0.0000324362,0.0000331965,0.0000339971,0.0000346720,0.0000353156,0.0000359896,0.0000363269,0.0000369624,0.0000374317,0.0000378772,0.0000382338,0.0000385889,0.0000390631,0.0000393613,0.0000396286,0.0000399191,0.0000401204,0.0000406299,0.0000407970,0.0000410552,0.0000413818,0.0000415418,0.0000418529,0.0000419614,0.0000421794,0.0000424648,0.0000426974,0.0000428165,0.0000429933,0.0000430795,0.0000433197,0.0000436195,0.0000436133,0.0000439164,0.0000439492,0.0000441667,0.0000443949,0.0000445112,0.0000444433,0.0000446775,0.0000447758,0.0000450002,0.0000450481,0.0000217350,0.0000241956,0.0000259780,0.0000274628,0.0000287095,0.0000297552,0.0000307834,0.0000315061,0.0000323090,0.0000329089,0.0000335498,0.0000341915,0.0000344982,0.0000350916,0.0000355409,0.0000359443,0.0000363202,0.0000366520,0.0000371384,0.0000373922,0.0000376326,0.0000379248,0.0000381065,0.0000385939,0.0000387485,0.0000390005,0.0000393257,0.0000394445,0.0000397412,0.0000398643,0.0000400613,0.0000403287,0.0000405767,0.0000406780,0.0000408198,0.0000409199,0.0000411695,0.0000414468,0.0000414530,0.0000417227,0.0000417583,0.0000419360,0.0000421705,0.0000422923,0.0000422265,0.0000424448,0.0000425252,0.0000427462,0.0000427956,0.0000206318,0.0000229685,0.0000246469,0.0000260447,0.0000272526,0.0000282445,0.0000292393,0.0000299201,0.0000306973,0.0000312486,0.0000318434,0.0000324783,0.0000327547,0.0000333230,0.0000337544,0.0000341440,0.0000345071,0.0000347996,0.0000352462,0.0000354938,0.0000357504,0.0000360140,0.0000361614,0.0000366636,0.0000368099,0.0000370416,0.0000373545,0.0000374886,0.0000377741,0.0000378666,0.0000380550,0.0000382939,0.0000385253,0.0000386647,0.0000387669,0.0000388768,0.0000391021,0.0000393676,0.0000393654,0.0000396382,0.0000396590,0.0000397961,0.0000400678,0.0000401829,0.0000401258,0.0000403107,0.0000404072,0.0000406014,0.0000406404,0.0000195997,0.0000217841,0.0000233959,0.0000247173,0.0000258612,0.0000268237,0.0000277766,0.0000283761,0.0000291476,0.0000296765,0.0000302186,0.0000308355,0.0000311008,0.0000316391,0.0000320542,0.0000324502,0.0000327919,0.0000330467,0.0000334745,0.0000337100,0.0000339260,0.0000342080,0.0000343555,0.0000348341,0.0000349585,0.0000351868,0.0000354796,0.0000356069,0.0000358730,0.0000359576,0.0000361401,0.0000363688,0.0000365848,0.0000367243,0.0000368392,0.0000369306,0.0000371320,0.0000373673,0.0000373829,0.0000376345,0.0000376570,0.0000378082,0.0000380692,0.0000381712,0.0000381190,0.0000383180,0.0000383845,0.0000385934,0.0000386035,0.0000186079,0.0000206753,0.0000222038,0.0000234506,0.0000245349,0.0000254592,0.0000263803,0.0000269292,0.0000276783,0.0000281850,0.0000287080,0.0000292664,0.0000295279,0.0000300423,0.0000304279,0.0000308131,0.0000311325,0.0000313805,0.0000317876,0.0000320231,0.0000322103,0.0000324983,0.0000326297,0.0000330860,0.0000332082,0.0000334207,0.0000337060,0.0000338117,0.0000340739,0.0000341487,0.0000343311,0.0000345552,0.0000347957,0.0000348713,0.0000349987,0.0000350871,0.0000352563,0.0000355039,0.0000355201,0.0000357585,0.0000357700,0.0000359116,0.0000361844,0.0000362664,0.0000362097,0.0000363969,0.0000364854,0.0000366615,0.0000366670,0.0000176472,0.0000196056,0.0000210625,0.0000222549,0.0000232873,0.0000241826,0.0000250408,0.0000255637,0.0000262810,0.0000267454,0.0000272516,0.0000277830,0.0000280331,0.0000285199,0.0000289144,0.0000292690,0.0000295747,0.0000298247,0.0000301981,0.0000303943,0.0000305920,0.0000308795,0.0000309864,0.0000314514,0.0000315356,0.0000317495,0.0000320419,0.0000321399,0.0000323748,0.0000324696,0.0000326062,0.0000328374,0.0000330656,0.0000331194,0.0000332383,0.0000333522,0.0000334879,0.0000337283,0.0000337408,0.0000339720,0.0000339669,0.0000341309,0.0000343717,0.0000344569,0.0000344076,0.0000345699,0.0000346741,0.0000348267,0.0000348349,0.0000167510,0.0000186134,0.0000199746,0.0000211481,0.0000221171,0.0000229788,0.0000237846,0.0000242903,0.0000249564,0.0000253952,0.0000258766,0.0000263842,0.0000266449,0.0000271046,0.0000274767,0.0000277803,0.0000280781,0.0000283087,0.0000286976,0.0000288741,0.0000290516,0.0000293427,0.0000294579,0.0000298688,0.0000299569,0.0000301428,0.0000304418,0.0000305383,0.0000307525,0.0000308282,0.0000309737,0.0000311893,0.0000314090,0.0000314616,0.0000315800,0.0000317065,0.0000318210,0.0000320129,0.0000320655,0.0000322878,0.0000322527,0.0000324053,0.0000326733,0.0000327354,0.0000326901,0.0000328591,0.0000329281,0.0000330885,0.0000331017,0.0000159023,0.0000176469,0.0000189642,0.0000200807,0.0000209886,0.0000218204,0.0000225822,0.0000230370,0.0000236916,0.0000241116,0.0000245843,0.0000250502,0.0000252950,0.0000257376,0.0000260974,0.0000263879,0.0000266701,0.0000268863,0.0000272701,0.0000274162,0.0000276220,0.0000278740,0.0000279698,0.0000283528,0.0000284566,0.0000286227,0.0000289056,0.0000290024,0.0000292154,0.0000292916,0.0000294175,0.0000296578,0.0000298495,0.0000298784,0.0000299854,0.0000301367,0.0000302223,0.0000303906,0.0000304678,0.0000306899,0.0000306522,0.0000307923,0.0000310401,0.0000311080,0.0000310335,0.0000311976,0.0000312897,0.0000314318,0.0000314719,0.0000151024,0.0000167436,0.0000179896,0.0000190812,0.0000199305,0.0000207314,0.0000214406,0.0000218633,0.0000225172,0.0000229141,0.0000233460,0.0000237674,0.0000240044,0.0000244514,0.0000247726,0.0000250507,0.0000253241,0.0000255421,0.0000258942,0.0000260548,0.0000262282,0.0000264648,0.0000265641,0.0000269187,0.0000270492,0.0000272039,0.0000274649,0.0000275448,0.0000277589,0.0000278260,0.0000279421,0.0000281709,0.0000283539,0.0000283732,0.0000284744,0.0000286234,0.0000287181,0.0000288795,0.0000289566,0.0000291631,0.0000291063,0.0000292523,0.0000295057,0.0000295399,0.0000294818,0.0000296565,0.0000297169,0.0000298681,0.0000298965,0.0000143285,0.0000159022,0.0000170788,0.0000181165,0.0000189159,0.0000196624,0.0000203791,0.0000207598,0.0000213822,0.0000217398,0.0000221702,0.0000225737,0.0000228064,0.0000232032,0.0000235358,0.0000237906,0.0000240450,0.0000242743,0.0000246134,0.0000247543,0.0000249128,0.0000251418,0.0000252480,0.0000255663,0.0000257135,0.0000258362,0.0000260981,0.0000261693,0.0000263533,0.0000264289,0.0000265520,0.0000267501,0.0000269312,0.0000269720,0.0000270502,0.0000271999,0.0000272856,0.0000274268,0.0000275180,0.0000276986,0.0000276712,0.0000278020,0.0000280261,0.0000280642,0.0000280078,0.0000281576,0.0000282314,0.0000283753,0.0000284060,0.0000135990,0.0000151031,0.0000162197,0.0000171980,0.0000179625,0.0000186664,0.0000193377,0.0000197209,0.0000203075,0.0000206411,0.0000210589,0.0000214310,0.0000216588,0.0000220442,0.0000223444,0.0000226128,0.0000228616,0.0000230607,0.0000233922,0.0000235160,0.0000236454,0.0000238689,0.0000239588,0.0000242822,0.0000244403,0.0000245272,0.0000247902,0.0000248544,0.0000250413,0.0000251070,0.0000252193,0.0000254093,0.0000255909,0.0000256125,0.0000256706,0.0000258629,0.0000258964,0.0000260617,0.0000261179,0.0000263061,0.0000262946,0.0000264112,0.0000266184,0.0000266553,0.0000266191,0.0000267407,0.0000268014,0.0000269627,0.0000269885,0.0000128939,0.0000143193,0.0000154012,0.0000163294,0.0000170467,0.0000177124,0.0000183489,0.0000187292,0.0000192807,0.0000195954,0.0000199833,0.0000203543,0.0000205673,0.0000209288,0.0000212087,0.0000214669,0.0000217063,0.0000219073,0.0000222251,0.0000223463,0.0000224683,0.0000226811,0.0000227578,0.0000230596,0.0000232154,0.0000232850,0.0000235142,0.0000235950,0.0000237810,0.0000238489,0.0000239631,0.0000241294,0.0000243011,0.0000243384,0.0000243871,0.0000245765,0.0000246116,0.0000247432,0.0000248192,0.0000250057,0.0000249664,0.0000250658,0.0000252952,0.0000253264,0.0000252974,0.0000254067,0.0000254445,0.0000256062,0.0000256447,0.0000122408,0.0000135913,0.0000146212,0.0000154920,0.0000162023,0.0000168273,0.0000174301,0.0000177901,0.0000183113,0.0000185955,0.0000189629,0.0000193379,0.0000195456,0.0000198736,0.0000201354,0.0000203706,0.0000205942,0.0000208091,0.0000210998,0.0000212481,0.0000213362,0.0000215558,0.0000216110,0.0000218793,0.0000220470,0.0000221131,0.0000223382,0.0000224237,0.0000225921,0.0000226488,0.0000227544,0.0000229383,0.0000230915,0.0000231121,0.0000231590,0.0000233473,0.0000233855,0.0000235118,0.0000235962,0.0000237780,0.0000237093,0.0000238382,0.0000240327,0.0000240319,0.0000240381,0.0000241434,0.0000241820,0.0000243326,0.0000243750,0.0000116303,0.0000129034,0.0000138789,0.0000147067,0.0000153681,0.0000159688,0.0000165381,0.0000168961,0.0000173735,0.0000176683,0.0000180134,0.0000183561,0.0000185591,0.0000188663,0.0000191044,0.0000193372,0.0000195693,0.0000197772,0.0000200196,0.0000201941,0.0000202678,0.0000204681,0.0000205303,0.0000207804,0.0000209276,0.0000209949,0.0000212170,0.0000213041,0.0000214717,0.0000215057,0.0000216145,0.0000217924,0.0000219308,0.0000219441,0.0000219783,0.0000221554,0.0000222116,0.0000223467,0.0000224102,0.0000225913,0.0000225114,0.0000226509,0.0000228391,0.0000228327,0.0000228512,0.0000229264,0.0000229647,0.0000231182,0.0000231734,0.0000110426,0.0000122583,0.0000131645,0.0000139544,0.0000146029,0.0000151709,0.0000156967,0.0000160399,0.0000165029,0.0000167694,0.0000171064,0.0000174375,0.0000176173,0.0000179206,0.0000181460,0.0000183690,0.0000185749,0.0000187890,0.0000190145,0.0000191894,0.0000192421,0.0000194367,0.0000194965,0.0000197436,0.0000198826,0.0000199438,0.0000201582,0.0000202382,0.0000204081,0.0000204218,0.0000205211,0.0000207026,0.0000208329,0.0000208520,0.0000208767,0.0000210464,0.0000211096,0.0000212328,0.0000212798,0.0000214786,0.0000213929,0.0000215054,0.0000217109,0.0000216888,0.0000217133,0.0000217893,0.0000218235,0.0000219857,0.0000219970,0.0000104803,0.0000116339,0.0000124979,0.0000132636,0.0000138557,0.0000143967,0.0000149114,0.0000152371,0.0000156656,0.0000159263,0.0000162558,0.0000165612,0.0000167315,0.0000170111,0.0000172422,0.0000174539,0.0000176330,0.0000178565,0.0000180507,0.0000182159,0.0000182843,0.0000184491,0.0000185260,0.0000187555,0.0000188690,0.0000189493,0.0000191476,0.0000192247,0.0000193769,0.0000194036,0.0000194965,0.0000196525,0.0000197922,0.0000198110,0.0000198427,0.0000199994,0.0000200581,0.0000201631,0.0000202127,0.0000203906,0.0000203175,0.0000204390,0.0000206246,0.0000205971,0.0000206337,0.0000207028,0.0000207302,0.0000208870,0.0000208996,0.0000099509,0.0000110432,0.0000118751,0.0000126030,0.0000131613,0.0000136759,0.0000141683,0.0000144695,0.0000148682,0.0000151356,0.0000154323,0.0000157201,0.0000158950,0.0000161504,0.0000163696,0.0000165823,0.0000167500,0.0000169509,0.0000171451,0.0000173029,0.0000173592,0.0000175310,0.0000176028,0.0000178076,0.0000179211,0.0000180051,0.0000181905,0.0000182844,0.0000184172,0.0000184132,0.0000185207,0.0000186686,0.0000187870,0.0000188214,0.0000188475,0.0000190130,0.0000190479,0.0000191490,0.0000191985,0.0000193800,0.0000193059,0.0000194458,0.0000195886,0.0000195634,0.0000196024,0.0000196594,0.0000197026,0.0000198511,0.0000198470,0.0000094477,0.0000104911,0.0000112771,0.0000119601,0.0000124874,0.0000129959,0.0000134460,0.0000137456,0.0000141218,0.0000143641,0.0000146516,0.0000149329,0.0000150909,0.0000153366,0.0000155449,0.0000157501,0.0000158957,0.0000160914,0.0000162567,0.0000164383,0.0000164889,0.0000166671,0.0000167151,0.0000169160,0.0000170125,0.0000171054,0.0000172715,0.0000173773,0.0000174928,0.0000174890,0.0000175867,0.0000177384,0.0000178458,0.0000178672,0.0000179041,0.0000180546,0.0000181106,0.0000181886,0.0000182338,0.0000183994,0.0000183240,0.0000184690,0.0000186221,0.0000185894,0.0000186374,0.0000186670,0.0000187303,0.0000188619,0.0000188694,0.0000089720,0.0000099589,0.0000107026,0.0000113566,0.0000118623,0.0000123487,0.0000127711,0.0000130536,0.0000134226,0.0000136432,0.0000139146,0.0000141812,0.0000143322,0.0000145611,0.0000147539,0.0000149682,0.0000150945,0.0000152795,0.0000154286,0.0000156216,0.0000156556,0.0000158419,0.0000158573,0.0000160750,0.0000161443,0.0000162482,0.0000164185,0.0000165066,0.0000166151,0.0000166170,0.0000167050,0.0000168626,0.0000169629,0.0000169860,0.0000170025,0.0000171446,0.0000172044,0.0000172703,0.0000173178,0.0000174879,0.0000174141,0.0000175268,0.0000176872,0.0000176711,0.0000176988,0.0000177456,0.0000178109,0.0000179168,0.0000179256,0.0000085275,0.0000094451,0.0000101628,0.0000107706,0.0000112737,0.0000117370,0.0000121316,0.0000123963,0.0000127358,0.0000129521,0.0000132117,0.0000134757,0.0000136117,0.0000138211,0.0000140096,0.0000142264,0.0000143371,0.0000145230,0.0000146569,0.0000148375,0.0000148764,0.0000150625,0.0000150524,0.0000152613,0.0000153305,0.0000154264,0.0000155865,0.0000157011,0.0000157908,0.0000157856,0.0000158628,0.0000160091,0.0000161027,0.0000161349,0.0000161473,0.0000162819,0.0000163595,0.0000164090,0.0000164410,0.0000166125,0.0000165365,0.0000166553,0.0000168136,0.0000167885,0.0000168275,0.0000168682,0.0000169094,0.0000170252,0.0000170369,0.0000081023,0.0000089612,0.0000096451,0.0000102265,0.0000107172,0.0000111421,0.0000115286,0.0000117675,0.0000120943,0.0000123012,0.0000125403,0.0000128061,0.0000129164,0.0000131250,0.0000132967,0.0000135211,0.0000136155,0.0000138033,0.0000139284,0.0000141042,0.0000141298,0.0000143058,0.0000143215,0.0000145038,0.0000145587,0.0000146477,0.0000148045,0.0000149182,0.0000150082,0.0000150003,0.0000150697,0.0000152035,0.0000153001,0.0000153343,0.0000153411,0.0000154607,0.0000155539,0.0000155998,0.0000156236,0.0000157821,0.0000157221,0.0000158274,0.0000159808,0.0000159514,0.0000159861,0.0000160247,0.0000160664,0.0000161798,0.0000161818,0.0000076972,0.0000085065,0.0000091551,0.0000097152,0.0000101677,0.0000105674,0.0000109423,0.0000111817,0.0000114870,0.0000116853,0.0000119162,0.0000121584,0.0000122584,0.0000124685,0.0000126257,0.0000128534,0.0000129426,0.0000131099,0.0000132349,0.0000134066,0.0000134324,0.0000135835,0.0000136163,0.0000137726,0.0000138239,0.0000139233,0.0000140655,0.0000141692,0.0000142645,0.0000142425,0.0000143126,0.0000144359,0.0000145465,0.0000145665,0.0000145713,0.0000146987,0.0000147904,0.0000148216,0.0000148447,0.0000150014,0.0000149212,0.0000150412,0.0000151837,0.0000151508,0.0000151905,0.0000152245,0.0000152537,0.0000153677,0.0000153698,0.0000073036,0.0000080984,0.0000086847,0.0000092291,0.0000096654,0.0000100420,0.0000104100,0.0000106165,0.0000109153,0.0000111120,0.0000113186,0.0000115475,0.0000116431,0.0000118577,0.0000119862,0.0000122240,0.0000122802,0.0000124457,0.0000125759,0.0000127385,0.0000127496,0.0000129104,0.0000129293,0.0000130759,0.0000131310,0.0000132353,0.0000133669,0.0000134655,0.0000135478,0.0000135134,0.0000136068,0.0000137139,0.0000138095,0.0000138294,0.0000138409,0.0000139606,0.0000140340,0.0000140820,0.0000141170,0.0000142653,0.0000141861,0.0000142954,0.0000144283,0.0000143886,0.0000144244,0.0000144720,0.0000144863,0.0000145905,0.0000146003,0.0000069225,0.0000076963,0.0000082462,0.0000087642,0.0000091762,0.0000095268,0.0000098858,0.0000100773,0.0000103598,0.0000105564,0.0000107564,0.0000109654,0.0000110552,0.0000112673,0.0000113891,0.0000116075,0.0000116683,0.0000118194,0.0000119394,0.0000120990,0.0000121055,0.0000122700,0.0000122896,0.0000124073,0.0000124679,0.0000125978,0.0000127146,0.0000127865,0.0000128849,0.0000128262,0.0000129186,0.0000130299,0.0000131239,0.0000131382,0.0000131467,0.0000132529,0.0000133430,0.0000133774,0.0000134036,0.0000135619,0.0000134811,0.0000135793,0.0000137112,0.0000136718,0.0000136919,0.0000137494,0.0000137722,0.0000138555,0.0000138721,0.0000065756,0.0000073049,0.0000078224,0.0000083150,0.0000087082,0.0000090392,0.0000093932,0.0000095713,0.0000098288,0.0000100287,0.0000102190,0.0000104217,0.0000105056,0.0000106953,0.0000108144,0.0000110177,0.0000110966,0.0000112248,0.0000113280,0.0000114744,0.0000114925,0.0000116562,0.0000116685,0.0000117816,0.0000118424,0.0000119685,0.0000120788,0.0000121368,0.0000122440,0.0000121812,0.0000122702,0.0000123751,0.0000124624,0.0000124772,0.0000124875,0.0000125954,0.0000126799,0.0000127078,0.0000127336,0.0000128879,0.0000128109,0.0000129025,0.0000130374,0.0000129706,0.0000130040,0.0000130799,0.0000130920,0.0000131769,0.0000131904,0.0000062360,0.0000069277,0.0000074373,0.0000078878,0.0000082738,0.0000085861,0.0000089250,0.0000090755,0.0000093365,0.0000095201,0.0000097038,0.0000098906,0.0000099726,0.0000101568,0.0000102511,0.0000104789,0.0000105402,0.0000106689,0.0000107682,0.0000108954,0.0000109074,0.0000110768,0.0000110884,0.0000112074,0.0000112540,0.0000113680,0.0000114859,0.0000115240,0.0000116452,0.0000115549,0.0000116640,0.0000117796,0.0000118446,0.0000118545,0.0000118686,0.0000119654,0.0000120582,0.0000120760,0.0000121025,0.0000122457,0.0000121690,0.0000122489,0.0000123908,0.0000123233,0.0000123521,0.0000124342,0.0000124462,0.0000125306,0.0000125324,0.0000059292,0.0000065732,0.0000070642,0.0000074889,0.0000078524,0.0000081676,0.0000084739,0.0000086250,0.0000088660,0.0000090497,0.0000092109,0.0000093917,0.0000094707,0.0000096444,0.0000097373,0.0000099570,0.0000100125,0.0000101420,0.0000102390,0.0000103555,0.0000103694,0.0000105236,0.0000105374,0.0000106467,0.0000107005,0.0000107915,0.0000109160,0.0000109499,0.0000110594,0.0000109724,0.0000110712,0.0000111862,0.0000112503,0.0000112687,0.0000112769,0.0000113744,0.0000114672,0.0000114713,0.0000114932,0.0000116354,0.0000115748,0.0000116390,0.0000117788,0.0000117105,0.0000117218,0.0000118123,0.0000118286,0.0000118939,0.0000119124,0.0000056281,0.0000062425,0.0000067027,0.0000071147,0.0000074504,0.0000077486,0.0000080501,0.0000081880,0.0000084123,0.0000085971,0.0000087418,0.0000089096,0.0000089919,0.0000091533,0.0000092397,0.0000094422,0.0000095142,0.0000096329,0.0000097194,0.0000098454,0.0000098496,0.0000099993,0.0000100053,0.0000101251,0.0000101703,0.0000102504,0.0000103655,0.0000103948,0.0000105120,0.0000104363,0.0000105189,0.0000106231,0.0000106911,0.0000107090,0.0000107231,0.0000108105,0.0000108936,0.0000109065,0.0000109174,0.0000110573,0.0000110132,0.0000110637,0.0000111925,0.0000111225,0.0000111439,0.0000112093,0.0000112270,0.0000112984,0.0000113354,0.0000053374,0.0000059303,0.0000063606,0.0000067541,0.0000070713,0.0000073557,0.0000076518,0.0000077654,0.0000079954,0.0000081628,0.0000083148,0.0000084630,0.0000085412,0.0000086839,0.0000087829,0.0000089682,0.0000090410,0.0000091467,0.0000092289,0.0000093531,0.0000093451,0.0000094953,0.0000095035,0.0000096219,0.0000096638,0.0000097268,0.0000098649,0.0000098729,0.0000099939,0.0000099108,0.0000099814,0.0000100952,0.0000101661,0.0000101841,0.0000101710,0.0000102697,0.0000103541,0.0000103532,0.0000103772,0.0000104901,0.0000104528,0.0000105102,0.0000106440,0.0000105695,0.0000105932,0.0000106472,0.0000106722,0.0000107366,0.0000107652,0.0000050663,0.0000056343,0.0000060427,0.0000064128,0.0000067177,0.0000069847,0.0000072680,0.0000073886,0.0000075876,0.0000077480,0.0000078921,0.0000080476,0.0000081090,0.0000082424,0.0000083423,0.0000085121,0.0000085810,0.0000086753,0.0000087509,0.0000088872,0.0000088720,0.0000090272,0.0000090252,0.0000091290,0.0000091823,0.0000092354,0.0000093838,0.0000093811,0.0000094887,0.0000094185,0.0000094837,0.0000095901,0.0000096469,0.0000096711,0.0000096607,0.0000097516,0.0000098467,0.0000098283,0.0000098626,0.0000099759,0.0000099318,0.0000099771,0.0000101137,0.0000100483,0.0000100728,0.0000101094,0.0000101544,0.0000102050,0.0000102330,0.0000048101,0.0000053401,0.0000057349,0.0000060921,0.0000063772,0.0000066321,0.0000069123,0.0000070158,0.0000072194,0.0000073625,0.0000074955,0.0000076389,0.0000077029,0.0000078299,0.0000079223,0.0000080871,0.0000081488,0.0000082312,0.0000083161,0.0000084321,0.0000084321,0.0000085806,0.0000085684,0.0000086757,0.0000087173,0.0000087710,0.0000089222,0.0000088955,0.0000090094,0.0000089384,0.0000090056,0.0000091072,0.0000091663,0.0000091769,0.0000091887,0.0000092770,0.0000093456,0.0000093385,0.0000093748,0.0000094822,0.0000094353,0.0000094846,0.0000096075,0.0000095587,0.0000095734,0.0000096015,0.0000096511,0.0000096983,0.0000097242,0.0000045656,0.0000050657,0.0000054523,0.0000057816,0.0000060618,0.0000063025,0.0000065652,0.0000066587,0.0000068593,0.0000069892,0.0000071317,0.0000072539,0.0000073167,0.0000074418,0.0000075133,0.0000076697,0.0000077377,0.0000078144,0.0000078992,0.0000080094,0.0000080087,0.0000081551,0.0000081350,0.0000082296,0.0000082830,0.0000083475,0.0000084826,0.0000084420,0.0000085506,0.0000084909,0.0000085510,0.0000086540,0.0000087111,0.0000086980,0.0000087352,0.0000088090,0.0000088924,0.0000088821,0.0000089109,0.0000090043,0.0000089507,0.0000090095,0.0000091155,0.0000090868,0.0000091015,0.0000091188,0.0000091686,0.0000092094,0.0000092384,0.0000043325,0.0000048113,0.0000051709,0.0000055032,0.0000057548,0.0000059803,0.0000062419,0.0000063222,0.0000065129,0.0000066263,0.0000067745,0.0000068875,0.0000069608,0.0000070727,0.0000071459,0.0000072776,0.0000073452,0.0000074158,0.0000074993,0.0000076098,0.0000076032,0.0000077419,0.0000077251,0.0000078197,0.0000078646,0.0000079257,0.0000080596,0.0000080071,0.0000081136,0.0000080664,0.0000081203,0.0000082277,0.0000082874,0.0000082636,0.0000083027,0.0000083730,0.0000084468,0.0000084325,0.0000084574,0.0000085585,0.0000085069,0.0000085560,0.0000086556,0.0000086433,0.0000086433,0.0000086545,0.0000087158,0.0000087531,0.0000087710,0.0000041088,0.0000045693,0.0000049138,0.0000052252,0.0000054692,0.0000056868,0.0000059195,0.0000060045,0.0000061822,0.0000062937,0.0000064381,0.0000065502,0.0000066070,0.0000067202,0.0000067844,0.0000069138,0.0000069790,0.0000070348,0.0000071260,0.0000072273,0.0000072295,0.0000073536,0.0000073344,0.0000074297,0.0000074711,0.0000075224,0.0000076589,0.0000076079,0.0000077161,0.0000076635,0.0000077175,0.0000078032,0.0000078704,0.0000078436,0.0000078906,0.0000079548,0.0000080134,0.0000080037,0.0000080313,0.0000081388,0.0000080813,0.0000081314,0.0000082190,0.0000082173,0.0000082053,0.0000082186,0.0000082801,0.0000083075,0.0000083435,0.0000039046,0.0000043454,0.0000046667,0.0000049600,0.0000052048,0.0000053998,0.0000056262,0.0000056965,0.0000058712,0.0000059743,0.0000061169,0.0000062150,0.0000062768,0.0000063837,0.0000064484,0.0000065617,0.0000066412,0.0000066822,0.0000067758,0.0000068771,0.0000068712,0.0000069911,0.0000069673,0.0000070518,0.0000071111,0.0000071479,0.0000072911,0.0000072299,0.0000073334,0.0000072628,0.0000073280,0.0000074078,0.0000074725,0.0000074464,0.0000075049,0.0000075566,0.0000076129,0.0000075999,0.0000076315,0.0000077331,0.0000076743,0.0000077307,0.0000078057,0.0000078131,0.0000077992,0.0000078020,0.0000078772,0.0000079089,0.0000079399,0.0000037102,0.0000041216,0.0000044252,0.0000047123,0.0000049363,0.0000051390,0.0000053459,0.0000054135,0.0000055838,0.0000056777,0.0000058034,0.0000059074,0.0000059627,0.0000060650,0.0000061171,0.0000062336,0.0000063039,0.0000063505,0.0000064379,0.0000065322,0.0000065193,0.0000066418,0.0000066234,0.0000067028,0.0000067641,0.0000067870,0.0000069214,0.0000068674,0.0000069693,0.0000068989,0.0000069609,0.0000070384,0.0000071057,0.0000070781,0.0000071294,0.0000071784,0.0000072316,0.0000072171,0.0000072545,0.0000073410,0.0000072907,0.0000073460,0.0000074150,0.0000074272,0.0000074092,0.0000074105,0.0000074822,0.0000075111,0.0000075411,0.0000035175,0.0000039113,0.0000041926,0.0000044692,0.0000046841,0.0000048767,0.0000050806,0.0000051399,0.0000052912,0.0000053991,0.0000055161,0.0000056041,0.0000056552,0.0000057632,0.0000058187,0.0000059251,0.0000059834,0.0000060276,0.0000061143,0.0000062095,0.0000061958,0.0000063125,0.0000062961,0.0000063759,0.0000064217,0.0000064490,0.0000065711,0.0000065319,0.0000066081,0.0000065505,0.0000066182,0.0000066958,0.0000067619,0.0000067166,0.0000067736,0.0000068173,0.0000068704,0.0000068627,0.0000068867,0.0000069749,0.0000069307,0.0000069772,0.0000070352,0.0000070525,0.0000070331,0.0000070378,0.0000071085,0.0000071507,0.0000071634,0.0000033434,0.0000037185,0.0000039891,0.0000042455,0.0000044481,0.0000046394,0.0000048186,0.0000048804,0.0000050272,0.0000051257,0.0000052425,0.0000053252,0.0000053735,0.0000054788,0.0000055351,0.0000056359,0.0000056839,0.0000057202,0.0000058057,0.0000058892,0.0000058874,0.0000059932,0.0000059859,0.0000060480,0.0000061010,0.0000061320,0.0000062462,0.0000062096,0.0000062772,0.0000062203,0.0000062890,0.0000063664,0.0000064188,0.0000063761,0.0000064286,0.0000064795,0.0000065198,0.0000065205,0.0000065410,0.0000066232,0.0000065836,0.0000066239,0.0000066884,0.0000067017,0.0000066780,0.0000066914,0.0000067506,0.0000068016,0.0000068113,0.0000031768,0.0000035315,0.0000037911,0.0000040282,0.0000042239,0.0000044048,0.0000045678,0.0000046349,0.0000047756,0.0000048716,0.0000049778,0.0000050654,0.0000050948,0.0000052098,0.0000052669,0.0000053573,0.0000054066,0.0000054317,0.0000055243,0.0000055899,0.0000055949,0.0000056849,0.0000056868,0.0000057449,0.0000057992,0.0000058220,0.0000059325,0.0000058896,0.0000059720,0.0000059014,0.0000059803,0.0000060510,0.0000060872,0.0000060597,0.0000061123,0.0000061579,0.0000061906,0.0000061942,0.0000062160,0.0000063013,0.0000062576,0.0000062913,0.0000063608,0.0000063670,0.0000063422,0.0000063600,0.0000064153,0.0000064759,0.0000064696,0.0000030123,0.0000033547,0.0000035931,0.0000038296,0.0000040052,0.0000041776,0.0000043390,0.0000043994,0.0000045255,0.0000046299,0.0000047287,0.0000048138,0.0000048442,0.0000049430,0.0000050046,0.0000050930,0.0000051348,0.0000051599,0.0000052476,0.0000053160,0.0000053164,0.0000053982,0.0000054065,0.0000054542,0.0000055069,0.0000055346,0.0000056361,0.0000055941,0.0000056682,0.0000056113,0.0000056764,0.0000057493,0.0000057915,0.0000057619,0.0000058140,0.0000058502,0.0000058788,0.0000058868,0.0000059020,0.0000059862,0.0000059504,0.0000059782,0.0000060503,0.0000060361,0.0000060355,0.0000060410,0.0000060997,0.0000061499,0.0000061532,0.0000028643,0.0000031848,0.0000034086,0.0000036357,0.0000037976,0.0000039765,0.0000041143,0.0000041719,0.0000043002,0.0000044026,0.0000044885,0.0000045740,0.0000046060,0.0000046885,0.0000047553,0.0000048339,0.0000048719,0.0000048989,0.0000049856,0.0000050477,0.0000050617,0.0000051327,0.0000051445,0.0000051812,0.0000052360,0.0000052657,0.0000053530,0.0000053153,0.0000053789,0.0000053250,0.0000053882,0.0000054674,0.0000054993,0.0000054757,0.0000055199,0.0000055610,0.0000055777,0.0000055955,0.0000056138,0.0000056899,0.0000056494,0.0000056813,0.0000057496,0.0000057373,0.0000057404,0.0000057441,0.0000057915,0.0000058461,0.0000058527,0.0000027140,0.0000030227,0.0000032367,0.0000034555,0.0000036058,0.0000037792,0.0000039144,0.0000039723,0.0000040899,0.0000041802,0.0000042609,0.0000043483,0.0000043782,0.0000044572,0.0000045167,0.0000045990,0.0000046291,0.0000046566,0.0000047378,0.0000047927,0.0000048087,0.0000048780,0.0000048895,0.0000049198,0.0000049786,0.0000050093,0.0000050936,0.0000050528,0.0000051106,0.0000050589,0.0000051240,0.0000051896,0.0000052241,0.0000052108,0.0000052417,0.0000053011,0.0000052889,0.0000053096,0.0000053343,0.0000054050,0.0000053699,0.0000054088,0.0000054666,0.0000054530,0.0000054568,0.0000054505,0.0000054959,0.0000055587,0.0000055581,0.0000025747,0.0000028711,0.0000030733,0.0000032743,0.0000034354,0.0000035886,0.0000037188,0.0000037732,0.0000038871,0.0000039645,0.0000040530,0.0000041299,0.0000041597,0.0000042368,0.0000042870,0.0000043668,0.0000043954,0.0000044279,0.0000045081,0.0000045488,0.0000045746,0.0000046344,0.0000046445,0.0000046693,0.0000047327,0.0000047577,0.0000048367,0.0000047966,0.0000048506,0.0000048025,0.0000048716,0.0000049344,0.0000049690,0.0000049523,0.0000049762,0.0000050431,0.0000050245,0.0000050496,0.0000050696,0.0000051305,0.0000050998,0.0000051533,0.0000051887,0.0000051722,0.0000051768,0.0000051697,0.0000052233,0.0000052811,0.0000052735,0.0000024495,0.0000027240,0.0000029149,0.0000031076,0.0000032572,0.0000034130,0.0000035378,0.0000035863,0.0000036959,0.0000037785,0.0000038444,0.0000039258,0.0000039476,0.0000040310,0.0000040694,0.0000041516,0.0000041703,0.0000042172,0.0000042877,0.0000043222,0.0000043501,0.0000044044,0.0000044124,0.0000044370,0.0000044934,0.0000045204,0.0000046019,0.0000045600,0.0000046142,0.0000045660,0.0000046333,0.0000046861,0.0000047246,0.0000046990,0.0000047235,0.0000047936,0.0000047777,0.0000048045,0.0000048223,0.0000048771,0.0000048504,0.0000048979,0.0000049353,0.0000049185,0.0000049178,0.0000049116,0.0000049555,0.0000050205,0.0000050112,0.0000023203,0.0000025873,0.0000027680,0.0000029513,0.0000030928,0.0000032417,0.0000033574,0.0000034125,0.0000035058,0.0000035870,0.0000036597,0.0000037251,0.0000037472,0.0000038324,0.0000038627,0.0000039339,0.0000039636,0.0000040037,0.0000040719,0.0000041106,0.0000041399,0.0000041885,0.0000041934,0.0000042182,0.0000042724,0.0000042965,0.0000043791,0.0000043300,0.0000043845,0.0000043489,0.0000043940,0.0000044497,0.0000044900,0.0000044534,0.0000044751,0.0000045527,0.0000045395,0.0000045688,0.0000045742,0.0000046331,0.0000046106,0.0000046586,0.0000046936,0.0000046734,0.0000046706,0.0000046607,0.0000047022,0.0000047697,0.0000047621,0.0000021997,0.0000024576,0.0000026311,0.0000028061,0.0000029338,0.0000030791,0.0000031833,0.0000032382,0.0000033358,0.0000034044,0.0000034779,0.0000035341,0.0000035567,0.0000036392,0.0000036730,0.0000037343,0.0000037586,0.0000038018,0.0000038715,0.0000039047,0.0000039308,0.0000039865,0.0000039829,0.0000040030,0.0000040533,0.0000040778,0.0000041554,0.0000041174,0.0000041658,0.0000041342,0.0000041733,0.0000042311,0.0000042668,0.0000042255,0.0000042515,0.0000043248,0.0000043151,0.0000043434,0.0000043411,0.0000044029,0.0000043889,0.0000044305,0.0000044488,0.0000044429,0.0000044346,0.0000044373,0.0000044643,0.0000045368,0.0000045232,0.0000020815,0.0000023376,0.0000025052,0.0000026659,0.0000027878,0.0000029153,0.0000030282,0.0000030700,0.0000031642,0.0000032307,0.0000033063,0.0000033565,0.0000033840,0.0000034637,0.0000034940,0.0000035465,0.0000035719,0.0000036094,0.0000036775,0.0000037105,0.0000037401,0.0000037921,0.0000037768,0.0000038071,0.0000038474,0.0000038726,0.0000039462,0.0000039031,0.0000039592,0.0000039241,0.0000039655,0.0000040390,0.0000040575,0.0000040134,0.0000040417,0.0000041073,0.0000040986,0.0000041320,0.0000041361,0.0000041943,0.0000041638,0.0000042059,0.0000042327,0.0000042222,0.0000042142,0.0000042211,0.0000042526,0.0000043136,0.0000043041,0.0000019780,0.0000022127,0.0000023752,0.0000025339,0.0000026399,0.0000027693,0.0000028766,0.0000029199,0.0000030060,0.0000030653,0.0000031458,0.0000031972,0.0000032195,0.0000032900,0.0000033161,0.0000033755,0.0000034030,0.0000034272,0.0000034941,0.0000035263,0.0000035513,0.0000036013,0.0000035856,0.0000036212,0.0000036594,0.0000036768,0.0000037443,0.0000037119,0.0000037510,0.0000037348,0.0000037698,0.0000038329,0.0000038580,0.0000038189,0.0000038408,0.0000039025,0.0000038999,0.0000039304,0.0000039340,0.0000039843,0.0000039564,0.0000040017,0.0000040131,0.0000040177,0.0000040086,0.0000040115,0.0000040328,0.0000040913,0.0000040869,0.0000018783,0.0000020988,0.0000022539,0.0000024098,0.0000025094,0.0000026247,0.0000027326,0.0000027756,0.0000028552,0.0000029150,0.0000029918,0.0000030394,0.0000030682,0.0000031202,0.0000031499,0.0000032092,0.0000032295,0.0000032533,0.0000033215,0.0000033512,0.0000033809,0.0000034207,0.0000034047,0.0000034379,0.0000034738,0.0000034942,0.0000035636,0.0000035231,0.0000035652,0.0000035429,0.0000035801,0.0000036461,0.0000036622,0.0000036244,0.0000036578,0.0000037114,0.0000037021,0.0000037276,0.0000037416,0.0000037872,0.0000037582,0.0000038066,0.0000038062,0.0000038109,0.0000038081,0.0000038124,0.0000038331,0.0000038803,0.0000038879,0.0000017856,0.0000019924,0.0000021455,0.0000022879,0.0000023770,0.0000024941,0.0000025940,0.0000026324,0.0000027013,0.0000027686,0.0000028439,0.0000028841,0.0000029164,0.0000029630,0.0000029800,0.0000030433,0.0000030640,0.0000030887,0.0000031551,0.0000031835,0.0000032092,0.0000032513,0.0000032359,0.0000032636,0.0000032958,0.0000033147,0.0000033793,0.0000033480,0.0000033904,0.0000033664,0.0000034054,0.0000034659,0.0000034845,0.0000034462,0.0000034711,0.0000035239,0.0000035244,0.0000035447,0.0000035618,0.0000036004,0.0000035770,0.0000036267,0.0000036129,0.0000036214,0.0000036189,0.0000036296,0.0000036366,0.0000036925,0.0000036897,0.0000016972,0.0000018905,0.0000020393,0.0000021751,0.0000022535,0.0000023741,0.0000024627,0.0000025002,0.0000025613,0.0000026306,0.0000027024,0.0000027389,0.0000027752,0.0000028187,0.0000028315,0.0000028896,0.0000029116,0.0000029393,0.0000029992,0.0000030172,0.0000030514,0.0000030929,0.0000030770,0.0000031047,0.0000031306,0.0000031528,0.0000032117,0.0000031755,0.0000032231,0.0000032038,0.0000032332,0.0000032964,0.0000033104,0.0000032806,0.0000032977,0.0000033446,0.0000033463,0.0000033610,0.0000033866,0.0000034190,0.0000034002,0.0000034506,0.0000034355,0.0000034396,0.0000034402,0.0000034583,0.0000034522,0.0000035052,0.0000035014,0.0000016090,0.0000017958,0.0000019344,0.0000020652,0.0000021348,0.0000022542,0.0000023317,0.0000023665,0.0000024337,0.0000024977,0.0000025634,0.0000026033,0.0000026320,0.0000026771,0.0000026949,0.0000027369,0.0000027664,0.0000027969,0.0000028478,0.0000028560,0.0000028951,0.0000029400,0.0000029204,0.0000029479,0.0000029750,0.0000029986,0.0000030436,0.0000030228,0.0000030615,0.0000030429,0.0000030797,0.0000031314,0.0000031408,0.0000031191,0.0000031314,0.0000031840,0.0000031784,0.0000031936,0.0000032185,0.0000032534,0.0000032311,0.0000032806,0.0000032605,0.0000032642,0.0000032703,0.0000032812,0.0000032806,0.0000033298,0.0000033178,0.0000015315,0.0000017063,0.0000018386,0.0000019641,0.0000020251,0.0000021395,0.0000022148,0.0000022492,0.0000023109,0.0000023715,0.0000024409,0.0000024695,0.0000024938,0.0000025452,0.0000025631,0.0000025938,0.0000026278,0.0000026577,0.0000027077,0.0000027159,0.0000027469,0.0000027948,0.0000027694,0.0000027972,0.0000028336,0.0000028562,0.0000028965,0.0000028736,0.0000029006,0.0000028803,0.0000029290,0.0000029764,0.0000029844,0.0000029638,0.0000029709,0.0000030328,0.0000030198,0.0000030385,0.0000030605,0.0000030894,0.0000030742,0.0000031267,0.0000030989,0.0000031057,0.0000031071,0.0000031177,0.0000031056,0.0000031642,0.0000031478,0.0000014557,0.0000016196,0.0000017454,0.0000018653,0.0000019203,0.0000020314,0.0000020959,0.0000021383,0.0000021997,0.0000022497,0.0000023192,0.0000023434,0.0000023689,0.0000024151,0.0000024360,0.0000024610,0.0000024983,0.0000025246,0.0000025715,0.0000025812,0.0000026024,0.0000026537,0.0000026299,0.0000026604,0.0000026924,0.0000027146,0.0000027503,0.0000027321,0.0000027566,0.0000027380,0.0000027857,0.0000028318,0.0000028322,0.0000028142,0.0000028234,0.0000028866,0.0000028678,0.0000028800,0.0000029121,0.0000029334,0.0000029250,0.0000029714,0.0000029414,0.0000029540,0.0000029588,0.0000029599,0.0000029531,0.0000030036,0.0000029871,0.0000013777,0.0000015363,0.0000016577,0.0000017704,0.0000018262,0.0000019318,0.0000019884,0.0000020272,0.0000020842,0.0000021315,0.0000022069,0.0000022220,0.0000022520,0.0000022915,0.0000023128,0.0000023421,0.0000023759,0.0000023985,0.0000024433,0.0000024541,0.0000024712,0.0000025105,0.0000024967,0.0000025240,0.0000025525,0.0000025762,0.0000026142,0.0000025942,0.0000026218,0.0000026018,0.0000026504,0.0000026941,0.0000026982,0.0000026740,0.0000026814,0.0000027469,0.0000027198,0.0000027402,0.0000027658,0.0000027914,0.0000027824,0.0000028237,0.0000027940,0.0000028071,0.0000028084,0.0000028156,0.0000028064,0.0000028563,0.0000028363,0.0000013087,0.0000014580,0.0000015674,0.0000016828,0.0000017319,0.0000018339,0.0000018857,0.0000019275,0.0000019879,0.0000020275,0.0000020938,0.0000021111,0.0000021421,0.0000021759,0.0000021952,0.0000022281,0.0000022625,0.0000022767,0.0000023244,0.0000023356,0.0000023472,0.0000023916,0.0000023692,0.0000024008,0.0000024269,0.0000024504,0.0000024814,0.0000024619,0.0000024961,0.0000024795,0.0000025147,0.0000025624,0.0000025626,0.0000025388,0.0000025440,0.0000026088,0.0000025855,0.0000026060,0.0000026262,0.0000026480,0.0000026445,0.0000026831,0.0000026570,0.0000026672,0.0000026626,0.0000026695,0.0000026656,0.0000027128,0.0000026978,0.0000012424,0.0000013847,0.0000014898,0.0000015967,0.0000016431,0.0000017403,0.0000017959,0.0000018342,0.0000018861,0.0000019273,0.0000019891,0.0000020037,0.0000020278,0.0000020667,0.0000020824,0.0000021129,0.0000021485,0.0000021590,0.0000022033,0.0000022167,0.0000022292,0.0000022666,0.0000022509,0.0000022807,0.0000023009,0.0000023272,0.0000023598,0.0000023369,0.0000023718,0.0000023513,0.0000023903,0.0000024336,0.0000024372,0.0000024122,0.0000024132,0.0000024774,0.0000024682,0.0000024803,0.0000024976,0.0000025114,0.0000025081,0.0000025540,0.0000025259,0.0000025316,0.0000025246,0.0000025351,0.0000025380,0.0000025647,0.0000025674,0.0000011818,0.0000013177,0.0000014186,0.0000015156,0.0000015639,0.0000016473,0.0000017060,0.0000017420,0.0000017880,0.0000018300,0.0000018884,0.0000018986,0.0000019242,0.0000019641,0.0000019766,0.0000020064,0.0000020406,0.0000020480,0.0000020999,0.0000021037,0.0000021184,0.0000021531,0.0000021369,0.0000021675,0.0000021925,0.0000022067,0.0000022417,0.0000022166,0.0000022539,0.0000022379,0.0000022639,0.0000023166,0.0000023174,0.0000022944,0.0000023058,0.0000023539,0.0000023448,0.0000023557,0.0000023755,0.0000023876,0.0000023857,0.0000024292,0.0000024037,0.0000024105,0.0000023978,0.0000024113,0.0000024141,0.0000024347,0.0000024374,0.0000011181,0.0000012490,0.0000013467,0.0000014416,0.0000014825,0.0000015661,0.0000016147,0.0000016504,0.0000016975,0.0000017378,0.0000017918,0.0000018006,0.0000018256,0.0000018663,0.0000018832,0.0000019057,0.0000019407,0.0000019452,0.0000019970,0.0000019933,0.0000020155,0.0000020495,0.0000020293,0.0000020610,0.0000020813,0.0000020944,0.0000021294,0.0000021068,0.0000021422,0.0000021279,0.0000021478,0.0000021957,0.0000022013,0.0000021781,0.0000021927,0.0000022351,0.0000022296,0.0000022370,0.0000022568,0.0000022741,0.0000022642,0.0000022989,0.0000022868,0.0000022917,0.0000022774,0.0000022938,0.0000022956,0.0000023127,0.0000023171,0.0000010598,0.0000011870,0.0000012766,0.0000013708,0.0000014066,0.0000014871,0.0000015298,0.0000015700,0.0000016082,0.0000016502,0.0000017057,0.0000017068,0.0000017394,0.0000017706,0.0000017934,0.0000018106,0.0000018409,0.0000018475,0.0000019000,0.0000018901,0.0000019123,0.0000019528,0.0000019272,0.0000019592,0.0000019749,0.0000019884,0.0000020245,0.0000020014,0.0000020270,0.0000020198,0.0000020370,0.0000020848,0.0000020917,0.0000020701,0.0000020860,0.0000021200,0.0000021182,0.0000021222,0.0000021476,0.0000021615,0.0000021448,0.0000021867,0.0000021726,0.0000021822,0.0000021651,0.0000021857,0.0000021822,0.0000021997,0.0000022054,0.0000010072,0.0000011268,0.0000012140,0.0000013004,0.0000013373,0.0000014083,0.0000014487,0.0000014901,0.0000015301,0.0000015684,0.0000016228,0.0000016177,0.0000016506,0.0000016887,0.0000017048,0.0000017201,0.0000017454,0.0000017497,0.0000018092,0.0000017982,0.0000018150,0.0000018557,0.0000018282,0.0000018638,0.0000018783,0.0000018817,0.0000019205,0.0000018985,0.0000019260,0.0000019177,0.0000019366,0.0000019831,0.0000019819,0.0000019666,0.0000019828,0.0000020108,0.0000020128,0.0000020196,0.0000020381,0.0000020571,0.0000020359,0.0000020767,0.0000020634,0.0000020759,0.0000020582,0.0000020765,0.0000020735,0.0000020914,0.0000020963,0.0000009537,0.0000010678,0.0000011511,0.0000012383,0.0000012681,0.0000013414,0.0000013757,0.0000014147,0.0000014570,0.0000014851,0.0000015379,0.0000015357,0.0000015652,0.0000016111,0.0000016206,0.0000016359,0.0000016567,0.0000016598,0.0000017152,0.0000017045,0.0000017256,0.0000017637,0.0000017356,0.0000017728,0.0000017849,0.0000017889,0.0000018244,0.0000018016,0.0000018275,0.0000018167,0.0000018376,0.0000018846,0.0000018885,0.0000018688,0.0000018793,0.0000019088,0.0000019166,0.0000019181,0.0000019385,0.0000019568,0.0000019364,0.0000019699,0.0000019591,0.0000019724,0.0000019545,0.0000019776,0.0000019701,0.0000019842,0.0000019935,0.0000009036,0.0000010183,0.0000010897,0.0000011766,0.0000012072,0.0000012745,0.0000013050,0.0000013408,0.0000013858,0.0000014140,0.0000014565,0.0000014608,0.0000014855,0.0000015303,0.0000015353,0.0000015518,0.0000015744,0.0000015756,0.0000016276,0.0000016193,0.0000016386,0.0000016737,0.0000016529,0.0000016849,0.0000016947,0.0000017022,0.0000017322,0.0000017058,0.0000017372,0.0000017201,0.0000017458,0.0000017972,0.0000017954,0.0000017789,0.0000017857,0.0000018152,0.0000018209,0.0000018254,0.0000018391,0.0000018554,0.0000018410,0.0000018688,0.0000018610,0.0000018750,0.0000018553,0.0000018833,0.0000018711,0.0000018848,0.0000018949,0.0000008562,0.0000009700,0.0000010344,0.0000011184,0.0000011448,0.0000012115,0.0000012412,0.0000012707,0.0000013163,0.0000013421,0.0000013822,0.0000013818,0.0000014117,0.0000014530,0.0000014597,0.0000014735,0.0000014936,0.0000014959,0.0000015473,0.0000015372,0.0000015550,0.0000015877,0.0000015725,0.0000015976,0.0000016119,0.0000016180,0.0000016462,0.0000016184,0.0000016521,0.0000016345,0.0000016567,0.0000017096,0.0000017048,0.0000016941,0.0000016971,0.0000017243,0.0000017269,0.0000017338,0.0000017464,0.0000017654,0.0000017499,0.0000017725,0.0000017726,0.0000017813,0.0000017654,0.0000017913,0.0000017832,0.0000017962,0.0000018005,0.0000008113,0.0000009201,0.0000009817,0.0000010595,0.0000010874,0.0000011491,0.0000011806,0.0000012122,0.0000012474,0.0000012723,0.0000013140,0.0000013110,0.0000013371,0.0000013811,0.0000013922,0.0000014031,0.0000014174,0.0000014187,0.0000014720,0.0000014586,0.0000014796,0.0000015105,0.0000014941,0.0000015188,0.0000015288,0.0000015406,0.0000015638,0.0000015366,0.0000015709,0.0000015548,0.0000015696,0.0000016259,0.0000016135,0.0000016122,0.0000016150,0.0000016366,0.0000016395,0.0000016495,0.0000016574,0.0000016774,0.0000016615,0.0000016861,0.0000016827,0.0000016996,0.0000016803,0.0000017039,0.0000016964,0.0000017071,0.0000017082,0.0000007697,0.0000008738,0.0000009315,0.0000010053,0.0000010313,0.0000010920,0.0000011219,0.0000011520,0.0000011852,0.0000012114,0.0000012543,0.0000012430,0.0000012670,0.0000013111,0.0000013254,0.0000013319,0.0000013460,0.0000013461,0.0000013959,0.0000013845,0.0000014081,0.0000014366,0.0000014174,0.0000014459,0.0000014495,0.0000014641,0.0000014855,0.0000014604,0.0000014947,0.0000014812,0.0000014864,0.0000015420,0.0000015320,0.0000015275,0.0000015394,0.0000015511,0.0000015613,0.0000015658,0.0000015775,0.0000015964,0.0000015795,0.0000015996,0.0000015975,0.0000016202,0.0000015950,0.0000016198,0.0000016140,0.0000016169,0.0000016235,0.0000007335,0.0000008295,0.0000008863,0.0000009557,0.0000009802,0.0000010363,0.0000010660,0.0000010942,0.0000011217,0.0000011504,0.0000011898,0.0000011789,0.0000012062,0.0000012451,0.0000012593,0.0000012654,0.0000012738,0.0000012829,0.0000013255,0.0000013130,0.0000013394,0.0000013581,0.0000013421,0.0000013665,0.0000013805,0.0000013870,0.0000014097,0.0000013875,0.0000014203,0.0000014077,0.0000014133,0.0000014627,0.0000014570,0.0000014505,0.0000014637,0.0000014755,0.0000014792,0.0000014837,0.0000014916,0.0000015176,0.0000014955,0.0000015203,0.0000015147,0.0000015439,0.0000015143,0.0000015402,0.0000015340,0.0000015377,0.0000015343,0.0000006978,0.0000007860,0.0000008369,0.0000009044,0.0000009354,0.0000009867,0.0000010144,0.0000010389,0.0000010648,0.0000010960,0.0000011332,0.0000011164,0.0000011454,0.0000011824,0.0000011961,0.0000012039,0.0000012132,0.0000012159,0.0000012597,0.0000012467,0.0000012733,0.0000012900,0.0000012775,0.0000013037,0.0000013094,0.0000013182,0.0000013451,0.0000013167,0.0000013450,0.0000013397,0.0000013421,0.0000013899,0.0000013849,0.0000013782,0.0000013943,0.0000014009,0.0000014088,0.0000014090,0.0000014197,0.0000014416,0.0000014165,0.0000014446,0.0000014397,0.0000014647,0.0000014415,0.0000014640,0.0000014610,0.0000014579,0.0000014562,0.0000006642,0.0000007500,0.0000007953,0.0000008598,0.0000008868,0.0000009442,0.0000009637,0.0000009859,0.0000010103,0.0000010409,0.0000010783,0.0000010621,0.0000010901,0.0000011232,0.0000011394,0.0000011448,0.0000011504,0.0000011579,0.0000011931,0.0000011836,0.0000012088,0.0000012309,0.0000012150,0.0000012373,0.0000012454,0.0000012533,0.0000012786,0.0000012484,0.0000012764,0.0000012714,0.0000012777,0.0000013163,0.0000013180,0.0000013141,0.0000013255,0.0000013363,0.0000013418,0.0000013448,0.0000013515,0.0000013700,0.0000013453,0.0000013693,0.0000013598,0.0000013897,0.0000013708,0.0000013969,0.0000013863,0.0000013822,0.0000013833,0.0000006310,0.0000007129,0.0000007560,0.0000008172,0.0000008410,0.0000008964,0.0000009123,0.0000009386,0.0000009583,0.0000009883,0.0000010244,0.0000010088,0.0000010338,0.0000010682,0.0000010845,0.0000010896,0.0000010923,0.0000011002,0.0000011331,0.0000011230,0.0000011490,0.0000011657,0.0000011564,0.0000011706,0.0000011862,0.0000011886,0.0000012158,0.0000011845,0.0000012154,0.0000012069,0.0000012141,0.0000012481,0.0000012551,0.0000012485,0.0000012545,0.0000012661,0.0000012743,0.0000012740,0.0000012861,0.0000013019,0.0000012802,0.0000013005,0.0000012925,0.0000013191,0.0000013041,0.0000013242,0.0000013130,0.0000013168,0.0000013178,0.0000005956,0.0000006758,0.0000007193,0.0000007771,0.0000007988,0.0000008518,0.0000008627,0.0000008982,0.0000009111,0.0000009375,0.0000009725,0.0000009615,0.0000009879,0.0000010140,0.0000010294,0.0000010323,0.0000010355,0.0000010483,0.0000010700,0.0000010666,0.0000010912,0.0000011065,0.0000010952,0.0000011139,0.0000011256,0.0000011296,0.0000011567,0.0000011230,0.0000011550,0.0000011464,0.0000011538,0.0000011886,0.0000011938,0.0000011862,0.0000011920,0.0000012069,0.0000012155,0.0000012123,0.0000012207,0.0000012368,0.0000012141,0.0000012345,0.0000012261,0.0000012532,0.0000012394,0.0000012645,0.0000012456,0.0000012434,0.0000012520,0.0000005664,0.0000006427,0.0000006811,0.0000007371,0.0000007564,0.0000008081,0.0000008172,0.0000008516,0.0000008630,0.0000008915,0.0000009213,0.0000009129,0.0000009360,0.0000009664,0.0000009791,0.0000009858,0.0000009812,0.0000009990,0.0000010160,0.0000010098,0.0000010394,0.0000010529,0.0000010349,0.0000010580,0.0000010692,0.0000010747,0.0000010982,0.0000010680,0.0000010982,0.0000010918,0.0000011006,0.0000011284,0.0000011320,0.0000011278,0.0000011284,0.0000011490,0.0000011555,0.0000011527,0.0000011529,0.0000011756,0.0000011555,0.0000011675,0.0000011647,0.0000011994,0.0000011743,0.0000012012,0.0000011792,0.0000011766,0.0000011876,0.0000005350,0.0000006094,0.0000006455,0.0000007031,0.0000007170,0.0000007685,0.0000007746,0.0000008094,0.0000008225,0.0000008489,0.0000008686,0.0000008631,0.0000008908,0.0000009190,0.0000009349,0.0000009356,0.0000009309,0.0000009442,0.0000009673,0.0000009610,0.0000009870,0.0000010009,0.0000009875,0.0000010076,0.0000010172,0.0000010187,0.0000010408,0.0000010168,0.0000010434,0.0000010355,0.0000010437,0.0000010713,0.0000010788,0.0000010728,0.0000010720,0.0000010963,0.0000010980,0.0000010954,0.0000010965,0.0000011163,0.0000011011,0.0000011047,0.0000011102,0.0000011432,0.0000011140,0.0000011466,0.0000011200,0.0000011199,0.0000011293,0.0000005081,0.0000005816,0.0000006136,0.0000006683,0.0000006811,0.0000007298,0.0000007383,0.0000007712,0.0000007809,0.0000008060,0.0000008265,0.0000008188,0.0000008431,0.0000008699,0.0000008874,0.0000008898,0.0000008862,0.0000008980,0.0000009225,0.0000009134,0.0000009364,0.0000009507,0.0000009392,0.0000009562,0.0000009712,0.0000009704,0.0000009874,0.0000009664,0.0000009931,0.0000009865,0.0000009925,0.0000010188,0.0000010224,0.0000010175,0.0000010188,0.0000010389,0.0000010436,0.0000010420,0.0000010429,0.0000010595,0.0000010454,0.0000010502,0.0000010584,0.0000010882,0.0000010573,0.0000010924,0.0000010656,0.0000010693,0.0000010709,0.0000004816,0.0000005529,0.0000005794,0.0000006343,0.0000006472,0.0000006931,0.0000007005,0.0000007342,0.0000007430,0.0000007635,0.0000007843,0.0000007805,0.0000007998,0.0000008262,0.0000008443,0.0000008458,0.0000008444,0.0000008551,0.0000008742,0.0000008679,0.0000008908,0.0000009048,0.0000008901,0.0000009118,0.0000009227,0.0000009237,0.0000009386,0.0000009161,0.0000009448,0.0000009324,0.0000009409,0.0000009712,0.0000009692,0.0000009655,0.0000009702,0.0000009867,0.0000009924,0.0000009893,0.0000009919,0.0000010027,0.0000009943,0.0000009972,0.0000010060,0.0000010370,0.0000010036,0.0000010384,0.0000010128,0.0000010192,0.0000010173,0.0000004563,0.0000005244,0.0000005499,0.0000006027,0.0000006133,0.0000006584,0.0000006673,0.0000006976,0.0000007053,0.0000007249,0.0000007469,0.0000007400,0.0000007607,0.0000007837,0.0000008028,0.0000008059,0.0000008033,0.0000008125,0.0000008298,0.0000008259,0.0000008491,0.0000008659,0.0000008441,0.0000008668,0.0000008785,0.0000008796,0.0000008918,0.0000008669,0.0000008997,0.0000008849,0.0000008883,0.0000009247,0.0000009233,0.0000009114,0.0000009216,0.0000009357,0.0000009412,0.0000009421,0.0000009407,0.0000009517,0.0000009409,0.0000009459,0.0000009525,0.0000009861,0.0000009562,0.0000009835,0.0000009630,0.0000009646,0.0000009661,0.0000004330,0.0000004963,0.0000005241,0.0000005730,0.0000005823,0.0000006269,0.0000006351,0.0000006617,0.0000006669,0.0000006884,0.0000007079,0.0000007015,0.0000007232,0.0000007454,0.0000007642,0.0000007627,0.0000007616,0.0000007718,0.0000007881,0.0000007846,0.0000008100,0.0000008210,0.0000008011,0.0000008238,0.0000008338,0.0000008347,0.0000008491,0.0000008238,0.0000008542,0.0000008381,0.0000008460,0.0000008754,0.0000008776,0.0000008688,0.0000008747,0.0000008893,0.0000008933,0.0000008979,0.0000008937,0.0000009038,0.0000008920,0.0000009007,0.0000009040,0.0000009400,0.0000009094,0.0000009326,0.0000009129,0.0000009198,0.0000009216,0.0000004141,0.0000004693,0.0000004966,0.0000005422,0.0000005543,0.0000005960,0.0000006035,0.0000006267,0.0000006313,0.0000006511,0.0000006712,0.0000006662,0.0000006867,0.0000007074,0.0000007274,0.0000007239,0.0000007205,0.0000007342,0.0000007504,0.0000007511,0.0000007725,0.0000007779,0.0000007611,0.0000007815,0.0000007897,0.0000007949,0.0000008053,0.0000007817,0.0000008100,0.0000007926,0.0000008020,0.0000008338,0.0000008327,0.0000008291,0.0000008317,0.0000008437,0.0000008514,0.0000008554,0.0000008492,0.0000008593,0.0000008471,0.0000008534,0.0000008588,0.0000008967,0.0000008676,0.0000008867,0.0000008664,0.0000008721,0.0000008762,0.0000003929,0.0000004464,0.0000004698,0.0000005106,0.0000005259,0.0000005648,0.0000005748,0.0000005949,0.0000005999,0.0000006197,0.0000006399,0.0000006337,0.0000006505,0.0000006726,0.0000006895,0.0000006909,0.0000006838,0.0000006974,0.0000007111,0.0000007130,0.0000007334,0.0000007383,0.0000007190,0.0000007415,0.0000007528,0.0000007582,0.0000007646,0.0000007438,0.0000007704,0.0000007493,0.0000007637,0.0000007925,0.0000007910,0.0000007854,0.0000007868,0.0000008005,0.0000008093,0.0000008133,0.0000008068,0.0000008140,0.0000008056,0.0000008118,0.0000008170,0.0000008537,0.0000008254,0.0000008402,0.0000008262,0.0000008286,0.0000008353,0.0000003731,0.0000004241,0.0000004474,0.0000004842,0.0000005018,0.0000005364,0.0000005497,0.0000005595,0.0000005706,0.0000005877,0.0000006098,0.0000006019,0.0000006184,0.0000006376,0.0000006559,0.0000006557,0.0000006503,0.0000006639,0.0000006759,0.0000006769,0.0000007001,0.0000006985,0.0000006818,0.0000007036,0.0000007170,0.0000007205,0.0000007261,0.0000007071,0.0000007308,0.0000007121,0.0000007216,0.0000007507,0.0000007500,0.0000007456,0.0000007438,0.0000007589,0.0000007702,0.0000007730,0.0000007640,0.0000007726,0.0000007647,0.0000007742,0.0000007781,0.0000008107,0.0000007803,0.0000007976,0.0000007824,0.0000007880,0.0000007953,0.0000003543,0.0000004018,0.0000004257,0.0000004590,0.0000004751,0.0000005115,0.0000005202,0.0000005302,0.0000005421,0.0000005575,0.0000005816,0.0000005719,0.0000005862,0.0000006070,0.0000006210,0.0000006215,0.0000006160,0.0000006295,0.0000006399,0.0000006442,0.0000006662,0.0000006664,0.0000006487,0.0000006660,0.0000006817,0.0000006849,0.0000006901,0.0000006737,0.0000006937,0.0000006775,0.0000006842,0.0000007127,0.0000007130,0.0000007082,0.0000007099,0.0000007155,0.0000007319,0.0000007365,0.0000007241,0.0000007317,0.0000007262,0.0000007367,0.0000007362,0.0000007706,0.0000007403,0.0000007571,0.0000007440,0.0000007451,0.0000007579,0.0000003360,0.0000003814,0.0000004040,0.0000004372,0.0000004486,0.0000004860,0.0000004981,0.0000005024,0.0000005143,0.0000005300,0.0000005498,0.0000005450,0.0000005569,0.0000005771,0.0000005897,0.0000005902,0.0000005816,0.0000005979,0.0000006093,0.0000006146,0.0000006307,0.0000006321,0.0000006161,0.0000006317,0.0000006454,0.0000006491,0.0000006528,0.0000006396,0.0000006613,0.0000006425,0.0000006510,0.0000006781,0.0000006742,0.0000006736,0.0000006697,0.0000006817,0.0000006961,0.0000007001,0.0000006897,0.0000006954,0.0000006910,0.0000006964,0.0000007031,0.0000007316,0.0000007023,0.0000007156,0.0000007040,0.0000007131,0.0000007196,0.0000003190,0.0000003630,0.0000003847,0.0000004171,0.0000004268,0.0000004612,0.0000004730,0.0000004785,0.0000004863,0.0000005030,0.0000005175,0.0000005149,0.0000005322,0.0000005454,0.0000005596,0.0000005602,0.0000005510,0.0000005680,0.0000005803,0.0000005848,0.0000006008,0.0000006016,0.0000005864,0.0000006030,0.0000006138,0.0000006173,0.0000006224,0.0000006065,0.0000006257,0.0000006118,0.0000006202,0.0000006444,0.0000006375,0.0000006361,0.0000006364,0.0000006453,0.0000006602,0.0000006664,0.0000006560,0.0000006634,0.0000006540,0.0000006631,0.0000006686,0.0000006977,0.0000006656,0.0000006807,0.0000006691,0.0000006753,0.0000006849,0.0000003032,0.0000003471,0.0000003659,0.0000003971,0.0000004066,0.0000004405,0.0000004507,0.0000004518,0.0000004592,0.0000004768,0.0000004906,0.0000004904,0.0000005080,0.0000005181,0.0000005296,0.0000005345,0.0000005240,0.0000005398,0.0000005516,0.0000005539,0.0000005717,0.0000005704,0.0000005577,0.0000005744,0.0000005847,0.0000005888,0.0000005933,0.0000005788,0.0000005934,0.0000005808,0.0000005889,0.0000006150,0.0000006050,0.0000006060,0.0000006031,0.0000006138,0.0000006259,0.0000006350,0.0000006240,0.0000006281,0.0000006193,0.0000006322,0.0000006344,0.0000006644,0.0000006304,0.0000006473,0.0000006371,0.0000006414,0.0000006460,0.0000002890,0.0000003293,0.0000003464,0.0000003788,0.0000003863,0.0000004163,0.0000004289,0.0000004266,0.0000004348,0.0000004535,0.0000004659,0.0000004646,0.0000004840,0.0000004927,0.0000005032,0.0000005086,0.0000004972,0.0000005104,0.0000005225,0.0000005244,0.0000005417,0.0000005404,0.0000005295,0.0000005476,0.0000005535,0.0000005596,0.0000005620,0.0000005495,0.0000005646,0.0000005519,0.0000005623,0.0000005829,0.0000005739,0.0000005759,0.0000005733,0.0000005823,0.0000005949,0.0000006026,0.0000005945,0.0000005957,0.0000005885,0.0000006023,0.0000006044,0.0000006328,0.0000005972,0.0000006142,0.0000006071,0.0000006113,0.0000006110,0.0000002753,0.0000003125,0.0000003292,0.0000003585,0.0000003666,0.0000003971,0.0000004092,0.0000004041,0.0000004143,0.0000004317,0.0000004444,0.0000004399,0.0000004567,0.0000004656,0.0000004785,0.0000004834,0.0000004725,0.0000004865,0.0000004964,0.0000005005,0.0000005177,0.0000005147,0.0000005029,0.0000005171,0.0000005256,0.0000005313,0.0000005346,0.0000005230,0.0000005351,0.0000005228,0.0000005344,0.0000005535,0.0000005449,0.0000005466,0.0000005433,0.0000005528,0.0000005681,0.0000005728,0.0000005650,0.0000005655,0.0000005616,0.0000005726,0.0000005755,0.0000006015,0.0000005700,0.0000005855,0.0000005769,0.0000005818,0.0000005783,0.0000002610,0.0000002962,0.0000003148,0.0000003418,0.0000003483,0.0000003777,0.0000003895,0.0000003874,0.0000003936,0.0000004107,0.0000004194,0.0000004161,0.0000004335,0.0000004427,0.0000004554,0.0000004595,0.0000004504,0.0000004615,0.0000004734,0.0000004755,0.0000004889,0.0000004891,0.0000004765,0.0000004918,0.0000004988,0.0000005072,0.0000005085,0.0000004978,0.0000005078,0.0000004978,0.0000005081,0.0000005282,0.0000005166,0.0000005206,0.0000005176,0.0000005241,0.0000005407,0.0000005448,0.0000005372,0.0000005374,0.0000005344,0.0000005438,0.0000005449,0.0000005732,0.0000005427,0.0000005573,0.0000005474,0.0000005530,0.0000005502,0.0000002475,0.0000002822,0.0000002974,0.0000003224,0.0000003285,0.0000003603,0.0000003699,0.0000003680,0.0000003750,0.0000003879,0.0000003982,0.0000003981,0.0000004115,0.0000004226,0.0000004311,0.0000004358,0.0000004268,0.0000004381,0.0000004515,0.0000004515,0.0000004640,0.0000004667,0.0000004505,0.0000004687,0.0000004734,0.0000004827,0.0000004811,0.0000004720,0.0000004828,0.0000004727,0.0000004815,0.0000005011,0.0000004900,0.0000004962,0.0000004929,0.0000004992,0.0000005135,0.0000005175,0.0000005123,0.0000005090,0.0000005063,0.0000005172,0.0000005153,0.0000005446,0.0000005145,0.0000005317,0.0000005216,0.0000005250,0.0000005230,0.0000002344,0.0000002682,0.0000002818,0.0000003058,0.0000003141,0.0000003418,0.0000003508,0.0000003492,0.0000003560,0.0000003676,0.0000003808,0.0000003778,0.0000003906,0.0000004029,0.0000004111,0.0000004156,0.0000004056,0.0000004192,0.0000004297,0.0000004294,0.0000004396,0.0000004424,0.0000004279,0.0000004457,0.0000004499,0.0000004596,0.0000004581,0.0000004475,0.0000004566,0.0000004481,0.0000004552,0.0000004766,0.0000004682,0.0000004697,0.0000004697,0.0000004753,0.0000004861,0.0000004919,0.0000004862,0.0000004855,0.0000004782,0.0000004912,0.0000004856,0.0000005177,0.0000004907,0.0000005078,0.0000004951,0.0000004999,0.0000004984,0.0000002230,0.0000002517,0.0000002666,0.0000002904,0.0000002971,0.0000003268,0.0000003337,0.0000003316,0.0000003366,0.0000003499,0.0000003590,0.0000003599,0.0000003710,0.0000003829,0.0000003907,0.0000003951,0.0000003840,0.0000003960,0.0000004092,0.0000004088,0.0000004164,0.0000004206,0.0000004069,0.0000004232,0.0000004299,0.0000004380,0.0000004367,0.0000004269,0.0000004327,0.0000004270,0.0000004322,0.0000004533,0.0000004458,0.0000004437,0.0000004485,0.0000004535,0.0000004648,0.0000004674,0.0000004620,0.0000004605,0.0000004516,0.0000004663,0.0000004588,0.0000004945,0.0000004699,0.0000004842,0.0000004689,0.0000004752,0.0000004741,0.0000002102,0.0000002383,0.0000002529,0.0000002758,0.0000002779,0.0000003110,0.0000003175,0.0000003135,0.0000003205,0.0000003318,0.0000003425,0.0000003408,0.0000003522,0.0000003617,0.0000003707,0.0000003766,0.0000003665,0.0000003751,0.0000003898,0.0000003876,0.0000003963,0.0000003996,0.0000003848,0.0000004034,0.0000004091,0.0000004166,0.0000004137,0.0000004051,0.0000004087,0.0000004070,0.0000004100,0.0000004292,0.0000004251,0.0000004227,0.0000004248,0.0000004318,0.0000004419,0.0000004420,0.0000004388,0.0000004359,0.0000004296,0.0000004419,0.0000004356,0.0000004692,0.0000004452,0.0000004608,0.0000004476,0.0000004499,0.0000004510,0.0000001991,0.0000002275,0.0000002411,0.0000002634,0.0000002630,0.0000002949,0.0000002999,0.0000002985,0.0000003025,0.0000003155,0.0000003235,0.0000003218,0.0000003347,0.0000003455,0.0000003531,0.0000003591,0.0000003470,0.0000003553,0.0000003696,0.0000003688,0.0000003772,0.0000003777,0.0000003654,0.0000003819,0.0000003901,0.0000003939,0.0000003939,0.0000003849,0.0000003896,0.0000003856,0.0000003886,0.0000004103,0.0000004047,0.0000004004,0.0000004059,0.0000004078,0.0000004206,0.0000004201,0.0000004169,0.0000004136,0.0000004081,0.0000004205,0.0000004153,0.0000004452,0.0000004216,0.0000004382,0.0000004242,0.0000004271,0.0000004300,0.0000001883,0.0000002155,0.0000002292,0.0000002500,0.0000002472,0.0000002796,0.0000002850,0.0000002851,0.0000002876,0.0000002995,0.0000003072,0.0000003068,0.0000003184,0.0000003303,0.0000003344,0.0000003414,0.0000003282,0.0000003372,0.0000003508,0.0000003504,0.0000003610,0.0000003586,0.0000003483,0.0000003621,0.0000003709,0.0000003755,0.0000003746,0.0000003651,0.0000003696,0.0000003664,0.0000003708,0.0000003920,0.0000003821,0.0000003820,0.0000003874,0.0000003878,0.0000003999,0.0000003996,0.0000003973,0.0000003925,0.0000003870,0.0000003987,0.0000003950,0.0000004223,0.0000003982,0.0000004172,0.0000004020,0.0000004060,0.0000004079,0.0000001793,0.0000002066,0.0000002164,0.0000002374,0.0000002365,0.0000002648,0.0000002729,0.0000002709,0.0000002749,0.0000002853,0.0000002912,0.0000002916,0.0000003033,0.0000003150,0.0000003177,0.0000003248,0.0000003119,0.0000003206,0.0000003330,0.0000003312,0.0000003428,0.0000003423,0.0000003302,0.0000003420,0.0000003553,0.0000003563,0.0000003582,0.0000003475,0.0000003529,0.0000003487,0.0000003520,0.0000003717,0.0000003646,0.0000003617,0.0000003681,0.0000003666,0.0000003809,0.0000003789,0.0000003774,0.0000003752,0.0000003704,0.0000003786,0.0000003746,0.0000004002,0.0000003760,0.0000003965,0.0000003795,0.0000003859,0.0000003850,0.0000001710,0.0000001966,0.0000002051,0.0000002290,0.0000002242,0.0000002519,0.0000002596,0.0000002575,0.0000002621,0.0000002715,0.0000002762,0.0000002778,0.0000002868,0.0000003010,0.0000003025,0.0000003082,0.0000002950,0.0000003048,0.0000003138,0.0000003148,0.0000003245,0.0000003262,0.0000003136,0.0000003262,0.0000003369,0.0000003388,0.0000003391,0.0000003297,0.0000003364,0.0000003315,0.0000003334,0.0000003550,0.0000003469,0.0000003416,0.0000003506,0.0000003477,0.0000003605,0.0000003587,0.0000003601,0.0000003543,0.0000003515,0.0000003627,0.0000003555,0.0000003813,0.0000003577,0.0000003786,0.0000003579,0.0000003673,0.0000003680,0.0000001634,0.0000001864,0.0000001946,0.0000002173,0.0000002137,0.0000002390,0.0000002475,0.0000002474,0.0000002487,0.0000002574,0.0000002623,0.0000002618,0.0000002719,0.0000002856,0.0000002866,0.0000002920,0.0000002810,0.0000002885,0.0000002994,0.0000002981,0.0000003079,0.0000003086,0.0000002994,0.0000003092,0.0000003204,0.0000003206,0.0000003222,0.0000003127,0.0000003179,0.0000003139,0.0000003184,0.0000003368,0.0000003305,0.0000003237,0.0000003352,0.0000003287,0.0000003443,0.0000003425,0.0000003414,0.0000003371,0.0000003340,0.0000003435,0.0000003394,0.0000003633,0.0000003407,0.0000003619,0.0000003409,0.0000003485,0.0000003482,0.0000001536,0.0000001788,0.0000001854,0.0000002072,0.0000002050,0.0000002271,0.0000002349,0.0000002349,0.0000002361,0.0000002442,0.0000002515,0.0000002501,0.0000002596,0.0000002720,0.0000002748,0.0000002771,0.0000002656,0.0000002738,0.0000002854,0.0000002821,0.0000002940,0.0000002934,0.0000002843,0.0000002933,0.0000003061,0.0000003022,0.0000003064,0.0000002955,0.0000003024,0.0000002981,0.0000003005,0.0000003209,0.0000003146,0.0000003084,0.0000003192,0.0000003136,0.0000003283,0.0000003239,0.0000003255,0.0000003224,0.0000003168,0.0000003243,0.0000003251,0.0000003459,0.0000003243,0.0000003447,0.0000003215,0.0000003305,0.0000003326,0.0000001464,0.0000001687,0.0000001772,0.0000001954,0.0000001949,0.0000002157,0.0000002218,0.0000002232,0.0000002235,0.0000002309,0.0000002385,0.0000002368,0.0000002472,0.0000002572,0.0000002595,0.0000002635,0.0000002532,0.0000002615,0.0000002701,0.0000002683,0.0000002793,0.0000002786,0.0000002668,0.0000002764,0.0000002916,0.0000002882,0.0000002930,0.0000002803,0.0000002863,0.0000002833,0.0000002839,0.0000003043,0.0000002984,0.0000002924,0.0000003046,0.0000002986,0.0000003138,0.0000003069,0.0000003088,0.0000003087,0.0000003012,0.0000003065,0.0000003065,0.0000003297,0.0000003070,0.0000003286,0.0000003066,0.0000003134,0.0000003170,0.0000001373,0.0000001597,0.0000001683,0.0000001852,0.0000001862,0.0000002064,0.0000002100,0.0000002118,0.0000002113,0.0000002200,0.0000002273,0.0000002261,0.0000002349,0.0000002447,0.0000002477,0.0000002509,0.0000002415,0.0000002471,0.0000002532,0.0000002568,0.0000002668,0.0000002637,0.0000002559,0.0000002642,0.0000002758,0.0000002728,0.0000002782,0.0000002639,0.0000002706,0.0000002694,0.0000002691,0.0000002871,0.0000002849,0.0000002776,0.0000002908,0.0000002830,0.0000003000,0.0000002914,0.0000002933,0.0000002937,0.0000002876,0.0000002924,0.0000002892,0.0000003135,0.0000002909,0.0000003109,0.0000002905,0.0000002977,0.0000003026,0.0000001302,0.0000001524,0.0000001609,0.0000001749,0.0000001759,0.0000001944,0.0000002001,0.0000002013,0.0000001995,0.0000002083,0.0000002161,0.0000002156,0.0000002238,0.0000002310,0.0000002380,0.0000002388,0.0000002283,0.0000002358,0.0000002402,0.0000002447,0.0000002549,0.0000002507,0.0000002433,0.0000002503,0.0000002608,0.0000002598,0.0000002652,0.0000002508,0.0000002575,0.0000002541,0.0000002574,0.0000002730,0.0000002726,0.0000002636,0.0000002767,0.0000002667,0.0000002862,0.0000002774,0.0000002798,0.0000002796,0.0000002721,0.0000002784,0.0000002755,0.0000002982,0.0000002738,0.0000002960,0.0000002755,0.0000002830,0.0000002882,0.0000001232,0.0000001438,0.0000001521,0.0000001664,0.0000001647,0.0000001851,0.0000001888,0.0000001908,0.0000001899,0.0000001968,0.0000002044,0.0000002027,0.0000002124,0.0000002198,0.0000002271,0.0000002268,0.0000002170,0.0000002257,0.0000002270,0.0000002333,0.0000002421,0.0000002369,0.0000002299,0.0000002389,0.0000002492,0.0000002477,0.0000002528,0.0000002399,0.0000002458,0.0000002414,0.0000002452,0.0000002596,0.0000002586,0.0000002508,0.0000002608,0.0000002538,0.0000002702,0.0000002652,0.0000002660,0.0000002652,0.0000002590,0.0000002648,0.0000002604,0.0000002829,0.0000002597,0.0000002825,0.0000002617,0.0000002684,0.0000002748,0.0000001162,0.0000001364,0.0000001457,0.0000001585,0.0000001573,0.0000001770,0.0000001809,0.0000001811,0.0000001794,0.0000001863,0.0000001930,0.0000001927,0.0000002025,0.0000002086,0.0000002147,0.0000002160,0.0000002065,0.0000002146,0.0000002160,0.0000002201,0.0000002309,0.0000002254,0.0000002190,0.0000002271,0.0000002358,0.0000002345,0.0000002402,0.0000002291,0.0000002345,0.0000002299,0.0000002327,0.0000002462,0.0000002446,0.0000002383,0.0000002479,0.0000002409,0.0000002564,0.0000002510,0.0000002529,0.0000002535,0.0000002464,0.0000002534,0.0000002472,0.0000002692,0.0000002460,0.0000002675,0.0000002492,0.0000002570,0.0000002609,0.0000001089,0.0000001298,0.0000001385,0.0000001508,0.0000001507,0.0000001683,0.0000001713,0.0000001713,0.0000001696,0.0000001760,0.0000001827,0.0000001832,0.0000001916,0.0000001993,0.0000002060,0.0000002053,0.0000001941,0.0000002045,0.0000002050,0.0000002092,0.0000002205,0.0000002149,0.0000002073,0.0000002162,0.0000002238,0.0000002242,0.0000002277,0.0000002183,0.0000002238,0.0000002185,0.0000002211,0.0000002350,0.0000002326,0.0000002256,0.0000002357,0.0000002304,0.0000002422,0.0000002404,0.0000002419,0.0000002416,0.0000002331,0.0000002408,0.0000002356,0.0000002549,0.0000002331,0.0000002536,0.0000002373,0.0000002440,0.0000002491,0.0000001034,0.0000001240,0.0000001315,0.0000001430,0.0000001445,0.0000001604,0.0000001641,0.0000001623,0.0000001626,0.0000001679,0.0000001745,0.0000001752,0.0000001810,0.0000001893,0.0000001949,0.0000001941,0.0000001850,0.0000001923,0.0000001959,0.0000002001,0.0000002097,0.0000002033,0.0000001972,0.0000002060,0.0000002126,0.0000002127,0.0000002155,0.0000002054,0.0000002145,0.0000002089,0.0000002103,0.0000002237,0.0000002187,0.0000002148,0.0000002251,0.0000002172,0.0000002306,0.0000002284,0.0000002314,0.0000002302,0.0000002227,0.0000002301,0.0000002249,0.0000002433,0.0000002201,0.0000002422,0.0000002234,0.0000002331,0.0000002346,0.0000000988,0.0000001179,0.0000001240,0.0000001356,0.0000001371,0.0000001536,0.0000001556,0.0000001531,0.0000001544,0.0000001584,0.0000001658,0.0000001670,0.0000001727,0.0000001820,0.0000001854,0.0000001835,0.0000001750,0.0000001831,0.0000001842,0.0000001909,0.0000001985,0.0000001940,0.0000001875,0.0000001970,0.0000002034,0.0000002028,0.0000002038,0.0000001950,0.0000002067,0.0000002000,0.0000001994,0.0000002128,0.0000002083,0.0000002040,0.0000002128,0.0000002077,0.0000002211,0.0000002167,0.0000002204,0.0000002179,0.0000002115,0.0000002199,0.0000002155,0.0000002327,0.0000002086,0.0000002287,0.0000002129,0.0000002221,0.0000002221,0.0000000943,0.0000001109,0.0000001179,0.0000001300,0.0000001306,0.0000001454,0.0000001479,0.0000001442,0.0000001476,0.0000001487,0.0000001568,0.0000001583,0.0000001633,0.0000001733,0.0000001766,0.0000001745,0.0000001658,0.0000001740,0.0000001740,0.0000001821,0.0000001878,0.0000001853,0.0000001764,0.0000001873,0.0000001933,0.0000001935,0.0000001933,0.0000001861,0.0000001960,0.0000001898,0.0000001889,0.0000002018,0.0000001963,0.0000001934,0.0000002019,0.0000001971,0.0000002104,0.0000002045,0.0000002097,0.0000002062,0.0000002025,0.0000002074,0.0000002057,0.0000002209,0.0000001981,0.0000002149,0.0000002012,0.0000002118,0.0000002134,0.0000000896,0.0000001040,0.0000001114,0.0000001239,0.0000001227,0.0000001369,0.0000001401,0.0000001368,0.0000001410,0.0000001411,0.0000001476,0.0000001505,0.0000001560,0.0000001636,0.0000001671,0.0000001660,0.0000001576,0.0000001652,0.0000001664,0.0000001718,0.0000001794,0.0000001767,0.0000001685,0.0000001774,0.0000001828,0.0000001835,0.0000001833,0.0000001778,0.0000001863,0.0000001818,0.0000001796,0.0000001916,0.0000001863,0.0000001833,0.0000001915,0.0000001856,0.0000002005,0.0000001963,0.0000001991,0.0000001969,0.0000001924,0.0000001971,0.0000001942,0.0000002090,0.0000001899,0.0000002029,0.0000001922,0.0000002011,0.0000002010,0.0000000847,0.0000000980,0.0000001056,0.0000001172,0.0000001185,0.0000001304,0.0000001320,0.0000001306,0.0000001318,0.0000001328,0.0000001406,0.0000001424,0.0000001472,0.0000001534,0.0000001577,0.0000001576,0.0000001493,0.0000001579,0.0000001571,0.0000001620,0.0000001710,0.0000001677,0.0000001605,0.0000001661,0.0000001744,0.0000001748,0.0000001746,0.0000001701,0.0000001768,0.0000001746,0.0000001703,0.0000001815,0.0000001779,0.0000001735,0.0000001819,0.0000001758,0.0000001921,0.0000001857,0.0000001886,0.0000001882,0.0000001830,0.0000001867,0.0000001847,0.0000001989,0.0000001806,0.0000001930,0.0000001815,0.0000001907,0.0000001914,0.0000000799,0.0000000934,0.0000001000,0.0000001122,0.0000001115,0.0000001234,0.0000001256,0.0000001246,0.0000001255,0.0000001257,0.0000001333,0.0000001357,0.0000001391,0.0000001444,0.0000001505,0.0000001498,0.0000001406,0.0000001497,0.0000001500,0.0000001520,0.0000001620,0.0000001594,0.0000001513,0.0000001559,0.0000001653,0.0000001664,0.0000001661,0.0000001624,0.0000001687,0.0000001668,0.0000001616,0.0000001740,0.0000001693,0.0000001642,0.0000001724,0.0000001663,0.0000001835,0.0000001752,0.0000001789,0.0000001809,0.0000001750,0.0000001768,0.0000001749,0.0000001883,0.0000001715,0.0000001837,0.0000001742,0.0000001793,0.0000001817,0.0000000747,0.0000000890,0.0000000958,0.0000001052,0.0000001056,0.0000001167,0.0000001205,0.0000001194,0.0000001201,0.0000001199,0.0000001272,0.0000001288,0.0000001329,0.0000001367,0.0000001439,0.0000001420,0.0000001338,0.0000001417,0.0000001428,0.0000001450,0.0000001529,0.0000001510,0.0000001435,0.0000001491,0.0000001568,0.0000001582,0.0000001594,0.0000001543,0.0000001604,0.0000001574,0.0000001534,0.0000001651,0.0000001617,0.0000001566,0.0000001643,0.0000001586,0.0000001751,0.0000001655,0.0000001692,0.0000001707,0.0000001658,0.0000001676,0.0000001654,0.0000001803,0.0000001638,0.0000001731,0.0000001674,0.0000001718,0.0000001731,0.0000000718,0.0000000844,0.0000000913,0.0000001017,0.0000001000,0.0000001109,0.0000001156,0.0000001132,0.0000001142,0.0000001138,0.0000001202,0.0000001225,0.0000001276,0.0000001299,0.0000001379,0.0000001336,0.0000001269,0.0000001343,0.0000001361,0.0000001383,0.0000001454,0.0000001417,0.0000001348,0.0000001411,0.0000001496,0.0000001484,0.0000001512,0.0000001464,0.0000001539,0.0000001504,0.0000001453,0.0000001555,0.0000001548,0.0000001499,0.0000001541,0.0000001501,0.0000001664,0.0000001575,0.0000001614,0.0000001613,0.0000001575,0.0000001596,0.0000001573,0.0000001707,0.0000001546,0.0000001637,0.0000001591,0.0000001619,0.0000001648,0.0000000683,0.0000000795,0.0000000863,0.0000000962,0.0000000965,0.0000001054,0.0000001072,0.0000001074,0.0000001088,0.0000001078,0.0000001133,0.0000001168,0.0000001211,0.0000001231,0.0000001320,0.0000001274,0.0000001210,0.0000001277,0.0000001300,0.0000001307,0.0000001381,0.0000001341,0.0000001285,0.0000001341,0.0000001427,0.0000001409,0.0000001438,0.0000001377,0.0000001468,0.0000001428,0.0000001380,0.0000001476,0.0000001483,0.0000001437,0.0000001471,0.0000001430,0.0000001576,0.0000001500,0.0000001521,0.0000001530,0.0000001492,0.0000001517,0.0000001501,0.0000001617,0.0000001475,0.0000001561,0.0000001513,0.0000001539,0.0000001554,0.0000000648,0.0000000758,0.0000000824,0.0000000917,0.0000000918,0.0000000993,0.0000001023,0.0000001030,0.0000001036,0.0000001012,0.0000001078,0.0000001112,0.0000001150,0.0000001149,0.0000001264,0.0000001209,0.0000001148,0.0000001216,0.0000001233,0.0000001232,0.0000001306,0.0000001270,0.0000001231,0.0000001270,0.0000001348,0.0000001336,0.0000001363,0.0000001314,0.0000001379,0.0000001362,0.0000001302,0.0000001388,0.0000001397,0.0000001358,0.0000001408,0.0000001342,0.0000001499,0.0000001437,0.0000001439,0.0000001449,0.0000001428,0.0000001424,0.0000001419,0.0000001528,0.0000001396,0.0000001474,0.0000001445,0.0000001481,0.0000001476,0.0000000620,0.0000000719,0.0000000784,0.0000000874,0.0000000876,0.0000000932,0.0000000966,0.0000000972,0.0000000994,0.0000000953,0.0000001029,0.0000001052,0.0000001085,0.0000001087,0.0000001201,0.0000001140,0.0000001088,0.0000001159,0.0000001162,0.0000001171,0.0000001254,0.0000001210,0.0000001161,0.0000001208,0.0000001281,0.0000001264,0.0000001306,0.0000001245,0.0000001322,0.0000001289,0.0000001231,0.0000001317,0.0000001329,0.0000001281,0.0000001334,0.0000001278,0.0000001426,0.0000001378,0.0000001369,0.0000001371,0.0000001369,0.0000001350,0.0000001349,0.0000001459,0.0000001328,0.0000001414,0.0000001380,0.0000001422,0.0000001402,0.0000000598,0.0000000686,0.0000000741,0.0000000838,0.0000000829,0.0000000893,0.0000000916,0.0000000915,0.0000000947,0.0000000913,0.0000000984,0.0000000994,0.0000001030,0.0000001040,0.0000001153,0.0000001068,0.0000001038,0.0000001096,0.0000001109,0.0000001120,0.0000001194,0.0000001157,0.0000001106,0.0000001144,0.0000001215,0.0000001208,0.0000001249,0.0000001192,0.0000001260,0.0000001233,0.0000001159,0.0000001247,0.0000001277,0.0000001222,0.0000001273,0.0000001217,0.0000001353,0.0000001304,0.0000001299,0.0000001306,0.0000001303,0.0000001273,0.0000001283,0.0000001370,0.0000001264,0.0000001341,0.0000001307,0.0000001349,0.0000001350,0.0000000572,0.0000000661,0.0000000708,0.0000000790,0.0000000801,0.0000000844,0.0000000880,0.0000000876,0.0000000898,0.0000000853,0.0000000945,0.0000000947,0.0000000979,0.0000000980,0.0000001086,0.0000000998,0.0000000992,0.0000001050,0.0000001059,0.0000001058,0.0000001134,0.0000001100,0.0000001055,0.0000001089,0.0000001161,0.0000001148,0.0000001178,0.0000001134,0.0000001192,0.0000001173,0.0000001082,0.0000001185,0.0000001225,0.0000001162,0.0000001218,0.0000001159,0.0000001277,0.0000001232,0.0000001236,0.0000001226,0.0000001244,0.0000001213,0.0000001200,0.0000001293,0.0000001197,0.0000001266,0.0000001232,0.0000001292,0.0000001270,0.0000000548,0.0000000634,0.0000000659,0.0000000756,0.0000000768,0.0000000803,0.0000000835,0.0000000832,0.0000000859,0.0000000812,0.0000000900,0.0000000902,0.0000000925,0.0000000924,0.0000001037,0.0000000940,0.0000000948,0.0000000999,0.0000001009,0.0000000992,0.0000001074,0.0000001053,0.0000000998,0.0000001033,0.0000001098,0.0000001087,0.0000001113,0.0000001086,0.0000001140,0.0000001124,0.0000001021,0.0000001115,0.0000001153,0.0000001097,0.0000001171,0.0000001091,0.0000001219,0.0000001144,0.0000001175,0.0000001169,0.0000001167,0.0000001159,0.0000001144,0.0000001232,0.0000001136,0.0000001197,0.0000001177,0.0000001236,0.0000001195,0.0000000519,0.0000000607,0.0000000612,0.0000000717,0.0000000735,0.0000000766,0.0000000796,0.0000000794,0.0000000820,0.0000000769,0.0000000844,0.0000000854,0.0000000893,0.0000000888,0.0000000978,0.0000000889,0.0000000897,0.0000000945,0.0000000965,0.0000000936,0.0000001028,0.0000001006,0.0000000947,0.0000000976,0.0000001038,0.0000001036,0.0000001063,0.0000001033,0.0000001090,0.0000001069,0.0000000967,0.0000001070,0.0000001096,0.0000001062,0.0000001119,0.0000001045,0.0000001152,0.0000001091,0.0000001112,0.0000001110,0.0000001117,0.0000001105,0.0000001093,0.0000001169,0.0000001075,0.0000001126,0.0000001128,0.0000001169,0.0000001143,0.0000000499,0.0000000582,0.0000000585,0.0000000677,0.0000000683,0.0000000735,0.0000000752,0.0000000751,0.0000000785,0.0000000732,0.0000000801,0.0000000807,0.0000000844,0.0000000850,0.0000000908,0.0000000855,0.0000000858,0.0000000901,0.0000000919,0.0000000890,0.0000000977,0.0000000951,0.0000000896,0.0000000927,0.0000000979,0.0000000993,0.0000001023,0.0000000997,0.0000001025,0.0000001022,0.0000000921,0.0000001031,0.0000001039,0.0000001012,0.0000001062,0.0000000994,0.0000001097,0.0000001038,0.0000001072,0.0000001060,0.0000001072,0.0000001059,0.0000001033,0.0000001104,0.0000001031,0.0000001063,0.0000001071,0.0000001108,0.0000001078,0.0000000483,0.0000000547,0.0000000555,0.0000000629,0.0000000651,0.0000000696,0.0000000723,0.0000000706,0.0000000745,0.0000000687,0.0000000749,0.0000000761,0.0000000808,0.0000000799,0.0000000860,0.0000000821,0.0000000812,0.0000000853,0.0000000881,0.0000000847,0.0000000928,0.0000000900,0.0000000858,0.0000000878,0.0000000941,0.0000000933,0.0000000963,0.0000000946,0.0000000978,0.0000000970,0.0000000874,0.0000000970,0.0000000993,0.0000000962,0.0000001005,0.0000000939,0.0000001030,0.0000000984,0.0000001030,0.0000001017,0.0000001006,0.0000001003,0.0000000988,0.0000001048,0.0000000983,0.0000001018,0.0000001011,0.0000001055,0.0000001028,0.0000000461,0.0000000507,0.0000000525,0.0000000584,0.0000000614,0.0000000659,0.0000000689,0.0000000674,0.0000000703,0.0000000653,0.0000000711,0.0000000727,0.0000000778,0.0000000762,0.0000000822,0.0000000770,0.0000000768,0.0000000812,0.0000000841,0.0000000796,0.0000000876,0.0000000849,0.0000000806,0.0000000835,0.0000000888,0.0000000894,0.0000000921,0.0000000899,0.0000000930,0.0000000918,0.0000000841,0.0000000923,0.0000000953,0.0000000916,0.0000000962,0.0000000896,0.0000000991,0.0000000927,0.0000000983,0.0000000951,0.0000000952,0.0000000953,0.0000000942,0.0000000997,0.0000000930,0.0000000951,0.0000000961,0.0000001006,0.0000000984,0.0000000439,0.0000000481,0.0000000499,0.0000000555,0.0000000585,0.0000000633,0.0000000637,0.0000000629,0.0000000667,0.0000000618,0.0000000674,0.0000000694,0.0000000729,0.0000000722,0.0000000772,0.0000000727,0.0000000734,0.0000000775,0.0000000796,0.0000000760,0.0000000828,0.0000000808,0.0000000760,0.0000000785,0.0000000850,0.0000000845,0.0000000877,0.0000000844,0.0000000880,0.0000000865,0.0000000813,0.0000000861,0.0000000897,0.0000000886,0.0000000909,0.0000000856,0.0000000948,0.0000000880,0.0000000940,0.0000000908,0.0000000887,0.0000000914,0.0000000891,0.0000000954,0.0000000893,0.0000000901,0.0000000919,0.0000000965,0.0000000942,0.0000000419,0.0000000453,0.0000000471,0.0000000531,0.0000000567,0.0000000601,0.0000000599,0.0000000600,0.0000000636,0.0000000582,0.0000000643,0.0000000666,0.0000000694,0.0000000694,0.0000000726,0.0000000689,0.0000000698,0.0000000733,0.0000000753,0.0000000720,0.0000000792,0.0000000766,0.0000000731,0.0000000741,0.0000000801,0.0000000818,0.0000000840,0.0000000791,0.0000000832,0.0000000829,0.0000000772,0.0000000813,0.0000000848,0.0000000835,0.0000000851,0.0000000813,0.0000000896,0.0000000838,0.0000000899,0.0000000860,0.0000000845,0.0000000863,0.0000000853,0.0000000908,0.0000000857,0.0000000860,0.0000000866,0.0000000909,0.0000000888,0.0000000402,0.0000000420,0.0000000445,0.0000000507,0.0000000540,0.0000000580,0.0000000560,0.0000000567,0.0000000604,0.0000000555,0.0000000615,0.0000000629,0.0000000659,0.0000000656,0.0000000693,0.0000000657,0.0000000660,0.0000000693,0.0000000717,0.0000000692,0.0000000753,0.0000000729,0.0000000689,0.0000000702,0.0000000748,0.0000000784,0.0000000808,0.0000000759,0.0000000801,0.0000000791,0.0000000729,0.0000000769,0.0000000814,0.0000000794,0.0000000817,0.0000000771,0.0000000848,0.0000000802,0.0000000855,0.0000000820,0.0000000808,0.0000000820,0.0000000818,0.0000000864,0.0000000818,0.0000000820,0.0000000816,0.0000000875,0.0000000849,0.0000000380,0.0000000397,0.0000000428,0.0000000486,0.0000000519,0.0000000542,0.0000000534,0.0000000537,0.0000000578,0.0000000527,0.0000000587,0.0000000595,0.0000000618,0.0000000621,0.0000000663,0.0000000622,0.0000000618,0.0000000653,0.0000000677,0.0000000654,0.0000000719,0.0000000694,0.0000000650,0.0000000663,0.0000000707,0.0000000750,0.0000000761,0.0000000727,0.0000000758,0.0000000744,0.0000000693,0.0000000728,0.0000000771,0.0000000755,0.0000000770,0.0000000725,0.0000000812,0.0000000768,0.0000000821,0.0000000771,0.0000000773,0.0000000783,0.0000000763,0.0000000819,0.0000000773,0.0000000766,0.0000000771,0.0000000833,0.0000000807,0.0000000368,0.0000000380,0.0000000409,0.0000000461,0.0000000492,0.0000000513,0.0000000514,0.0000000507,0.0000000550,0.0000000496,0.0000000557,0.0000000565,0.0000000586,0.0000000579,0.0000000631,0.0000000595,0.0000000581,0.0000000613,0.0000000636,0.0000000615,0.0000000677,0.0000000654,0.0000000620,0.0000000626,0.0000000669,0.0000000709,0.0000000725,0.0000000687,0.0000000720,0.0000000710,0.0000000655,0.0000000690,0.0000000730,0.0000000722,0.0000000727,0.0000000688,0.0000000774,0.0000000731,0.0000000782,0.0000000730,0.0000000730,0.0000000741,0.0000000722,0.0000000783,0.0000000731,0.0000000731,0.0000000726,0.0000000790,0.0000000763,0.0000000350,0.0000000365,0.0000000388,0.0000000432,0.0000000464,0.0000000488,0.0000000488,0.0000000483,0.0000000524,0.0000000473,0.0000000514,0.0000000540,0.0000000553,0.0000000552,0.0000000596,0.0000000554,0.0000000553,0.0000000580,0.0000000604,0.0000000595,0.0000000647,0.0000000626,0.0000000599,0.0000000595,0.0000000639,0.0000000674,0.0000000695,0.0000000661,0.0000000674,0.0000000666,0.0000000615,0.0000000654,0.0000000705,0.0000000683,0.0000000689,0.0000000654,0.0000000737,0.0000000703,0.0000000745,0.0000000698,0.0000000696,0.0000000694,0.0000000692,0.0000000735,0.0000000695,0.0000000695,0.0000000694,0.0000000749,0.0000000719,0.0000000343,0.0000000349,0.0000000369,0.0000000409,0.0000000443,0.0000000461,0.0000000468,0.0000000460,0.0000000498,0.0000000452,0.0000000492,0.0000000509,0.0000000524,0.0000000527,0.0000000577,0.0000000525,0.0000000523,0.0000000553,0.0000000578,0.0000000565,0.0000000617,0.0000000586,0.0000000570,0.0000000565,0.0000000609,0.0000000641,0.0000000657,0.0000000630,0.0000000637,0.0000000635,0.0000000590,0.0000000628,0.0000000673,0.0000000646,0.0000000654,0.0000000612,0.0000000694,0.0000000669,0.0000000693,0.0000000669,0.0000000649,0.0000000656,0.0000000650,0.0000000701,0.0000000667,0.0000000667,0.0000000649,0.0000000706,0.0000000671,0.0000000326,0.0000000336,0.0000000352,0.0000000380,0.0000000417,0.0000000431,0.0000000442,0.0000000440,0.0000000473,0.0000000423,0.0000000465,0.0000000485,0.0000000489,0.0000000509,0.0000000545,0.0000000489,0.0000000492,0.0000000527,0.0000000548,0.0000000533,0.0000000589,0.0000000551,0.0000000540,0.0000000532,0.0000000577,0.0000000604,0.0000000635,0.0000000593,0.0000000612,0.0000000604,0.0000000554,0.0000000600,0.0000000647,0.0000000614,0.0000000616,0.0000000584,0.0000000650,0.0000000635,0.0000000644,0.0000000643,0.0000000625,0.0000000623,0.0000000612,0.0000000671,0.0000000640,0.0000000627,0.0000000620,0.0000000671,0.0000000632,0.0000000314,0.0000000322,0.0000000331,0.0000000366,0.0000000394,0.0000000408,0.0000000417,0.0000000425,0.0000000456,0.0000000397,0.0000000444,0.0000000458,0.0000000465,0.0000000476,0.0000000513,0.0000000467,0.0000000468,0.0000000502,0.0000000521,0.0000000506,0.0000000565,0.0000000530,0.0000000506,0.0000000500,0.0000000536,0.0000000575,0.0000000611,0.0000000565,0.0000000579,0.0000000572,0.0000000524,0.0000000574,0.0000000608,0.0000000572,0.0000000585,0.0000000554,0.0000000613,0.0000000606,0.0000000616,0.0000000614,0.0000000597,0.0000000589,0.0000000581,0.0000000640,0.0000000618,0.0000000586,0.0000000582,0.0000000632,0.0000000611,0.0000000297,0.0000000307,0.0000000321,0.0000000340,0.0000000372,0.0000000391,0.0000000391,0.0000000408,0.0000000433,0.0000000373,0.0000000422,0.0000000442,0.0000000448,0.0000000448,0.0000000492,0.0000000440,0.0000000445,0.0000000485,0.0000000495,0.0000000474,0.0000000538,0.0000000511,0.0000000479,0.0000000464,0.0000000515,0.0000000539,0.0000000582,0.0000000540,0.0000000552,0.0000000542,0.0000000494,0.0000000540,0.0000000591,0.0000000545,0.0000000558,0.0000000527,0.0000000585,0.0000000579,0.0000000591,0.0000000579,0.0000000574,0.0000000558,0.0000000551,0.0000000613,0.0000000586,0.0000000559,0.0000000556,0.0000000605,0.0000000569,0.0000000283,0.0000000301,0.0000000307,0.0000000326,0.0000000351,0.0000000374,0.0000000382,0.0000000389,0.0000000416,0.0000000360,0.0000000400,0.0000000424,0.0000000428,0.0000000434,0.0000000477,0.0000000418,0.0000000420,0.0000000454,0.0000000473,0.0000000447,0.0000000519,0.0000000487,0.0000000456,0.0000000439,0.0000000489,0.0000000513,0.0000000548,0.0000000515,0.0000000522,0.0000000525,0.0000000467,0.0000000516,0.0000000552,0.0000000513,0.0000000534,0.0000000501,0.0000000560,0.0000000552,0.0000000567,0.0000000549,0.0000000543,0.0000000532,0.0000000520,0.0000000583,0.0000000570,0.0000000533,0.0000000530,0.0000000569,0.0000000539,0.0000000267,0.0000000286,0.0000000289,0.0000000310,0.0000000333,0.0000000358,0.0000000362,0.0000000369,0.0000000399,0.0000000351,0.0000000387,0.0000000398,0.0000000405,0.0000000408,0.0000000461,0.0000000397,0.0000000397,0.0000000431,0.0000000447,0.0000000430,0.0000000493,0.0000000460,0.0000000431,0.0000000419,0.0000000456,0.0000000487,0.0000000522,0.0000000496,0.0000000494,0.0000000495,0.0000000443,0.0000000484,0.0000000524,0.0000000488,0.0000000506,0.0000000479,0.0000000530,0.0000000516,0.0000000533,0.0000000528,0.0000000517,0.0000000507,0.0000000501,0.0000000551,0.0000000538,0.0000000508,0.0000000502,0.0000000548,0.0000000509,0.0000000261,0.0000000279,0.0000000273,0.0000000296,0.0000000322,0.0000000338,0.0000000346,0.0000000350,0.0000000383,0.0000000338,0.0000000362,0.0000000385,0.0000000383,0.0000000390,0.0000000439,0.0000000376,0.0000000380,0.0000000410,0.0000000420,0.0000000409,0.0000000469,0.0000000439,0.0000000404,0.0000000401,0.0000000443,0.0000000471,0.0000000485,0.0000000463,0.0000000464,0.0000000460,0.0000000422,0.0000000455,0.0000000498,0.0000000464,0.0000000477,0.0000000453,0.0000000504,0.0000000488,0.0000000513,0.0000000509,0.0000000492,0.0000000487,0.0000000463,0.0000000517,0.0000000514,0.0000000486,0.0000000471,0.0000000526,0.0000000488,0.0000000246,0.0000000267,0.0000000258,0.0000000283,0.0000000302,0.0000000321,0.0000000329,0.0000000331,0.0000000364,0.0000000324,0.0000000340,0.0000000363,0.0000000363,0.0000000368,0.0000000418,0.0000000358,0.0000000362,0.0000000393,0.0000000402,0.0000000387,0.0000000446,0.0000000417,0.0000000382,0.0000000382,0.0000000418,0.0000000449,0.0000000469,0.0000000442,0.0000000438,0.0000000432,0.0000000399,0.0000000436,0.0000000475,0.0000000434,0.0000000446,0.0000000437,0.0000000472,0.0000000465,0.0000000483,0.0000000489,0.0000000475,0.0000000458,0.0000000444,0.0000000496,0.0000000486,0.0000000455,0.0000000453,0.0000000496,0.0000000462,0.0000000239,0.0000000254,0.0000000242,0.0000000270,0.0000000293,0.0000000311,0.0000000308,0.0000000311,0.0000000342,0.0000000314,0.0000000319,0.0000000347,0.0000000346,0.0000000349,0.0000000403,0.0000000340,0.0000000342,0.0000000378,0.0000000383,0.0000000361,0.0000000419,0.0000000393,0.0000000365,0.0000000362,0.0000000400,0.0000000425,0.0000000449,0.0000000421,0.0000000411,0.0000000407,0.0000000382,0.0000000417,0.0000000459,0.0000000412,0.0000000425,0.0000000421,0.0000000443,0.0000000443,0.0000000452,0.0000000460,0.0000000453,0.0000000443,0.0000000419,0.0000000466,0.0000000456,0.0000000431,0.0000000426,0.0000000469,0.0000000443,0.0000000228,0.0000000237,0.0000000226,0.0000000260,0.0000000275,0.0000000301,0.0000000292,0.0000000302,0.0000000324,0.0000000300,0.0000000304,0.0000000332,0.0000000334,0.0000000332,0.0000000381,0.0000000323,0.0000000327,0.0000000362,0.0000000367,0.0000000344,0.0000000394,0.0000000372,0.0000000346,0.0000000346,0.0000000377,0.0000000407,0.0000000435,0.0000000400,0.0000000394,0.0000000387,0.0000000365,0.0000000389,0.0000000440,0.0000000388,0.0000000407,0.0000000395,0.0000000418,0.0000000422,0.0000000425,0.0000000432,0.0000000434,0.0000000427,0.0000000390,0.0000000443,0.0000000436,0.0000000402,0.0000000401,0.0000000446,0.0000000420,0.0000000214,0.0000000226,0.0000000214,0.0000000244,0.0000000258,0.0000000291,0.0000000278,0.0000000286,0.0000000305,0.0000000283,0.0000000294,0.0000000320,0.0000000320,0.0000000320,0.0000000365,0.0000000304,0.0000000307,0.0000000339,0.0000000342,0.0000000331,0.0000000372,0.0000000356,0.0000000330,0.0000000327,0.0000000360,0.0000000383,0.0000000413,0.0000000379,0.0000000372,0.0000000372,0.0000000352,0.0000000364,0.0000000417,0.0000000371,0.0000000383,0.0000000369,0.0000000396,0.0000000401,0.0000000407,0.0000000402,0.0000000416,0.0000000407,0.0000000363,0.0000000421,0.0000000412,0.0000000382,0.0000000378,0.0000000423,0.0000000400,0.0000000207,0.0000000213,0.0000000208,0.0000000233,0.0000000245,0.0000000282,0.0000000267,0.0000000273,0.0000000294,0.0000000265,0.0000000275,0.0000000300,0.0000000299,0.0000000308,0.0000000353,0.0000000288,0.0000000293,0.0000000318,0.0000000320,0.0000000311,0.0000000360,0.0000000337,0.0000000313,0.0000000313,0.0000000344,0.0000000364,0.0000000395,0.0000000352,0.0000000352,0.0000000355,0.0000000324,0.0000000351,0.0000000396,0.0000000354,0.0000000365,0.0000000356,0.0000000379,0.0000000382,0.0000000390,0.0000000387,0.0000000407,0.0000000395,0.0000000350,0.0000000400,0.0000000394,0.0000000362,0.0000000361,0.0000000405,0.0000000377,0.0000000201,0.0000000201,0.0000000202,0.0000000222,0.0000000231,0.0000000267,0.0000000252,0.0000000257,0.0000000276,0.0000000250,0.0000000257,0.0000000282,0.0000000280,0.0000000294,0.0000000333,0.0000000273,0.0000000272,0.0000000299,0.0000000304,0.0000000292,0.0000000344,0.0000000321,0.0000000291,0.0000000292,0.0000000330,0.0000000343,0.0000000378,0.0000000337,0.0000000339,0.0000000332,0.0000000305,0.0000000324,0.0000000369,0.0000000337,0.0000000343,0.0000000344,0.0000000356,0.0000000368,0.0000000373,0.0000000362,0.0000000382,0.0000000380,0.0000000333,0.0000000387,0.0000000376,0.0000000347,0.0000000345,0.0000000387,0.0000000355,0.0000000190,0.0000000187,0.0000000190,0.0000000206,0.0000000222,0.0000000260,0.0000000244,0.0000000242,0.0000000267,0.0000000235,0.0000000243,0.0000000271,0.0000000267,0.0000000273,0.0000000325,0.0000000257,0.0000000263,0.0000000288,0.0000000291,0.0000000276,0.0000000327,0.0000000300,0.0000000273,0.0000000281,0.0000000304,0.0000000323,0.0000000359,0.0000000328,0.0000000318,0.0000000313,0.0000000288,0.0000000309,0.0000000347,0.0000000321,0.0000000325,0.0000000332,0.0000000339,0.0000000350,0.0000000363,0.0000000341,0.0000000360,0.0000000359,0.0000000319,0.0000000371,0.0000000354,0.0000000320,0.0000000324,0.0000000359,0.0000000330,0.0000000183,0.0000000176,0.0000000180,0.0000000193,0.0000000216,0.0000000247,0.0000000234,0.0000000226,0.0000000253,0.0000000228,0.0000000222,0.0000000261,0.0000000248,0.0000000258,0.0000000309,0.0000000243,0.0000000247,0.0000000276,0.0000000269,0.0000000268,0.0000000311,0.0000000286,0.0000000257,0.0000000266,0.0000000292,0.0000000311,0.0000000344,0.0000000310,0.0000000298,0.0000000299,0.0000000273,0.0000000298,0.0000000333,0.0000000303,0.0000000313,0.0000000318,0.0000000323,0.0000000327,0.0000000351,0.0000000320,0.0000000343,0.0000000343,0.0000000304,0.0000000353,0.0000000342,0.0000000311,0.0000000317,0.0000000340,0.0000000310,0.0000000174,0.0000000167,0.0000000173,0.0000000188,0.0000000205,0.0000000238,0.0000000226,0.0000000217,0.0000000240,0.0000000210,0.0000000209,0.0000000238,0.0000000229,0.0000000238,0.0000000299,0.0000000229,0.0000000237,0.0000000257,0.0000000254,0.0000000256,0.0000000296,0.0000000273,0.0000000241,0.0000000256,0.0000000274,0.0000000294,0.0000000324,0.0000000296,0.0000000288,0.0000000286,0.0000000265,0.0000000284,0.0000000314,0.0000000288,0.0000000296,0.0000000302,0.0000000302,0.0000000313,0.0000000334,0.0000000315,0.0000000316,0.0000000328,0.0000000290,0.0000000332,0.0000000329,0.0000000298,0.0000000298,0.0000000317,0.0000000300,0.0000000168,0.0000000157,0.0000000164,0.0000000181,0.0000000197,0.0000000222,0.0000000215,0.0000000210,0.0000000227,0.0000000204,0.0000000201,0.0000000216,0.0000000220,0.0000000223,0.0000000281,0.0000000215,0.0000000222,0.0000000238,0.0000000239,0.0000000243,0.0000000287,0.0000000256,0.0000000235,0.0000000249,0.0000000262,0.0000000278,0.0000000307,0.0000000281,0.0000000278,0.0000000275,0.0000000250,0.0000000270,0.0000000300,0.0000000274,0.0000000285,0.0000000283,0.0000000287,0.0000000299,0.0000000323,0.0000000301,0.0000000299,0.0000000312,0.0000000274,0.0000000322,0.0000000303,0.0000000284,0.0000000288,0.0000000298,0.0000000287,0.0000000159,0.0000000152,0.0000000156,0.0000000173,0.0000000185,0.0000000209,0.0000000201,0.0000000203,0.0000000216,0.0000000196,0.0000000190,0.0000000206,0.0000000209,0.0000000208,0.0000000268,0.0000000206,0.0000000209,0.0000000223,0.0000000227,0.0000000229,0.0000000268,0.0000000246,0.0000000222,0.0000000242,0.0000000245,0.0000000267,0.0000000296,0.0000000271,0.0000000265,0.0000000260,0.0000000241,0.0000000250,0.0000000285,0.0000000256,0.0000000269,0.0000000269,0.0000000277,0.0000000287,0.0000000310,0.0000000292,0.0000000288,0.0000000294,0.0000000261,0.0000000302,0.0000000279,0.0000000272,0.0000000270,0.0000000283,0.0000000269,0.0000000152,0.0000000145,0.0000000147,0.0000000167,0.0000000171,0.0000000204,0.0000000184,0.0000000187,0.0000000210,0.0000000186,0.0000000180,0.0000000195,0.0000000201,0.0000000195,0.0000000252,0.0000000192,0.0000000201,0.0000000213,0.0000000219,0.0000000219,0.0000000256,0.0000000235,0.0000000209,0.0000000225,0.0000000234,0.0000000255,0.0000000290,0.0000000262,0.0000000250,0.0000000249,0.0000000227,0.0000000234,0.0000000277,0.0000000245,0.0000000254,0.0000000257,0.0000000272,0.0000000268,0.0000000296,0.0000000271,0.0000000273,0.0000000289,0.0000000251,0.0000000283,0.0000000266,0.0000000255,0.0000000258,0.0000000270,0.0000000257,0.0000000143,0.0000000139,0.0000000141,0.0000000161,0.0000000159,0.0000000193,0.0000000177,0.0000000174,0.0000000204,0.0000000174,0.0000000172,0.0000000178,0.0000000187,0.0000000185,0.0000000240,0.0000000185,0.0000000187,0.0000000198,0.0000000210,0.0000000208,0.0000000246,0.0000000223,0.0000000203,0.0000000213,0.0000000223,0.0000000243,0.0000000279,0.0000000251,0.0000000238,0.0000000238,0.0000000213,0.0000000221,0.0000000257,0.0000000232,0.0000000245,0.0000000249,0.0000000258,0.0000000258,0.0000000280,0.0000000254,0.0000000260,0.0000000266,0.0000000237,0.0000000267,0.0000000256,0.0000000243,0.0000000240,0.0000000257,0.0000000248,0.0000000130,0.0000000130,0.0000000134,0.0000000154,0.0000000148,0.0000000180,0.0000000167,0.0000000161,0.0000000194,0.0000000171,0.0000000162,0.0000000167,0.0000000179,0.0000000176,0.0000000228,0.0000000176,0.0000000174,0.0000000193,0.0000000205,0.0000000192,0.0000000238,0.0000000209,0.0000000197,0.0000000198,0.0000000212,0.0000000233,0.0000000261,0.0000000239,0.0000000225,0.0000000226,0.0000000196,0.0000000206,0.0000000241,0.0000000225,0.0000000232,0.0000000239,0.0000000245,0.0000000241,0.0000000271,0.0000000246,0.0000000245,0.0000000253,0.0000000224,0.0000000253,0.0000000247,0.0000000229,0.0000000230,0.0000000241,0.0000000235};
vec Han_data2={
0.9999880344,0.9999881993,0.9999883509,0.9999883960,0.9999885269,0.9999885428,0.9999886450,0.9999886693,0.9999887157,0.9999886761,0.9999887518,0.9999887831,0.9999888026,0.9999887599,0.9999887482,0.9999886869,0.9999887357,0.9999886871,0.9999886885,0.9999886386,0.9999886042,0.9999885726,0.9999886180,0.9999884695,0.9999884548,0.9999883337,0.9999883376,0.9999882864,0.9999881766,0.9999881352,0.9999880846,0.9999880095,0.9999879152,0.9999878824,0.9999878174,0.9999877175,0.9999876696,0.9999875615,0.9999874977,0.9999873760,0.9999872960,0.9999872130,0.9999871510,0.9999870507,0.9999869279,0.9999868486,0.9999867509,0.9999867365,0.9999864644,0.9999863838,0.9999862460,0.9999861992,0.9999860146,0.9999858661,0.9999857658,0.9999855408,0.9999853663,0.9999852238,0.9999850967,0.9999849437,0.9999847259,0.9999845626,0.9999843817,0.9999841309,0.9999839802,0.9999837699,0.9999835044,0.9999832558,0.9999829723,0.9999827425,0.9999824850,0.9999821300,0.9999818935,0.9999815265,0.9999812677,0.9999809225,0.9999804318,0.9999799707,0.9999796017,0.9999791286,0.9999786900,0.9999780911,0.9999775082,0.9999768673,0.9999761371,0.9999754330,0.9999744151,0.9999734754,0.9999723701,0.9999711967,0.9999698092,0.9999682455,0.9999662364,0.9999641325,0.9999611242,0.9999573881,0.9999524576,0.9999450682,0.9999327917,0.9999051338,0.8177377147,0.8192624447,0.8205953929,0.8217746482,0.8227975198,0.8237245972,0.8245363914,0.8252849514,0.8259551621,0.8265704013,0.8271234096,0.8276272316,0.8280983365,0.8285294622,0.8289348098,0.8293149794,0.8296534608,0.8299865316,0.8302863158,0.8305641897,0.8308427519,0.8310921070,0.8313270901,0.8315521929,0.8317727642,0.8319563736,0.8321508368,0.8323367291,0.8325088263,0.8326682873,0.8328334069,0.8329741620,0.8331163155,0.8332488761,0.8333855137,0.8335089902,0.8336309884,0.8337488065,0.8338516168,0.8339655768,0.8340650229,0.8341673948,0.8342584673,0.8343556089,0.8344379445,0.8345282585,0.8346122116,0.8346929978,0.8347646197,0.8348396648,0.8349123202,0.8349820510,0.8350531601,0.8351251091,0.8351820456,0.8352507437,0.8353094439,0.8353640780,0.8354120601,0.8354776345,0.8355277037,0.8355804221,0.8356326315,0.8356853782,0.8357316293,0.8357789539,0.8358214633,0.8358687043,0.8359129972,0.8359569171,0.8359898732,0.8360416685,0.8360802121,0.8361098742,0.8361563506,0.8361886749,0.8362320947,0.8362568334,0.8362915737,0.8363287513,0.8363628015,0.8363877522,0.8364367263,0.8364556574,0.8364846710,0.8365121531,0.8365446323,0.8365796768,0.8365984171,0.8366331498,0.8366559576,0.8366793749,0.8367138837,0.8367392026,0.8367640294,0.8367863307,0.8368125337,0.8368297823,0.8368560792,0.8368746703,0.7354003433,0.7373313974,0.7390202043,0.7405267054,0.7418310423,0.7430208964,0.7440665924,0.7450271059,0.7458917565,0.7466886997,0.7474119524,0.7480676393,0.7486793514,0.7492476032,0.7497709451,0.7502628929,0.7507121769,0.7511446071,0.7515416496,0.7519098155,0.7522703632,0.7526023512,0.7529087697,0.7532074503,0.7534951954,0.7537405948,0.7539976156,0.7542364402,0.7544689897,0.7546807417,0.7548997243,0.7550857503,0.7552773115,0.7554551069,0.7556348902,0.7557983777,0.7559646520,0.7561180662,0.7562549703,0.7564072371,0.7565357358,0.7566715221,0.7567924330,0.7569247903,0.7570348454,0.7571531593,0.7572631466,0.7573676910,0.7574662786,0.7575700408,0.7576648634,0.7577538581,0.7578532630,0.7579470222,0.7580257697,0.7581166106,0.7581923411,0.7582689085,0.7583346827,0.7584220461,0.7584855000,0.7585587203,0.7586343564,0.7586916471,0.7587610249,0.7588221202,0.7588810849,0.7589397443,0.7590075815,0.7590617062,0.7591112345,0.7591719198,0.7592274919,0.7592684911,0.7593341665,0.7593760967,0.7594248574,0.7594712694,0.7595165729,0.7595615801,0.7596032096,0.7596440179,0.7597029152,0.7597321691,0.7597717699,0.7598111844,0.7598532263,0.7599012776,0.7599244686,0.7599658297,0.7599976943,0.7600381836,0.7600775717,0.7601081093,0.7601467932,0.7601736929,0.7602088346,0.7602355537,0.7602683912,0.7602945059,0.6715222789,0.6736721904,0.6755523036,0.6772358611,0.6787028174,0.6800386238,0.6812156271,0.6822994643,0.6832758477,0.6841757968,0.6849925125,0.6857380931,0.6864360388,0.6870840833,0.6876779505,0.6882373680,0.6887489672,0.6892408898,0.6896969567,0.6901219357,0.6905333120,0.6909121037,0.6912598350,0.6916043885,0.6919322167,0.6922230861,0.6925114651,0.6927879551,0.6930566973,0.6932993627,0.6935513111,0.6937665402,0.6939882076,0.6941927670,0.6943941507,0.6945852188,0.6947750406,0.6949537211,0.6951155361,0.6952887381,0.6954389737,0.6955925749,0.6957367639,0.6958843423,0.6960164093,0.6961507593,0.6962790130,0.6963962105,0.6965179804,0.6966377000,0.6967475986,0.6968485502,0.6969644732,0.6970702375,0.6971564965,0.6972617187,0.6973548128,0.6974418306,0.6975262585,0.6976227655,0.6976984114,0.6977770228,0.6978675064,0.6979340466,0.6980148445,0.6980870497,0.6981560068,0.6982286354,0.6983045511,0.6983653522,0.6984247072,0.6984889820,0.6985533538,0.6986062551,0.6986772952,0.6987279747,0.6987837811,0.6988381046,0.6988894963,0.6989445159,0.6989945529,0.6990437085,0.6991075097,0.6991442761,0.6991927293,0.6992396606,0.6992851458,0.6993355858,0.6993712096,0.6994157765,0.6994530169,0.6995035983,0.6995467077,0.6995818517,0.6996222923,0.6996548838,0.6997003625,0.6997266587,0.6997681641,0.6998002187,0.6179779036,0.6202531521,0.6222450820,0.6240273507,0.6255885932,0.6270058137,0.6282684312,0.6294265707,0.6304703218,0.6314382598,0.6323157384,0.6331129885,0.6338643497,0.6345617682,0.6352027938,0.6358065653,0.6363587455,0.6368917055,0.6373792860,0.6378382633,0.6382838443,0.6386980842,0.6390702789,0.6394453856,0.6397978597,0.6401240412,0.6404299024,0.6407282047,0.6410251380,0.6412834059,0.6415553346,0.6417945272,0.6420347851,0.6422593713,0.6424808763,0.6426839838,0.6428876981,0.6430832520,0.6432635070,0.6434499429,0.6436146849,0.6437857490,0.6439375925,0.6440970036,0.6442455786,0.6443899334,0.6445289021,0.6446551394,0.6447904476,0.6449207784,0.6450401010,0.6451520282,0.6452766042,0.6453945356,0.6454906686,0.6456008444,0.6457075073,0.6458003899,0.6458942284,0.6459988212,0.6460798613,0.6461712360,0.6462678675,0.6463424584,0.6464290129,0.6465106092,0.6465848016,0.6466589383,0.6467465778,0.6468110692,0.6468811118,0.6469467996,0.6470180109,0.6470775176,0.6471503314,0.6472081557,0.6472716088,0.6473345901,0.6473897338,0.6474519221,0.6475060380,0.6475591348,0.6476256736,0.6476642083,0.6477191626,0.6477739946,0.6478203920,0.6478741844,0.6479172348,0.6479680448,0.6480058402,0.6480627998,0.6481068543,0.6481490158,0.6481927784,0.6482263148,0.6482770000,0.6483075279,0.6483541048,0.6483845568,0.5714930598,0.5738306619,0.5758860513,0.5777265550,0.5793415982,0.5808072096,0.5821142305,0.5833223616,0.5844032988,0.5854114039,0.5863265319,0.5871582226,0.5879445017,0.5886691250,0.5893445762,0.5899698787,0.5905503999,0.5911107257,0.5916227452,0.5920989515,0.5925703431,0.5930010087,0.5933941131,0.5937890877,0.5941617920,0.5945042202,0.5948248755,0.5951370537,0.5954499968,0.5957235401,0.5960095954,0.5962650827,0.5965158565,0.5967531113,0.5969861412,0.5972002573,0.5974170236,0.5976222907,0.5978116664,0.5980107602,0.5981852723,0.5983674206,0.5985244327,0.5986933750,0.5988482322,0.5990003371,0.5991460397,0.5992839291,0.5994287211,0.5995674782,0.5996903673,0.5998131260,0.5999373125,0.6000662998,0.6001666237,0.6002863304,0.6004006045,0.6004954399,0.6005955343,0.6007100670,0.6007918942,0.6008914082,0.6009884317,0.6010723746,0.6011617027,0.6012461892,0.6013256324,0.6014027875,0.6014982728,0.6015678196,0.6016430680,0.6017128323,0.6017853988,0.6018529339,0.6019310151,0.6019878921,0.6020538053,0.6021207173,0.6021793880,0.6022503394,0.6023029611,0.6023668701,0.6024291885,0.6024759868,0.6025298504,0.6025902444,0.6026399199,0.6026956956,0.6027431559,0.6027925510,0.6028354028,0.6028995250,0.6029445506,0.6029879499,0.6030346450,0.6030711483,0.6031237389,0.6031597774,0.6032036730,0.6032413940,0.5303248954,0.5326868278,0.5347662398,0.5366342245,0.5382776984,0.5397674355,0.5410985833,0.5423301041,0.5434352853,0.5444619647,0.5453995221,0.5462531328,0.5470579282,0.5477984715,0.5484925602,0.5491376110,0.5497297873,0.5503046144,0.5508335896,0.5513242005,0.5518081894,0.5522516635,0.5526564214,0.5530652623,0.5534492269,0.5538022974,0.5541339635,0.5544527249,0.5547803001,0.5550621852,0.5553578864,0.5556224059,0.5558810671,0.5561235347,0.5563648991,0.5565918799,0.5568112275,0.5570284389,0.5572215756,0.5574262422,0.5576087280,0.5577951512,0.5579592824,0.5581339458,0.5582928782,0.5584523448,0.5586036715,0.5587466213,0.5588976940,0.5590373841,0.5591683012,0.5592974055,0.5594218432,0.5595568652,0.5596628263,0.5597871538,0.5599038117,0.5600023267,0.5601036004,0.5602262119,0.5603105398,0.5604141968,0.5605153534,0.5606029919,0.5606969749,0.5607857144,0.5608672369,0.5609470377,0.5610448869,0.5611186785,0.5612000849,0.5612673339,0.5613446978,0.5614172885,0.5614946626,0.5615562601,0.5616260562,0.5616931384,0.5617555080,0.5618248552,0.5618835155,0.5619497834,0.5620136531,0.5620623357,0.5621211546,0.5621775901,0.5622330445,0.5622947233,0.5623417705,0.5623913074,0.5624403220,0.5625045472,0.5625481487,0.5625958094,0.5626460698,0.5626834404,0.5627364225,0.5627761407,0.5628207947,0.5628590257,0.4933856571,0.4957504522,0.4978289995,0.4997027698,0.5013520076,0.5028482551,0.5041895713,0.5054289774,0.5065431768,0.5075770834,0.5085242338,0.5093876946,0.5102033372,0.5109552690,0.5116551468,0.5123110856,0.5129111148,0.5134922744,0.5140282464,0.5145275408,0.5150201559,0.5154702148,0.5158830599,0.5162973881,0.5166874195,0.5170482806,0.5173870422,0.5177106859,0.5180454880,0.5183320683,0.5186352895,0.5189039039,0.5191667323,0.5194109637,0.5196621156,0.5198917547,0.5201158962,0.5203387598,0.5205320225,0.5207435354,0.5209279651,0.5211213241,0.5212915350,0.5214718353,0.5216309016,0.5217933467,0.5219501836,0.5220958601,0.5222501382,0.5223916515,0.5225244092,0.5226571770,0.5227887070,0.5229234057,0.5230317003,0.5231618385,0.5232782037,0.5233777290,0.5234831228,0.5236094996,0.5236962899,0.5237994366,0.5239052897,0.5239933926,0.5240930963,0.5241846822,0.5242651826,0.5243514309,0.5244462245,0.5245254067,0.5246069930,0.5246761904,0.5247537811,0.5248288980,0.5249073245,0.5249715810,0.5250468638,0.5251110273,0.5251798160,0.5252501504,0.5253073068,0.5253756028,0.5254414000,0.5254899923,0.5255519108,0.5256125633,0.5256634720,0.5257298841,0.5257739099,0.5258305477,0.5258829988,0.5259481590,0.5259905337,0.5260411502,0.5260920389,0.5261274405,0.5261863826,0.5262242146,0.5262665807,0.5263131159,0.4599413600,0.4622872035,0.4643513789,0.4662141647,0.4678550027,0.4693488639,0.4706846311,0.4719212491,0.4730348649,0.4740662384,0.4750159809,0.4758821693,0.4767006681,0.4774503132,0.4781550851,0.4788126318,0.4794189296,0.4800019023,0.4805402851,0.4810433079,0.4815327687,0.4819919100,0.4824056922,0.4828278305,0.4832157849,0.4835785595,0.4839241372,0.4842489895,0.4845852452,0.4848761702,0.4851818311,0.4854508324,0.4857189913,0.4859652780,0.4862215294,0.4864538684,0.4866811731,0.4869015869,0.4870993396,0.4873125479,0.4875012177,0.4876939738,0.4878660787,0.4880550615,0.4882098104,0.4883739192,0.4885331882,0.4886822339,0.4888365264,0.4889805288,0.4891115732,0.4892512048,0.4893817167,0.4895192777,0.4896286151,0.4897619461,0.4898753226,0.4899793085,0.4900876657,0.4902131072,0.4903033427,0.4904071372,0.4905167518,0.4906029688,0.4907061916,0.4907971839,0.4908807655,0.4909691087,0.4910636300,0.4911414797,0.4912285952,0.4912987094,0.4913764283,0.4914556778,0.4915314838,0.4915982050,0.4916757743,0.4917426586,0.4918062229,0.4918823479,0.4919362443,0.4920085003,0.4920759941,0.4921228167,0.4921888709,0.4922499753,0.4923005585,0.4923648034,0.4924083690,0.4924689354,0.4925215611,0.4925896508,0.4926319072,0.4926827211,0.4927313966,0.4927721282,0.4928319278,0.4928724598,0.4929172973,0.4929577966,0.4294626959,0.4317754394,0.4338072541,0.4356486627,0.4372722658,0.4387488019,0.4400749114,0.4413022574,0.4424069178,0.4434325700,0.4443756686,0.4452382393,0.4460505140,0.4467968287,0.4474993627,0.4481559654,0.4487592690,0.4493414893,0.4498760908,0.4503794336,0.4508700321,0.4513309405,0.4517431236,0.4521625328,0.4525521768,0.4529163761,0.4532636343,0.4535887271,0.4539221537,0.4542167553,0.4545234768,0.4547893967,0.4550608633,0.4553069579,0.4555625384,0.4557972649,0.4560266137,0.4562460041,0.4564485260,0.4566602676,0.4568505987,0.4570416225,0.4572168105,0.4574039741,0.4575579818,0.4577241295,0.4578864293,0.4580363323,0.4581888978,0.4583353034,0.4584662651,0.4586066261,0.4587379710,0.4588772790,0.4589851128,0.4591198027,0.4592376567,0.4593412223,0.4594472975,0.4595720766,0.4596665495,0.4597739599,0.4598769784,0.4599676535,0.4600713593,0.4601607676,0.4602484982,0.4603353067,0.4604322191,0.4605129218,0.4605973189,0.4606663158,0.4607465634,0.4608241671,0.4609001092,0.4609684909,0.4610452801,0.4611167857,0.4611837741,0.4612537546,0.4613088499,0.4613814538,0.4614503040,0.4614960880,0.4615648048,0.4616286823,0.4616780097,0.4617440515,0.4617851890,0.4618456671,0.4618994385,0.4619701846,0.4620097335,0.4620637787,0.4621148877,0.4621553863,0.4622129661,0.4622537026,0.4622987462,0.4623401189,0.4015432140,0.4038100043,0.4058104542,0.4076147457,0.4092113324,0.4106682492,0.4119743405,0.4131863635,0.4142749771,0.4152897731,0.4162217975,0.4170763717,0.4178767063,0.4186182164,0.4193132452,0.4199636257,0.4205647267,0.4211386938,0.4216713701,0.4221738559,0.4226607192,0.4231146130,0.4235271857,0.4239408361,0.4243292151,0.4246914290,0.4250388217,0.4253645252,0.4256955828,0.4259895375,0.4262957264,0.4265582976,0.4268297803,0.4270749084,0.4273300350,0.4275640796,0.4277911447,0.4280108168,0.4282119013,0.4284240423,0.4286146714,0.4288045042,0.4289828318,0.4291653378,0.4293185069,0.4294862478,0.4296479883,0.4297980063,0.4299514302,0.4300979879,0.4302286376,0.4303706296,0.4305003753,0.4306358652,0.4307486758,0.4308828294,0.4309988155,0.4311054126,0.4312091262,0.4313345303,0.4314262070,0.4315342326,0.4316390283,0.4317294641,0.4318304606,0.4319195884,0.4320094642,0.4320929736,0.4321920902,0.4322766842,0.4323608980,0.4324296095,0.4325101012,0.4325852825,0.4326618635,0.4327320315,0.4328085564,0.4328796970,0.4329456971,0.4330197116,0.4330750595,0.4331434682,0.4332121457,0.4332587351,0.4333292932,0.4333893167,0.4334414708,0.4335099745,0.4335496290,0.4336132133,0.4336649168,0.4337317463,0.4337725226,0.4338273644,0.4338761369,0.4339220184,0.4339764273,0.4340196147,0.4340653022,0.4341043918,0.3758683584,0.3780811825,0.3800382999,0.3818011281,0.3833653321,0.3847936408,0.3860805241,0.3872696969,0.3883387265,0.3893367741,0.3902536439,0.3910983752,0.3918855074,0.3926161603,0.3933028067,0.3939444372,0.3945382127,0.3951035680,0.3956324589,0.3961298912,0.3966064915,0.3970587910,0.3974659379,0.3978747701,0.3982587563,0.3986157420,0.3989625610,0.3992864456,0.3996129231,0.3999060013,0.4002109472,0.4004694730,0.4007368920,0.4009796531,0.4012341052,0.4014703931,0.4016919008,0.4019092134,0.4021097816,0.4023183445,0.4025068501,0.4027011111,0.4028752901,0.4030578032,0.4032084596,0.4033734853,0.4035381168,0.4036873269,0.4038385041,0.4039845177,0.4041126385,0.4042557445,0.4043841750,0.4045190414,0.4046314536,0.4047613279,0.4048816623,0.4049837713,0.4050903749,0.4052146337,0.4053065331,0.4054138787,0.4055150769,0.4056079800,0.4057085678,0.4057986469,0.4058876290,0.4059694275,0.4060693641,0.4061543808,0.4062352168,0.4063037175,0.4063896113,0.4064601612,0.4065323760,0.4066060379,0.4066808456,0.4067546299,0.4068185887,0.4068915494,0.4069462249,0.4070156376,0.4070840718,0.4071314924,0.4072000070,0.4072600977,0.4073138403,0.4073839478,0.4074198478,0.4074843114,0.4075346342,0.4076034753,0.4076423023,0.4076973086,0.4077460670,0.4077892519,0.4078455002,0.4078867120,0.4079374966,0.4079763753,0.3521761145,0.3543278958,0.3562387335,0.3579609347,0.3594902762,0.3608826995,0.3621414929,0.3633040185,0.3643554250,0.3653317752,0.3662300280,0.3670591555,0.3678303343,0.3685488112,0.3692288551,0.3698529329,0.3704415113,0.3709971671,0.3715191318,0.3720084976,0.3724773653,0.3729223398,0.3733265890,0.3737267375,0.3741069135,0.3744574530,0.3748001495,0.3751195094,0.3754431028,0.3757333158,0.3760332569,0.3762873226,0.3765526920,0.3767933445,0.3770448401,0.3772750195,0.3774951247,0.3777098483,0.3779076256,0.3781153283,0.3783019261,0.3784943477,0.3786656480,0.3788448509,0.3789963543,0.3791571285,0.3793199540,0.3794698955,0.3796184157,0.3797646484,0.3798934273,0.3800350283,0.3801599429,0.3802950698,0.3804069397,0.3805346387,0.3806516527,0.3807555797,0.3808604099,0.3809827055,0.3810736841,0.3811788726,0.3812814519,0.3813723450,0.3814711273,0.3815617074,0.3816504679,0.3817306843,0.3818320869,0.3819161322,0.3819928100,0.3820646326,0.3821491792,0.3822175225,0.3822885197,0.3823651893,0.3824385302,0.3825113278,0.3825746702,0.3826445686,0.3827028330,0.3827709687,0.3828405899,0.3828871194,0.3829547491,0.3830137561,0.3830669241,0.3831334118,0.3831728952,0.3832367487,0.3832866181,0.3833527910,0.3833906482,0.3834497216,0.3834939595,0.3835379135,0.3835932913,0.3836345250,0.3836852368,0.3837239184,0.3302568194,0.3323468237,0.3342048826,0.3358791988,0.3373697945,0.3387291250,0.3399581293,0.3410907287,0.3421151850,0.3430712778,0.3439505148,0.3447629839,0.3455166498,0.3462215821,0.3468866448,0.3475010053,0.3480747766,0.3486198825,0.3491328681,0.3496123723,0.3500735937,0.3505140028,0.3509106557,0.3513006048,0.3516733286,0.3520215076,0.3523557785,0.3526711194,0.3529901258,0.3532738696,0.3535673158,0.3538184633,0.3540815673,0.3543190932,0.3545669872,0.3547944289,0.3550081073,0.3552235345,0.3554178322,0.3556220781,0.3558049825,0.3559949722,0.3561643848,0.3563415162,0.3564913874,0.3566489032,0.3568113871,0.3569566219,0.3571021643,0.3572479903,0.3573745332,0.3575142384,0.3576374532,0.3577710485,0.3578831407,0.3580080549,0.3581229256,0.3582255672,0.3583304947,0.3584508453,0.3585414424,0.3586442206,0.3587478648,0.3588366768,0.3589322176,0.3590232135,0.3591124423,0.3591894030,0.3592905501,0.3593746710,0.3594497278,0.3595187877,0.3596070417,0.3596707692,0.3597403674,0.3598196498,0.3598894074,0.3599637904,0.3600265281,0.3600975952,0.3601557101,0.3602198581,0.3602863210,0.3603326576,0.3603994300,0.3604569526,0.3605106222,0.3605802894,0.3606163052,0.3606777695,0.3607307894,0.3607933713,0.3608339045,0.3608881934,0.3609356031,0.3609806493,0.3610327891,0.3610732850,0.3611236458,0.3611630380,0.3099308897,0.3119579316,0.3137632550,0.3153841121,0.3168328068,0.3181571170,0.3193522008,0.3204564503,0.3214527497,0.3223887251,0.3232447609,0.3240384067,0.3247738622,0.3254626763,0.3261154928,0.3267130986,0.3272736535,0.3278066201,0.3283097306,0.3287800742,0.3292319933,0.3296628582,0.3300479641,0.3304316873,0.3307969371,0.3311397715,0.3314703054,0.3317773318,0.3320910692,0.3323677415,0.3326544171,0.3329033817,0.3331603693,0.3333945999,0.3336414355,0.3338618247,0.3340719959,0.3342794907,0.3344755554,0.3346740133,0.3348546257,0.3350398627,0.3352073981,0.3353827997,0.3355324100,0.3356849062,0.3358440699,0.3359886545,0.3361302153,0.3362736676,0.3363998317,0.3365353306,0.3366606633,0.3367903113,0.3368996138,0.3370207085,0.3371379050,0.3372372261,0.3373434927,0.3374583976,0.3375482180,0.3376482831,0.3377487432,0.3378373516,0.3379326098,0.3380226653,0.3381119645,0.3381877622,0.3382833735,0.3383699357,0.3384433372,0.3385127499,0.3385971665,0.3386618612,0.3387315630,0.3388060716,0.3388755506,0.3389501069,0.3390097511,0.3390827852,0.3391389031,0.3391996909,0.3392666182,0.3393150998,0.3393816055,0.3394357791,0.3394872542,0.3395566003,0.3395963291,0.3396558384,0.3397087194,0.3397675209,0.3398087983,0.3398611954,0.3399120001,0.3399534746,0.3400034758,0.3400440221,0.3400932987,0.3401326347,0.2910470372,0.2930090313,0.2947576086,0.2963303906,0.2977325388,0.2990198108,0.3001794138,0.3012550079,0.3022223777,0.3031347896,0.3039657860,0.3047393489,0.3054588439,0.3061276310,0.3067640392,0.3073464281,0.3078937637,0.3084114668,0.3089058204,0.3093641254,0.3098066414,0.3102274481,0.3106027555,0.3109782717,0.3113363633,0.3116721354,0.3119962102,0.3122977571,0.3126013982,0.3128723147,0.3131531488,0.3133987026,0.3136499966,0.3138803081,0.3141220887,0.3143378722,0.3145414170,0.3147483145,0.3149400507,0.3151334941,0.3153114804,0.3154937994,0.3156578176,0.3158298771,0.3159759020,0.3161236935,0.3162827757,0.3164246933,0.3165649032,0.3167039824,0.3168261464,0.3169604461,0.3170841008,0.3172068641,0.3173177287,0.3174379174,0.3175514105,0.3176530889,0.3177551341,0.3178680505,0.3179572330,0.3180560951,0.3181533015,0.3182386663,0.3183321534,0.3184239617,0.3185119547,0.3185843130,0.3186791546,0.3187624657,0.3188362929,0.3189037078,0.3189837676,0.3190511993,0.3191227656,0.3191900722,0.3192634998,0.3193358845,0.3193913886,0.3194668511,0.3195183353,0.3195805347,0.3196440994,0.3196930298,0.3197593466,0.3198095860,0.3198643645,0.3199307749,0.3199683712,0.3200274606,0.3200798528,0.3201383212,0.3201803323,0.3202289446,0.3202818952,0.3203224936,0.3203720818,0.3204125898,0.3204597894,0.3204991776,0.2734745505,0.2753698310,0.2770614966,0.2785833643,0.2799404444,0.2811878650,0.2823137091,0.2833562046,0.2842971983,0.2851807186,0.2859888077,0.2867425143,0.2874401739,0.2880933114,0.2887121578,0.2892799807,0.2898109019,0.2903152882,0.2907981961,0.2912444858,0.2916786762,0.2920869109,0.2924533864,0.2928198294,0.2931672405,0.2934954490,0.2938125092,0.2941072971,0.2944033644,0.2946677751,0.2949438581,0.2951853853,0.2954303815,0.2956556154,0.2958921222,0.2961014632,0.2963015277,0.2965019880,0.2966922698,0.2968793803,0.2970533022,0.2972311692,0.2973929196,0.2975608562,0.2977043424,0.2978488839,0.2980041345,0.2981444841,0.2982826431,0.2984182381,0.2985376397,0.2986693902,0.2987902427,0.2989116503,0.2990207313,0.2991364701,0.2992483748,0.2993486711,0.2994491398,0.2995586069,0.2996454638,0.2997446772,0.2998399509,0.2999233528,0.3000133968,0.3001057873,0.3001900249,0.3002626817,0.3003538793,0.3004354791,0.3005071125,0.3005743581,0.3006519772,0.3007214846,0.3007900986,0.3008541074,0.3009288527,0.3009977756,0.3010525436,0.3011253592,0.3011755312,0.3012393390,0.3013009664,0.3013498292,0.3014127055,0.3014630979,0.3015164432,0.3015835179,0.3016205009,0.3016779930,0.3017284693,0.3017858711,0.3018290385,0.3018757347,0.3019263242,0.3019661940,0.3020174564,0.3020578982,0.3021021983,0.3021425766,0.2571029591,0.2589295290,0.2605634064,0.2620333978,0.2633439950,0.2645546042,0.2656436690,0.2666518821,0.2675639348,0.2684232347,0.2692040384,0.2699368163,0.2706132972,0.2712471927,0.2718512109,0.2724016695,0.2729172087,0.2734079523,0.2738787507,0.2743102203,0.2747357178,0.2751317039,0.2754892245,0.2758467158,0.2761871037,0.2765039739,0.2768135024,0.2771020258,0.2773909564,0.2776466021,0.2779179455,0.2781541714,0.2783941287,0.2786133279,0.2788410080,0.2790491251,0.2792434682,0.2794382888,0.2796228749,0.2798063725,0.2799745452,0.2801502572,0.2803095045,0.2804717813,0.2806120077,0.2807553464,0.2809053984,0.2810452638,0.2811766300,0.2813126214,0.2814278921,0.2815563142,0.2816738913,0.2817960249,0.2819002013,0.2820117229,0.2821230924,0.2822193623,0.2823190077,0.2824239632,0.2825124488,0.2826080030,0.2827026201,0.2827828546,0.2828732002,0.2829607925,0.2830425889,0.2831144198,0.2832034718,0.2832863070,0.2833545921,0.2834217229,0.2834952898,0.2835646607,0.2836316078,0.2836957442,0.2837696997,0.2838337374,0.2838894958,0.2839592046,0.2840102000,0.2840722904,0.2841309050,0.2841820986,0.2842424398,0.2842937500,0.2843454852,0.2844111421,0.2844447998,0.2845028834,0.2845528437,0.2846061401,0.2846474298,0.2846943614,0.2847470571,0.2847832570,0.2848354666,0.2848755245,0.2849185059,0.2849601489,0.2418244417,0.2435852810,0.2451632005,0.2465808414,0.2478478159,0.2490138659,0.2500700482,0.2510443089,0.2519257197,0.2527591329,0.2535157644,0.2542247287,0.2548799072,0.2554949897,0.2560831065,0.2566169472,0.2571181831,0.2575941644,0.2580516767,0.2584715511,0.2588836219,0.2592670798,0.2596178812,0.2599650812,0.2602946145,0.2606050758,0.2609020272,0.2611865271,0.2614657644,0.2617142187,0.2619799313,0.2622091232,0.2624451189,0.2626572407,0.2628786848,0.2630838724,0.2632717364,0.2634606344,0.2636395838,0.2638202484,0.2639854441,0.2641549702,0.2643104026,0.2644688623,0.2646071859,0.2647488045,0.2648945392,0.2650297596,0.2651571024,0.2652893892,0.2654022198,0.2655285995,0.2656440874,0.2657620211,0.2658634915,0.2659730534,0.2660828079,0.2661776800,0.2662726082,0.2663784952,0.2664642936,0.2665557795,0.2666464492,0.2667279938,0.2668137419,0.2669006305,0.2669787625,0.2670525243,0.2671366380,0.2672193321,0.2672841312,0.2673529157,0.2674223445,0.2674931155,0.2675540715,0.2676169323,0.2676926147,0.2677520930,0.2678082241,0.2678773635,0.2679269234,0.2679869430,0.2680446523,0.2680960449,0.2681530659,0.2682067654,0.2682559218,0.2683170376,0.2683526146,0.2684065204,0.2684575814,0.2685100059,0.2685499693,0.2685978161,0.2686465302,0.2686838083,0.2687338447,0.2687748312,0.2688149961,0.2688580013,0.2275563826,0.2292507632,0.2307705546,0.2321362288,0.2333585861,0.2344816797,0.2355041092,0.2364451109,0.2372956609,0.2381053011,0.2388373934,0.2395227955,0.2401558982,0.2407519790,0.2413221247,0.2418387064,0.2423258706,0.2427875605,0.2432304297,0.2436394523,0.2440374438,0.2444097877,0.2447512926,0.2450867205,0.2454081980,0.2457095253,0.2459985754,0.2462755089,0.2465461048,0.2467882666,0.2470469955,0.2472685781,0.2474995097,0.2477045280,0.2479194315,0.2481220531,0.2483054519,0.2484882184,0.2486611731,0.2488384143,0.2489988934,0.2491651929,0.2493171860,0.2494690291,0.2496044805,0.2497413999,0.2498837979,0.2500154396,0.2501425361,0.2502677473,0.2503801447,0.2505016054,0.2506130397,0.2507307454,0.2508295486,0.2509379996,0.2510424878,0.2511352188,0.2512280414,0.2513309224,0.2514139655,0.2515036371,0.2515910533,0.2516717076,0.2517545739,0.2518397506,0.2519160984,0.2519886644,0.2520710595,0.2521499903,0.2522133424,0.2522810993,0.2523501884,0.2524186646,0.2524798400,0.2525377437,0.2526150857,0.2526706091,0.2527270471,0.2527934501,0.2528424582,0.2529032825,0.2529577946,0.2530072874,0.2530635007,0.2531161935,0.2531639373,0.2532220195,0.2532580607,0.2533124143,0.2533617897,0.2534130605,0.2534516412,0.2534982261,0.2535459318,0.2535813115,0.2536311984,0.2536705479,0.2537100439,0.2537515676,0.2142145485,0.2158450609,0.2173062491,0.2186203083,0.2198003597,0.2208825959,0.2218676525,0.2227775238,0.2235962301,0.2243802931,0.2250875826,0.2257517148,0.2263633629,0.2269402441,0.2274915554,0.2279915525,0.2284630290,0.2289108586,0.2293384169,0.2297371032,0.2301229543,0.2304839224,0.2308124513,0.2311423212,0.2314511394,0.2317439136,0.2320247315,0.2322950093,0.2325563265,0.2327906230,0.2330433769,0.2332594383,0.2334829677,0.2336848260,0.2338896008,0.2340869505,0.2342666489,0.2344438805,0.2346129203,0.2347810294,0.2349402611,0.2351014596,0.2352495911,0.2353952673,0.2355277206,0.2356626241,0.2358011243,0.2359292751,0.2360517731,0.2361747091,0.2362845929,0.2364028705,0.2365105286,0.2366270430,0.2367227095,0.2368283793,0.2369313400,0.2370198127,0.2371102639,0.2372102364,0.2372927449,0.2373775070,0.2374648061,0.2375410704,0.2376231135,0.2377055461,0.2377789996,0.2378504630,0.2379307254,0.2380074605,0.2380688755,0.2381346430,0.2382041255,0.2382695465,0.2383289327,0.2383855945,0.2384608860,0.2385154771,0.2385702840,0.2386344123,0.2386836168,0.2387420152,0.2387948063,0.2388454510,0.2388971050,0.2389501716,0.2389957437,0.2390551866,0.2390915580,0.2391417064,0.2391891605,0.2392400758,0.2392775125,0.2393238881,0.2393696646,0.2394053205,0.2394554781,0.2394912946,0.2395305494,0.2395702939,0.2017293057,0.2032956108,0.2047005123,0.2059640787,0.2071028147,0.2081425979,0.2090954582,0.2099720169,0.2107619732,0.2115151983,0.2121998767,0.2128422535,0.2134312203,0.2139908261,0.2145220553,0.2150069921,0.2154639917,0.2158961607,0.2163109263,0.2166958147,0.2170706584,0.2174206399,0.2177390270,0.2180581148,0.2183564171,0.2186411957,0.2189146920,0.2191725147,0.2194298727,0.2196547980,0.2199013350,0.2201109918,0.2203264971,0.2205242914,0.2207217334,0.2209128817,0.2210871583,0.2212591968,0.2214254989,0.2215891807,0.2217431678,0.2219000012,0.2220426442,0.2221858437,0.2223133938,0.2224454143,0.2225798661,0.2227053725,0.2228225260,0.2229423837,0.2230483399,0.2231646254,0.2232678211,0.2233829234,0.2234753871,0.2235774726,0.2236771156,0.2237630135,0.2238516644,0.2239496253,0.2240293057,0.2241114223,0.2241960306,0.2242723502,0.2243517754,0.2244315886,0.2245054718,0.2245729953,0.2246528131,0.2247237640,0.2247836225,0.2248483663,0.2249165638,0.2249779164,0.2250370832,0.2250926560,0.2251672545,0.2252200630,0.2252733204,0.2253350427,0.2253820141,0.2254426483,0.2254921478,0.2255416673,0.2255932381,0.2256446180,0.2256886350,0.2257466935,0.2257822594,0.2258289349,0.2258768027,0.2259280911,0.2259610148,0.2260082923,0.2260517601,0.2260866874,0.2261360641,0.2261701613,0.2262100249,0.2262475166,0.1900381021,0.1915392768,0.1928910447,0.1941056899,0.1952009567,0.1962019265,0.1971184274,0.1979622779,0.1987245376,0.1994505217,0.2001114301,0.2007309693,0.2013006415,0.2018399523,0.2023514901,0.2028243860,0.2032663968,0.2036792337,0.2040820449,0.2044568456,0.2048170929,0.2051541765,0.2054647699,0.2057706263,0.2060621585,0.2063390609,0.2066027648,0.2068530717,0.2071022547,0.2073201051,0.2075588183,0.2077588173,0.2079721898,0.2081631067,0.2083542953,0.2085398634,0.2087081838,0.2088760706,0.2090364116,0.2091961641,0.2093450527,0.2094960488,0.2096354165,0.2097748441,0.2098991371,0.2100277578,0.2101553639,0.2102791119,0.2103931360,0.2105074163,0.2106120293,0.2107252119,0.2108270333,0.2109374217,0.2110275368,0.2111249443,0.2112225073,0.2113077643,0.2113927849,0.2114876645,0.2115639030,0.2116461085,0.2117273178,0.2118003439,0.2118774160,0.2119557965,0.2120271780,0.2120932126,0.2121702683,0.2122400004,0.2123006986,0.2123611302,0.2124290188,0.2124884242,0.2125435202,0.2125994508,0.2126695012,0.2127224957,0.2127754701,0.2128357397,0.2128814014,0.2129395386,0.2129872530,0.2130351771,0.2130845940,0.2131343039,0.2131792783,0.2132362789,0.2132719078,0.2133146955,0.2133607112,0.2134104999,0.2134463735,0.2134895260,0.2135314238,0.2135648572,0.2136143908,0.2136485226,0.2136865740,0.2137219816,0.1790791038,0.1805192694,0.1818182047,0.1829829017,0.1840373480,0.1850008826,0.1858825919,0.1866947812,0.1874264196,0.1881262613,0.1887661442,0.1893622805,0.1899114215,0.1904317275,0.1909251226,0.1913824466,0.1918086580,0.1922063546,0.1925973769,0.1929586440,0.1933073292,0.1936336163,0.1939334804,0.1942301561,0.1945122395,0.1947790188,0.1950362005,0.1952761550,0.1955174902,0.1957273073,0.1959599498,0.1961521161,0.1963604935,0.1965452546,0.1967305194,0.1969102219,0.1970722921,0.1972356107,0.1973924041,0.1975449794,0.1976914565,0.1978399261,0.1979725840,0.1981090248,0.1982290992,0.1983523073,0.1984779358,0.1985954268,0.1987087258,0.1988172842,0.1989188086,0.1990287288,0.1991289918,0.1992360341,0.1993237392,0.1994178372,0.1995118671,0.1995942831,0.1996785517,0.1997712108,0.1998435460,0.1999231653,0.2000029474,0.2000716630,0.2001471189,0.2002253829,0.2002947366,0.2003563244,0.2004323494,0.2004992031,0.2005575135,0.2006190198,0.2006839609,0.2007434494,0.2007947379,0.2008506028,0.2009167555,0.2009680726,0.2010212457,0.2010781198,0.2011252535,0.2011799327,0.2012261681,0.2012735655,0.2013180001,0.2013689451,0.2014129888,0.2014681114,0.2015031358,0.2015452486,0.2015917537,0.2016350503,0.2016734063,0.2017139297,0.2017542433,0.2017891660,0.2018364676,0.2018692106,0.2019059980,0.2019420765,0.1688008962,0.1701826649,0.1714299463,0.1725478213,0.1735618028,0.1744858665,0.1753358799,0.1761144760,0.1768193829,0.1774937903,0.1781095677,0.1786837476,0.1792130468,0.1797168646,0.1801886036,0.1806306392,0.1810427344,0.1814279734,0.1818038088,0.1821532785,0.1824895526,0.1828050710,0.1830956258,0.1833805125,0.1836554013,0.1839137587,0.1841614365,0.1843925570,0.1846270337,0.1848278926,0.1850531201,0.1852390114,0.1854415399,0.1856209031,0.1857992256,0.1859734330,0.1861298173,0.1862891724,0.1864420304,0.1865896018,0.1867308421,0.1868736910,0.1870036629,0.1871349282,0.1872532004,0.1873702022,0.1874944446,0.1876048128,0.1877168218,0.1878214874,0.1879203921,0.1880267614,0.1881239553,0.1882264663,0.1883126798,0.1884031524,0.1884952368,0.1885754423,0.1886565586,0.1887457204,0.1888186068,0.1888959809,0.1889710440,0.1890368926,0.1891101304,0.1891866675,0.1892540413,0.1893143596,0.1893871653,0.1894519988,0.1895095966,0.1895681386,0.1896286600,0.1896894999,0.1897392409,0.1897957985,0.1898572200,0.1899088585,0.1899589817,0.1900149502,0.1900619348,0.1901134753,0.1901577359,0.1902048897,0.1902487827,0.1902953660,0.1903405275,0.1903937007,0.1904296678,0.1904696956,0.1905130549,0.1905565389,0.1905955505,0.1906328546,0.1906726387,0.1907063900,0.1907505954,0.1907822529,0.1908194949,0.1908532934,0.1591560469,0.1604804337,0.1616766103,0.1627496662,0.1637221855,0.1646112506,0.1654275263,0.1661769679,0.1668558808,0.1675065587,0.1680959336,0.1686479251,0.1691576274,0.1696415117,0.1700990386,0.1705242882,0.1709213851,0.1712926714,0.1716553477,0.1719925295,0.1723168360,0.1726219695,0.1729030491,0.1731782293,0.1734428330,0.1736923286,0.1739316370,0.1741572673,0.1743829150,0.1745756300,0.1747947966,0.1749747321,0.1751690799,0.1753437808,0.1755143088,0.1756848723,0.1758344725,0.1759888498,0.1761363524,0.1762804338,0.1764180081,0.1765552156,0.1766810583,0.1768085974,0.1769226259,0.1770363497,0.1771563823,0.1772642281,0.1773727643,0.1774713476,0.1775690608,0.1776727246,0.1777656376,0.1778632489,0.1779498417,0.1780353623,0.1781233287,0.1782049073,0.1782825737,0.1783687267,0.1784385708,0.1785131116,0.1785866764,0.1786505262,0.1787214641,0.1787945522,0.1788609891,0.1789191879,0.1789897650,0.1790526212,0.1791089766,0.1791650378,0.1792231019,0.1792827153,0.1793327006,0.1793876670,0.1794447082,0.1794957541,0.1795455475,0.1795984564,0.1796430673,0.1796933323,0.1797360281,0.1797821212,0.1798257723,0.1798719410,0.1799146256,0.1799656925,0.1800013635,0.1800402580,0.1800841964,0.1801236473,0.1801595136,0.1801978869,0.1802358973,0.1802684061,0.1803138739,0.1803416065,0.1803789914,0.1804138729,0.1500994996,0.1513698887,0.1525164088,0.1535455264,0.1544798276,0.1553318577,0.1561172995,0.1568345095,0.1574888432,0.1581139210,0.1586829637,0.1592130410,0.1597026187,0.1601676089,0.1606081650,0.1610199112,0.1614028038,0.1617594191,0.1621095562,0.1624350078,0.1627459748,0.1630414553,0.1633134465,0.1635789774,0.1638323606,0.1640731561,0.1643062813,0.1645242087,0.1647392363,0.1649274634,0.1651386542,0.1653113121,0.1655015516,0.1656693838,0.1658344236,0.1659979246,0.1661436041,0.1662951760,0.1664368725,0.1665736255,0.1667062133,0.1668402796,0.1669592961,0.1670843870,0.1671966818,0.1673062591,0.1674214334,0.1675269353,0.1676319639,0.1677271756,0.1678208928,0.1679218899,0.1680127928,0.1681059880,0.1681888366,0.1682725474,0.1683571239,0.1684373323,0.1685094585,0.1685938204,0.1686636146,0.1687347821,0.1688048592,0.1688689250,0.1689388478,0.1690072787,0.1690705316,0.1691288669,0.1691973933,0.1692579276,0.1693133589,0.1693679102,0.1694243204,0.1694819380,0.1695295258,0.1695814623,0.1696394102,0.1696881556,0.1697360731,0.1697876418,0.1698308708,0.1698795672,0.1699206031,0.1699628641,0.1700070994,0.1700502033,0.1700935443,0.1701414401,0.1701781533,0.1702142423,0.1702583623,0.1702958163,0.1703327685,0.1703690475,0.1704048937,0.1704362606,0.1704801382,0.1705045637,0.1705441566,0.1705788208,0.1415916431,0.1428085651,0.1439086121,0.1448959358,0.1457901146,0.1466085294,0.1473629672,0.1480515138,0.1486805009,0.1492808052,0.1498268652,0.1503378432,0.1508069021,0.1512557566,0.1516809710,0.1520755259,0.1524433597,0.1527897574,0.1531272090,0.1534390684,0.1537391218,0.1540231113,0.1542855691,0.1545409327,0.1547866977,0.1550190417,0.1552435530,0.1554545356,0.1556600373,0.1558430842,0.1560462319,0.1562131964,0.1563976852,0.1565598818,0.1567189972,0.1568780679,0.1570169263,0.1571629455,0.1573008999,0.1574323054,0.1575623238,0.1576902577,0.1578059604,0.1579278032,0.1580355462,0.1581417075,0.1582518562,0.1583542705,0.1584545265,0.1585480609,0.1586379809,0.1587369873,0.1588248560,0.1589134335,0.1589955908,0.1590748700,0.1591590831,0.1592338990,0.1593042249,0.1593842359,0.1594529106,0.1595222619,0.1595918226,0.1596536244,0.1597202174,0.1597880349,0.1598483225,0.1599035992,0.1599703848,0.1600291331,0.1600826157,0.1601351999,0.1601902635,0.1602465316,0.1602931024,0.1603421245,0.1603970321,0.1604476369,0.1604928679,0.1605421170,0.1605849722,0.1606306986,0.1606692135,0.1607128864,0.1607543448,0.1607980249,0.1608384368,0.1608863242,0.1609220556,0.1609538691,0.1609980135,0.1610355953,0.1610715178,0.1611062951,0.1611400029,0.1611700226,0.1612130409,0.1612357996,0.1612756063,0.1613081327,0.1335969783,0.1347641778,0.1358156288,0.1367598967,0.1376203337,0.1384044170,0.1391304501,0.1397883784,0.1403929540,0.1409691203,0.1414943208,0.1419843276,0.1424360228,0.1428684023,0.1432774721,0.1436565211,0.1440102910,0.1443438896,0.1446689752,0.1449688979,0.1452566157,0.1455303447,0.1457834244,0.1460315947,0.1462666171,0.1464919658,0.1467066682,0.1469102903,0.1471101778,0.1472868196,0.1474808697,0.1476439655,0.1478214574,0.1479783709,0.1481316089,0.1482850650,0.1484194665,0.1485591185,0.1486917890,0.1488206644,0.1489443152,0.1490700747,0.1491799557,0.1492982532,0.1494019825,0.1495044173,0.1496125492,0.1497095514,0.1498073883,0.1498982908,0.1499838985,0.1500787627,0.1501662981,0.1502505644,0.1503314622,0.1504065973,0.1504884277,0.1505626531,0.1506291685,0.1507051878,0.1507725396,0.1508388283,0.1509069960,0.1509689657,0.1510307340,0.1510969401,0.1511553175,0.1512097933,0.1512733137,0.1513295445,0.1513806208,0.1514339470,0.1514866919,0.1515417751,0.1515855878,0.1516316746,0.1516862162,0.1517359391,0.1517785539,0.1518271375,0.1518688467,0.1519132208,0.1519494913,0.1519923124,0.1520307800,0.1520739559,0.1521147146,0.1521590258,0.1521952276,0.1522267038,0.1522682362,0.1523037372,0.1523395541,0.1523745102,0.1524046263,0.1524350178,0.1524765234,0.1524973122,0.1525378549,0.1525692551,0.1260811836,0.1271963948,0.1282029841,0.1291084474,0.1299334665,0.1306841073,0.1313803495,0.1320113306,0.1325920718,0.1331438294,0.1336499339,0.1341202179,0.1345543510,0.1349712149,0.1353634437,0.1357273431,0.1360702437,0.1363907380,0.1367032927,0.1369920915,0.1372678107,0.1375313980,0.1377757087,0.1380159763,0.1382402613,0.1384574031,0.1386649410,0.1388609012,0.1390525623,0.1392234681,0.1394100313,0.1395682910,0.1397396641,0.1398897014,0.1400361991,0.1401873227,0.1403174711,0.1404513520,0.1405798411,0.1407016968,0.1408221339,0.1409448875,0.1410506027,0.1411627558,0.1412646517,0.1413642104,0.1414680218,0.1415625834,0.1416566203,0.1417444568,0.1418258829,0.1419176371,0.1420018037,0.1420842174,0.1421616791,0.1422347496,0.1423161474,0.1423874933,0.1424512776,0.1425231028,0.1425873242,0.1426527312,0.1427186092,0.1427778100,0.1428381382,0.1429009615,0.1429603332,0.1430114527,0.1430734943,0.1431280743,0.1431769803,0.1432299400,0.1432781857,0.1433316332,0.1433744712,0.1434175681,0.1434718874,0.1435178790,0.1435608069,0.1436084020,0.1436487942,0.1436916049,0.1437252867,0.1437680281,0.1438052732,0.1438488069,0.1438861057,0.1439294797,0.1439667751,0.1439962271,0.1440338516,0.1440691853,0.1441043978,0.1441382663,0.1441679322,0.1441965443,0.1442368461,0.1442570098,0.1442955388,0.1443276755,0.1190122224,0.1200780163,0.1210417532,0.1219086339,0.1226981156,0.1234174889,0.1240843285,0.1246907569,0.1252473936,0.1257773023,0.1262636226,0.1267125187,0.1271311513,0.1275326396,0.1279073970,0.1282582245,0.1285880958,0.1288976246,0.1291943394,0.1294747141,0.1297391556,0.1299936182,0.1302297354,0.1304613085,0.1306763602,0.1308869082,0.1310858053,0.1312738997,0.1314580104,0.1316240187,0.1318020107,0.1319565537,0.1321206477,0.1322656161,0.1324075850,0.1325531021,0.1326785495,0.1328076807,0.1329311250,0.1330485189,0.1331651297,0.1332816626,0.1333862953,0.1334937963,0.1335918720,0.1336876869,0.1337871605,0.1338793707,0.1339692286,0.1340550018,0.1341331691,0.1342210562,0.1343015606,0.1343829991,0.1344589947,0.1345278100,0.1346051332,0.1346750065,0.1347373870,0.1348060351,0.1348669309,0.1349323112,0.1349962545,0.1350518706,0.1351109978,0.1351710630,0.1352293343,0.1352786280,0.1353387096,0.1353896774,0.1354390601,0.1354899631,0.1355368665,0.1355867257,0.1356292966,0.1356714249,0.1357250007,0.1357675226,0.1358096076,0.1358578791,0.1358963514,0.1359375770,0.1359686884,0.1360093279,0.1360450053,0.1360874953,0.1361229745,0.1361650368,0.1362021774,0.1362311221,0.1362679974,0.1363022741,0.1363369454,0.1363684100,0.1363986827,0.1364246862,0.1364633088,0.1364827832,0.1365188747,0.1365505676,0.1123594636,0.1133795543,0.1143005823,0.1151322502,0.1158890818,0.1165740819,0.1172154252,0.1177978419,0.1183313188,0.1188377590,0.1193054801,0.1197350378,0.1201369810,0.1205227347,0.1208819011,0.1212197544,0.1215367434,0.1218346147,0.1221203476,0.1223886992,0.1226433031,0.1228895966,0.1231156478,0.1233369876,0.1235450157,0.1237488328,0.1239388686,0.1241210793,0.1242992293,0.1244589030,0.1246293553,0.1247773541,0.1249367845,0.1250765195,0.1252119613,0.1253535993,0.1254734715,0.1255976419,0.1257174220,0.1258294788,0.1259421579,0.1260543086,0.1261558776,0.1262597608,0.1263547572,0.1264474188,0.1265424650,0.1266327970,0.1267195024,0.1268015290,0.1268783244,0.1269616714,0.1270396054,0.1271173848,0.1271902334,0.1272573317,0.1273309169,0.1273986901,0.1274594119,0.1275248418,0.1275863371,0.1276486385,0.1277091912,0.1277636775,0.1278198610,0.1278775886,0.1279350964,0.1279817261,0.1280403504,0.1280890511,0.1281380666,0.1281876207,0.1282322977,0.1282812764,0.1283204560,0.1283617498,0.1284153198,0.1284551934,0.1284956991,0.1285417724,0.1285780001,0.1286182454,0.1286478210,0.1286899838,0.1287222225,0.1287634908,0.1287975001,0.1288387407,0.1288761384,0.1289027463,0.1289373395,0.1289724040,0.1290035629,0.1290364589,0.1290650048,0.1290902264,0.1291264090,0.1291463545,0.1291804792,0.1292123697,0.1060978509,0.1070728474,0.1079535657,0.1087487642,0.1094738045,0.1101296530,0.1107438186,0.1113014424,0.1118138021,0.1123006185,0.1127474913,0.1131596008,0.1135471914,0.1139160331,0.1142603687,0.1145864750,0.1148908370,0.1151757737,0.1154509982,0.1157092311,0.1159538404,0.1161912175,0.1164075392,0.1166205541,0.1168211724,0.1170153390,0.1171989937,0.1173741024,0.1175478423,0.1176999895,0.1178650712,0.1180068466,0.1181596571,0.1182955576,0.1184251498,0.1185612096,0.1186785286,0.1187972580,0.1189133483,0.1190185638,0.1191283120,0.1192365667,0.1193353571,0.1194352337,0.1195281224,0.1196154430,0.1197085458,0.1197947967,0.1198771013,0.1199565851,0.1200310454,0.1201111918,0.1201868347,0.1202622426,0.1203315992,0.1203965981,0.1204683749,0.1205341394,0.1205916270,0.1206535687,0.1207136633,0.1207745325,0.1208333231,0.1208848970,0.1209413115,0.1209924568,0.1210514374,0.1210946575,0.1211530124,0.1212001392,0.1212470021,0.1212928638,0.1213367041,0.1213863408,0.1214225418,0.1214627645,0.1215148065,0.1215520802,0.1215930329,0.1216359885,0.1216723821,0.1217105407,0.1217403206,0.1217805447,0.1218107901,0.1218507959,0.1218836957,0.1219210263,0.1219598664,0.1219854506,0.1220191634,0.1220519460,0.1220818973,0.1221159495,0.1221399473,0.1221679204,0.1222017254,0.1222202760,0.1222542681,0.1222848600,0.1002017364,0.1011342340,0.1019751113,0.1027359307,0.1034300597,0.1040573378,0.1046464209,0.1051804678,0.1056714804,0.1061375043,0.1065660142,0.1069622207,0.1073338609,0.1076877252,0.1080168169,0.1083299905,0.1086217760,0.1088964767,0.1091609344,0.1094083683,0.1096446188,0.1098735667,0.1100805109,0.1102865278,0.1104773094,0.1106645733,0.1108418417,0.1110082633,0.1111767386,0.1113231256,0.1114819867,0.1116185896,0.1117669387,0.1118972860,0.1120211810,0.1121518659,0.1122652055,0.1123807456,0.1124903882,0.1125926073,0.1126986690,0.1128018184,0.1128969156,0.1129947826,0.1130843324,0.1131690078,0.1132573258,0.1133416706,0.1134200455,0.1134974051,0.1135687061,0.1136458501,0.1137193658,0.1137915357,0.1138592421,0.1139209654,0.1139895716,0.1140543858,0.1141086885,0.1141687546,0.1142275693,0.1142864753,0.1143412694,0.1143927348,0.1144474427,0.1144950017,0.1145530867,0.1145950103,0.1146521112,0.1146969930,0.1147399611,0.1147861231,0.1148280157,0.1148740335,0.1149102181,0.1149508262,0.1150001845,0.1150357179,0.1150756966,0.1151163617,0.1151540526,0.1151891404,0.1152179504,0.1152562204,0.1152856551,0.1153231710,0.1153557076,0.1153916733,0.1154288536,0.1154542324,0.1154856705,0.1155183855,0.1155478281,0.1155804891,0.1156032147,0.1156303613,0.1156640394,0.1156826882,0.1157137508,0.1157434130,0.0946496320,0.0955399590,0.0963429196,0.0970717161,0.0977344885,0.0983350439,0.0988990071,0.0994104388,0.0998801103,0.1003273056,0.1007386104,0.1011179872,0.1014745141,0.1018133617,0.1021290251,0.1024292149,0.1027093127,0.1029742557,0.1032288656,0.1034647664,0.1036936310,0.1039121429,0.1041119429,0.1043092150,0.1044925670,0.1046732060,0.1048432059,0.1050034463,0.1051652827,0.1053059129,0.1054580508,0.1055896376,0.1057336202,0.1058590323,0.1059774717,0.1061021086,0.1062126327,0.1063236687,0.1064277442,0.1065279580,0.1066303656,0.1067276883,0.1068214627,0.1069139849,0.1070007660,0.1070828079,0.1071674458,0.1072496008,0.1073233941,0.1073983936,0.1074693328,0.1075433344,0.1076118469,0.1076815415,0.1077473215,0.1078068865,0.1078717481,0.1079357845,0.1079872530,0.1080457435,0.1081025536,0.1081601980,0.1082131431,0.1082634767,0.1083160694,0.1083602532,0.1084169593,0.1084564709,0.1085106668,0.1085558268,0.1085963892,0.1086415588,0.1086823435,0.1087248852,0.1087615252,0.1088009399,0.1088461547,0.1088819784,0.1089206830,0.1089597778,0.1089961493,0.1090309605,0.1090595325,0.1090962108,0.1091237284,0.1091596826,0.1091918574,0.1092268731,0.1092607060,0.1092871409,0.1093168328,0.1093480417,0.1093772109,0.1094082877,0.1094292672,0.1094566771,0.1094876585,0.1095061053,0.1095377602,0.1095651054,0.0894183372,0.0902693178,0.0910363304,0.0917335324,0.0923664180,0.0929403036,0.0934818408,0.0939702130,0.0944202489,0.0948481447,0.0952424850,0.0956064797,0.0959479085,0.0962721136,0.0965754026,0.0968635041,0.0971317402,0.0973855474,0.0976319639,0.0978578452,0.0980770791,0.0982859452,0.0984787285,0.0986679595,0.0988446141,0.0990177123,0.0991810082,0.0993350527,0.0994913210,0.0996268790,0.0997731016,0.0998999501,0.1000366044,0.1001593952,0.1002722222,0.1003908520,0.1004974609,0.1006066425,0.1007051366,0.1008016676,0.1009009441,0.1009944133,0.1010842904,0.1011735674,0.1012574153,0.1013352151,0.1014170126,0.1014974095,0.1015675148,0.1016398996,0.1017072796,0.1017798069,0.1018453921,0.1019139849,0.1019758045,0.1020335972,0.1020964996,0.1021570169,0.1022072719,0.1022627198,0.1023183110,0.1023741854,0.1024260000,0.1024722153,0.1025239590,0.1025668551,0.1026222793,0.1026590590,0.1027117425,0.1027551326,0.1027939061,0.1028376000,0.1028770142,0.1029186793,0.1029534296,0.1029927576,0.1030361160,0.1030693977,0.1031079098,0.1031445298,0.1031802347,0.1032129444,0.1032403541,0.1032765984,0.1033043232,0.1033371151,0.1033686503,0.1034019800,0.1034355719,0.1034607611,0.1034893153,0.1035197108,0.1035475564,0.1035785905,0.1035972429,0.1036262810,0.1036542396,0.1036711760,0.1037026858,0.1037275751,0.0844888920,0.0853010986,0.0860337261,0.0867001615,0.0873048720,0.0878543846,0.0883722496,0.0888399745,0.0892704619,0.0896808762,0.0900573869,0.0904068304,0.0907346921,0.0910448436,0.0913351997,0.0916109392,0.0918684883,0.0921125957,0.0923485268,0.0925659529,0.0927760050,0.0929767759,0.0931621205,0.0933430938,0.0935132690,0.0936809137,0.0938361611,0.0939837609,0.0941340738,0.0942647461,0.0944044886,0.0945270798,0.0946589342,0.0947744491,0.0948862337,0.0949996503,0.0951011727,0.0952061823,0.0953020119,0.0953954901,0.0954897397,0.0955805562,0.0956664721,0.0957519475,0.0958324896,0.0959071027,0.0959872041,0.0960634213,0.0961320716,0.0962010728,0.0962650868,0.0963351361,0.0963990635,0.0964650894,0.0965245371,0.0965789284,0.0966404072,0.0967004575,0.0967475438,0.0967997557,0.0968538658,0.0969069487,0.0969579295,0.0970025955,0.0970527734,0.0970939913,0.0971473504,0.0971828415,0.0972333082,0.0972748947,0.0973115800,0.0973561973,0.0973926247,0.0974322330,0.0974658663,0.0975050577,0.0975460137,0.0975787434,0.0976155846,0.0976516698,0.0976846014,0.0977165692,0.0977432275,0.0977778077,0.0978043690,0.0978355459,0.0978661147,0.0978997608,0.0979303196,0.0979559176,0.0979827097,0.0980128147,0.0980394817,0.0980683411,0.0980866322,0.0981145299,0.0981420164,0.0981597883,0.0981890377,0.0982127042,0.0798412258,0.0806175776,0.0813176660,0.0819525846,0.0825319875,0.0830569110,0.0835521824,0.0840005769,0.0844146122,0.0848047942,0.0851654560,0.0855004192,0.0858135955,0.0861113466,0.0863897900,0.0866537917,0.0869005198,0.0871348081,0.0873616588,0.0875701585,0.0877721131,0.0879638763,0.0881421802,0.0883166558,0.0884790973,0.0886410705,0.0887890027,0.0889304438,0.0890761406,0.0892027269,0.0893369347,0.0894524815,0.0895785631,0.0896908424,0.0897994603,0.0899073093,0.0900053923,0.0901052077,0.0901976640,0.0902885396,0.0903777835,0.0904654492,0.0905478896,0.0906310623,0.0907083592,0.0907797234,0.0908568505,0.0909302952,0.0909972563,0.0910622024,0.0911248175,0.0911931562,0.0912528613,0.0913179166,0.0913765430,0.0914253996,0.0914854702,0.0915435465,0.0915889722,0.0916400153,0.0916911007,0.0917431659,0.0917914535,0.0918361192,0.0918831907,0.0919224747,0.0919740597,0.0920091305,0.0920576365,0.0920961962,0.0921325149,0.0921744624,0.0922104405,0.0922487715,0.0922809518,0.0923192834,0.0923573525,0.0923899484,0.0924261923,0.0924600989,0.0924920849,0.0925223592,0.0925473918,0.0925803620,0.0926062801,0.0926371285,0.0926670726,0.0926982644,0.0927296972,0.0927534657,0.0927793794,0.0928062468,0.0928341219,0.0928610336,0.0928789502,0.0929073710,0.0929331334,0.0929515591,0.0929798853,0.0930016299,0.0754597380,0.0761999815,0.0768693237,0.0774758042,0.0780299725,0.0785313830,0.0790049395,0.0794336191,0.0798306193,0.0802031106,0.0805486520,0.0808705817,0.0811696757,0.0814546504,0.0817213771,0.0819756145,0.0822103380,0.0824355336,0.0826531699,0.0828538873,0.0830470711,0.0832305023,0.0834016392,0.0835698480,0.0837255911,0.0838800345,0.0840225329,0.0841595185,0.0842978639,0.0844195574,0.0845500850,0.0846612220,0.0847805358,0.0848902323,0.0849935619,0.0850972840,0.0851915912,0.0852886360,0.0853753138,0.0854632086,0.0855488001,0.0856344319,0.0857123719,0.0857934544,0.0858668606,0.0859368118,0.0860098208,0.0860792802,0.0861459805,0.0862089489,0.0862692046,0.0863342592,0.0863911108,0.0864536272,0.0865085239,0.0865571603,0.0866146327,0.0866706596,0.0867150304,0.0867639211,0.0868122930,0.0868630153,0.0869091204,0.0869521340,0.0869981461,0.0870348512,0.0870856993,0.0871192977,0.0871668672,0.0872025282,0.0872364444,0.0872788064,0.0873123810,0.0873500532,0.0873811252,0.0874172719,0.0874524500,0.0874859409,0.0875203760,0.0875540669,0.0875843621,0.0876139517,0.0876374941,0.0876696001,0.0876940009,0.0877240610,0.0877528427,0.0877833563,0.0878132689,0.0878361316,0.0878606644,0.0878880909,0.0879141627,0.0879394907,0.0879567717,0.0879833558,0.0880089873,0.0880266059,0.0880527172,0.0880749537,0.0713276887,0.0720342683,0.0726730332,0.0732509800,0.0737812966,0.0742604192,0.0747148335,0.0751240088,0.0755027122,0.0758598961,0.0761897112,0.0764983582,0.0767847795,0.0770578986,0.0773138078,0.0775572139,0.0777829920,0.0779982550,0.0782068625,0.0783997740,0.0785837424,0.0787598342,0.0789254184,0.0790861606,0.0792348208,0.0793847123,0.0795205134,0.0796516651,0.0797847343,0.0799020070,0.0800261486,0.0801331114,0.0802473016,0.0803533548,0.0804521252,0.0805527682,0.0806425628,0.0807358381,0.0808213889,0.0809051971,0.0809859612,0.0810677439,0.0811428378,0.0812210995,0.0812926909,0.0813593569,0.0814301795,0.0814976823,0.0815603910,0.0816220160,0.0816787685,0.0817419687,0.0817949290,0.0818565875,0.0819093168,0.0819562753,0.0820110359,0.0820661598,0.0821091596,0.0821552967,0.0822015660,0.0822495867,0.0822946058,0.0823362617,0.0823815824,0.0824154709,0.0824650809,0.0824972477,0.0825448275,0.0825767706,0.0826102479,0.0826499012,0.0826825473,0.0827203583,0.0827485694,0.0827832598,0.0828171135,0.0828500261,0.0828837249,0.0829153247,0.0829442263,0.0829747711,0.0829962602,0.0830278438,0.0830498695,0.0830802441,0.0831076337,0.0831361898,0.0831645056,0.0831873229,0.0832098831,0.0832379291,0.0832623937,0.0832868454,0.0833044081,0.0833287373,0.0833540541,0.0833712476,0.0833943858,0.0834173096,0.0674287923,0.0681039668,0.0687143995,0.0692651678,0.0697717649,0.0702303180,0.0706641663,0.0710560081,0.0714187402,0.0717581422,0.0720749367,0.0723708365,0.0726440984,0.0729061868,0.0731526404,0.0733845635,0.0736006806,0.0738079812,0.0740066062,0.0741923997,0.0743694777,0.0745376702,0.0746961429,0.0748496662,0.0749929006,0.0751363845,0.0752664668,0.0753938135,0.0755203478,0.0756331364,0.0757520739,0.0758552902,0.0759636714,0.0760659057,0.0761619279,0.0762577522,0.0763429944,0.0764338063,0.0765170411,0.0765968773,0.0766735504,0.0767523648,0.0768246572,0.0768990427,0.0769678386,0.0770329107,0.0771006118,0.0771657772,0.0772269274,0.0772857583,0.0773416170,0.0774004126,0.0774511344,0.0775112759,0.0775623333,0.0776062465,0.0776604254,0.0777127883,0.0777532400,0.0777981381,0.0778424679,0.0778894159,0.0779314588,0.0779711263,0.0780168815,0.0780489446,0.0780973115,0.0781263274,0.0781735056,0.0782036664,0.0782367593,0.0782753590,0.0783057565,0.0783421086,0.0783678430,0.0784032653,0.0784348477,0.0784662925,0.0785002109,0.0785305379,0.0785571932,0.0785873149,0.0786081243,0.0786386070,0.0786597889,0.0786883653,0.0787142736,0.0787430166,0.0787691667,0.0787907081,0.0788140084,0.0788401034,0.0788650908,0.0788877723,0.0789058888,0.0789272392,0.0789521020,0.0789694029,0.0789911863,0.0790130531,0.0637517560,0.0643964319,0.0649776103,0.0655032155,0.0659867694,0.0664263587,0.0668412671,0.0672149451,0.0675614348,0.0678868510,0.0681888169,0.0684730373,0.0687351721,0.0689855103,0.0692208689,0.0694434119,0.0696497554,0.0698475483,0.0700390592,0.0702180090,0.0703872760,0.0705487135,0.0707007013,0.0708470086,0.0709860582,0.0711236210,0.0712469041,0.0713686771,0.0714912495,0.0716004229,0.0717146059,0.0718127810,0.0719168684,0.0720146064,0.0721078545,0.0721989002,0.0722800998,0.0723675955,0.0724477226,0.0725238717,0.0725979176,0.0726737577,0.0727446389,0.0728144662,0.0728800022,0.0729428997,0.0730087974,0.0730702587,0.0731290127,0.0731864013,0.0732387856,0.0732962601,0.0733446776,0.0734032591,0.0734525732,0.0734946621,0.0735462305,0.0735980747,0.0736342607,0.0736791457,0.0737213352,0.0737668017,0.0738054865,0.0738457442,0.0738897756,0.0739198643,0.0739673339,0.0739944511,0.0740404717,0.0740682097,0.0741015101,0.0741372471,0.0741664756,0.0742019275,0.0742271242,0.0742630568,0.0742911613,0.0743216043,0.0743541254,0.0743841226,0.0744083509,0.0744382738,0.0744580914,0.0744878116,0.0745089059,0.0745360398,0.0745607038,0.0745879158,0.0746131238,0.0746338176,0.0746564033,0.0746804201,0.0747059866,0.0747280959,0.0747449813,0.0747652888,0.0747900218,0.0748068479,0.0748261289,0.0748482304,0.0602819825,0.0608962293,0.0614508939,0.0619516892,0.0624148629,0.0628328081,0.0632313830,0.0635869426,0.0639204762,0.0642305676,0.0645193163,0.0647905772,0.0650424704,0.0652814958,0.0655067303,0.0657205770,0.0659183202,0.0661071331,0.0662902489,0.0664619282,0.0666246999,0.0667793606,0.0669250995,0.0670653583,0.0671983242,0.0673302872,0.0674493926,0.0675647181,0.0676828208,0.0677873816,0.0678974033,0.0679910569,0.0680912373,0.0681849366,0.0682741479,0.0683619476,0.0684409534,0.0685226109,0.0686006198,0.0686747307,0.0687445796,0.0688181577,0.0688856427,0.0689530999,0.0690161974,0.0690754742,0.0691408418,0.0691981978,0.0692548065,0.0693099229,0.0693609431,0.0694168625,0.0694622910,0.0695190303,0.0695663653,0.0696063942,0.0696561796,0.0697064088,0.0697413041,0.0697849982,0.0698247969,0.0698688149,0.0699057841,0.0699433421,0.0699863600,0.0700162915,0.0700628306,0.0700873547,0.0701313908,0.0701578797,0.0701896532,0.0702244386,0.0702524228,0.0702873436,0.0703107334,0.0703453177,0.0703717633,0.0704031227,0.0704331654,0.0704614012,0.0704851939,0.0705142363,0.0705320861,0.0705624719,0.0705820685,0.0706079818,0.0706322739,0.0706585133,0.0706819245,0.0707022513,0.0707236480,0.0707476138,0.0707726397,0.0707925006,0.0708099581,0.0708292384,0.0708528445,0.0708690689,0.0708879571,0.0709087566,0.0570067153,0.0575942206,0.0581217133,0.0585987872,0.0590419997,0.0594408581,0.0598209829,0.0601626500,0.0604807027,0.0607770061,0.0610526027,0.0613119674,0.0615537425,0.0617821281,0.0619985883,0.0622029744,0.0623911647,0.0625721427,0.0627473009,0.0629128060,0.0630692727,0.0632171852,0.0633560208,0.0634901758,0.0636185507,0.0637446216,0.0638581173,0.0639693182,0.0640829571,0.0641829227,0.0642887783,0.0643778806,0.0644737439,0.0645645283,0.0646503034,0.0647354092,0.0648096579,0.0648883280,0.0649641863,0.0650355344,0.0651016318,0.0651718252,0.0652368413,0.0653018448,0.0653628622,0.0654177930,0.0654839387,0.0655368160,0.0655915950,0.0656444509,0.0656940720,0.0657481572,0.0657911180,0.0658451358,0.0658905246,0.0659298731,0.0659767014,0.0660250372,0.0660600476,0.0661014402,0.0661394231,0.0661818199,0.0662164896,0.0662532440,0.0662950656,0.0663249445,0.0663681004,0.0663921424,0.0664329051,0.0664598640,0.0664894034,0.0665229160,0.0665504735,0.0665843668,0.0666073438,0.0666395686,0.0666650967,0.0666960514,0.0667238139,0.0667520190,0.0667732959,0.0668020099,0.0668198596,0.0668490239,0.0668679698,0.0668926555,0.0669158269,0.0669414208,0.0669636822,0.0669832977,0.0670036099,0.0670263270,0.0670498408,0.0670693668,0.0670868881,0.0671048943,0.0671274498,0.0671442204,0.0671614522,0.0671813456,0.0539142107,0.0544743690,0.0549789419,0.0554336537,0.0558565074,0.0562372801,0.0566004746,0.0569263665,0.0572307727,0.0575141465,0.0577768568,0.0580270592,0.0582570762,0.0584757234,0.0586833633,0.0588787012,0.0590584882,0.0592317518,0.0593995552,0.0595579455,0.0597078988,0.0598489831,0.0599833382,0.0601111866,0.0602347050,0.0603557813,0.0604641141,0.0605717632,0.0606789736,0.0607747131,0.0608759054,0.0609633247,0.0610540450,0.0611409432,0.0612227071,0.0613053316,0.0613767573,0.0614521060,0.0615237135,0.0615939118,0.0616567186,0.0617239981,0.0617862905,0.0618499920,0.0619084267,0.0619615627,0.0620237787,0.0620749986,0.0621271023,0.0621780551,0.0622244811,0.0622777146,0.0623183756,0.0623702640,0.0624141945,0.0624527506,0.0624984803,0.0625431862,0.0625777132,0.0626161079,0.0626522888,0.0626948140,0.0627271383,0.0627639641,0.0628026559,0.0628316894,0.0628720422,0.0628967099,0.0629351915,0.0629618145,0.0629905024,0.0630218575,0.0630482827,0.0630808306,0.0631035249,0.0631351097,0.0631583266,0.0631878147,0.0632152950,0.0632418212,0.0632622034,0.0632895193,0.0633069630,0.0633351983,0.0633535635,0.0633779726,0.0633995020,0.0634236404,0.0634454478,0.0634644176,0.0634839658,0.0635057013,0.0635284013,0.0635471742,0.0635645952,0.0635816869,0.0636036769,0.0636199881,0.0636359845,0.0636552486,0.0509936821,0.0515291854,0.0520110458,0.0524441445,0.0528486167,0.0532113971,0.0535578123,0.0538705816,0.0541602552,0.0544307076,0.0546823898,0.0549223477,0.0551418715,0.0553505670,0.0555499395,0.0557366036,0.0559083491,0.0560750395,0.0562354370,0.0563863319,0.0565300737,0.0566652760,0.0567940125,0.0569155920,0.0570340213,0.0571503586,0.0572549605,0.0573576792,0.0574605502,0.0575527434,0.0576485648,0.0577332572,0.0578190824,0.0579037664,0.0579825220,0.0580614561,0.0581298441,0.0582025992,0.0582710469,0.0583385847,0.0583985592,0.0584626218,0.0585226566,0.0585833607,0.0586406721,0.0586913181,0.0587504469,0.0588003179,0.0588502152,0.0588989393,0.0589430618,0.0589950720,0.0590345740,0.0590832688,0.0591257608,0.0591625022,0.0592064871,0.0592497754,0.0592837125,0.0593189175,0.0593544913,0.0593958044,0.0594266347,0.0594631657,0.0594979737,0.0595277436,0.0595658405,0.0595895637,0.0596264068,0.0596531654,0.0596790558,0.0597101899,0.0597349373,0.0597654411,0.0597872151,0.0598193105,0.0598407422,0.0598691210,0.0598955791,0.0599209600,0.0599394756,0.0599668524,0.0599841359,0.0600098578,0.0600292258,0.0600520392,0.0600729383,0.0600956219,0.0601171216,0.0601348891,0.0601546126,0.0601749488,0.0601970450,0.0602136259,0.0602312443,0.0602478304,0.0602684156,0.0602858034,0.0602998563,0.0603185656,0.0482376334,0.0487473877,0.0492063150,0.0496201698,0.0500068583,0.0503534648,0.0506837772,0.0509821883,0.0512594315,0.0515174510,0.0517578958,0.0519886270,0.0521971342,0.0523972934,0.0525873626,0.0527666752,0.0529309385,0.0530901058,0.0532440604,0.0533880745,0.0535259138,0.0536553863,0.0537787139,0.0538954600,0.0540089977,0.0541198493,0.0542208614,0.0543196779,0.0544165186,0.0545050523,0.0545974332,0.0546785431,0.0547609132,0.0548430754,0.0549171056,0.0549930588,0.0550589205,0.0551282443,0.0551941771,0.0552593846,0.0553171244,0.0553782531,0.0554356097,0.0554937571,0.0555487703,0.0555976728,0.0556548425,0.0557016785,0.0557502257,0.0557974482,0.0558392654,0.0558893607,0.0559268164,0.0559748720,0.0560146378,0.0560504997,0.0560919491,0.0561331581,0.0561664464,0.0562006230,0.0562343740,0.0562736419,0.0563025189,0.0563383936,0.0563721057,0.0564007101,0.0564378650,0.0564604566,0.0564938441,0.0565216953,0.0565467809,0.0565757640,0.0565992734,0.0566291363,0.0566495539,0.0566816766,0.0567012315,0.0567294877,0.0567540946,0.0567788136,0.0567966253,0.0568228464,0.0568382063,0.0568639742,0.0568829535,0.0569043561,0.0569245603,0.0569462052,0.0569673628,0.0569855796,0.0570019324,0.0570231858,0.0570437631,0.0570594432,0.0570773587,0.0570936444,0.0571126564,0.0571293986,0.0571429394,0.0571598109,0.0456348884,0.0461202547,0.0465582966,0.0469530989,0.0473223087,0.0476530765,0.0479675428,0.0482529349,0.0485176716,0.0487635916,0.0489937076,0.0492152506,0.0494146242,0.0496047250,0.0497870938,0.0499578702,0.0501153830,0.0502678475,0.0504152573,0.0505531454,0.0506860544,0.0508083869,0.0509269871,0.0510392056,0.0511479321,0.0512536314,0.0513513635,0.0514450712,0.0515383685,0.0516225159,0.0517125225,0.0517899004,0.0518681657,0.0519472655,0.0520195284,0.0520905055,0.0521542599,0.0522205235,0.0522837801,0.0523464687,0.0524014387,0.0524594146,0.0525154351,0.0525714293,0.0526238565,0.0526707781,0.0527255816,0.0527705890,0.0528163348,0.0528619855,0.0529029215,0.0529496790,0.0529867377,0.0530326205,0.0530721800,0.0531055170,0.0531456411,0.0531849839,0.0532166632,0.0532500579,0.0532827104,0.0533194301,0.0533477622,0.0533807296,0.0534147447,0.0534413672,0.0534777386,0.0534989398,0.0535296722,0.0535576443,0.0535825827,0.0536102419,0.0536330334,0.0536599779,0.0536806432,0.0537116258,0.0537300776,0.0537574301,0.0537820972,0.0538044700,0.0538223827,0.0538475584,0.0538621076,0.0538863531,0.0539058830,0.0539250474,0.0539449605,0.0539652297,0.0539859608,0.0540031073,0.0540193238,0.0540395832,0.0540588300,0.0540748371,0.0540913972,0.0541080168,0.0541257740,0.0541417598,0.0541544508,0.0541710983,0.0431761979,0.0436378511,0.0440561298,0.0444327443,0.0447857180,0.0451011570,0.0454015702,0.0456733851,0.0459262095,0.0461601624,0.0463805671,0.0465939109,0.0467842071,0.0469654238,0.0471395815,0.0473030107,0.0474537263,0.0475988492,0.0477397304,0.0478723241,0.0479993649,0.0481167032,0.0482305117,0.0483384040,0.0484407488,0.0485433211,0.0486360504,0.0487269627,0.0488163820,0.0488962700,0.0489819121,0.0490572001,0.0491324330,0.0492070002,0.0492768291,0.0493446146,0.0494058476,0.0494687145,0.0495302657,0.0495913916,0.0496429412,0.0496989637,0.0497523242,0.0498073026,0.0498552553,0.0499013062,0.0499541027,0.0499969418,0.0500416787,0.0500857615,0.0501236848,0.0501682783,0.0502039826,0.0502485729,0.0502871066,0.0503183144,0.0503582369,0.0503955143,0.0504248220,0.0504572984,0.0504884758,0.0505237538,0.0505509120,0.0505825410,0.0506151590,0.0506411656,0.0506750257,0.0506959388,0.0507261242,0.0507526958,0.0507771152,0.0508028769,0.0508245737,0.0508512724,0.0508709765,0.0508998158,0.0509182066,0.0509446658,0.0509680592,0.0509895817,0.0510072102,0.0510301500,0.0510463300,0.0510688286,0.0510859963,0.0511059847,0.0511232494,0.0511452520,0.0511630852,0.0511809893,0.0511959104,0.0512155569,0.0512340043,0.0512494406,0.0512658756,0.0512811442,0.0512990814,0.0513146489,0.0513262837,0.0513430084,0.0408540136,0.0412925510,0.0416915170,0.0420509086,0.0423881531,0.0426887452,0.0429749853,0.0432353707,0.0434770689,0.0437007775,0.0439109176,0.0441141001,0.0442962482,0.0444692330,0.0446365907,0.0447915412,0.0449372290,0.0450752464,0.0452096950,0.0453366661,0.0454588375,0.0455713341,0.0456795320,0.0457822426,0.0458808966,0.0459789743,0.0460675868,0.0461545896,0.0462403218,0.0463176857,0.0463995668,0.0464712475,0.0465429617,0.0466149758,0.0466823488,0.0467468368,0.0468054162,0.0468663514,0.0469253168,0.0469835890,0.0470321476,0.0470858705,0.0471369221,0.0471900352,0.0472365792,0.0472806595,0.0473307739,0.0473720445,0.0474154849,0.0474575845,0.0474950380,0.0475366039,0.0475711095,0.0476138897,0.0476514698,0.0476811950,0.0477184294,0.0477544733,0.0477828887,0.0478147979,0.0478441796,0.0478771679,0.0479037211,0.0479338377,0.0479652592,0.0479904967,0.0480223885,0.0480427357,0.0480726244,0.0480971441,0.0481211051,0.0481458879,0.0481669867,0.0481922869,0.0482107300,0.0482381586,0.0482564349,0.0482810646,0.0483048355,0.0483254637,0.0483423362,0.0483638800,0.0483808085,0.0484008846,0.0484180252,0.0484364633,0.0484536147,0.0484753814,0.0484924256,0.0485095670,0.0485228112,0.0485420700,0.0485606936,0.0485746211,0.0485909020,0.0486050754,0.0486222998,0.0486379902,0.0486491440,0.0486652928,0.0386590727,0.0390765014,0.0394580153,0.0398004153,0.0401214381,0.0404087711,0.0406809539,0.0409299948,0.0411620884,0.0413752665,0.0415755894,0.0417697654,0.0419421720,0.0421088234,0.0422683701,0.0424171264,0.0425568677,0.0426883985,0.0428165117,0.0429373333,0.0430559558,0.0431641698,0.0432677510,0.0433642432,0.0434592843,0.0435529470,0.0436378466,0.0437208024,0.0438033353,0.0438781525,0.0439561927,0.0440246627,0.0440934616,0.0441627289,0.0442263317,0.0442889982,0.0443444134,0.0444038461,0.0444600250,0.0445162355,0.0445627157,0.0446130553,0.0446620392,0.0447137497,0.0447586028,0.0448004245,0.0448484698,0.0448885975,0.0449297092,0.0449694134,0.0450062489,0.0450455922,0.0450788464,0.0451200817,0.0451559854,0.0451838317,0.0452199631,0.0452555890,0.0452815325,0.0453131817,0.0453412743,0.0453726761,0.0453990510,0.0454265033,0.0454563157,0.0454815519,0.0455121436,0.0455308673,0.0455602412,0.0455836831,0.0456066613,0.0456305059,0.0456506753,0.0456747712,0.0456929808,0.0457194878,0.0457373256,0.0457605391,0.0457837355,0.0458033983,0.0458189770,0.0458396263,0.0458567353,0.0458757695,0.0458913567,0.0459088502,0.0459267086,0.0459469556,0.0459642895,0.0459803480,0.0459936644,0.0460110304,0.0460291671,0.0460423034,0.0460569881,0.0460711020,0.0460890420,0.0461028592,0.0461134193,0.0461299873,0.0365852168,0.0369829861,0.0373467261,0.0376733077,0.0379787203,0.0382541310,0.0385119437,0.0387503532,0.0389722976,0.0391765792,0.0393671682,0.0395531183,0.0397169352,0.0398763742,0.0400294244,0.0401719433,0.0403049509,0.0404310368,0.0405526408,0.0406684628,0.0407816024,0.0408853733,0.0409843848,0.0410771442,0.0411685382,0.0412572424,0.0413389508,0.0414187354,0.0414976428,0.0415690275,0.0416441030,0.0417088745,0.0417746895,0.0418416924,0.0419029061,0.0419625763,0.0420161701,0.0420735726,0.0421266738,0.0421806033,0.0422258529,0.0422735010,0.0423205003,0.0423699218,0.0424126317,0.0424527825,0.0424981137,0.0425373702,0.0425765046,0.0426153322,0.0426492665,0.0426879845,0.0427199622,0.0427587616,0.0427940880,0.0428216834,0.0428550712,0.0428894771,0.0429143160,0.0429440732,0.0429719894,0.0430024641,0.0430263829,0.0430537202,0.0430823703,0.0431065228,0.0431360064,0.0431532161,0.0431809121,0.0432035321,0.0432262310,0.0432502274,0.0432683464,0.0432926830,0.0433096794,0.0433345148,0.0433512832,0.0433739912,0.0433967761,0.0434140311,0.0434304466,0.0434496548,0.0434657073,0.0434845356,0.0434997809,0.0435160341,0.0435345330,0.0435534287,0.0435697747,0.0435843967,0.0435978787,0.0436147365,0.0436328511,0.0436441830,0.0436584764,0.0436727792,0.0436893376,0.0437028427,0.0437127800,0.0437294339,0.0346249986,0.0350041406,0.0353511917,0.0356628696,0.0359530747,0.0362160876,0.0364622132,0.0366893981,0.0369014938,0.0370969572,0.0372786955,0.0374566356,0.0376117087,0.0377653516,0.0379114643,0.0380468265,0.0381746526,0.0382952276,0.0384109412,0.0385224588,0.0386300541,0.0387294381,0.0388241969,0.0389133754,0.0390005357,0.0390856271,0.0391629514,0.0392391568,0.0393159715,0.0393838186,0.0394562954,0.0395179348,0.0395809021,0.0396449672,0.0397037097,0.0397612749,0.0398127855,0.0398669607,0.0399181469,0.0399697381,0.0400128647,0.0400590998,0.0401045173,0.0401510531,0.0401932721,0.0402296209,0.0402735719,0.0403114776,0.0403493138,0.0403865418,0.0404193059,0.0404559000,0.0404866598,0.0405238398,0.0405577609,0.0405842927,0.0406163475,0.0406489910,0.0406732969,0.0407013920,0.0407288423,0.0407568731,0.0407805293,0.0408072959,0.0408350013,0.0408570572,0.0408865445,0.0409020306,0.0409285379,0.0409500433,0.0409727597,0.0409957607,0.0410137882,0.0410364411,0.0410523771,0.0410761237,0.0410920260,0.0411141991,0.0411361843,0.0411520632,0.0411686153,0.0411867278,0.0412025316,0.0412202740,0.0412339982,0.0412510774,0.0412681859,0.0412864358,0.0413022427,0.0413158987,0.0413291412,0.0413456274,0.0413631343,0.0413745146,0.0413873072,0.0414017255,0.0414167659,0.0414296591,0.0414401566,0.0414558025,0.0327727570,0.0331335454,0.0334644283,0.0337614657,0.0340383317,0.0342896035,0.0345230368,0.0347409994,0.0349425119,0.0351284571,0.0353032926,0.0354726746,0.0356220836,0.0357679203,0.0359069150,0.0360377932,0.0361591934,0.0362743367,0.0363854337,0.0364919801,0.0365940950,0.0366894773,0.0367801483,0.0368654590,0.0369503810,0.0370302591,0.0371045178,0.0371775551,0.0372510367,0.0373162635,0.0373856303,0.0374443596,0.0375047753,0.0375662712,0.0376230091,0.0376764324,0.0377268176,0.0377780653,0.0378282467,0.0378764651,0.0379182499,0.0379610710,0.0380063454,0.0380512658,0.0380905449,0.0381260781,0.0381678598,0.0382042445,0.0382407714,0.0382763724,0.0383078946,0.0383423460,0.0383725286,0.0384078304,0.0384395219,0.0384660746,0.0384971926,0.0385280100,0.0385516360,0.0385782222,0.0386052100,0.0386308823,0.0386548834,0.0386793546,0.0387066170,0.0387275042,0.0387564370,0.0387705767,0.0387958469,0.0388171307,0.0388388081,0.0388604252,0.0388773344,0.0388992835,0.0389151911,0.0389373417,0.0389533882,0.0389733692,0.0389955443,0.0390109510,0.0390263798,0.0390430090,0.0390592607,0.0390769157,0.0390896304,0.0391056477,0.0391220894,0.0391401647,0.0391545789,0.0391681986,0.0391803569,0.0391967147,0.0392127402,0.0392235975,0.0392362508,0.0392497035,0.0392653931,0.0392773637,0.0392873388,0.0393020013,0.0310218115,0.0313654026,0.0316803481,0.0319628125,0.0322275848,0.0324672112,0.0326891754,0.0328978753,0.0330901549,0.0332676630,0.0334344232,0.0335960669,0.0337388774,0.0338788356,0.0340116846,0.0341361450,0.0342525503,0.0343622285,0.0344691272,0.0345701412,0.0346675934,0.0347587685,0.0348467107,0.0349274426,0.0350084682,0.0350856133,0.0351557033,0.0352260958,0.0352956081,0.0353588115,0.0354254124,0.0354825745,0.0355395735,0.0355977442,0.0356527485,0.0357034917,0.0357519896,0.0358014871,0.0358497703,0.0358957816,0.0359365378,0.0359757082,0.0360193901,0.0360620433,0.0361000041,0.0361346278,0.0361751391,0.0362088573,0.0362438013,0.0362780097,0.0363091579,0.0363410599,0.0363702857,0.0364040672,0.0364345631,0.0364603150,0.0364901046,0.0365194490,0.0365428238,0.0365667357,0.0365927532,0.0366180899,0.0366405994,0.0366643266,0.0366905712,0.0367111245,0.0367386780,0.0367524349,0.0367771909,0.0367966626,0.0368172475,0.0368383240,0.0368543259,0.0368756517,0.0368909764,0.0369121192,0.0369277145,0.0369467607,0.0369682500,0.0369828997,0.0369978378,0.0370136437,0.0370291608,0.0370468566,0.0370587603,0.0370734024,0.0370902944,0.0371072207,0.0371202241,0.0371336383,0.0371455637,0.0371612544,0.0371767328,0.0371868266,0.0371991170,0.0372119292,0.0372270707,0.0372394352,0.0372479883,0.0372621366,0.0293666487,0.0296932172,0.0299937786,0.0302625471,0.0305147523,0.0307435504,0.0309553352,0.0311540809,0.0313374963,0.0315079132,0.0316662578,0.0318215683,0.0319569567,0.0320912195,0.0322179403,0.0323366751,0.0324478636,0.0325521325,0.0326543637,0.0327511537,0.0328447826,0.0329325220,0.0330162204,0.0330932320,0.0331707312,0.0332446434,0.0333114585,0.0333788611,0.0334451763,0.0335057950,0.0335699643,0.0336247844,0.0336793924,0.0337344594,0.0337871203,0.0338364185,0.0338826352,0.0339286974,0.0339761517,0.0340199077,0.0340597088,0.0340963865,0.0341387476,0.0341787149,0.0342151209,0.0342492673,0.0342876633,0.0343202306,0.0343527699,0.0343867442,0.0344158028,0.0344464653,0.0344740580,0.0345070001,0.0345369784,0.0345605545,0.0345905090,0.0346175886,0.0346402021,0.0346626282,0.0346879951,0.0347124585,0.0347341824,0.0347569962,0.0347816822,0.0348011401,0.0348274530,0.0348412993,0.0348653556,0.0348835136,0.0349040611,0.0349230084,0.0349391994,0.0349596967,0.0349731528,0.0349942223,0.0350086984,0.0350266877,0.0350482336,0.0350617559,0.0350766410,0.0350916612,0.0351059837,0.0351236651,0.0351341986,0.0351486606,0.0351659807,0.0351817629,0.0351934473,0.0352071907,0.0352191768,0.0352337293,0.0352474846,0.0352577053,0.0352696780,0.0352820964,0.0352960679,0.0353094785,0.0353176022,0.0353301588,0.0278012510,0.0281113724,0.0283986909,0.0286551334,0.0288949892,0.0291138833,0.0293146704,0.0295057542,0.0296800508,0.0298425281,0.0299939785,0.0301420191,0.0302708767,0.0303991053,0.0305208590,0.0306342436,0.0307397818,0.0308398976,0.0309379171,0.0310308015,0.0311200234,0.0312036635,0.0312825545,0.0313572905,0.0314312666,0.0315019629,0.0315660585,0.0316298874,0.0316942673,0.0317516228,0.0318136796,0.0318653345,0.0319178239,0.0319711417,0.0320198856,0.0320679246,0.0321118974,0.0321568364,0.0322018507,0.0322433481,0.0322821540,0.0323167117,0.0323574109,0.0323958703,0.0324309158,0.0324633801,0.0324999588,0.0325315875,0.0325630285,0.0325952690,0.0326229379,0.0326530338,0.0326787921,0.0327105332,0.0327393541,0.0327617753,0.0327907104,0.0328165041,0.0328380275,0.0328601757,0.0328840740,0.0329081283,0.0329280782,0.0329497726,0.0329737483,0.0329920659,0.0330182857,0.0330313859,0.0330540117,0.0330714131,0.0330912704,0.0331091973,0.0331239552,0.0331442532,0.0331579587,0.0331783175,0.0331914654,0.0332081921,0.0332294245,0.0332421742,0.0332566606,0.0332706892,0.0332847810,0.0333007562,0.0333116747,0.0333248420,0.0333431535,0.0333574403,0.0333687483,0.0333820297,0.0333930915,0.0334075704,0.0334208609,0.0334308818,0.0334419408,0.0334539694,0.0334673208,0.0334797793,0.0334876740,0.0334995204,0.0263204560,0.0266160824,0.0268901452,0.0271342716,0.0273630538,0.0275712108,0.0277639942,0.0279455905,0.0281115086,0.0282663019,0.0284116079,0.0285521413,0.0286750825,0.0287988675,0.0289139313,0.0290225210,0.0291233663,0.0292192220,0.0293124720,0.0294016300,0.0294870865,0.0295666937,0.0296417609,0.0297136282,0.0297845125,0.0298521465,0.0299135092,0.0299741673,0.0300357362,0.0300912392,0.0301506505,0.0301991141,0.0302503294,0.0303000163,0.0303475787,0.0303931696,0.0304354678,0.0304784074,0.0305224250,0.0305614746,0.0305985514,0.0306322226,0.0306712684,0.0307079146,0.0307408405,0.0307721089,0.0308066801,0.0308376329,0.0308682867,0.0308984276,0.0309253517,0.0309543606,0.0309778171,0.0310091021,0.0310364986,0.0310583272,0.0310855753,0.0311107436,0.0311311307,0.0311526710,0.0311750054,0.0311978402,0.0312179370,0.0312388595,0.0312615261,0.0312791371,0.0313041295,0.0313169689,0.0313379423,0.0313553120,0.0313741263,0.0313901868,0.0314052171,0.0314251828,0.0314375052,0.0314569634,0.0314699980,0.0314857327,0.0315058190,0.0315180081,0.0315320645,0.0315459966,0.0315594932,0.0315749405,0.0315853223,0.0315980897,0.0316159241,0.0316292504,0.0316404753,0.0316531743,0.0316633006,0.0316780191,0.0316905304,0.0316995122,0.0317110271,0.0317218764,0.0317349624,0.0317471062,0.0317539795,0.0317657075,0.0249208783,0.0252017278,0.0254628445,0.0256960765,0.0259135344,0.0261126037,0.0262958686,0.0264693084,0.0266267639,0.0267752066,0.0269143187,0.0270482582,0.0271651056,0.0272832139,0.0273932871,0.0274971868,0.0275940448,0.0276848649,0.0277740918,0.0278595577,0.0279404521,0.0280176467,0.0280891014,0.0281578420,0.0282259266,0.0282904016,0.0283485658,0.0284068214,0.0284655838,0.0285183984,0.0285759405,0.0286219438,0.0286714937,0.0287185065,0.0287639793,0.0288074213,0.0288479738,0.0288888887,0.0289314280,0.0289685802,0.0290040423,0.0290364363,0.0290739706,0.0291082603,0.0291404295,0.0291700636,0.0292038978,0.0292338943,0.0292618630,0.0292919997,0.0293175025,0.0293453211,0.0293672607,0.0293975755,0.0294248370,0.0294450575,0.0294702898,0.0294945657,0.0295145784,0.0295345885,0.0295566726,0.0295784987,0.0295974371,0.0296180025,0.0296387328,0.0296564182,0.0296807816,0.0296928239,0.0297131540,0.0297293038,0.0297475361,0.0297632354,0.0297783177,0.0297965006,0.0298085671,0.0298273809,0.0298389664,0.0298534623,0.0298737052,0.0298852176,0.0298989668,0.0299116750,0.0299248155,0.0299393146,0.0299498031,0.0299626529,0.0299797715,0.0299918644,0.0300026103,0.0300156569,0.0300252522,0.0300387949,0.0300509983,0.0300593946,0.0300703995,0.0300805393,0.0300933578,0.0301051417,0.0301116445,0.0301231891,0.0235967741,0.0238642918,0.0241125970,0.0243347847,0.0245424220,0.0247330597,0.0249076812,0.0250726843,0.0252229098,0.0253640207,0.0254974277,0.0256253454,0.0257364391,0.0258488224,0.0259547951,0.0260536778,0.0261465993,0.0262327107,0.0263174805,0.0263989537,0.0264772001,0.0265500026,0.0266185128,0.0266844502,0.0267499525,0.0268111717,0.0268670218,0.0269227415,0.0269790022,0.0270294192,0.0270839216,0.0271284483,0.0271759180,0.0272207173,0.0272645339,0.0273058993,0.0273446886,0.0273828858,0.0274250464,0.0274604073,0.0274947923,0.0275251978,0.0275610748,0.0275943153,0.0276243607,0.0276534962,0.0276856175,0.0277140183,0.0277421872,0.0277691293,0.0277947318,0.0278209503,0.0278415678,0.0278714572,0.0278969720,0.0279163144,0.0279411580,0.0279641004,0.0279831124,0.0280016464,0.0280230870,0.0280445153,0.0280630250,0.0280822282,0.0281025372,0.0281194614,0.0281429695,0.0281536994,0.0281732272,0.0281888386,0.0282063516,0.0282217169,0.0282358720,0.0282538646,0.0282649179,0.0282831904,0.0282938151,0.0283082958,0.0283265716,0.0283377777,0.0283520738,0.0283637366,0.0283766791,0.0283905552,0.0284005818,0.0284130533,0.0284286501,0.0284405564,0.0284514838,0.0284634783,0.0284725658,0.0284864977,0.0284974976,0.0285053143,0.0285164419,0.0285251185,0.0285386659,0.0285490960,0.0285561969,0.0285661838,0.0223444518,0.0225995052,0.0228355958,0.0230473500,0.0232452709,0.0234269572,0.0235932858,0.0237507592,0.0238936878,0.0240282401,0.0241557584,0.0242783684,0.0243845296,0.0244916177,0.0245924864,0.0246873689,0.0247754631,0.0248577962,0.0249385513,0.0250171136,0.0250920553,0.0251604135,0.0252263267,0.0252896718,0.0253522147,0.0254109017,0.0254642246,0.0255176783,0.0255713909,0.0256194782,0.0256717979,0.0257145470,0.0257598028,0.0258021804,0.0258442777,0.0258838859,0.0259209836,0.0259573275,0.0259978684,0.0260314116,0.0260645479,0.0260939368,0.0261278821,0.0261594431,0.0261891686,0.0262169569,0.0262480143,0.0262749355,0.0263011974,0.0263271833,0.0263518471,0.0263768942,0.0263967844,0.0264254050,0.0264497462,0.0264684139,0.0264918486,0.0265141639,0.0265322558,0.0265506036,0.0265704007,0.0265913482,0.0266093796,0.0266272273,0.0266469719,0.0266628077,0.0266847953,0.0266958517,0.0267139776,0.0267294742,0.0267457059,0.0267609000,0.0267747869,0.0267917673,0.0268025709,0.0268195371,0.0268306962,0.0268438269,0.0268613958,0.0268719273,0.0268864402,0.0268966042,0.0269091114,0.0269224095,0.0269320738,0.0269446222,0.0269595806,0.0269711906,0.0269812037,0.0269920722,0.0270014151,0.0270146212,0.0270245010,0.0270326756,0.0270433206,0.0270518652,0.0270655319,0.0270745314,0.0270814459,0.0270911613,0.0211596005,0.0214030639,0.0216278443,0.0218292777,0.0220180340,0.0221910614,0.0223500968,0.0224995830,0.0226362125,0.0227649438,0.0228859584,0.0230030721,0.0231046955,0.0232062461,0.0233025482,0.0233935238,0.0234778565,0.0235562376,0.0236325965,0.0237079416,0.0237800901,0.0238459233,0.0239085737,0.0239684782,0.0240286251,0.0240850917,0.0241357430,0.0241872733,0.0242382148,0.0242849843,0.0243342890,0.0243752790,0.0244185743,0.0244586574,0.0244990238,0.0245367403,0.0245725522,0.0246067488,0.0246458921,0.0246783013,0.0247094535,0.0247388816,0.0247703574,0.0248012250,0.0248298991,0.0248557025,0.0248850066,0.0249112673,0.0249366303,0.0249614194,0.0249851050,0.0250087513,0.0250284294,0.0250555157,0.0250790782,0.0250967486,0.0251194473,0.0251396041,0.0251573994,0.0251749473,0.0251942133,0.0252143772,0.0252319618,0.0252487967,0.0252681301,0.0252826066,0.0253034032,0.0253140204,0.0253318905,0.0253467258,0.0253623279,0.0253768182,0.0253900868,0.0254062384,0.0254167108,0.0254333141,0.0254436961,0.0254561831,0.0254733271,0.0254835410,0.0254969772,0.0255076393,0.0255192085,0.0255312659,0.0255412530,0.0255530636,0.0255674791,0.0255774351,0.0255878229,0.0255985020,0.0256080975,0.0256199484,0.0256294403,0.0256371849,0.0256480097,0.0256557954,0.0256689059,0.0256771865,0.0256839986,0.0256937378,0.0200389766,0.0202703642,0.0204846289,0.0206767205,0.0208564391,0.0210214562,0.0211729546,0.0213154995,0.0214452464,0.0215681602,0.0216843136,0.0217959575,0.0218936239,0.0219896067,0.0220818184,0.0221685213,0.0222492340,0.0223233953,0.0223967199,0.0224694198,0.0225377936,0.0226006883,0.0226605581,0.0227177013,0.0227754331,0.0228298337,0.0228777857,0.0229264028,0.0229758316,0.0230201632,0.0230676983,0.0231069368,0.0231482083,0.0231859210,0.0232257569,0.0232617464,0.0232953534,0.0233281049,0.0233659901,0.0233961987,0.0234258696,0.0234546298,0.0234842524,0.0235141332,0.0235418868,0.0235660610,0.0235947423,0.0236194006,0.0236431956,0.0236671959,0.0236900335,0.0237134234,0.0237310316,0.0237573512,0.0237801123,0.0237975256,0.0238187202,0.0238383157,0.0238554758,0.0238720254,0.0238896839,0.0239096815,0.0239264556,0.0239431051,0.0239606986,0.0239752410,0.0239951828,0.0240051478,0.0240223968,0.0240372938,0.0240513353,0.0240653238,0.0240783068,0.0240937552,0.0241040933,0.0241192131,0.0241297742,0.0241410341,0.0241577929,0.0241679960,0.0241806797,0.0241910442,0.0242022765,0.0242137964,0.0242230637,0.0242345134,0.0242476129,0.0242572536,0.0242673138,0.0242770046,0.0242866383,0.0242987877,0.0243074094,0.0243147465,0.0243249619,0.0243330422,0.0243455672,0.0243535458,0.0243601724,0.0243688240,0.0189787334,0.0191991261,0.0194029871,0.0195863530,0.0197569881,0.0199147648,0.0200588645,0.0201946741,0.0203189737,0.0204356721,0.0205465129,0.0206528584,0.0207466568,0.0208373249,0.0209263836,0.0210086660,0.0210853326,0.0211560672,0.0212266616,0.0212956720,0.0213611650,0.0214214625,0.0214783630,0.0215336961,0.0215877200,0.0216403792,0.0216858435,0.0217331278,0.0217798181,0.0218219367,0.0218672543,0.0219056656,0.0219442124,0.0219803943,0.0220187664,0.0220531374,0.0220850803,0.0221168164,0.0221533029,0.0221818868,0.0222109772,0.0222384919,0.0222666786,0.0222950517,0.0223218235,0.0223446878,0.0223719509,0.0223957519,0.0224179325,0.0224411397,0.0224630616,0.0224854270,0.0225022277,0.0225277652,0.0225493564,0.0225663189,0.0225860820,0.0226050823,0.0226214238,0.0226370184,0.0226542252,0.0226735179,0.0226887744,0.0227054227,0.0227230880,0.0227367501,0.0227552303,0.0227651535,0.0227808963,0.0227956906,0.0228091686,0.0228223643,0.0228349811,0.0228500827,0.0228594635,0.0228742688,0.0228837086,0.0228952203,0.0229109811,0.0229212159,0.0229326589,0.0229430899,0.0229542061,0.0229649381,0.0229738555,0.0229845327,0.0229972482,0.0230061939,0.0230163177,0.0230252335,0.0230350637,0.0230463685,0.0230546956,0.0230618555,0.0230707036,0.0230795096,0.0230908020,0.0230986325,0.0231046578,0.0231133601,0.0179751405,0.0181853499,0.0183793823,0.0185546545,0.0187166444,0.0188668529,0.0190044804,0.0191344228,0.0192526675,0.0193635451,0.0194692817,0.0195705656,0.0196607428,0.0197471613,0.0198322207,0.0199109245,0.0199832936,0.0200514804,0.0201189804,0.0201845600,0.0202467040,0.0203042057,0.0203591168,0.0204119358,0.0204635199,0.0205141243,0.0205572550,0.0206023584,0.0206469975,0.0206869683,0.0207306508,0.0207674314,0.0208043294,0.0208387461,0.0208764324,0.0209081999,0.0209391252,0.0209692386,0.0210043752,0.0210308689,0.0210594504,0.0210860985,0.0211127627,0.0211395397,0.0211652659,0.0211873993,0.0212136369,0.0212358790,0.0212572908,0.0212796977,0.0213001644,0.0213216066,0.0213381633,0.0213622628,0.0213830236,0.0213996702,0.0214184971,0.0214358964,0.0214525167,0.0214674246,0.0214839020,0.0215024423,0.0215160295,0.0215330612,0.0215495184,0.0215633656,0.0215803309,0.0215897108,0.0216058091,0.0216185835,0.0216320242,0.0216442057,0.0216566779,0.0216714373,0.0216805925,0.0216942101,0.0217033479,0.0217148902,0.0217298271,0.0217392729,0.0217495507,0.0217606094,0.0217708154,0.0217814346,0.0217895080,0.0218003235,0.0218127607,0.0218205282,0.0218300900,0.0218387789,0.0218489189,0.0218591268,0.0218676129,0.0218742590,0.0218828457,0.0218913044,0.0219016410,0.0219093833,0.0219146747,0.0219235441,0.0170253158,0.0172262190,0.0174108050,0.0175773164,0.0177316001,0.0178754909,0.0180061678,0.0181305506,0.0182432639,0.0183486370,0.0184490205,0.0185461391,0.0186317316,0.0187148410,0.0187961716,0.0188704398,0.0189401470,0.0190055209,0.0190694898,0.0191320948,0.0191917316,0.0192465565,0.0192988263,0.0193493157,0.0193987182,0.0194469756,0.0194886894,0.0195306932,0.0195744914,0.0196126175,0.0196540175,0.0196885085,0.0197242288,0.0197576277,0.0197930286,0.0198238714,0.0198530468,0.0198821093,0.0199156823,0.0199411570,0.0199686114,0.0199935561,0.0200190055,0.0200442787,0.0200696855,0.0200907235,0.0201161817,0.0201370237,0.0201577581,0.0201794717,0.0201991177,0.0202190197,0.0202353872,0.0202584808,0.0202780484,0.0202933093,0.0203121254,0.0203288130,0.0203447616,0.0203586371,0.0203747439,0.0203925233,0.0204053945,0.0204216553,0.0204370036,0.0204509731,0.0204673230,0.0204758051,0.0204914815,0.0205031170,0.0205160602,0.0205285810,0.0205399108,0.0205540349,0.0205625573,0.0205752278,0.0205847136,0.0205959687,0.0206104410,0.0206193936,0.0206286455,0.0206396718,0.0206493289,0.0206589890,0.0206670877,0.0206770529,0.0206895717,0.0206970174,0.0207061278,0.0207146806,0.0207240863,0.0207344141,0.0207422541,0.0207486735,0.0207561810,0.0207645907,0.0207744811,0.0207818788,0.0207872679,0.0207962893,0.0161272139,0.0163179618,0.0164943462,0.0166524099,0.0167993690,0.0169370198,0.0170609141,0.0171796764,0.0172869248,0.0173874738,0.0174835299,0.0175759193,0.0176573242,0.0177371799,0.0178145721,0.0178852964,0.0179518988,0.0180143879,0.0180753350,0.0181352448,0.0181924254,0.0182447540,0.0182946345,0.0183427105,0.0183889546,0.0184361078,0.0184759825,0.0185159653,0.0185574561,0.0185937639,0.0186332463,0.0186668440,0.0187009516,0.0187329718,0.0187669580,0.0187964165,0.0188239987,0.0188519597,0.0188838688,0.0189080684,0.0189338971,0.0189587969,0.0189831564,0.0190065594,0.0190317266,0.0190519120,0.0190759048,0.0190959846,0.0191154402,0.0191368157,0.0191551994,0.0191738658,0.0191894341,0.0192121996,0.0192307041,0.0192452232,0.0192630150,0.0192792762,0.0192947267,0.0193077529,0.0193233206,0.0193400256,0.0193525056,0.0193679597,0.0193830803,0.0193962984,0.0194117862,0.0194197896,0.0194349603,0.0194466751,0.0194585639,0.0194697880,0.0194817063,0.0194946050,0.0195033192,0.0195153024,0.0195246698,0.0195356693,0.0195491530,0.0195574932,0.0195664604,0.0195770398,0.0195863271,0.0195955627,0.0196032535,0.0196125253,0.0196250496,0.0196322714,0.0196401403,0.0196488076,0.0196584043,0.0196680636,0.0196756083,0.0196818683,0.0196880968,0.0196969476,0.0197065799,0.0197132294,0.0197186075,0.0197273521,0.0152770046,0.0154582766,0.0156269276,0.0157769332,0.0159170535,0.0160482878,0.0161660459,0.0162795369,0.0163813656,0.0164776004,0.0165694496,0.0166575845,0.0167349245,0.0168109026,0.0168851280,0.0169522737,0.0170160447,0.0170759360,0.0171337950,0.0171910583,0.0172453188,0.0172956618,0.0173433525,0.0173894364,0.0174334097,0.0174781495,0.0175162118,0.0175544507,0.0175946674,0.0176287510,0.0176667799,0.0176990583,0.0177314168,0.0177623957,0.0177938747,0.0178228749,0.0178492585,0.0178757429,0.0179061320,0.0179295969,0.0179540161,0.0179788613,0.0180018085,0.0180234130,0.0180476967,0.0180675187,0.0180902665,0.0181093976,0.0181284962,0.0181482323,0.0181659104,0.0181837785,0.0181989544,0.0182207599,0.0182383863,0.0182519959,0.0182691272,0.0182847694,0.0182993082,0.0183115974,0.0183263466,0.0183428997,0.0183552470,0.0183694304,0.0183843079,0.0183967686,0.0184116339,0.0184190069,0.0184335342,0.0184441972,0.0184565453,0.0184668089,0.0184789396,0.0184904554,0.0184987130,0.0185104973,0.0185197662,0.0185303050,0.0185426454,0.0185509282,0.0185596689,0.0185698571,0.0185787473,0.0185876015,0.0185940864,0.0186036343,0.0186162824,0.0186229336,0.0186295218,0.0186385031,0.0186477181,0.0186567253,0.0186643965,0.0186709096,0.0186759302,0.0186846706,0.0186934965,0.0187003879,0.0187057119,0.0187129835,0.0144722316,0.0146449869,0.0148058789,0.0149486972,0.0150821829,0.0152066697,0.0153191653,0.0154271178,0.0155240889,0.0156156821,0.0157035286,0.0157872501,0.0158616078,0.0159338065,0.0160047666,0.0160687553,0.0161300933,0.0161871445,0.0162417870,0.0162967138,0.0163488386,0.0163966282,0.0164421192,0.0164865669,0.0165280859,0.0165709030,0.0166070618,0.0166439948,0.0166823619,0.0167146787,0.0167512621,0.0167814558,0.0168130843,0.0168428652,0.0168728983,0.0169001989,0.0169256969,0.0169507354,0.0169801116,0.0170025004,0.0170260581,0.0170494383,0.0170713730,0.0170918113,0.0171157163,0.0171340722,0.0171556601,0.0171739750,0.0171926844,0.0172116856,0.0172284744,0.0172451077,0.0172595445,0.0172808961,0.0172973763,0.0173104731,0.0173274833,0.0173418700,0.0173564739,0.0173685856,0.0173822783,0.0173977134,0.0174102105,0.0174232399,0.0174372303,0.0174491446,0.0174637428,0.0174707203,0.0174847702,0.0174946268,0.0175061901,0.0175167747,0.0175276814,0.0175389274,0.0175472849,0.0175582147,0.0175675235,0.0175768092,0.0175892192,0.0175966793,0.0176057205,0.0176149860,0.0176235397,0.0176317670,0.0176385441,0.0176479111,0.0176596093,0.0176657878,0.0176723156,0.0176809657,0.0176893129,0.0176978917,0.0177054941,0.0177118454,0.0177168899,0.0177249302,0.0177338115,0.0177401070,0.0177448313,0.0177520971,0.0137108155,0.0138753318,0.0140281131,0.0141635683,0.0142911080,0.0144100877,0.0145168337,0.0146197784,0.0147124274,0.0147991707,0.0148829552,0.0149630256,0.0150341756,0.0151027509,0.0151705794,0.0152321022,0.0152897091,0.0153451921,0.0153974784,0.0154496911,0.0154991097,0.0155443533,0.0155879053,0.0156310541,0.0156705716,0.0157111857,0.0157459908,0.0157813505,0.0158178630,0.0158492657,0.0158832682,0.0159125290,0.0159424990,0.0159708256,0.0159999008,0.0160260391,0.0160502949,0.0160740624,0.0161021367,0.0161233487,0.0161460932,0.0161682874,0.0161889693,0.0162090641,0.0162318647,0.0162492377,0.0162699258,0.0162873391,0.0163057739,0.0163238267,0.0163391507,0.0163563399,0.0163693797,0.0163896455,0.0164056491,0.0164183825,0.0164343944,0.0164482241,0.0164621534,0.0164736960,0.0164873315,0.0165016839,0.0165133470,0.0165267624,0.0165395606,0.0165510357,0.0165647257,0.0165708151,0.0165848942,0.0165947478,0.0166061272,0.0166156796,0.0166263509,0.0166372518,0.0166448437,0.0166555308,0.0166643962,0.0166732438,0.0166851024,0.0166926989,0.0167007773,0.0167102538,0.0167179464,0.0167262005,0.0167326981,0.0167414133,0.0167522990,0.0167586687,0.0167653926,0.0167724616,0.0167812354,0.0167894506,0.0167965513,0.0168025360,0.0168074978,0.0168152191,0.0168237464,0.0168296361,0.0168340658,0.0168412446,0.0129900717,0.0131463018,0.0132924517,0.0134206960,0.0135425757,0.0136551962,0.0137571499,0.0138549545,0.0139439342,0.0140264426,0.0141063558,0.0141822590,0.0142508497,0.0143157150,0.0143808714,0.0144398044,0.0144947802,0.0145476258,0.0145968887,0.0146474068,0.0146937255,0.0147372119,0.0147791609,0.0148196369,0.0148578348,0.0148972398,0.0149297244,0.0149636103,0.0149984701,0.0150288432,0.0150604354,0.0150893102,0.0151177528,0.0151449317,0.0151725406,0.0151974554,0.0152202778,0.0152439104,0.0152701216,0.0152909132,0.0153125099,0.0153339643,0.0153535979,0.0153725240,0.0153944060,0.0154111880,0.0154303375,0.0154478426,0.0154650186,0.0154824262,0.0154969182,0.0155137335,0.0155253424,0.0155453135,0.0155600221,0.0155732090,0.0155880350,0.0156009968,0.0156146872,0.0156258192,0.0156388395,0.0156522966,0.0156645296,0.0156767618,0.0156891333,0.0157000825,0.0157130482,0.0157185638,0.0157321075,0.0157418168,0.0157524133,0.0157616871,0.0157718544,0.0157821438,0.0157898961,0.0157999209,0.0158086684,0.0158166412,0.0158277484,0.0158353500,0.0158432765,0.0158518401,0.0158591593,0.0158673882,0.0158738924,0.0158817798,0.0158921725,0.0158981010,0.0159053514,0.0159115081,0.0159202268,0.0159283719,0.0159347783,0.0159401721,0.0159451928,0.0159524778,0.0159610583,0.0159666500,0.0159706987,0.0159778312,0.0123068993,0.0124564624,0.0125953462,0.0127178697,0.0128332959,0.0129407950,0.0130376629,0.0131310806,0.0132163128,0.0132945636,0.0133708143,0.0134432078,0.0135093460,0.0135703353,0.0136326044,0.0136888599,0.0137414822,0.0137916315,0.0138387648,0.0138869653,0.0139312406,0.0139729614,0.0140132196,0.0140511594,0.0140875358,0.0141255766,0.0141563019,0.0141888119,0.0142217362,0.0142507387,0.0142814050,0.0143092236,0.0143358718,0.0143617871,0.0143886058,0.0144120934,0.0144344258,0.0144564461,0.0144814146,0.0145018257,0.0145222885,0.0145426358,0.0145615797,0.0145794604,0.0146006253,0.0146166277,0.0146350047,0.0146513676,0.0146679665,0.0146846248,0.0146983385,0.0147147421,0.0147258493,0.0147446035,0.0147590102,0.0147715389,0.0147860641,0.0147979887,0.0148112825,0.0148218283,0.0148343045,0.0148472960,0.0148590506,0.0148706081,0.0148828833,0.0148934298,0.0149056828,0.0149105315,0.0149233640,0.0149326203,0.0149432543,0.0149521390,0.0149620300,0.0149723216,0.0149791889,0.0149884721,0.0149966733,0.0150043965,0.0150149879,0.0150219843,0.0150298788,0.0150384532,0.0150448004,0.0150530842,0.0150596833,0.0150670534,0.0150771323,0.0150825180,0.0150893514,0.0150953440,0.0151037937,0.0151116549,0.0151168873,0.0151227472,0.0151281375,0.0151343583,0.0151430104,0.0151482555,0.0151526978,0.0151588053,0.0116607847,0.0118031374,0.0119355129,0.0120518903,0.0121616194,0.0122641537,0.0123560147,0.0124451568,0.0125264591,0.0126013841,0.0126737768,0.0127436060,0.0128063774,0.0128643161,0.0129238926,0.0129774980,0.0130278166,0.0130755193,0.0131203748,0.0131661662,0.0132086854,0.0132483966,0.0132869337,0.0133227539,0.0133582067,0.0133943844,0.0134227837,0.0134543278,0.0134863170,0.0135135670,0.0135433805,0.0135695945,0.0135952165,0.0136201173,0.0136451376,0.0136681838,0.0136888914,0.0137104717,0.0137338379,0.0137538357,0.0137732456,0.0137921930,0.0138111251,0.0138277191,0.0138478792,0.0138637507,0.0138810894,0.0138965374,0.0139127799,0.0139282312,0.0139414283,0.0139569774,0.0139680369,0.0139856710,0.0139994704,0.0140117500,0.0140259927,0.0140368208,0.0140494117,0.0140596669,0.0140717496,0.0140841002,0.0140954717,0.0141063958,0.0141179457,0.0141282694,0.0141399846,0.0141445332,0.0141568209,0.0141655843,0.0141756166,0.0141843481,0.0141943294,0.0142039272,0.0142103631,0.0142192054,0.0142274544,0.0142343857,0.0142445242,0.0142511923,0.0142589755,0.0142672414,0.0142726939,0.0142814858,0.0142873528,0.0142945104,0.0143041104,0.0143091134,0.0143160019,0.0143214485,0.0143295779,0.0143373619,0.0143417017,0.0143479476,0.0143533405,0.0143593058,0.0143671907,0.0143726896,0.0143765732,0.0143822967,0.0110494023,0.0111842103,0.0113101407,0.0114212905,0.0115259733,0.0116238708,0.0117108420,0.0117961389,0.0118730585,0.0119447498,0.0120140521,0.0120804434,0.0121404508,0.0121949246,0.0122521798,0.0123037346,0.0123518159,0.0123969152,0.0124397952,0.0124835823,0.0125241618,0.0125615622,0.0125988515,0.0126324451,0.0126670378,0.0127014480,0.0127284120,0.0127587428,0.0127894879,0.0128155907,0.0128439190,0.0128686519,0.0128928844,0.0129167809,0.0129412537,0.0129627582,0.0129825669,0.0130033653,0.0130256029,0.0130449327,0.0130630646,0.0130811905,0.0130997256,0.0131155762,0.0131343754,0.0131502526,0.0131661363,0.0131810159,0.0131964305,0.0132119308,0.0132237426,0.0132391995,0.0132498709,0.0132666420,0.0132796324,0.0132910603,0.0133049888,0.0133150014,0.0133269352,0.0133370486,0.0133484239,0.0133609266,0.0133716181,0.0133815477,0.0133927956,0.0134030042,0.0134135954,0.0134181473,0.0134302398,0.0134382899,0.0134482030,0.0134563579,0.0134660247,0.0134753489,0.0134810499,0.0134896004,0.0134976763,0.0135049472,0.0135143196,0.0135208542,0.0135281303,0.0135358270,0.0135408817,0.0135497250,0.0135547131,0.0135620243,0.0135712483,0.0135752063,0.0135825906,0.0135878805,0.0135952281,0.0136032917,0.0136071139,0.0136129803,0.0136182100,0.0136243045,0.0136313382,0.0136369614,0.0136403773,0.0136457110,0.0104701530,0.0105985051,0.0107181383,0.0108240264,0.0109240788,0.0110170024,0.0110998142,0.0111812016,0.0112543450,0.0113230619,0.0113886643,0.0114520033,0.0115094973,0.0115615507,0.0116155679,0.0116649716,0.0117108396,0.0117541964,0.0117947115,0.0118370937,0.0118750194,0.0119114034,0.0119468413,0.0119791026,0.0120118664,0.0120446477,0.0120707303,0.0120993792,0.0121284553,0.0121540785,0.0121809309,0.0122037743,0.0122278401,0.0122507189,0.0122737921,0.0122943751,0.0123133038,0.0123327154,0.0123542904,0.0123729279,0.0123900689,0.0124075789,0.0124254346,0.0124403984,0.0124583724,0.0124732202,0.0124887695,0.0125030039,0.0125181329,0.0125324502,0.0125433023,0.0125584649,0.0125685284,0.0125847398,0.0125969597,0.0126079001,0.0126211072,0.0126308469,0.0126425215,0.0126520283,0.0126633742,0.0126753092,0.0126851232,0.0126950775,0.0127054037,0.0127153149,0.0127256335,0.0127305005,0.0127412102,0.0127487951,0.0127583462,0.0127663958,0.0127757909,0.0127844613,0.0127907541,0.0127975462,0.0128059182,0.0128125104,0.0128215648,0.0128280417,0.0128346411,0.0128424489,0.0128471946,0.0128560053,0.0128602807,0.0128676185,0.0128766048,0.0128802867,0.0128869198,0.0128923072,0.0128989287,0.0129067837,0.0129108816,0.0129157442,0.0129212271,0.0129268133,0.0129339644,0.0129391348,0.0129423720,0.0129479916,0.0099217728,0.0100437583,0.0101571777,0.0102582145,0.0103540242,0.0104424405,0.0105207508,0.0105986601,0.0106683432,0.0107339898,0.0107963205,0.0108565697,0.0109113506,0.0109614827,0.0110128489,0.0110596543,0.0111031741,0.0111450928,0.0111834365,0.0112242254,0.0112601846,0.0112953873,0.0113288340,0.0113596818,0.0113909258,0.0114223096,0.0114473319,0.0114745856,0.0115028029,0.0115267059,0.0115526042,0.0115741518,0.0115974592,0.0116187537,0.0116410162,0.0116612790,0.0116790511,0.0116970000,0.0117177824,0.0117353917,0.0117522429,0.0117690670,0.0117860047,0.0118004218,0.0118172438,0.0118315859,0.0118465787,0.0118601294,0.0118747239,0.0118882873,0.0118986990,0.0119134405,0.0119228935,0.0119380890,0.0119499147,0.0119600209,0.0119729885,0.0119824442,0.0119935020,0.0120025456,0.0120129047,0.0120252326,0.0120347227,0.0120439370,0.0120538415,0.0120631852,0.0120724814,0.0120779550,0.0120879758,0.0120954347,0.0121042194,0.0121122897,0.0121208151,0.0121290592,0.0121352053,0.0121421622,0.0121500405,0.0121563618,0.0121643241,0.0121714208,0.0121771399,0.0121842825,0.0121891878,0.0121981421,0.0122018759,0.0122088195,0.0122173532,0.0122211797,0.0122275201,0.0122321915,0.0122391858,0.0122461745,0.0122506057,0.0122550035,0.0122605333,0.0122654707,0.0122726551,0.0122775037,0.0122807919,0.0122856436,0.0094020560,0.0095186917,0.0096267109,0.0097224103,0.0098138087,0.0098984719,0.0099726499,0.0100467074,0.0101132451,0.0101761267,0.0102353014,0.0102926264,0.0103446805,0.0103926762,0.0104415702,0.0104862160,0.0105278037,0.0105674672,0.0106046795,0.0106435595,0.0106779443,0.0107116510,0.0107429308,0.0107722558,0.0108015631,0.0108320459,0.0108559742,0.0108821846,0.0109090068,0.0109319915,0.0109568186,0.0109772091,0.0109997119,0.0110198414,0.0110415530,0.0110606518,0.0110773269,0.0110948189,0.0111144543,0.0111311015,0.0111475529,0.0111638104,0.0111791058,0.0111932298,0.0112096373,0.0112231938,0.0112374985,0.0112505526,0.0112648910,0.0112781390,0.0112877493,0.0113013945,0.0113107959,0.0113251980,0.0113362638,0.0113463320,0.0113584358,0.0113672974,0.0113780139,0.0113869188,0.0113968358,0.0114081085,0.0114177099,0.0114259260,0.0114349944,0.0114448630,0.0114536528,0.0114590196,0.0114686990,0.0114752314,0.0114840372,0.0114917638,0.0115001440,0.0115078487,0.0115139683,0.0115204274,0.0115279722,0.0115337355,0.0115416422,0.0115484421,0.0115538639,0.0115603054,0.0115653925,0.0115742855,0.0115773961,0.0115840647,0.0115922562,0.0115962119,0.0116020637,0.0116067239,0.0116129030,0.0116193202,0.0116237161,0.0116285712,0.0116340679,0.0116384256,0.0116457474,0.0116498939,0.0116529854,0.0116575469,0.0089098950,0.0090209078,0.0091243836,0.0092147328,0.0093020554,0.0093827488,0.0094535964,0.0095235522,0.0095873971,0.0096473270,0.0097039377,0.0097584346,0.0098077855,0.0098536270,0.0099003304,0.0099431896,0.0099825608,0.0100202753,0.0100564729,0.0100928876,0.0101261614,0.0101578429,0.0101876895,0.0102156616,0.0102433507,0.0102729851,0.0102956812,0.0103209540,0.0103461976,0.0103682272,0.0103918319,0.0104111017,0.0104330322,0.0104519005,0.0104730362,0.0104912114,0.0105073454,0.0105238207,0.0105425157,0.0105588154,0.0105737750,0.0105895716,0.0106039838,0.0106175523,0.0106339632,0.0106463799,0.0106601666,0.0106724557,0.0106864418,0.0106989190,0.0107085149,0.0107207915,0.0107299675,0.0107441394,0.0107550179,0.0107640601,0.0107759239,0.0107840742,0.0107945306,0.0108027823,0.0108120828,0.0108231972,0.0108324257,0.0108402314,0.0108490589,0.0108584655,0.0108665577,0.0108719562,0.0108813340,0.0108871014,0.0108957654,0.0109033370,0.0109114773,0.0109186260,0.0109243523,0.0109305207,0.0109378169,0.0109436571,0.0109509614,0.0109573235,0.0109628147,0.0109688854,0.0109741187,0.0109826019,0.0109850415,0.0109916263,0.0109991155,0.0110029913,0.0110090635,0.0110136200,0.0110196778,0.0110252478,0.0110296760,0.0110345175,0.0110395263,0.0110434205,0.0110509646,0.0110542695,0.0110578009,0.0110622030,0.0084441991,0.0085495447,0.0086485834,0.0087339829,0.0088170010,0.0088944588,0.0089617969,0.0090285522,0.0090890106,0.0091464777,0.0091999157,0.0092515762,0.0092990346,0.0093427655,0.0093878051,0.0094282616,0.0094658300,0.0095018455,0.0095362878,0.0095709174,0.0096025107,0.0096330497,0.0096613479,0.0096886020,0.0097144685,0.0097430097,0.0097645751,0.0097889629,0.0098128978,0.0098340656,0.0098566087,0.0098746798,0.0098955991,0.0099139749,0.0099337897,0.0099514800,0.0099668952,0.0099826988,0.0100003438,0.0100162067,0.0100302666,0.0100453818,0.0100592425,0.0100724576,0.0100873519,0.0100997452,0.0101123906,0.0101242255,0.0101381034,0.0101498279,0.0101589926,0.0101706312,0.0101797475,0.0101929160,0.0102038540,0.0102122146,0.0102235770,0.0102311542,0.0102412287,0.0102490697,0.0102581622,0.0102686290,0.0102776334,0.0102850735,0.0102934442,0.0103024786,0.0103099274,0.0103153380,0.0103245766,0.0103297870,0.0103378658,0.0103452309,0.0103527754,0.0103598422,0.0103652317,0.0103714026,0.0103783585,0.0103842309,0.0103908135,0.0103971975,0.0104023674,0.0104082983,0.0104131006,0.0104206847,0.0104236189,0.0104300040,0.0104369792,0.0104410931,0.0104467957,0.0104508080,0.0104564507,0.0104616953,0.0104657042,0.0104707557,0.0104753833,0.0104786980,0.0104864489,0.0104897440,0.0104932080,0.0104975126,0.0080029313,0.0081033177,0.0081977554,0.0082785961,0.0083579261,0.0084318123,0.0084956418,0.0085594117,0.0086169704,0.0086719564,0.0087225856,0.0087718827,0.0088168883,0.0088586425,0.0089018126,0.0089404878,0.0089758190,0.0090105013,0.0090429917,0.0090763832,0.0091066612,0.0091354419,0.0091625600,0.0091885908,0.0092131853,0.0092404350,0.0092614363,0.0092847294,0.0093072507,0.0093276452,0.0093487999,0.0093660091,0.0093861327,0.0094036935,0.0094230561,0.0094393869,0.0094544025,0.0094697057,0.0094860750,0.0095017382,0.0095151877,0.0095295787,0.0095429248,0.0095553729,0.0095697676,0.0095814244,0.0095934673,0.0096047022,0.0096180849,0.0096289402,0.0096381422,0.0096491211,0.0096577733,0.0096701667,0.0096807434,0.0096895048,0.0096992597,0.0097068959,0.0097164204,0.0097240906,0.0097326554,0.0097428181,0.0097512695,0.0097585791,0.0097667891,0.0097753400,0.0097821376,0.0097875705,0.0097963571,0.0098013148,0.0098085253,0.0098160861,0.0098236081,0.0098297941,0.0098351101,0.0098411735,0.0098480986,0.0098531721,0.0098598325,0.0098657139,0.0098705247,0.0098766193,0.0098810618,0.0098880368,0.0098909597,0.0098974535,0.0099036464,0.0099080486,0.0099131722,0.0099167894,0.0099225432,0.0099275376,0.0099314156,0.0099363140,0.0099408149,0.0099438750,0.0099516319,0.0099545125,0.0099571692,0.0099619462,0.0075848732,0.0076805538,0.0077704304,0.0078475851,0.0079230514,0.0079933231,0.0080540740,0.0081153410,0.0081697587,0.0082221564,0.0082700409,0.0083171668,0.0083605207,0.0084002220,0.0084410922,0.0084779988,0.0085115188,0.0085448899,0.0085758423,0.0086076476,0.0086365603,0.0086641707,0.0086897744,0.0087148415,0.0087382661,0.0087642294,0.0087840276,0.0088067542,0.0088278762,0.0088477137,0.0088673447,0.0088844308,0.0089028538,0.0089200818,0.0089388315,0.0089540507,0.0089683758,0.0089830297,0.0089983534,0.0090137813,0.0090269194,0.0090403040,0.0090531002,0.0090651056,0.0090783194,0.0090899174,0.0091013510,0.0091118642,0.0091249992,0.0091353536,0.0091441515,0.0091546290,0.0091626316,0.0091747926,0.0091847385,0.0091932277,0.0092027374,0.0092099745,0.0092186616,0.0092260085,0.0092344559,0.0092442966,0.0092522762,0.0092591469,0.0092671439,0.0092747661,0.0092815952,0.0092866982,0.0092953790,0.0093000831,0.0093070251,0.0093143015,0.0093217209,0.0093268701,0.0093326070,0.0093386723,0.0093445807,0.0093498167,0.0093560163,0.0093617611,0.0093663768,0.0093722040,0.0093762367,0.0093834188,0.0093857809,0.0093918979,0.0093981936,0.0094022504,0.0094073042,0.0094106374,0.0094159758,0.0094210350,0.0094245153,0.0094290267,0.0094335982,0.0094363556,0.0094439165,0.0094462429,0.0094492524,0.0094538545,0.0071891527,0.0072799518,0.0073656708,0.0074390743,0.0075113587,0.0075780610,0.0076357432,0.0076939364,0.0077458115,0.0077962767,0.0078416978,0.0078862478,0.0079274250,0.0079657015,0.0080044938,0.0080397758,0.0080720654,0.0081035316,0.0081329749,0.0081631215,0.0081913449,0.0082168902,0.0082416484,0.0082655441,0.0082879073,0.0083128951,0.0083311274,0.0083533383,0.0083731758,0.0083924864,0.0084110612,0.0084276728,0.0084449379,0.0084615718,0.0084793497,0.0084937505,0.0085079214,0.0085215901,0.0085361772,0.0085509679,0.0085638387,0.0085765319,0.0085888423,0.0086003313,0.0086127142,0.0086240811,0.0086350321,0.0086449630,0.0086572863,0.0086669244,0.0086756165,0.0086856027,0.0086931636,0.0087047905,0.0087145799,0.0087224993,0.0087319240,0.0087382249,0.0087469252,0.0087535807,0.0087618510,0.0087712508,0.0087790396,0.0087857568,0.0087928269,0.0088001118,0.0088069581,0.0088115878,0.0088207962,0.0088249716,0.0088314953,0.0088385324,0.0088454556,0.0088502951,0.0088560220,0.0088616782,0.0088672084,0.0088721774,0.0088779920,0.0088837380,0.0088881170,0.0088934539,0.0088973471,0.0089047491,0.0089065203,0.0089124383,0.0089186253,0.0089227121,0.0089273900,0.0089304301,0.0089357708,0.0089404722,0.0089437621,0.0089479153,0.0089524641,0.0089557079,0.0089621501,0.0089645948,0.0089672312,0.0089720694,0.0068144030,0.0069007751,0.0069822146,0.0070522300,0.0071212730,0.0071842169,0.0072391867,0.0072946917,0.0073443478,0.0073924347,0.0074353075,0.0074779911,0.0075176240,0.0075538025,0.0075905300,0.0076245340,0.0076551067,0.0076853430,0.0077132852,0.0077419836,0.0077692605,0.0077933620,0.0078169956,0.0078397627,0.0078612238,0.0078848835,0.0079024381,0.0079233304,0.0079422266,0.0079607964,0.0079789225,0.0079947464,0.0080110778,0.0080266139,0.0080436229,0.0080576766,0.0080710794,0.0080845030,0.0080979361,0.0081122482,0.0081245188,0.0081366979,0.0081485175,0.0081591587,0.0081708954,0.0081819804,0.0081921171,0.0082017064,0.0082133985,0.0082233862,0.0082315801,0.0082407087,0.0082479907,0.0082585400,0.0082687844,0.0082763521,0.0082853026,0.0082913778,0.0082994917,0.0083059284,0.0083137320,0.0083227900,0.0083298402,0.0083368997,0.0083433226,0.0083501905,0.0083569105,0.0083612263,0.0083698467,0.0083743242,0.0083799249,0.0083870444,0.0083938366,0.0083980772,0.0084038249,0.0084089110,0.0084146038,0.0084195254,0.0084248240,0.0084302791,0.0084342145,0.0084397491,0.0084431778,0.0084504416,0.0084522203,0.0084576535,0.0084640716,0.0084676645,0.0084724962,0.0084750049,0.0084799191,0.0084844483,0.0084873717,0.0084915126,0.0084960421,0.0084994878,0.0085052113,0.0085079125,0.0085100376,0.0085148975,0.0064596682,0.0065414535,0.0066190852,0.0066856134,0.0067516313,0.0068113327,0.0068635611,0.0069170456,0.0069639470,0.0070098754,0.0070507287,0.0070911716,0.0071292023,0.0071632202,0.0071984473,0.0072308218,0.0072597077,0.0072890866,0.0073154988,0.0073429440,0.0073686913,0.0073919708,0.0074142947,0.0074358052,0.0074566509,0.0074792296,0.0074957041,0.0075159579,0.0075337788,0.0075517415,0.0075689452,0.0075840635,0.0075999275,0.0076145365,0.0076308210,0.0076441042,0.0076571519,0.0076698409,0.0076827536,0.0076961319,0.0077079432,0.0077196254,0.0077309906,0.0077409560,0.0077523109,0.0077629488,0.0077723901,0.0077815398,0.0077929824,0.0078024183,0.0078101397,0.0078186309,0.0078257287,0.0078358451,0.0078456106,0.0078528735,0.0078611935,0.0078677225,0.0078753210,0.0078811231,0.0078891079,0.0078971271,0.0079036101,0.0079106669,0.0079171736,0.0079237898,0.0079300829,0.0079341590,0.0079423286,0.0079467359,0.0079519296,0.0079590419,0.0079653097,0.0079695248,0.0079750262,0.0079798191,0.0079850706,0.0079893820,0.0079944179,0.0080001328,0.0080041028,0.0080090503,0.0080124546,0.0080197007,0.0080212570,0.0080261845,0.0080326312,0.0080359030,0.0080404106,0.0080429036,0.0080473660,0.0080518844,0.0080545329,0.0080586986,0.0080631174,0.0080664860,0.0080720821,0.0080743055,0.0080766190,0.0080812495,0.0061230720,0.0062011938,0.0062750454,0.0063384008,0.0064014323,0.0064581928,0.0065077724,0.0065587243,0.0066033286,0.0066471161,0.0066860961,0.0067242772,0.0067607125,0.0067930896,0.0068264834,0.0068572192,0.0068848660,0.0069129683,0.0069387039,0.0069643191,0.0069888946,0.0070114321,0.0070326375,0.0070531725,0.0070727716,0.0070941082,0.0071100356,0.0071298393,0.0071466421,0.0071640161,0.0071802906,0.0071943938,0.0072096522,0.0072236639,0.0072393949,0.0072521025,0.0072646342,0.0072761678,0.0072888656,0.0073018005,0.0073133181,0.0073240268,0.0073346548,0.0073442708,0.0073555448,0.0073658413,0.0073744642,0.0073833629,0.0073935877,0.0074028812,0.0074104796,0.0074188813,0.0074256685,0.0074347071,0.0074446971,0.0074511644,0.0074591309,0.0074653652,0.0074723965,0.0074782641,0.0074860820,0.0074935838,0.0074995969,0.0075066936,0.0075126381,0.0075197031,0.0075254997,0.0075291823,0.0075369577,0.0075410285,0.0075460840,0.0075525444,0.0075589271,0.0075628686,0.0075681968,0.0075728806,0.0075777771,0.0075818833,0.0075870830,0.0075922509,0.0075963600,0.0076011316,0.0076040129,0.0076108170,0.0076123453,0.0076169988,0.0076231892,0.0076264430,0.0076307990,0.0076330411,0.0076372431,0.0076416438,0.0076441180,0.0076483059,0.0076524675,0.0076558491,0.0076608747,0.0076630046,0.0076650713,0.0076699092,0.0058043333,0.0058785946,0.0059489482,0.0060092339,0.0060690517,0.0061231087,0.0061705537,0.0062189688,0.0062615038,0.0063030855,0.0063404134,0.0063769387,0.0064117148,0.0064419876,0.0064740759,0.0065034455,0.0065300926,0.0065566034,0.0065812853,0.0066054013,0.0066289035,0.0066507753,0.0066709310,0.0066905257,0.0067092698,0.0067299339,0.0067444596,0.0067636008,0.0067793903,0.0067960769,0.0068118086,0.0068247805,0.0068399945,0.0068531711,0.0068679525,0.0068800861,0.0068921927,0.0069031833,0.0069150159,0.0069277484,0.0069387997,0.0069490290,0.0069591512,0.0069678876,0.0069792265,0.0069887558,0.0069971313,0.0070051483,0.0070153557,0.0070244850,0.0070310121,0.0070397951,0.0070460880,0.0070546307,0.0070645093,0.0070702160,0.0070778408,0.0070839228,0.0070904842,0.0070962541,0.0071038037,0.0071108219,0.0071168536,0.0071234528,0.0071289290,0.0071359006,0.0071413760,0.0071451206,0.0071524543,0.0071565353,0.0071609679,0.0071673644,0.0071732326,0.0071775326,0.0071823159,0.0071865405,0.0071915474,0.0071953496,0.0072001623,0.0072051640,0.0072092251,0.0072135878,0.0072165262,0.0072228192,0.0072245119,0.0072288216,0.0072345504,0.0072377843,0.0072421844,0.0072445350,0.0072479390,0.0072524627,0.0072546776,0.0072588276,0.0072626676,0.0072664972,0.0072708856,0.0072731035,0.0072750543,0.0072794743,0.0055024739,0.0055730471,0.0056398759,0.0056973951,0.0057543859,0.0058058331,0.0058505957,0.0058972841,0.0059376118,0.0059770247,0.0060127794,0.0060477085,0.0060806647,0.0061098885,0.0061401404,0.0061679991,0.0061932281,0.0062189272,0.0062420540,0.0062654173,0.0062878956,0.0063085741,0.0063278046,0.0063465960,0.0063646041,0.0063842477,0.0063982965,0.0064162894,0.0064313974,0.0064472219,0.0064624329,0.0064747673,0.0064892221,0.0065018248,0.0065158263,0.0065278200,0.0065388201,0.0065491987,0.0065604964,0.0065728112,0.0065833323,0.0065933514,0.0066030690,0.0066112260,0.0066219428,0.0066310438,0.0066392663,0.0066469445,0.0066567545,0.0066652161,0.0066716585,0.0066800641,0.0066858377,0.0066939705,0.0067034211,0.0067092194,0.0067162239,0.0067221013,0.0067287212,0.0067338115,0.0067413647,0.0067479116,0.0067535458,0.0067601706,0.0067650481,0.0067718884,0.0067769443,0.0067808283,0.0067878571,0.0067914443,0.0067955918,0.0068018741,0.0068076635,0.0068115220,0.0068161382,0.0068201125,0.0068251472,0.0068285232,0.0068332138,0.0068380212,0.0068421852,0.0068461053,0.0068491628,0.0068549137,0.0068565307,0.0068609749,0.0068659795,0.0068694613,0.0068733739,0.0068758729,0.0068792636,0.0068832434,0.0068854171,0.0068893473,0.0068927900,0.0068965453,0.0069008992,0.0069030477,0.0069048825,0.0069091976,0.0052166375,0.0052833413,0.0053468878,0.0054019178,0.0054561030,0.0055050633,0.0055476834,0.0055925635,0.0056306077,0.0056682006,0.0057020584,0.0057354742,0.0057668673,0.0057947082,0.0058237383,0.0058499366,0.0058743322,0.0058987048,0.0059207627,0.0059434130,0.0059642931,0.0059842714,0.0060026200,0.0060199665,0.0060378283,0.0060562945,0.0060696823,0.0060869686,0.0061012743,0.0061164297,0.0061309174,0.0061425122,0.0061568799,0.0061682392,0.0061817804,0.0061931561,0.0062038269,0.0062136955,0.0062245232,0.0062365067,0.0062461802,0.0062558108,0.0062654696,0.0062730333,0.0062832326,0.0062919642,0.0063000049,0.0063071529,0.0063163732,0.0063248341,0.0063310297,0.0063390824,0.0063441139,0.0063516875,0.0063608958,0.0063668762,0.0063730645,0.0063792071,0.0063852404,0.0063903327,0.0063974024,0.0064035361,0.0064090221,0.0064153198,0.0064196891,0.0064263658,0.0064313238,0.0064351738,0.0064417509,0.0064454223,0.0064493406,0.0064550349,0.0064611352,0.0064642793,0.0064688200,0.0064726665,0.0064774993,0.0064807743,0.0064852652,0.0064895830,0.0064938167,0.0064971899,0.0065004484,0.0065060227,0.0065076846,0.0065119176,0.0065167162,0.0065197039,0.0065238167,0.0065257096,0.0065295174,0.0065329684,0.0065350298,0.0065388929,0.0065421616,0.0065459637,0.0065499181,0.0065517666,0.0065538140,0.0065580773,0.0049452658,0.0050091427,0.0050696593,0.0051216553,0.0051734865,0.0052203376,0.0052604246,0.0053038800,0.0053397121,0.0053751919,0.0054078496,0.0054393213,0.0054694410,0.0054956986,0.0055235838,0.0055483124,0.0055716304,0.0055950908,0.0056161777,0.0056377290,0.0056576700,0.0056766325,0.0056939577,0.0057106878,0.0057275727,0.0057451494,0.0057582145,0.0057745548,0.0057882470,0.0058028280,0.0058166084,0.0058276880,0.0058413891,0.0058521865,0.0058647690,0.0058762060,0.0058861127,0.0058958167,0.0059058755,0.0059175008,0.0059265568,0.0059357787,0.0059449253,0.0059523556,0.0059619958,0.0059706703,0.0059777870,0.0059850281,0.0059935049,0.0060017837,0.0060076697,0.0060152297,0.0060201938,0.0060270729,0.0060361257,0.0060420684,0.0060474146,0.0060533925,0.0060594048,0.0060643478,0.0060710526,0.0060767915,0.0060818018,0.0060881428,0.0060922037,0.0060985727,0.0061035271,0.0061070015,0.0061131471,0.0061169413,0.0061205984,0.0061260906,0.0061320065,0.0061349920,0.0061392171,0.0061429988,0.0061476079,0.0061507926,0.0061550966,0.0061591526,0.0061629802,0.0061663429,0.0061698394,0.0061750075,0.0061762806,0.0061806717,0.0061851257,0.0061882666,0.0061921716,0.0061936431,0.0061973322,0.0062008550,0.0062028534,0.0062064190,0.0062092993,0.0062131344,0.0062167539,0.0062187324,0.0062206555,0.0062245508,0.0046885819,0.0047492988,0.0048067556,0.0048564103,0.0049058858,0.0049501841,0.0049884232,0.0050299272,0.0050637086,0.0050974592,0.0051285825,0.0051588958,0.0051870556,0.0052124772,0.0052391985,0.0052627375,0.0052848861,0.0053073528,0.0053273916,0.0053479083,0.0053668577,0.0053850115,0.0054013944,0.0054174727,0.0054338193,0.0054503105,0.0054626609,0.0054783232,0.0054914847,0.0055053884,0.0055186249,0.0055292720,0.0055418946,0.0055526210,0.0055642423,0.0055754744,0.0055847776,0.0055940944,0.0056037018,0.0056148833,0.0056237244,0.0056321444,0.0056408690,0.0056475819,0.0056573424,0.0056657027,0.0056726079,0.0056792661,0.0056874330,0.0056954814,0.0057009803,0.0057078904,0.0057127148,0.0057195212,0.0057279848,0.0057336639,0.0057389020,0.0057443532,0.0057506263,0.0057553410,0.0057610858,0.0057670726,0.0057716756,0.0057777462,0.0057816526,0.0057877421,0.0057923660,0.0057958187,0.0058018753,0.0058053613,0.0058090153,0.0058138395,0.0058197265,0.0058227017,0.0058265770,0.0058306472,0.0058348430,0.0058376113,0.0058420114,0.0058454562,0.0058495637,0.0058527654,0.0058559557,0.0058610958,0.0058618710,0.0058660600,0.0058706934,0.0058739713,0.0058773143,0.0058786694,0.0058821882,0.0058854070,0.0058877522,0.0058911980,0.0058935375,0.0058972461,0.0059006378,0.0059026216,0.0059044301,0.0059078877,0.0044451536,0.0045027534,0.0045578487,0.0046047866,0.0046521685,0.0046940810,0.0047309636,0.0047701027,0.0048022975,0.0048342478,0.0048640984,0.0048930690,0.0049196843,0.0049440363,0.0049692572,0.0049920807,0.0050127933,0.0050341799,0.0050535573,0.0050729906,0.0050911326,0.0051085302,0.0051242105,0.0051395288,0.0051550239,0.0051708645,0.0051825149,0.0051974629,0.0052098745,0.0052231781,0.0052360451,0.0052456984,0.0052582068,0.0052684486,0.0052794333,0.0052899643,0.0052988690,0.0053079953,0.0053170465,0.0053278660,0.0053362699,0.0053440623,0.0053525677,0.0053590904,0.0053682764,0.0053762257,0.0053824515,0.0053889167,0.0053972075,0.0054047729,0.0054101291,0.0054165352,0.0054213084,0.0054279541,0.0054358386,0.0054413980,0.0054466417,0.0054515866,0.0054571357,0.0054622282,0.0054670505,0.0054731885,0.0054774485,0.0054834872,0.0054872964,0.0054930069,0.0054973028,0.0055004714,0.0055065264,0.0055098130,0.0055134122,0.0055178612,0.0055236606,0.0055263836,0.0055300234,0.0055338527,0.0055376227,0.0055404189,0.0055448693,0.0055481470,0.0055516739,0.0055548799,0.0055582464,0.0055631281,0.0055640423,0.0055675958,0.0055722324,0.0055751216,0.0055784359,0.0055799032,0.0055833536,0.0055864891,0.0055888913,0.0055919267,0.0055939900,0.0055976093,0.0056007968,0.0056027154,0.0056046095,0.0056075334,0.0042145238,0.0042696287,0.0043220348,0.0043663834,0.0044113798,0.0044514155,0.0044863936,0.0045238288,0.0045545548,0.0045849539,0.0046134363,0.0046408628,0.0046663156,0.0046896400,0.0047133868,0.0047351099,0.0047549144,0.0047757024,0.0047938932,0.0048123450,0.0048296052,0.0048459116,0.0048610108,0.0048756555,0.0048905533,0.0049056135,0.0049168325,0.0049312674,0.0049429855,0.0049556902,0.0049679698,0.0049769223,0.0049888790,0.0049988025,0.0050094846,0.0050191829,0.0050277497,0.0050363966,0.0050451696,0.0050553434,0.0050636819,0.0050711922,0.0050792166,0.0050854272,0.0050942705,0.0051019368,0.0051078272,0.0051138678,0.0051218535,0.0051289183,0.0051341083,0.0051406470,0.0051450831,0.0051512190,0.0051583792,0.0051638654,0.0051690670,0.0051735175,0.0051791619,0.0051839231,0.0051885440,0.0051943602,0.0051983182,0.0052043011,0.0052079213,0.0052129377,0.0052176843,0.0052205055,0.0052259387,0.0052294953,0.0052329614,0.0052368693,0.0052426540,0.0052451092,0.0052485302,0.0052523094,0.0052556948,0.0052584977,0.0052626441,0.0052660732,0.0052693018,0.0052726590,0.0052757684,0.0052799936,0.0052812145,0.0052848777,0.0052892264,0.0052919036,0.0052949931,0.0052960405,0.0052999117,0.0053028166,0.0053053183,0.0053080540,0.0053097027,0.0053130867,0.0053162427,0.0053181244,0.0053201682,0.0053228425,0.0039959262,0.0040483941,0.0040983980,0.0041402520,0.0041834331,0.0042213417,0.0042549960,0.0042904484,0.0043196354,0.0043485504,0.0043757560,0.0044017094,0.0044258680,0.0044484549,0.0044708413,0.0044917491,0.0045103844,0.0045305699,0.0045475529,0.0045651372,0.0045817570,0.0045973055,0.0046116433,0.0046255374,0.0046395660,0.0046541122,0.0046648289,0.0046787603,0.0046897502,0.0047020418,0.0047136470,0.0047223203,0.0047338133,0.0047430686,0.0047533094,0.0047624406,0.0047705954,0.0047787500,0.0047872830,0.0047972239,0.0048050457,0.0048122970,0.0048198118,0.0048258873,0.0048343650,0.0048415392,0.0048472475,0.0048530397,0.0048604134,0.0048675464,0.0048721792,0.0048784162,0.0048829215,0.0048888822,0.0048955966,0.0049005113,0.0049056903,0.0049097044,0.0049150651,0.0049195529,0.0049244010,0.0049299685,0.0049333059,0.0049394144,0.0049426806,0.0049475981,0.0049518641,0.0049548437,0.0049601731,0.0049633910,0.0049669538,0.0049704872,0.0049758069,0.0049784169,0.0049819075,0.0049853425,0.0049883043,0.0049908198,0.0049947260,0.0049983170,0.0050017190,0.0050046758,0.0050075323,0.0050117361,0.0050128424,0.0050161730,0.0050202084,0.0050231837,0.0050261322,0.0050271531,0.0050310869,0.0050334910,0.0050359226,0.0050385582,0.0050400244,0.0050430951,0.0050462016,0.0050481572,0.0050499460,0.0050525587,0.0037888428,0.0038388651,0.0038861656,0.0039260134,0.0039673073,0.0040033829,0.0040354796,0.0040691064,0.0040967462,0.0041242971,0.0041504250,0.0041752520,0.0041980097,0.0042196346,0.0042411153,0.0042609806,0.0042789564,0.0042980032,0.0043140836,0.0043308828,0.0043466397,0.0043614162,0.0043751820,0.0043885225,0.0044017707,0.0044152967,0.0044256557,0.0044393216,0.0044496642,0.0044610815,0.0044726542,0.0044809779,0.0044919756,0.0045005253,0.0045102054,0.0045188789,0.0045271979,0.0045345311,0.0045429109,0.0045522355,0.0045598032,0.0045666225,0.0045737196,0.0045797226,0.0045876921,0.0045945236,0.0045998435,0.0046053217,0.0046128283,0.0046193968,0.0046237225,0.0046296284,0.0046338041,0.0046397129,0.0046460895,0.0046506700,0.0046559162,0.0046594902,0.0046648292,0.0046689332,0.0046736857,0.0046789055,0.0046820014,0.0046880718,0.0046912501,0.0046956056,0.0046999724,0.0047026557,0.0047078038,0.0047111599,0.0047144204,0.0047175520,0.0047227103,0.0047253523,0.0047286292,0.0047319583,0.0047347075,0.0047370576,0.0047409458,0.0047443723,0.0047473634,0.0047503381,0.0047530544,0.0047575792,0.0047582719,0.0047613816,0.0047654219,0.0047681336,0.0047709228,0.0047719135,0.0047757752,0.0047777860,0.0047805324,0.0047830920,0.0047844955,0.0047871108,0.0047903905,0.0047920952,0.0047936039,0.0047961203,0.0035925035,0.0036402431,0.0036850623,0.0037230041,0.0037627052,0.0037968958,0.0038271717,0.0038592927,0.0038856461,0.0039118069,0.0039369147,0.0039605410,0.0039819116,0.0040027082,0.0040230031,0.0040419151,0.0040590835,0.0040773886,0.0040924824,0.0041087081,0.0041237022,0.0041377013,0.0041509177,0.0041637541,0.0041760422,0.0041889797,0.0041991239,0.0042118940,0.0042220297,0.0042328539,0.0042438514,0.0042520444,0.0042625756,0.0042704397,0.0042796210,0.0042877873,0.0042955449,0.0043028722,0.0043108662,0.0043198069,0.0043268716,0.0043332052,0.0043403619,0.0043461686,0.0043537090,0.0043600461,0.0043653359,0.0043703945,0.0043777972,0.0043840879,0.0043882884,0.0043938962,0.0043974509,0.0044029829,0.0044095092,0.0044138281,0.0044190147,0.0044226827,0.0044275782,0.0044312594,0.0044357498,0.0044407974,0.0044439887,0.0044497201,0.0044524787,0.0044567008,0.0044609361,0.0044634096,0.0044683562,0.0044714643,0.0044748147,0.0044778356,0.0044825054,0.0044851518,0.0044881622,0.0044913947,0.0044942469,0.0044964792,0.0045003735,0.0045034905,0.0045063881,0.0045091783,0.0045112410,0.0045158871,0.0045165270,0.0045196888,0.0045236004,0.0045262489,0.0045288184,0.0045297847,0.0045333168,0.0045354137,0.0045379248,0.0045407131,0.0045417416,0.0045445935,0.0045469875,0.0045491192,0.0045505182,0.0045527122,0.0034063178,0.0034520602,0.0034948321,0.0035307712,0.0035684567,0.0036009668,0.0036298173,0.0036607106,0.0036854319,0.0037104861,0.0037344887,0.0037567558,0.0037772991,0.0037968615,0.0038163364,0.0038345091,0.0038508476,0.0038683108,0.0038824925,0.0038981658,0.0039123346,0.0039256543,0.0039381445,0.0039506790,0.0039623011,0.0039745173,0.0039842074,0.0039961854,0.0040060943,0.0040163550,0.0040266718,0.0040345388,0.0040447542,0.0040521440,0.0040610298,0.0040686291,0.0040761184,0.0040829445,0.0040908656,0.0040993072,0.0041060483,0.0041121369,0.0041189600,0.0041243776,0.0041319276,0.0041379525,0.0041428917,0.0041477153,0.0041547950,0.0041606610,0.0041646405,0.0041703593,0.0041734095,0.0041788049,0.0041851935,0.0041893214,0.0041938660,0.0041975567,0.0042025231,0.0042056862,0.0042101049,0.0042150578,0.0042182621,0.0042234672,0.0042259273,0.0042302974,0.0042339369,0.0042367204,0.0042412467,0.0042442433,0.0042475909,0.0042503284,0.0042546799,0.0042573072,0.0042601747,0.0042631705,0.0042659377,0.0042680939,0.0042718521,0.0042746334,0.0042775350,0.0042799410,0.0042821877,0.0042868390,0.0042872847,0.0042900474,0.0042939069,0.0042966058,0.0042990228,0.0042999198,0.0043033387,0.0043054750,0.0043076286,0.0043104396,0.0043113525,0.0043141526,0.0043162268,0.0043184373,0.0043196456,0.0043217364,0.0032301482,0.0032735368,0.0033143516,0.0033482358,0.0033844679,0.0034153861,0.0034427551,0.0034720661,0.0034955920,0.0035195171,0.0035422986,0.0035634265,0.0035832855,0.0036018634,0.0036201043,0.0036378471,0.0036530086,0.0036700695,0.0036830966,0.0036983631,0.0037118253,0.0037245696,0.0037365679,0.0037483234,0.0037596437,0.0037709235,0.0037802310,0.0037916410,0.0038011761,0.0038109860,0.0038209713,0.0038284875,0.0038381133,0.0038453182,0.0038537861,0.0038609303,0.0038683341,0.0038745503,0.0038820787,0.0038901293,0.0038964553,0.0039024814,0.0039089211,0.0039144167,0.0039212710,0.0039271017,0.0039315499,0.0039361496,0.0039432974,0.0039490381,0.0039524931,0.0039581382,0.0039610125,0.0039660654,0.0039718724,0.0039760250,0.0039803176,0.0039839604,0.0039886796,0.0039918140,0.0039959794,0.0040005181,0.0040038854,0.0040089903,0.0040113139,0.0040152796,0.0040188526,0.0040212931,0.0040259572,0.0040286412,0.0040321848,0.0040343297,0.0040384726,0.0040412489,0.0040437428,0.0040467280,0.0040493740,0.0040515162,0.0040553254,0.0040576809,0.0040603050,0.0040627130,0.0040652529,0.0040693045,0.0040698046,0.0040720850,0.0040761916,0.0040787012,0.0040807556,0.0040819359,0.0040850578,0.0040870996,0.0040892274,0.0040918592,0.0040927228,0.0040955691,0.0040973066,0.0040993364,0.0041007078,0.0041026299,0.0030632635,0.0031042688,0.0031429569,0.0031754675,0.0032099329,0.0032391314,0.0032654315,0.0032933370,0.0033158636,0.0033384911,0.0033603387,0.0033802671,0.0033992129,0.0034172006,0.0034343543,0.0034510145,0.0034654937,0.0034820994,0.0034944383,0.0035087887,0.0035217454,0.0035338656,0.0035453992,0.0035565876,0.0035671389,0.0035781482,0.0035871849,0.0035977497,0.0036068482,0.0036160307,0.0036257558,0.0036327205,0.0036421373,0.0036487466,0.0036572947,0.0036637775,0.0036708479,0.0036766745,0.0036838220,0.0036915746,0.0036977438,0.0037036371,0.0037097168,0.0037149150,0.0037217053,0.0037270012,0.0037313960,0.0037356411,0.0037424737,0.0037480114,0.0037514313,0.0037567939,0.0037592622,0.0037643650,0.0037696824,0.0037738765,0.0037777801,0.0037810851,0.0037857349,0.0037889185,0.0037927812,0.0037972307,0.0038004022,0.0038055972,0.0038076733,0.0038111782,0.0038147533,0.0038168556,0.0038215339,0.0038239770,0.0038274127,0.0038295718,0.0038333473,0.0038362291,0.0038385102,0.0038415013,0.0038439516,0.0038457770,0.0038497428,0.0038516477,0.0038542722,0.0038564720,0.0038589925,0.0038628052,0.0038633224,0.0038653998,0.0038694704,0.0038717101,0.0038737507,0.0038748193,0.0038780014,0.0038799297,0.0038819848,0.0038841236,0.0038852875,0.0038879067,0.0038896886,0.0038917932,0.0038928026,0.0038946071,0.0029051061,0.0029439545,0.0029806708,0.0030115065,0.0030447516,0.0030721404,0.0030971992,0.0031240058,0.0031452036,0.0031670336,0.0031876022,0.0032064377,0.0032244309,0.0032416976,0.0032579672,0.0032737310,0.0032878980,0.0033035582,0.0033154193,0.0033290812,0.0033415316,0.0033530433,0.0033641951,0.0033745843,0.0033848084,0.0033951191,0.0034036160,0.0034140685,0.0034224552,0.0034312187,0.0034403981,0.0034473444,0.0034561374,0.0034625971,0.0034707900,0.0034765483,0.0034834544,0.0034893660,0.0034959498,0.0035034325,0.0035092356,0.0035149982,0.0035206941,0.0035257987,0.0035321605,0.0035371405,0.0035416010,0.0035455680,0.0035519161,0.0035570531,0.0035604533,0.0035659440,0.0035680611,0.0035729371,0.0035780122,0.0035820400,0.0035857534,0.0035886629,0.0035932896,0.0035964193,0.0035996703,0.0036043634,0.0036072543,0.0036122260,0.0036141041,0.0036173813,0.0036210515,0.0036228982,0.0036276775,0.0036298832,0.0036330783,0.0036351308,0.0036386482,0.0036414738,0.0036435550,0.0036464262,0.0036488190,0.0036505059,0.0036543375,0.0036564391,0.0036589375,0.0036608996,0.0036631291,0.0036671170,0.0036675721,0.0036696801,0.0036732873,0.0036752919,0.0036776504,0.0036783747,0.0036814955,0.0036829751,0.0036853159,0.0036872207,0.0036883051,0.0036909794,0.0036926131,0.0036945573,0.0036954655,0.0036972911,0.0027550219,0.0027919908,0.0028268894,0.0028562101,0.0028878247,0.0029140856,0.0029375849,0.0029633483,0.0029835076,0.0030041656,0.0030236942,0.0030416388,0.0030589610,0.0030752239,0.0030910061,0.0031057095,0.0031193468,0.0031343806,0.0031455400,0.0031586955,0.0031704711,0.0031816038,0.0031920120,0.0032020899,0.0032116405,0.0032216517,0.0032297978,0.0032397176,0.0032477434,0.0032560340,0.0032647338,0.0032712668,0.0032796814,0.0032859502,0.0032938114,0.0032989547,0.0033058852,0.0033113523,0.0033176861,0.0033248811,0.0033303606,0.0033358312,0.0033412412,0.0033461291,0.0033521612,0.0033570177,0.0033614795,0.0033651442,0.0033708960,0.0033760953,0.0033792687,0.0033844622,0.0033864567,0.0033912594,0.0033960987,0.0033997641,0.0034035485,0.0034064306,0.0034106003,0.0034134979,0.0034166684,0.0034211270,0.0034238916,0.0034290076,0.0034304988,0.0034334958,0.0034371760,0.0034391042,0.0034436006,0.0034455447,0.0034488208,0.0034504807,0.0034537434,0.0034566031,0.0034585491,0.0034616545,0.0034636491,0.0034653633,0.0034687421,0.0034710065,0.0034735565,0.0034752048,0.0034774911,0.0034811319,0.0034816511,0.0034837446,0.0034869906,0.0034887588,0.0034912968,0.0034922212,0.0034951192,0.0034965504,0.0034985161,0.0035003602,0.0035013031,0.0035042361,0.0035054603,0.0035076195,0.0035085637,0.0035100122,0.0026125527,0.0026477909,0.0026810840,0.0027089769,0.0027388218,0.0027641981,0.0027862266,0.0028109595,0.0028303338,0.0028498117,0.0028683169,0.0028855611,0.0029019353,0.0029175236,0.0029324270,0.0029466423,0.0029594466,0.0029737164,0.0029842557,0.0029971601,0.0030081753,0.0030188438,0.0030288885,0.0030384455,0.0030476653,0.0030568757,0.0030647840,0.0030743040,0.0030819867,0.0030899436,0.0030979911,0.0031043377,0.0031123817,0.0031181875,0.0031258404,0.0031308718,0.0031374137,0.0031424478,0.0031488603,0.0031553299,0.0031607735,0.0031658746,0.0031710933,0.0031757706,0.0031814968,0.0031861593,0.0031904336,0.0031939495,0.0031996439,0.0032044103,0.0032073326,0.0032124893,0.0032141693,0.0032189452,0.0032235106,0.0032272450,0.0032305905,0.0032333826,0.0032373648,0.0032400698,0.0032431870,0.0032474421,0.0032498967,0.0032549753,0.0032563764,0.0032591678,0.0032625880,0.0032645982,0.0032689338,0.0032706539,0.0032736498,0.0032752469,0.0032782937,0.0032811097,0.0032832091,0.0032861388,0.0032880077,0.0032895578,0.0032930701,0.0032949065,0.0032976195,0.0032990930,0.0033013986,0.0033048053,0.0033050777,0.0033072396,0.0033104602,0.0033120624,0.0033142864,0.0033154325,0.0033182529,0.0033195588,0.0033212578,0.0033231228,0.0033241341,0.0033267108,0.0033281100,0.0033299216,0.0033311240,0.0033324028,0.0024777127,0.0025113073,0.0025429168,0.0025693512,0.0025980241,0.0026218847,0.0026428331,0.0026663692,0.0026848256,0.0027032190,0.0027211167,0.0027375199,0.0027532451,0.0027676155,0.0027820165,0.0027954230,0.0028077706,0.0028215281,0.0028316169,0.0028437074,0.0028542209,0.0028643065,0.0028739753,0.0028829126,0.0028919731,0.0029007437,0.0029082824,0.0029172225,0.0029248036,0.0029320123,0.0029400164,0.0029458450,0.0029537253,0.0029592122,0.0029666326,0.0029712396,0.0029774070,0.0029824616,0.0029883966,0.0029946329,0.0029995988,0.0030046050,0.0030095165,0.0030140771,0.0030194024,0.0030240025,0.0030283100,0.0030315968,0.0030367614,0.0030416062,0.0030442583,0.0030494155,0.0030510662,0.0030554933,0.0030597975,0.0030632457,0.0030665185,0.0030691929,0.0030729284,0.0030755925,0.0030785585,0.0030824573,0.0030848424,0.0030895588,0.0030911952,0.0030937733,0.0030969709,0.0030988992,0.0031030805,0.0031047923,0.0031074199,0.0031089473,0.0031121445,0.0031148358,0.0031167835,0.0031196632,0.0031213377,0.0031229882,0.0031263757,0.0031279569,0.0031304896,0.0031317707,0.0031340829,0.0031372670,0.0031376166,0.0031397442,0.0031429564,0.0031443632,0.0031463932,0.0031477235,0.0031500710,0.0031513718,0.0031530630,0.0031548330,0.0031557572,0.0031580784,0.0031596910,0.0031613865,0.0031626610,0.0031638953,0.0023499009,0.0023818377,0.0024118074,0.0024370093,0.0024645319,0.0024869302,0.0025068949,0.0025293766,0.0025468963,0.0025644072,0.0025815519,0.0025970619,0.0026120417,0.0026255363,0.0026393278,0.0026525193,0.0026640349,0.0026772639,0.0026867997,0.0026982332,0.0027081341,0.0027179786,0.0027270152,0.0027356697,0.0027442501,0.0027526394,0.0027597151,0.0027683531,0.0027757129,0.0027825904,0.0027902531,0.0027958361,0.0028028898,0.0028082305,0.0028154232,0.0028198720,0.0028255733,0.0028306344,0.0028362792,0.0028421922,0.0028470654,0.0028516993,0.0028562238,0.0028607010,0.0028656987,0.0028701031,0.0028742992,0.0028774905,0.0028821873,0.0028868364,0.0028896803,0.0028946324,0.0028961230,0.0029003054,0.0029044331,0.0029079282,0.0029108635,0.0029133154,0.0029167981,0.0029195385,0.0029223132,0.0029260181,0.0029281204,0.0029326139,0.0029342812,0.0029368754,0.0029396792,0.0029417278,0.0029457639,0.0029475986,0.0029497235,0.0029512050,0.0029545632,0.0029567878,0.0029586940,0.0029616344,0.0029629789,0.0029646461,0.0029679244,0.0029693301,0.0029718927,0.0029732819,0.0029752177,0.0029784039,0.0029787082,0.0029806405,0.0029838999,0.0029853603,0.0029871826,0.0029884155,0.0029904025,0.0029919766,0.0029935314,0.0029951269,0.0029961579,0.0029982363,0.0029997784,0.0030012203,0.0030027821,0.0030037889,0.0022285798,0.0022588849,0.0022876122,0.0023114047,0.0023379051,0.0023589393,0.0023779042,0.0023995022,0.0024161907,0.0024328289,0.0024490252,0.0024638392,0.0024781186,0.0024910413,0.0025042315,0.0025165663,0.0025277339,0.0025403162,0.0025494050,0.0025602027,0.0025694961,0.0025789745,0.0025876349,0.0025959778,0.0026040366,0.0026121694,0.0026187525,0.0026271860,0.0026341742,0.0026407437,0.0026481304,0.0026531419,0.0026599434,0.0026650725,0.0026720609,0.0026761541,0.0026818686,0.0026865092,0.0026918640,0.0026974878,0.0027022437,0.0027064667,0.0027109451,0.0027153899,0.0027200548,0.0027241350,0.0027284219,0.0027314484,0.0027356281,0.0027400778,0.0027427532,0.0027477576,0.0027490104,0.0027529943,0.0027569077,0.0027603462,0.0027630548,0.0027653906,0.0027687754,0.0027713039,0.0027741313,0.0027774037,0.0027795314,0.0027836390,0.0027855160,0.0027879105,0.0027907575,0.0027924545,0.0027963793,0.0027982180,0.0028000595,0.0028016404,0.0028048726,0.0028069586,0.0028087242,0.0028116779,0.0028127584,0.0028146565,0.0028174701,0.0028187830,0.0028214252,0.0028226868,0.0028248618,0.0028278006,0.0028278761,0.0028298490,0.0028327944,0.0028342629,0.0028359318,0.0028372264,0.0028390757,0.0028405456,0.0028422001,0.0028435647,0.0028446273,0.0028466231,0.0028481187,0.0028494046,0.0028509908,0.0028517045,0.0021136888,0.0021424771,0.0021699090,0.0021925605,0.0022176213,0.0022374972,0.0022557294,0.0022762939,0.0022920790,0.0023079502,0.0023232525,0.0023379206,0.0023510725,0.0023634357,0.0023760809,0.0023877941,0.0023983999,0.0024105495,0.0024192694,0.0024293896,0.0024380216,0.0024472291,0.0024555759,0.0024634199,0.0024711166,0.0024786785,0.0024850706,0.0024931058,0.0024998282,0.0025059975,0.0025130626,0.0025179106,0.0025244207,0.0025295064,0.0025356664,0.0025398822,0.0025453543,0.0025499340,0.0025547973,0.0025603205,0.0025646432,0.0025689859,0.0025730013,0.0025773293,0.0025816451,0.0025857481,0.0025895273,0.0025927778,0.0025966999,0.0026006123,0.0026034416,0.0026083373,0.0026093392,0.0026133289,0.0026169861,0.0026204058,0.0026225523,0.0026250980,0.0026282315,0.0026308019,0.0026335236,0.0026367764,0.0026387067,0.0026425690,0.0026444086,0.0026466067,0.0026493263,0.0026509570,0.0026547322,0.0026563473,0.0026580545,0.0026594935,0.0026627049,0.0026648550,0.0026662303,0.0026691896,0.0026703406,0.0026721107,0.0026747794,0.0026762667,0.0026785632,0.0026796854,0.0026818804,0.0026846976,0.0026848873,0.0026865281,0.0026895108,0.0026910240,0.0026924827,0.0026936692,0.0026954681,0.0026965971,0.0026982886,0.0026997759,0.0027008939,0.0027027515,0.0027040281,0.0027052501,0.0027068616,0.0027075780,0.0020047243,0.0020323270,0.0020582247,0.0020795953,0.0021035665,0.0021225550,0.0021397806,0.0021595409,0.0021746343,0.0021896614,0.0022041233,0.0022182574,0.0022306332,0.0022426702,0.0022545157,0.0022655954,0.0022759590,0.0022872548,0.0022956791,0.0023052819,0.0023134274,0.0023221264,0.0023302450,0.0023377543,0.0023448360,0.0023522611,0.0023582043,0.0023660966,0.0023721550,0.0023781784,0.0023847942,0.0023894555,0.0023960534,0.0024008798,0.0024064298,0.0024104950,0.0024158878,0.0024201900,0.0024249136,0.0024300903,0.0024342176,0.0024384079,0.0024421409,0.0024463689,0.0024504408,0.0024545513,0.0024579472,0.0024609704,0.0024647334,0.0024684998,0.0024712182,0.0024758887,0.0024769066,0.0024806886,0.0024841651,0.0024874354,0.0024895508,0.0024921435,0.0024950771,0.0024973090,0.0024998488,0.0025030536,0.0025047828,0.0025084591,0.0025103046,0.0025123274,0.0025150024,0.0025165702,0.0025201292,0.0025216182,0.0025236396,0.0025249028,0.0025278482,0.0025299296,0.0025313769,0.0025340745,0.0025351880,0.0025368622,0.0025393083,0.0025407530,0.0025429021,0.0025440612,0.0025461558,0.0025488977,0.0025490442,0.0025506535,0.0025532469,0.0025550538,0.0025561896,0.0025572868,0.0025591204,0.0025602684,0.0025617035,0.0025632954,0.0025642723,0.0025663248,0.0025672684,0.0025685725,0.0025699472,0.0025706913,0.0019014792,0.0019278808,0.0019522413,0.0019726560,0.0019955213,0.0020135715,0.0020298124,0.0020487923,0.0020629451,0.0020773371,0.0020910809,0.0021044627,0.0021163591,0.0021278491,0.0021390972,0.0021496591,0.0021596891,0.0021701388,0.0021783605,0.0021874223,0.0021951964,0.0022035350,0.0022114228,0.0022184910,0.0022251875,0.0022323472,0.0022380340,0.0022454930,0.0022512570,0.0022571859,0.0022632998,0.0022675754,0.0022742356,0.0022786664,0.0022838549,0.0022879883,0.0022928377,0.0022971226,0.0023015189,0.0023064628,0.0023105528,0.0023144791,0.0023180871,0.0023220648,0.0023262608,0.0023296615,0.0023332221,0.0023359593,0.0023395925,0.0023431078,0.0023458540,0.0023504089,0.0023510462,0.0023549975,0.0023580612,0.0023612564,0.0023630932,0.0023658775,0.0023688478,0.0023708483,0.0023731578,0.0023761958,0.0023779693,0.0023814162,0.0023831443,0.0023850419,0.0023874993,0.0023891373,0.0023925742,0.0023939663,0.0023959357,0.0023970316,0.0023999324,0.0024017669,0.0024033116,0.0024058772,0.0024067371,0.0024083602,0.0024106020,0.0024122309,0.0024141269,0.0024151260,0.0024172741,0.0024197874,0.0024203744,0.0024215713,0.0024242163,0.0024257824,0.0024270482,0.0024280469,0.0024297430,0.0024306341,0.0024321235,0.0024335457,0.0024346132,0.0024365810,0.0024373219,0.0024388286,0.0024403100,0.0024407960,0.0018036889,0.0018285743,0.0018517794,0.0018714426,0.0018929010,0.0019102294,0.0019257781,0.0019437718,0.0019571534,0.0019709524,0.0019840204,0.0019967804,0.0020080543,0.0020188816,0.0020297052,0.0020398043,0.0020493780,0.0020591884,0.0020668357,0.0020758256,0.0020831040,0.0020910703,0.0020986237,0.0021053959,0.0021116488,0.0021186594,0.0021237865,0.0021309324,0.0021366265,0.0021421763,0.0021481027,0.0021520924,0.0021584701,0.0021626100,0.0021673557,0.0021715419,0.0021761141,0.0021802419,0.0021846478,0.0021892404,0.0021932415,0.0021967371,0.0022004376,0.0022042227,0.0022081452,0.0022112654,0.0022145805,0.0022174533,0.0022208342,0.0022243288,0.0022268672,0.0022311711,0.0022316702,0.0022355592,0.0022385096,0.0022412405,0.0022431857,0.0022459955,0.0022488316,0.0022507140,0.0022530089,0.0022558328,0.0022574024,0.0022607807,0.0022625360,0.0022641496,0.0022664924,0.0022683717,0.0022715781,0.0022725537,0.0022746497,0.0022757075,0.0022784348,0.0022802331,0.0022815675,0.0022840789,0.0022848155,0.0022863204,0.0022885624,0.0022901888,0.0022921078,0.0022930505,0.0022951328,0.0022973122,0.0022982018,0.0022992701,0.0023016920,0.0023031418,0.0023044651,0.0023053777,0.0023068618,0.0023076670,0.0023093128,0.0023105037,0.0023114725,0.0023133530,0.0023141147,0.0023155625,0.0023169128,0.0023174145,0.0017108707,0.0017346565,0.0017566017,0.0017753645,0.0017956783,0.0018122877,0.0018270233,0.0018441592,0.0018568458,0.0018698582,0.0018823482,0.0018943751,0.0019052906,0.0019156509,0.0019260257,0.0019353800,0.0019445399,0.0019538741,0.0019610860,0.0019700167,0.0019766550,0.0019843637,0.0019915410,0.0019979731,0.0020040259,0.0020106391,0.0020156933,0.0020224288,0.0020277867,0.0020331757,0.0020386860,0.0020426013,0.0020487595,0.0020524604,0.0020571234,0.0020609916,0.0020652594,0.0020693311,0.0020736686,0.0020779429,0.0020818144,0.0020852215,0.0020887285,0.0020923764,0.0020960761,0.0020991528,0.0021021786,0.0021050922,0.0021081330,0.0021114282,0.0021138127,0.0021181400,0.0021186540,0.0021221916,0.0021249847,0.0021276910,0.0021295010,0.0021321872,0.0021347750,0.0021364976,0.0021388857,0.0021416743,0.0021429232,0.0021464214,0.0021480139,0.0021494370,0.0021517532,0.0021534969,0.0021565497,0.0021574779,0.0021595604,0.0021606430,0.0021630782,0.0021648126,0.0021662062,0.0021685735,0.0021690814,0.0021707065,0.0021729396,0.0021742532,0.0021761462,0.0021771362,0.0021791146,0.0021812950,0.0021820966,0.0021829750,0.0021851774,0.0021867735,0.0021879778,0.0021888508,0.0021904528,0.0021910565,0.0021926673,0.0021937012,0.0021945744,0.0021965991,0.0021971474,0.0021987245,0.0021999033,0.0022004439,0.0016227695,0.0016453744,0.0016662951,0.0016841867,0.0017033310,0.0017192433,0.0017331486,0.0017496956,0.0017614739,0.0017741612,0.0017859724,0.0017973809,0.0018079174,0.0018177642,0.0018275633,0.0018365445,0.0018452451,0.0018540369,0.0018609924,0.0018693737,0.0018759424,0.0018831789,0.0018898666,0.0018960206,0.0019018714,0.0019080810,0.0019129188,0.0019194007,0.0019243797,0.0019297525,0.0019349946,0.0019387902,0.0019443901,0.0019480851,0.0019524599,0.0019562474,0.0019603478,0.0019640812,0.0019684782,0.0019725048,0.0019761257,0.0019792406,0.0019827986,0.0019861454,0.0019896430,0.0019925934,0.0019954891,0.0019983079,0.0020011319,0.0020043229,0.0020064514,0.0020109203,0.0020111352,0.0020145515,0.0020171649,0.0020198440,0.0020216344,0.0020243219,0.0020265724,0.0020281199,0.0020304005,0.0020332114,0.0020343537,0.0020377598,0.0020392462,0.0020407228,0.0020430543,0.0020443674,0.0020474561,0.0020483612,0.0020501636,0.0020510487,0.0020535613,0.0020552681,0.0020566687,0.0020588842,0.0020593637,0.0020609901,0.0020628243,0.0020643664,0.0020661844,0.0020670107,0.0020689326,0.0020709560,0.0020718114,0.0020727537,0.0020746458,0.0020763156,0.0020772530,0.0020782364,0.0020798159,0.0020805212,0.0020818920,0.0020829235,0.0020836175,0.0020854913,0.0020860698,0.0020877844,0.0020887735,0.0020893931,0.0015392709,0.0015607515,0.0015808731,0.0015977027,0.0016159056,0.0016309892,0.0016442959,0.0016601085,0.0016713452,0.0016834067,0.0016945868,0.0017054597,0.0017154841,0.0017247940,0.0017342822,0.0017427415,0.0017511204,0.0017593799,0.0017660021,0.0017740360,0.0017801525,0.0017869167,0.0017936421,0.0017993520,0.0018051049,0.0018107663,0.0018155997,0.0018218043,0.0018264557,0.0018317336,0.0018364780,0.0018402120,0.0018453416,0.0018490776,0.0018531805,0.0018568354,0.0018608124,0.0018644883,0.0018684756,0.0018724248,0.0018756202,0.0018787391,0.0018822105,0.0018851752,0.0018888754,0.0018913898,0.0018943085,0.0018967713,0.0018994203,0.0019027348,0.0019047363,0.0019091140,0.0019091291,0.0019123409,0.0019150496,0.0019173861,0.0019191330,0.0019216649,0.0019240211,0.0019253989,0.0019275366,0.0019302854,0.0019312479,0.0019346452,0.0019361353,0.0019374940,0.0019397215,0.0019410796,0.0019440055,0.0019448739,0.0019466123,0.0019473920,0.0019497614,0.0019515125,0.0019528588,0.0019547294,0.0019552610,0.0019567670,0.0019582580,0.0019599543,0.0019618352,0.0019626317,0.0019643768,0.0019663532,0.0019671780,0.0019680589,0.0019699183,0.0019713868,0.0019725915,0.0019732898,0.0019747123,0.0019752134,0.0019766427,0.0019779264,0.0019783697,0.0019800624,0.0019807272,0.0019824757,0.0019834531,0.0019839620,0.0014600299,0.0014803429,0.0014998796,0.0015158403,0.0015331078,0.0015473812,0.0015600529,0.0015750694,0.0015857066,0.0015972742,0.0016079070,0.0016182790,0.0016278743,0.0016367780,0.0016457334,0.0016538175,0.0016616672,0.0016696873,0.0016760045,0.0016835656,0.0016894627,0.0016958619,0.0017022483,0.0017076826,0.0017130755,0.0017185386,0.0017230675,0.0017290496,0.0017334743,0.0017386009,0.0017429617,0.0017466743,0.0017515192,0.0017550244,0.0017589285,0.0017625603,0.0017663581,0.0017699238,0.0017737501,0.0017773206,0.0017802950,0.0017834729,0.0017866253,0.0017895118,0.0017929100,0.0017953907,0.0017981577,0.0018006562,0.0018032642,0.0018061905,0.0018081088,0.0018123370,0.0018123418,0.0018155920,0.0018181238,0.0018202400,0.0018217791,0.0018243932,0.0018264959,0.0018278484,0.0018298196,0.0018326450,0.0018334892,0.0018367007,0.0018381639,0.0018395546,0.0018415169,0.0018429119,0.0018455738,0.0018465217,0.0018481032,0.0018489284,0.0018512596,0.0018527947,0.0018540573,0.0018558793,0.0018563994,0.0018577850,0.0018591842,0.0018610103,0.0018627132,0.0018636732,0.0018652373,0.0018670753,0.0018680564,0.0018686764,0.0018704440,0.0018717709,0.0018728471,0.0018735621,0.0018750478,0.0018754567,0.0018769278,0.0018782388,0.0018784684,0.0018800909,0.0018805878,0.0018823467,0.0018834138,0.0018836918,0.0013849497,0.0014041814,0.0014228591,0.0014381367,0.0014545378,0.0014681252,0.0014801609,0.0014943550,0.0015045919,0.0015156071,0.0015256856,0.0015355684,0.0015446752,0.0015532729,0.0015616218,0.0015694637,0.0015770567,0.0015845470,0.0015904222,0.0015977404,0.0016033449,0.0016093985,0.0016154937,0.0016206637,0.0016259642,0.0016310192,0.0016354982,0.0016410713,0.0016453525,0.0016502861,0.0016544284,0.0016578829,0.0016624410,0.0016659410,0.0016695563,0.0016730800,0.0016768014,0.0016801319,0.0016837145,0.0016870708,0.0016899183,0.0016929136,0.0016959647,0.0016988255,0.0017018696,0.0017043613,0.0017069869,0.0017096017,0.0017117688,0.0017146167,0.0017165011,0.0017206328,0.0017202602,0.0017235459,0.0017262914,0.0017280157,0.0017293409,0.0017320175,0.0017339146,0.0017353523,0.0017370798,0.0017397551,0.0017406744,0.0017436847,0.0017453062,0.0017464696,0.0017483846,0.0017498228,0.0017521621,0.0017530610,0.0017544554,0.0017555775,0.0017575318,0.0017591602,0.0017601926,0.0017619389,0.0017626482,0.0017640104,0.0017652164,0.0017670644,0.0017685761,0.0017697996,0.0017709854,0.0017727841,0.0017736927,0.0017745317,0.0017760870,0.0017771011,0.0017783582,0.0017789193,0.0017803508,0.0017809239,0.0017821804,0.0017834312,0.0017834999,0.0017852384,0.0017856062,0.0017874229,0.0017883484,0.0017884608,0.0013138324,0.0013320624,0.0013498099,0.0013644162,0.0013797818,0.0013929019,0.0014044294,0.0014179447,0.0014274865,0.0014380304,0.0014476015,0.0014570465,0.0014656910,0.0014738886,0.0014820116,0.0014894540,0.0014967118,0.0015036859,0.0015092509,0.0015162719,0.0015214896,0.0015273127,0.0015332888,0.0015381689,0.0015430340,0.0015480128,0.0015521872,0.0015574122,0.0015618265,0.0015663373,0.0015702497,0.0015737290,0.0015780225,0.0015814352,0.0015848202,0.0015881441,0.0015917505,0.0015949876,0.0015983832,0.0016015519,0.0016041366,0.0016070142,0.0016098108,0.0016127127,0.0016155882,0.0016180678,0.0016205611,0.0016229493,0.0016249243,0.0016277785,0.0016294834,0.0016333937,0.0016332933,0.0016361803,0.0016387847,0.0016405116,0.0016418533,0.0016442303,0.0016462203,0.0016475090,0.0016493007,0.0016516411,0.0016525566,0.0016554771,0.0016568822,0.0016581719,0.0016598810,0.0016613828,0.0016635032,0.0016645134,0.0016657344,0.0016668117,0.0016687416,0.0016702379,0.0016711727,0.0016729290,0.0016735480,0.0016751714,0.0016760976,0.0016777133,0.0016792267,0.0016804735,0.0016814408,0.0016832137,0.0016841586,0.0016849589,0.0016864983,0.0016875753,0.0016884852,0.0016892286,0.0016904979,0.0016910972,0.0016921788,0.0016933768,0.0016936251,0.0016950433,0.0016956318,0.0016971625,0.0016980968,0.0016982647,0.0012463993,0.0012636920,0.0012804765,0.0012944444,0.0013090500,0.0013215370,0.0013326094,0.0013453831,0.0013544954,0.0013643090,0.0013735718,0.0013825751,0.0013908594,0.0013987044,0.0014063942,0.0014134784,0.0014203382,0.0014271111,0.0014322249,0.0014389876,0.0014440322,0.0014497075,0.0014552095,0.0014597043,0.0014645657,0.0014692606,0.0014732046,0.0014782630,0.0014823506,0.0014867527,0.0014906313,0.0014939783,0.0014978647,0.0015011052,0.0015043709,0.0015076773,0.0015109794,0.0015140754,0.0015174274,0.0015204021,0.0015228952,0.0015254875,0.0015282714,0.0015310219,0.0015337278,0.0015360423,0.0015384443,0.0015407620,0.0015426137,0.0015452264,0.0015469389,0.0015507698,0.0015505186,0.0015533086,0.0015557949,0.0015573954,0.0015587590,0.0015609641,0.0015626949,0.0015642011,0.0015658714,0.0015681392,0.0015688991,0.0015717093,0.0015731337,0.0015742858,0.0015758760,0.0015774194,0.0015795548,0.0015802507,0.0015815221,0.0015823768,0.0015846144,0.0015860119,0.0015867759,0.0015885401,0.0015889341,0.0015904601,0.0015914639,0.0015929803,0.0015944485,0.0015957020,0.0015965513,0.0015982226,0.0015991755,0.0015998238,0.0016014232,0.0016023279,0.0016032507,0.0016039971,0.0016051469,0.0016057618,0.0016067195,0.0016079532,0.0016082078,0.0016096342,0.0016101199,0.0016116692,0.0016125782,0.0016127005,0.0011824737,0.0011988133,0.0012147524,0.0012280938,0.0012419705,0.0012538492,0.0012644612,0.0012765988,0.0012852418,0.0012946574,0.0013034736,0.0013119530,0.0013198048,0.0013273707,0.0013346353,0.0013413270,0.0013478996,0.0013543928,0.0013591449,0.0013658110,0.0013704180,0.0013758470,0.0013812461,0.0013853918,0.0013898197,0.0013945547,0.0013984052,0.0014033424,0.0014071954,0.0014111017,0.0014150488,0.0014180193,0.0014217707,0.0014249065,0.0014280762,0.0014312216,0.0014343950,0.0014372209,0.0014405327,0.0014433119,0.0014457268,0.0014483272,0.0014508734,0.0014533996,0.0014559486,0.0014580133,0.0014604208,0.0014626018,0.0014645429,0.0014670454,0.0014686054,0.0014723682,0.0014719882,0.0014746887,0.0014771057,0.0014785149,0.0014799040,0.0014820611,0.0014836738,0.0014851495,0.0014866067,0.0014888637,0.0014893616,0.0014921834,0.0014936468,0.0014948423,0.0014962375,0.0014976824,0.0014996586,0.0015004919,0.0015016984,0.0015025462,0.0015045690,0.0015059528,0.0015066363,0.0015084198,0.0015086628,0.0015100166,0.0015111062,0.0015125563,0.0015138635,0.0015152971,0.0015159563,0.0015176582,0.0015185274,0.0015189857,0.0015206885,0.0015215397,0.0015224791,0.0015232118,0.0015243530,0.0015249110,0.0015255910,0.0015268228,0.0015269675,0.0015286040,0.0015288745,0.0015305189,0.0015311952,0.0015313571,0.0011218298,0.0011374120,0.0011524507,0.0011651394,0.0011782308,0.0011896517,0.0011997125,0.0012111660,0.0012195052,0.0012284293,0.0012368723,0.0012449078,0.0012525498,0.0012596492,0.0012665248,0.0012730102,0.0012791177,0.0012853744,0.0012900089,0.0012963281,0.0013006929,0.0013057512,0.0013110766,0.0013148566,0.0013190869,0.0013236182,0.0013273222,0.0013320005,0.0013356573,0.0013393648,0.0013432252,0.0013460890,0.0013496738,0.0013525922,0.0013554767,0.0013585564,0.0013617029,0.0013643079,0.0013675569,0.0013702117,0.0013723646,0.0013748951,0.0013774024,0.0013797893,0.0013821429,0.0013842111,0.0013865279,0.0013885090,0.0013903788,0.0013925651,0.0013941573,0.0013977403,0.0013976028,0.0013999265,0.0014022762,0.0014036716,0.0014049238,0.0014072353,0.0014086028,0.0014099960,0.0014115392,0.0014134627,0.0014140219,0.0014167659,0.0014181387,0.0014191908,0.0014207659,0.0014218532,0.0014239325,0.0014247208,0.0014259308,0.0014265548,0.0014287102,0.0014297651,0.0014306086,0.0014321020,0.0014324781,0.0014338440,0.0014348404,0.0014360894,0.0014374344,0.0014388466,0.0014393866,0.0014410642,0.0014418469,0.0014422471,0.0014439186,0.0014448401,0.0014457031,0.0014463432,0.0014476134,0.0014479355,0.0014485931,0.0014498889,0.0014500672,0.0014515917,0.0014518583,0.0014533500,0.0014540006,0.0014542865,0.0010642002,0.0010790994,0.0010932614,0.0011054293,0.0011178110,0.0011287806,0.0011383224,0.0011492742,0.0011571413,0.0011655993,0.0011736216,0.0011813337,0.0011886196,0.0011955020,0.0012019460,0.0012081944,0.0012139508,0.0012198120,0.0012242698,0.0012302941,0.0012345275,0.0012394016,0.0012442627,0.0012479598,0.0012518914,0.0012563421,0.0012598548,0.0012642999,0.0012678693,0.0012712800,0.0012749461,0.0012775807,0.0012811246,0.0012840201,0.0012866811,0.0012896035,0.0012925231,0.0012952241,0.0012981102,0.0013008287,0.0013028605,0.0013051498,0.0013076787,0.0013098722,0.0013122434,0.0013140513,0.0013162033,0.0013182636,0.0013199140,0.0013221941,0.0013236676,0.0013270129,0.0013269376,0.0013291567,0.0013313773,0.0013325722,0.0013339175,0.0013361162,0.0013373007,0.0013388210,0.0013399220,0.0013419421,0.0013426001,0.0013452335,0.0013464526,0.0013474605,0.0013489792,0.0013501481,0.0013520462,0.0013528061,0.0013538833,0.0013547076,0.0013564448,0.0013575127,0.0013584149,0.0013598116,0.0013602121,0.0013614426,0.0013624992,0.0013636714,0.0013648691,0.0013661781,0.0013668383,0.0013683858,0.0013690357,0.0013695529,0.0013711936,0.0013720942,0.0013728488,0.0013735468,0.0013746644,0.0013750698,0.0013755384,0.0013766958,0.0013769482,0.0013784367,0.0013787229,0.0013802298,0.0013806075,0.0013809582,0.0010096056,0.0010238058,0.0010372849,0.0010489316,0.0010605720,0.0010709833,0.0010801382,0.0010905396,0.0010981715,0.0011060253,0.0011137054,0.0011208610,0.0011278793,0.0011344742,0.0011405716,0.0011466289,0.0011521528,0.0011577025,0.0011619646,0.0011675627,0.0011715597,0.0011763397,0.0011810904,0.0011845480,0.0011881934,0.0011924547,0.0011958474,0.0012001311,0.0012034411,0.0012066452,0.0012102730,0.0012126842,0.0012161124,0.0012188578,0.0012213601,0.0012240963,0.0012269901,0.0012295156,0.0012323047,0.0012348264,0.0012367663,0.0012389467,0.0012415106,0.0012434514,0.0012457434,0.0012475988,0.0012494734,0.0012516272,0.0012529651,0.0012552353,0.0012567439,0.0012599517,0.0012599233,0.0012621439,0.0012639193,0.0012649874,0.0012664285,0.0012685765,0.0012697260,0.0012712724,0.0012723782,0.0012741527,0.0012747876,0.0012772882,0.0012784363,0.0012794318,0.0012810278,0.0012818908,0.0012838612,0.0012844897,0.0012856367,0.0012865424,0.0012879044,0.0012890019,0.0012896383,0.0012913162,0.0012915018,0.0012928563,0.0012937753,0.0012948565,0.0012962096,0.0012973526,0.0012980172,0.0012993230,0.0013001594,0.0013006128,0.0013021403,0.0013030784,0.0013037059,0.0013042157,0.0013053680,0.0013057132,0.0013062189,0.0013072842,0.0013076265,0.0013089907,0.0013091755,0.0013106883,0.0013110557,0.0013114895,0.0009578562,0.0009713916,0.0009842913,0.0009952666,0.0010063549,0.0010162193,0.0010248613,0.0010348107,0.0010421081,0.0010495682,0.0010568625,0.0010637238,0.0010703724,0.0010765736,0.0010825085,0.0010883280,0.0010934643,0.0010988903,0.0011029143,0.0011082296,0.0011120010,0.0011165111,0.0011211401,0.0011242428,0.0011279005,0.0011320745,0.0011350911,0.0011391815,0.0011424088,0.0011454639,0.0011488741,0.0011511872,0.0011543625,0.0011570160,0.0011594078,0.0011620575,0.0011649384,0.0011670751,0.0011697832,0.0011722907,0.0011740866,0.0011762516,0.0011785486,0.0011806482,0.0011825026,0.0011844213,0.0011862736,0.0011883604,0.0011895662,0.0011916747,0.0011931598,0.0011962367,0.0011961823,0.0011983230,0.0012001093,0.0012009652,0.0012022669,0.0012043412,0.0012056868,0.0012070606,0.0012081086,0.0012099545,0.0012104176,0.0012127107,0.0012139895,0.0012147436,0.0012165024,0.0012171157,0.0012189351,0.0012196375,0.0012206486,0.0012215198,0.0012229124,0.0012240572,0.0012244917,0.0012261055,0.0012263090,0.0012275635,0.0012284503,0.0012294854,0.0012308290,0.0012320202,0.0012324385,0.0012338609,0.0012345606,0.0012350462,0.0012364739,0.0012375017,0.0012381393,0.0012384495,0.0012395547,0.0012398623,0.0012404057,0.0012413516,0.0012417723,0.0012430391,0.0012431574,0.0012446292,0.0012451121,0.0012454741,0.0009087160,0.0009215517,0.0009338859,0.0009442353,0.0009549457,0.0009642493,0.0009724919,0.0009817550,0.0009889063,0.0009960241,0.0010029071,0.0010094851,0.0010158973,0.0010217237,0.0010273781,0.0010328589,0.0010379390,0.0010428930,0.0010466658,0.0010518308,0.0010554034,0.0010596480,0.0010641143,0.0010671124,0.0010705813,0.0010745781,0.0010773317,0.0010813496,0.0010844178,0.0010875063,0.0010904871,0.0010926286,0.0010958095,0.0010983629,0.0011005864,0.0011031214,0.0011060134,0.0011079204,0.0011104543,0.0011129098,0.0011145641,0.0011167945,0.0011188570,0.0011208986,0.0011226921,0.0011243503,0.0011263944,0.0011284318,0.0011295941,0.0011313203,0.0011327882,0.0011355474,0.0011357131,0.0011377382,0.0011394005,0.0011403356,0.0011414585,0.0011435240,0.0011445914,0.0011462113,0.0011470692,0.0011488582,0.0011493062,0.0011514987,0.0011526810,0.0011533308,0.0011551254,0.0011556191,0.0011576012,0.0011580969,0.0011589669,0.0011598137,0.0011611622,0.0011622304,0.0011626271,0.0011642265,0.0011645760,0.0011656204,0.0011664947,0.0011675566,0.0011688545,0.0011699914,0.0011703702,0.0011717061,0.0011723201,0.0011729730,0.0011741076,0.0011749998,0.0011756728,0.0011760191,0.0011770718,0.0011774010,0.0011779423,0.0011788027,0.0011791297,0.0011803966,0.0011804685,0.0011819172,0.0011823256,0.0011828804,0.0008621306,0.0008744354,0.0008860789,0.0008959433,0.0009059784,0.0009150109,0.0009228418,0.0009317018,0.0009383739,0.0009451859,0.0009517484,0.0009581406,0.0009640476,0.0009697403,0.0009751216,0.0009803769,0.0009850156,0.0009898943,0.0009933928,0.0009983586,0.0010017820,0.0010058739,0.0010101044,0.0010128123,0.0010162353,0.0010200054,0.0010226598,0.0010264084,0.0010293741,0.0010323028,0.0010351107,0.0010372134,0.0010402556,0.0010425759,0.0010448068,0.0010472694,0.0010500285,0.0010518006,0.0010542288,0.0010565496,0.0010579627,0.0010602050,0.0010622645,0.0010641857,0.0010659292,0.0010675027,0.0010694655,0.0010713500,0.0010725016,0.0010739903,0.0010755948,0.0010781591,0.0010783243,0.0010802148,0.0010818104,0.0010828805,0.0010838146,0.0010858695,0.0010867830,0.0010881898,0.0010891195,0.0010908965,0.0010912376,0.0010934387,0.0010945271,0.0010950755,0.0010968700,0.0010972297,0.0010991970,0.0010997011,0.0011005829,0.0011013836,0.0011025996,0.0011035947,0.0011039648,0.0011056739,0.0011059897,0.0011068374,0.0011077393,0.0011086672,0.0011099542,0.0011110705,0.0011114116,0.0011127475,0.0011131384,0.0011138701,0.0011149883,0.0011156771,0.0011164866,0.0011167121,0.0011178615,0.0011180656,0.0011185052,0.0011194152,0.0011196743,0.0011209286,0.0011210747,0.0011224550,0.0011227605,0.0011233176,0.0008179535,0.0008296279,0.0008408681,0.0008502041,0.0008597056,0.0008682508,0.0008756176,0.0008841433,0.0008904668,0.0008969426,0.0009031471,0.0009092634,0.0009147700,0.0009203381,0.0009254996,0.0009305306,0.0009349310,0.0009395649,0.0009428844,0.0009476390,0.0009508544,0.0009549195,0.0009588331,0.0009613964,0.0009648245,0.0009681827,0.0009707874,0.0009743156,0.0009772805,0.0009799333,0.0009825200,0.0009847730,0.0009875173,0.0009896714,0.0009918017,0.0009942345,0.0009968384,0.0009985591,0.0010008385,0.0010029816,0.0010044151,0.0010065787,0.0010085150,0.0010103490,0.0010119343,0.0010136080,0.0010153765,0.0010171177,0.0010182662,0.0010196866,0.0010213561,0.0010237136,0.0010238114,0.0010257134,0.0010271216,0.0010282592,0.0010291230,0.0010310573,0.0010318625,0.0010332177,0.0010340749,0.0010357797,0.0010360129,0.0010382328,0.0010393431,0.0010398247,0.0010415718,0.0010418861,0.0010437662,0.0010442402,0.0010451424,0.0010459440,0.0010470937,0.0010480583,0.0010482657,0.0010499751,0.0010500526,0.0010510000,0.0010519328,0.0010527179,0.0010539037,0.0010549613,0.0010554180,0.0010567257,0.0010569682,0.0010577369,0.0010588276,0.0010593753,0.0010602760,0.0010605041,0.0010615436,0.0010619184,0.0010622235,0.0010629242,0.0010633504,0.0010644543,0.0010646068,0.0010659804,0.0010663885,0.0010667973,0.0007760535,0.0007871615,0.0007977644,0.0008067567,0.0008157385,0.0008240272,0.0008308967,0.0008390630,0.0008450839,0.0008512976,0.0008571399,0.0008629159,0.0008681308,0.0008735544,0.0008783365,0.0008833321,0.0008872978,0.0008917774,0.0008950307,0.0008995630,0.0009024733,0.0009063127,0.0009101582,0.0009126060,0.0009159094,0.0009191051,0.0009215271,0.0009247886,0.0009277289,0.0009302884,0.0009327315,0.0009348225,0.0009373914,0.0009393456,0.0009415449,0.0009438793,0.0009464205,0.0009479473,0.0009502426,0.0009522843,0.0009535377,0.0009556257,0.0009575339,0.0009593029,0.0009608709,0.0009622597,0.0009640530,0.0009656876,0.0009666681,0.0009681098,0.0009697067,0.0009718001,0.0009720079,0.0009738898,0.0009752159,0.0009762624,0.0009772131,0.0009789608,0.0009797772,0.0009810156,0.0009818291,0.0009834953,0.0009837252,0.0009858504,0.0009868178,0.0009873708,0.0009889840,0.0009892606,0.0009911280,0.0009914894,0.0009923853,0.0009930768,0.0009942924,0.0009951543,0.0009954635,0.0009971122,0.0009969958,0.0009978270,0.0009989514,0.0009998280,0.0010008429,0.0010018143,0.0010021684,0.0010034731,0.0010036501,0.0010044359,0.0010053459,0.0010060743,0.0010069348,0.0010070205,0.0010081269,0.0010084507,0.0010087371,0.0010094531,0.0010099108,0.0010108451,0.0010110520,0.0010123882,0.0010127134,0.0010130621,0.0007363183,0.0007468913,0.0007570471,0.0007655473,0.0007741062,0.0007819741,0.0007883788,0.0007961017,0.0008020474,0.0008078571,0.0008134141,0.0008189453,0.0008240284,0.0008290170,0.0008337270,0.0008383352,0.0008421749,0.0008465059,0.0008495685,0.0008537528,0.0008566309,0.0008603249,0.0008638867,0.0008662400,0.0008694274,0.0008725262,0.0008747975,0.0008778829,0.0008807579,0.0008831548,0.0008855171,0.0008874125,0.0008899580,0.0008916838,0.0008937830,0.0008960244,0.0008983462,0.0008998474,0.0009021204,0.0009040977,0.0009053069,0.0009073177,0.0009091602,0.0009107229,0.0009122712,0.0009137585,0.0009152964,0.0009167761,0.0009178041,0.0009191751,0.0009206927,0.0009228248,0.0009230266,0.0009246972,0.0009260700,0.0009270467,0.0009278685,0.0009295477,0.0009304321,0.0009315516,0.0009321621,0.0009338317,0.0009340394,0.0009359709,0.0009369175,0.0009375063,0.0009391364,0.0009393172,0.0009410592,0.0009415154,0.0009423126,0.0009429670,0.0009441140,0.0009449650,0.0009452527,0.0009468092,0.0009468004,0.0009475342,0.0009485993,0.0009494046,0.0009503271,0.0009512668,0.0009517361,0.0009528730,0.0009531408,0.0009539090,0.0009546486,0.0009554351,0.0009563019,0.0009563743,0.0009574197,0.0009577496,0.0009579128,0.0009586667,0.0009591118,0.0009598837,0.0009600246,0.0009614608,0.0009618299,0.0009620571,0.0006985786,0.0007086712,0.0007183984,0.0007263685,0.0007346937,0.0007420546,0.0007481329,0.0007554751,0.0007611298,0.0007667870,0.0007719070,0.0007772290,0.0007820099,0.0007867588,0.0007912569,0.0007956780,0.0007993026,0.0008034816,0.0008063403,0.0008103610,0.0008131719,0.0008166080,0.0008199349,0.0008221957,0.0008253194,0.0008282835,0.0008304219,0.0008333717,0.0008361836,0.0008383156,0.0008405936,0.0008425533,0.0008447740,0.0008465333,0.0008486028,0.0008506737,0.0008529084,0.0008544263,0.0008564180,0.0008583767,0.0008595753,0.0008613559,0.0008632160,0.0008646309,0.0008661356,0.0008676267,0.0008690475,0.0008703027,0.0008714010,0.0008727563,0.0008742199,0.0008761194,0.0008763923,0.0008779628,0.0008792992,0.0008802712,0.0008810120,0.0008826619,0.0008834699,0.0008845045,0.0008851621,0.0008867828,0.0008868152,0.0008887982,0.0008897516,0.0008902101,0.0008918375,0.0008920012,0.0008936363,0.0008939587,0.0008948001,0.0008954806,0.0008966170,0.0008972868,0.0008976607,0.0008992290,0.0008991144,0.0008999186,0.0009008066,0.0009016073,0.0009022956,0.0009033517,0.0009037957,0.0009049001,0.0009050541,0.0009059609,0.0009066176,0.0009073448,0.0009081503,0.0009081440,0.0009091989,0.0009095115,0.0009097125,0.0009104300,0.0009108608,0.0009114816,0.0009117085,0.0009131804,0.0009136098,0.0009136938,0.0006628476,0.0006723457,0.0006817156,0.0006892794,0.0006972099,0.0007042194,0.0007099957,0.0007169169,0.0007223940,0.0007276553,0.0007325166,0.0007376879,0.0007421702,0.0007468121,0.0007510275,0.0007552178,0.0007586262,0.0007626584,0.0007653764,0.0007692804,0.0007718740,0.0007749694,0.0007783896,0.0007805166,0.0007834087,0.0007861977,0.0007884329,0.0007912176,0.0007937997,0.0007958626,0.0007979460,0.0007998886,0.0008019503,0.0008036988,0.0008057025,0.0008077014,0.0008097879,0.0008112051,0.0008129881,0.0008149932,0.0008162055,0.0008178357,0.0008196559,0.0008209948,0.0008222609,0.0008237017,0.0008251436,0.0008262242,0.0008273017,0.0008287384,0.0008301020,0.0008318582,0.0008320994,0.0008336035,0.0008349856,0.0008359615,0.0008365361,0.0008381428,0.0008389296,0.0008399564,0.0008405505,0.0008420363,0.0008420214,0.0008439463,0.0008449286,0.0008452762,0.0008469491,0.0008470491,0.0008485475,0.0008488916,0.0008496021,0.0008504545,0.0008514066,0.0008519726,0.0008523921,0.0008539788,0.0008537823,0.0008546269,0.0008554933,0.0008562932,0.0008567647,0.0008577675,0.0008583146,0.0008593780,0.0008594984,0.0008604111,0.0008608188,0.0008616320,0.0008624732,0.0008625194,0.0008634294,0.0008637888,0.0008638971,0.0008647198,0.0008649389,0.0008656416,0.0008658409,0.0008672843,0.0008677569,0.0008677310,0.0006289306,0.0006379977,0.0006469518,0.0006540149,0.0006616307,0.0006683973,0.0006738553,0.0006803880,0.0006855394,0.0006906800,0.0006951989,0.0007001423,0.0007044268,0.0007087727,0.0007129182,0.0007169354,0.0007201287,0.0007240339,0.0007264949,0.0007303249,0.0007327205,0.0007355210,0.0007389206,0.0007410069,0.0007436143,0.0007462947,0.0007484718,0.0007510380,0.0007535323,0.0007554829,0.0007576399,0.0007594004,0.0007614475,0.0007630461,0.0007650140,0.0007667804,0.0007687706,0.0007702056,0.0007718564,0.0007738196,0.0007749200,0.0007763972,0.0007781793,0.0007795880,0.0007806570,0.0007821516,0.0007833838,0.0007843872,0.0007853917,0.0007868321,0.0007881184,0.0007899453,0.0007900719,0.0007914449,0.0007929580,0.0007937081,0.0007943259,0.0007957560,0.0007966059,0.0007976525,0.0007981681,0.0007994917,0.0007995396,0.0008014314,0.0008023216,0.0008026505,0.0008041494,0.0008042919,0.0008058292,0.0008060899,0.0008068323,0.0008075813,0.0008085263,0.0008090597,0.0008094567,0.0008110028,0.0008107496,0.0008116641,0.0008124102,0.0008131839,0.0008136093,0.0008145681,0.0008149989,0.0008160822,0.0008161375,0.0008170155,0.0008175371,0.0008182544,0.0008191121,0.0008191112,0.0008199316,0.0008203358,0.0008203655,0.0008212504,0.0008214233,0.0008220417,0.0008222998,0.0008235933,0.0008240812,0.0008240641,0.0005967491,0.0006054210,0.0006138781,0.0006205223,0.0006278888,0.0006342875,0.0006395004,0.0006456865,0.0006505486,0.0006555093,0.0006597828,0.0006645291,0.0006685815,0.0006728152,0.0006766471,0.0006805103,0.0006834523,0.0006873379,0.0006896243,0.0006932448,0.0006955437,0.0006982179,0.0007013774,0.0007033943,0.0007058473,0.0007085253,0.0007105036,0.0007130863,0.0007153469,0.0007173059,0.0007192880,0.0007210150,0.0007229682,0.0007243181,0.0007263043,0.0007280917,0.0007298719,0.0007313206,0.0007328123,0.0007347223,0.0007357383,0.0007371737,0.0007388944,0.0007402551,0.0007411425,0.0007425600,0.0007438566,0.0007448348,0.0007457456,0.0007470516,0.0007482715,0.0007499720,0.0007500850,0.0007514640,0.0007529107,0.0007535999,0.0007543308,0.0007556175,0.0007564038,0.0007573638,0.0007578628,0.0007591796,0.0007593177,0.0007610465,0.0007618937,0.0007621860,0.0007636918,0.0007636979,0.0007652579,0.0007654680,0.0007661157,0.0007669062,0.0007677826,0.0007683669,0.0007688019,0.0007701076,0.0007699097,0.0007708633,0.0007714515,0.0007722603,0.0007726238,0.0007735230,0.0007739797,0.0007749574,0.0007750620,0.0007759396,0.0007764249,0.0007770182,0.0007779445,0.0007779463,0.0007787642,0.0007790876,0.0007791506,0.0007800529,0.0007801707,0.0007807174,0.0007809997,0.0007821179,0.0007826852,0.0007826813,0.0005663012,0.0005745209,0.0005824806,0.0005888696,0.0005958184,0.0006019046,0.0006069674,0.0006128247,0.0006174165,0.0006220822,0.0006261544,0.0006305781,0.0006346509,0.0006385904,0.0006422614,0.0006459703,0.0006486842,0.0006524502,0.0006545978,0.0006580682,0.0006602803,0.0006628647,0.0006657053,0.0006676775,0.0006700693,0.0006726108,0.0006744026,0.0006769681,0.0006790785,0.0006809560,0.0006828530,0.0006845437,0.0006864543,0.0006875962,0.0006895513,0.0006912749,0.0006929529,0.0006943612,0.0006957138,0.0006975145,0.0006984809,0.0006999870,0.0007015219,0.0007028298,0.0007037559,0.0007050361,0.0007062575,0.0007071439,0.0007081075,0.0007092784,0.0007103274,0.0007121397,0.0007121379,0.0007135868,0.0007148884,0.0007156006,0.0007162077,0.0007175988,0.0007181959,0.0007192012,0.0007196385,0.0007209541,0.0007210472,0.0007226617,0.0007235329,0.0007238238,0.0007252196,0.0007251793,0.0007267005,0.0007269070,0.0007274884,0.0007281824,0.0007291584,0.0007296901,0.0007300794,0.0007313738,0.0007312274,0.0007320537,0.0007325791,0.0007334188,0.0007337930,0.0007345461,0.0007350293,0.0007359836,0.0007360711,0.0007369150,0.0007373180,0.0007378538,0.0007388400,0.0007387287,0.0007395739,0.0007398671,0.0007400045,0.0007407689,0.0007409622,0.0007413688,0.0007417063,0.0007428546,0.0007434026,0.0007433393,0.0005374208,0.0005451443,0.0005528612,0.0005588812,0.0005653266,0.0005712135,0.0005760321,0.0005816414,0.0005859206,0.0005904405,0.0005943105,0.0005985733,0.0006023339,0.0006061301,0.0006096876,0.0006131600,0.0006156864,0.0006193792,0.0006214073,0.0006246243,0.0006268037,0.0006292854,0.0006320051,0.0006337596,0.0006360805,0.0006385004,0.0006401727,0.0006426795,0.0006446512,0.0006464398,0.0006483510,0.0006498679,0.0006517564,0.0006527509,0.0006546331,0.0006562677,0.0006579889,0.0006592680,0.0006606309,0.0006622186,0.0006632760,0.0006646330,0.0006660762,0.0006673505,0.0006681691,0.0006693674,0.0006705481,0.0006714202,0.0006723247,0.0006735021,0.0006744384,0.0006762392,0.0006761127,0.0006776510,0.0006788213,0.0006794563,0.0006801135,0.0006813787,0.0006820283,0.0006828897,0.0006833525,0.0006846518,0.0006848384,0.0006862421,0.0006871090,0.0006873580,0.0006887438,0.0006886183,0.0006901648,0.0006902341,0.0006908828,0.0006915446,0.0006924313,0.0006930216,0.0006933389,0.0006945785,0.0006944478,0.0006952706,0.0006956060,0.0006964842,0.0006968429,0.0006975106,0.0006980309,0.0006988713,0.0006991078,0.0006999058,0.0007002351,0.0007006866,0.0007016271,0.0007015451,0.0007023364,0.0007027078,0.0007027243,0.0007035078,0.0007037557,0.0007041520,0.0007044338,0.0007056355,0.0007060814,0.0007060340,0.0005099822,0.0005173039,0.0005246970,0.0005303683,0.0005363778,0.0005421701,0.0005466930,0.0005520330,0.0005560994,0.0005604235,0.0005640754,0.0005681149,0.0005717032,0.0005754023,0.0005786322,0.0005820526,0.0005844698,0.0005878706,0.0005899329,0.0005929009,0.0005950193,0.0005974454,0.0006000231,0.0006016991,0.0006038507,0.0006061349,0.0006077743,0.0006102216,0.0006120036,0.0006137925,0.0006154730,0.0006170205,0.0006187920,0.0006197153,0.0006214487,0.0006231294,0.0006247748,0.0006259965,0.0006272168,0.0006288287,0.0006297354,0.0006310642,0.0006323584,0.0006336347,0.0006344820,0.0006356503,0.0006366337,0.0006375171,0.0006383604,0.0006394492,0.0006403791,0.0006421439,0.0006420792,0.0006435821,0.0006445509,0.0006451890,0.0006457812,0.0006470796,0.0006476175,0.0006484270,0.0006488975,0.0006500707,0.0006503025,0.0006516595,0.0006525926,0.0006527718,0.0006540161,0.0006539163,0.0006553673,0.0006554023,0.0006560681,0.0006567221,0.0006575998,0.0006581657,0.0006584249,0.0006596445,0.0006594522,0.0006602210,0.0006606191,0.0006613796,0.0006617733,0.0006623939,0.0006629781,0.0006637286,0.0006640588,0.0006647097,0.0006649743,0.0006653901,0.0006663269,0.0006662847,0.0006670538,0.0006674279,0.0006674973,0.0006681795,0.0006684038,0.0006687625,0.0006690528,0.0006701714,0.0006705418,0.0006705434,0.0004839957,0.0004909498,0.0004979205,0.0005033319,0.0005090412,0.0005145476,0.0005188382,0.0005238790,0.0005278258,0.0005319764,0.0005353904,0.0005392443,0.0005426820,0.0005460882,0.0005492365,0.0005524550,0.0005548655,0.0005579893,0.0005599888,0.0005629166,0.0005648570,0.0005671487,0.0005696695,0.0005712129,0.0005731594,0.0005755095,0.0005770693,0.0005793305,0.0005809554,0.0005827224,0.0005842745,0.0005857651,0.0005874623,0.0005884674,0.0005899753,0.0005915414,0.0005931979,0.0005943440,0.0005955021,0.0005969958,0.0005979214,0.0005991430,0.0006003564,0.0006015754,0.0006024360,0.0006035243,0.0006044097,0.0006053724,0.0006060670,0.0006072588,0.0006081247,0.0006098033,0.0006096821,0.0006112019,0.0006120142,0.0006127043,0.0006132283,0.0006145098,0.0006149502,0.0006155919,0.0006161734,0.0006173662,0.0006175666,0.0006187957,0.0006197356,0.0006199382,0.0006210541,0.0006209948,0.0006223193,0.0006223696,0.0006230017,0.0006237071,0.0006245412,0.0006250320,0.0006252448,0.0006264278,0.0006263166,0.0006270002,0.0006273592,0.0006281716,0.0006284805,0.0006290891,0.0006296244,0.0006302918,0.0006306637,0.0006312930,0.0006316029,0.0006318992,0.0006329228,0.0006328397,0.0006335474,0.0006338797,0.0006339240,0.0006345938,0.0006348499,0.0006351027,0.0006354267,0.0006364752,0.0006369356,0.0006368786,0.0004593414,0.0004658877,0.0004725854,0.0004776547,0.0004830466,0.0004884108,0.0004924275,0.0004972116,0.0005010069,0.0005049878,0.0005082614,0.0005117701,0.0005150789,0.0005183185,0.0005213071,0.0005244167,0.0005266309,0.0005297092,0.0005315959,0.0005343211,0.0005362028,0.0005384547,0.0005407082,0.0005422852,0.0005441343,0.0005464329,0.0005478834,0.0005499115,0.0005514983,0.0005532235,0.0005547905,0.0005561177,0.0005577913,0.0005587378,0.0005601680,0.0005616025,0.0005632060,0.0005643002,0.0005654911,0.0005668881,0.0005676449,0.0005690764,0.0005700675,0.0005711773,0.0005720150,0.0005731074,0.0005738886,0.0005748311,0.0005754845,0.0005765575,0.0005773784,0.0005790464,0.0005789323,0.0005803758,0.0005810965,0.0005818700,0.0005822817,0.0005835465,0.0005839361,0.0005844513,0.0005851663,0.0005861983,0.0005864811,0.0005877041,0.0005884952,0.0005886793,0.0005897739,0.0005897191,0.0005910130,0.0005910517,0.0005916592,0.0005923181,0.0005930686,0.0005935556,0.0005937817,0.0005948396,0.0005948100,0.0005953953,0.0005958227,0.0005965605,0.0005968839,0.0005974784,0.0005980549,0.0005986539,0.0005990193,0.0005996302,0.0005998037,0.0006001562,0.0006010800,0.0006010350,0.0006017074,0.0006019516,0.0006021108,0.0006026238,0.0006029002,0.0006033249,0.0006035448,0.0006044653,0.0006049160,0.0006049697,0.0004358673,0.0004420490,0.0004484982,0.0004532865,0.0004583688,0.0004634815,0.0004673503,0.0004719851,0.0004755679,0.0004793192,0.0004823842,0.0004857804,0.0004888972,0.0004919570,0.0004948429,0.0004978362,0.0004999767,0.0005028049,0.0005045898,0.0005071793,0.0005089436,0.0005111378,0.0005132640,0.0005147001,0.0005165401,0.0005187794,0.0005201667,0.0005219948,0.0005234980,0.0005252030,0.0005267685,0.0005280354,0.0005295410,0.0005304264,0.0005319140,0.0005332551,0.0005347378,0.0005357670,0.0005369845,0.0005382642,0.0005391049,0.0005403277,0.0005411489,0.0005424057,0.0005431495,0.0005441237,0.0005449455,0.0005458019,0.0005464658,0.0005474900,0.0005483023,0.0005498822,0.0005496318,0.0005511031,0.0005518054,0.0005525377,0.0005529302,0.0005541533,0.0005545483,0.0005549657,0.0005557077,0.0005566213,0.0005569241,0.0005580559,0.0005588441,0.0005590139,0.0005600523,0.0005600308,0.0005613548,0.0005613990,0.0005618007,0.0005625388,0.0005632218,0.0005637105,0.0005638265,0.0005649144,0.0005649298,0.0005654125,0.0005658900,0.0005665265,0.0005668743,0.0005674479,0.0005679629,0.0005685964,0.0005688892,0.0005694994,0.0005697156,0.0005699965,0.0005708742,0.0005708539,0.0005714839,0.0005717874,0.0005718700,0.0005723539,0.0005725372,0.0005731088,0.0005732442,0.0005740902,0.0005745648,0.0005746506,0.0004136272,0.0004194837,0.0004256344,0.0004301265,0.0004349674,0.0004399591,0.0004435059,0.0004479167,0.0004514056,0.0004549835,0.0004579614,0.0004609240,0.0004640544,0.0004669842,0.0004697761,0.0004727096,0.0004745480,0.0004772738,0.0004789317,0.0004815122,0.0004831197,0.0004852060,0.0004872690,0.0004886127,0.0004904398,0.0004925315,0.0004938429,0.0004955627,0.0004969055,0.0004986956,0.0005001762,0.0005014521,0.0005028548,0.0005035587,0.0005050782,0.0005063683,0.0005077347,0.0005086149,0.0005098047,0.0005110293,0.0005118693,0.0005130446,0.0005138702,0.0005150597,0.0005156975,0.0005166194,0.0005174523,0.0005182419,0.0005189462,0.0005199488,0.0005206188,0.0005222963,0.0005219693,0.0005232870,0.0005240066,0.0005247233,0.0005250647,0.0005262429,0.0005265984,0.0005269804,0.0005277369,0.0005285824,0.0005288414,0.0005299673,0.0005307491,0.0005308250,0.0005318599,0.0005318163,0.0005330552,0.0005331099,0.0005335486,0.0005341540,0.0005349021,0.0005354038,0.0005354834,0.0005365025,0.0005365207,0.0005368976,0.0005374327,0.0005380374,0.0005384397,0.0005389694,0.0005393594,0.0005400374,0.0005403499,0.0005408529,0.0005409732,0.0005413106,0.0005421954,0.0005421639,0.0005427613,0.0005429936,0.0005431526,0.0005436368,0.0005436862,0.0005443440,0.0005444487,0.0005452826,0.0005456969,0.0005458427,0.0003924837,0.0003980706,0.0004040324,0.0004081999,0.0004128385,0.0004175214,0.0004208388,0.0004251767,0.0004283537,0.0004318466,0.0004347221,0.0004374237,0.0004404783,0.0004433700,0.0004458477,0.0004487196,0.0004505227,0.0004530474,0.0004547399,0.0004570670,0.0004587041,0.0004606670,0.0004625931,0.0004639779,0.0004657464,0.0004675565,0.0004689470,0.0004704525,0.0004718288,0.0004734407,0.0004748921,0.0004762212,0.0004774633,0.0004781476,0.0004795503,0.0004807324,0.0004820374,0.0004829577,0.0004841702,0.0004852207,0.0004859862,0.0004871938,0.0004879384,0.0004890759,0.0004897188,0.0004904835,0.0004913637,0.0004920856,0.0004927414,0.0004936908,0.0004943966,0.0004958507,0.0004956940,0.0004968972,0.0004975052,0.0004982954,0.0004985692,0.0004997317,0.0005000719,0.0005004549,0.0005011766,0.0005018928,0.0005022240,0.0005032932,0.0005040205,0.0005041660,0.0005050367,0.0005051358,0.0005062391,0.0005062980,0.0005069099,0.0005072702,0.0005080617,0.0005084831,0.0005085284,0.0005095045,0.0005095834,0.0005098699,0.0005103663,0.0005109926,0.0005113711,0.0005119650,0.0005123145,0.0005129273,0.0005133082,0.0005137306,0.0005137683,0.0005140906,0.0005150318,0.0005149202,0.0005155316,0.0005158222,0.0005158106,0.0005162614,0.0005163470,0.0005169823,0.0005171520,0.0005180357,0.0005183891,0.0005184495,0.0003724272,0.0003778023,0.0003834154,0.0003874241,0.0003918372,0.0003962133,0.0003994457,0.0004035064,0.0004066114,0.0004098853,0.0004126460,0.0004152302,0.0004181915,0.0004208316,0.0004232183,0.0004259745,0.0004276265,0.0004301015,0.0004316920,0.0004338719,0.0004354714,0.0004373870,0.0004391706,0.0004404097,0.0004421540,0.0004438085,0.0004452336,0.0004466484,0.0004479693,0.0004494657,0.0004508846,0.0004521296,0.0004533547,0.0004540078,0.0004553239,0.0004564513,0.0004576488,0.0004585467,0.0004597654,0.0004607278,0.0004614870,0.0004626500,0.0004633277,0.0004644287,0.0004650328,0.0004657309,0.0004665421,0.0004672811,0.0004679270,0.0004688096,0.0004694345,0.0004708920,0.0004706901,0.0004719055,0.0004724353,0.0004731911,0.0004734685,0.0004745993,0.0004749028,0.0004753360,0.0004759864,0.0004766063,0.0004769145,0.0004780144,0.0004786839,0.0004788627,0.0004796449,0.0004796686,0.0004807356,0.0004809115,0.0004813622,0.0004818206,0.0004824411,0.0004828968,0.0004829236,0.0004839570,0.0004839327,0.0004842548,0.0004847648,0.0004853290,0.0004857322,0.0004862298,0.0004865984,0.0004871477,0.0004875694,0.0004879196,0.0004879603,0.0004882573,0.0004891524,0.0004890984,0.0004896655,0.0004898938,0.0004898090,0.0004903815,0.0004903399,0.0004909597,0.0004912632,0.0004920301,0.0004923943,0.0004924112,0.0003534513,0.0003584928,0.0003638198,0.0003676477,0.0003718988,0.0003760104,0.0003791448,0.0003829775,0.0003859090,0.0003889853,0.0003916923,0.0003942334,0.0003969143,0.0003994737,0.0004016971,0.0004043790,0.0004059539,0.0004083091,0.0004098172,0.0004119711,0.0004134094,0.0004152153,0.0004169016,0.0004181114,0.0004198471,0.0004214144,0.0004226966,0.0004240790,0.0004253228,0.0004268539,0.0004281920,0.0004292682,0.0004303773,0.0004311939,0.0004322654,0.0004334110,0.0004344866,0.0004353650,0.0004366233,0.0004375081,0.0004381234,0.0004392254,0.0004399005,0.0004410270,0.0004416097,0.0004422140,0.0004430523,0.0004438370,0.0004443366,0.0004452084,0.0004458684,0.0004472058,0.0004469637,0.0004481895,0.0004486127,0.0004492982,0.0004496426,0.0004506667,0.0004510476,0.0004514478,0.0004519845,0.0004525454,0.0004529174,0.0004539844,0.0004546593,0.0004547930,0.0004554738,0.0004554899,0.0004565097,0.0004567742,0.0004571576,0.0004576003,0.0004581754,0.0004586262,0.0004587726,0.0004596264,0.0004596081,0.0004598891,0.0004604400,0.0004609897,0.0004613014,0.0004619393,0.0004621827,0.0004627089,0.0004630728,0.0004633673,0.0004634141,0.0004637263,0.0004645574,0.0004646201,0.0004650112,0.0004652568,0.0004652032,0.0004658350,0.0004656853,0.0004662426,0.0004665669,0.0004672456,0.0004676365,0.0004676345,0.0003354526,0.0003402042,0.0003452650,0.0003489742,0.0003529182,0.0003568543,0.0003599076,0.0003635093,0.0003662544,0.0003691945,0.0003718687,0.0003741750,0.0003768138,0.0003792080,0.0003813607,0.0003838909,0.0003853017,0.0003876508,0.0003890495,0.0003911394,0.0003924781,0.0003942081,0.0003957249,0.0003969405,0.0003986070,0.0004000812,0.0004013471,0.0004026175,0.0004038511,0.0004052104,0.0004065980,0.0004075899,0.0004086226,0.0004094130,0.0004104368,0.0004115869,0.0004125271,0.0004133409,0.0004145946,0.0004154185,0.0004160140,0.0004170584,0.0004176925,0.0004186932,0.0004193277,0.0004198624,0.0004206746,0.0004214154,0.0004219093,0.0004227446,0.0004235051,0.0004247282,0.0004243928,0.0004257104,0.0004260172,0.0004267675,0.0004269865,0.0004279878,0.0004283272,0.0004287168,0.0004292079,0.0004297850,0.0004302314,0.0004310865,0.0004318016,0.0004319330,0.0004326126,0.0004325682,0.0004335383,0.0004337842,0.0004341677,0.0004346394,0.0004351023,0.0004356165,0.0004357522,0.0004365048,0.0004365775,0.0004369133,0.0004373622,0.0004378829,0.0004381160,0.0004387621,0.0004389231,0.0004394756,0.0004397518,0.0004400949,0.0004401696,0.0004404049,0.0004412246,0.0004412687,0.0004415833,0.0004418832,0.0004418703,0.0004424939,0.0004422771,0.0004428422,0.0004432012,0.0004437558,0.0004441583,0.0004440856,0.0003183190,0.0003228269,0.0003277372,0.0003311670,0.0003348977,0.0003387239,0.0003415670,0.0003449910,0.0003476739,0.0003504471,0.0003529870,0.0003552016,0.0003577761,0.0003600633,0.0003620626,0.0003644419,0.0003657682,0.0003679600,0.0003692926,0.0003713748,0.0003725768,0.0003742808,0.0003756857,0.0003768893,0.0003783685,0.0003798197,0.0003810768,0.0003822898,0.0003834814,0.0003847376,0.0003860802,0.0003870325,0.0003880209,0.0003887363,0.0003897037,0.0003908036,0.0003917204,0.0003925299,0.0003936670,0.0003943383,0.0003950301,0.0003960269,0.0003965680,0.0003975894,0.0003981816,0.0003986390,0.0003994773,0.0004002164,0.0004006670,0.0004014628,0.0004020882,0.0004033267,0.0004029681,0.0004041475,0.0004045373,0.0004052457,0.0004055304,0.0004064527,0.0004068157,0.0004071867,0.0004076456,0.0004081529,0.0004086688,0.0004094043,0.0004101109,0.0004102292,0.0004108501,0.0004108788,0.0004117722,0.0004119802,0.0004123955,0.0004128320,0.0004132109,0.0004136453,0.0004139276,0.0004145752,0.0004146762,0.0004150123,0.0004153250,0.0004158996,0.0004161049,0.0004168493,0.0004169034,0.0004173804,0.0004177000,0.0004180568,0.0004181459,0.0004183142,0.0004190816,0.0004191357,0.0004194182,0.0004198255,0.0004197057,0.0004202013,0.0004201541,0.0004206312,0.0004209512,0.0004214694,0.0004218614,0.0004218284,0.0003020620,0.0003063556,0.0003110346,0.0003143474,0.0003178260,0.0003214862,0.0003242081,0.0003275341,0.0003300271,0.0003326528,0.0003350437,0.0003372149,0.0003396282,0.0003418482,0.0003437371,0.0003459824,0.0003472455,0.0003493259,0.0003506617,0.0003525600,0.0003537470,0.0003554293,0.0003567410,0.0003577224,0.0003593097,0.0003606331,0.0003618043,0.0003629808,0.0003640937,0.0003653084,0.0003665892,0.0003674924,0.0003685042,0.0003690874,0.0003700447,0.0003710383,0.0003719514,0.0003726935,0.0003738827,0.0003744623,0.0003751148,0.0003760960,0.0003766144,0.0003775565,0.0003781060,0.0003786102,0.0003793094,0.0003800382,0.0003805798,0.0003812797,0.0003818645,0.0003829824,0.0003826526,0.0003838085,0.0003841085,0.0003849011,0.0003851971,0.0003859355,0.0003863678,0.0003866708,0.0003871738,0.0003876389,0.0003881572,0.0003888116,0.0003894935,0.0003896158,0.0003901075,0.0003902155,0.0003910265,0.0003912554,0.0003915996,0.0003920234,0.0003924677,0.0003928553,0.0003932117,0.0003936998,0.0003938400,0.0003941993,0.0003944991,0.0003950509,0.0003952654,0.0003958983,0.0003959502,0.0003963838,0.0003967429,0.0003970508,0.0003971811,0.0003973101,0.0003981053,0.0003981401,0.0003982986,0.0003987552,0.0003986403,0.0003991445,0.0003991019,0.0003995584,0.0003997815,0.0004002916,0.0004006468,0.0004006217,0.0002865862,0.0002907279,0.0002952069,0.0002983655,0.0003016853,0.0003051067,0.0003076709,0.0003108930,0.0003132512,0.0003158050,0.0003180066,0.0003201583,0.0003224255,0.0003245129,0.0003263308,0.0003285340,0.0003296975,0.0003318010,0.0003330401,0.0003347377,0.0003358690,0.0003374311,0.0003387103,0.0003396494,0.0003410787,0.0003423111,0.0003435310,0.0003446086,0.0003457090,0.0003468395,0.0003480625,0.0003489626,0.0003499407,0.0003504434,0.0003513908,0.0003523318,0.0003531660,0.0003539534,0.0003550141,0.0003556098,0.0003562324,0.0003572039,0.0003575916,0.0003586314,0.0003590604,0.0003594973,0.0003602997,0.0003608993,0.0003613828,0.0003621020,0.0003626260,0.0003636499,0.0003633306,0.0003645804,0.0003648140,0.0003655313,0.0003658023,0.0003664461,0.0003669233,0.0003672141,0.0003676480,0.0003680487,0.0003686361,0.0003692577,0.0003698426,0.0003699737,0.0003704869,0.0003706171,0.0003713217,0.0003715927,0.0003718873,0.0003723499,0.0003726711,0.0003730878,0.0003735429,0.0003739426,0.0003740381,0.0003743843,0.0003746854,0.0003751935,0.0003753921,0.0003760266,0.0003760508,0.0003764794,0.0003768242,0.0003771306,0.0003772050,0.0003773601,0.0003781804,0.0003780985,0.0003783087,0.0003787573,0.0003786313,0.0003791228,0.0003791270,0.0003795425,0.0003797181,0.0003802041,0.0003806236,0.0003805255,0.0002719684,0.0002759893,0.0002801631,0.0002831980,0.0002863328,0.0002896233,0.0002919816,0.0002950874,0.0002973222,0.0002997332,0.0003018846,0.0003039485,0.0003060190,0.0003080680,0.0003098472,0.0003118580,0.0003130321,0.0003149093,0.0003162195,0.0003178604,0.0003188653,0.0003203371,0.0003215577,0.0003224292,0.0003238405,0.0003250417,0.0003262318,0.0003271859,0.0003282744,0.0003293396,0.0003305553,0.0003313214,0.0003323395,0.0003328086,0.0003336529,0.0003345556,0.0003353495,0.0003360465,0.0003371974,0.0003376762,0.0003382533,0.0003391752,0.0003395615,0.0003405243,0.0003409648,0.0003413269,0.0003421759,0.0003426958,0.0003431474,0.0003438783,0.0003443909,0.0003453430,0.0003450187,0.0003461554,0.0003464596,0.0003471564,0.0003473980,0.0003480183,0.0003484416,0.0003487168,0.0003491279,0.0003494986,0.0003500914,0.0003506253,0.0003512774,0.0003513901,0.0003519125,0.0003519994,0.0003527039,0.0003528880,0.0003531938,0.0003536621,0.0003539555,0.0003544164,0.0003547867,0.0003552102,0.0003552804,0.0003556232,0.0003558574,0.0003564002,0.0003566295,0.0003571275,0.0003571370,0.0003575592,0.0003579112,0.0003581792,0.0003582889,0.0003584035,0.0003592566,0.0003591667,0.0003592592,0.0003597389,0.0003596483,0.0003601147,0.0003601391,0.0003605002,0.0003606634,0.0003611553,0.0003615463,0.0003614488,0.0002581007,0.0002619019,0.0002659620,0.0002687923,0.0002717312,0.0002748992,0.0002771587,0.0002801050,0.0002822012,0.0002845852,0.0002866041,0.0002884857,0.0002905634,0.0002924784,0.0002941928,0.0002960722,0.0002971817,0.0002989596,0.0003002224,0.0003018238,0.0003027819,0.0003041562,0.0003052664,0.0003061351,0.0003074561,0.0003086572,0.0003097591,0.0003106259,0.0003116676,0.0003126945,0.0003138695,0.0003145751,0.0003156091,0.0003160163,0.0003167710,0.0003176989,0.0003184113,0.0003190743,0.0003202296,0.0003206243,0.0003211881,0.0003220587,0.0003224788,0.0003233612,0.0003237950,0.0003240715,0.0003249600,0.0003253723,0.0003258619,0.0003265671,0.0003270520,0.0003280393,0.0003276361,0.0003287571,0.0003290794,0.0003297028,0.0003299272,0.0003305034,0.0003309458,0.0003311733,0.0003315955,0.0003319245,0.0003325503,0.0003330040,0.0003336204,0.0003337023,0.0003343005,0.0003342433,0.0003349887,0.0003352182,0.0003355320,0.0003358945,0.0003362591,0.0003366550,0.0003369597,0.0003373748,0.0003374409,0.0003377238,0.0003379600,0.0003385317,0.0003386963,0.0003391417,0.0003391692,0.0003396673,0.0003399628,0.0003402520,0.0003402855,0.0003404262,0.0003412501,0.0003412145,0.0003413086,0.0003417079,0.0003415863,0.0003420362,0.0003421401,0.0003424462,0.0003426157,0.0003430140,0.0003434341,0.0003433395,0.0002449700,0.0002485981,0.0002524638,0.0002550623,0.0002579164,0.0002609102,0.0002630609,0.0002658406,0.0002679154,0.0002700902,0.0002720310,0.0002738914,0.0002758269,0.0002776964,0.0002793333,0.0002811112,0.0002821125,0.0002838279,0.0002850010,0.0002866202,0.0002874788,0.0002888068,0.0002898118,0.0002906887,0.0002919081,0.0002930965,0.0002940520,0.0002949033,0.0002958732,0.0002968753,0.0002980069,0.0002987600,0.0002996790,0.0003000976,0.0003007632,0.0003016734,0.0003022700,0.0003030088,0.0003040567,0.0003045085,0.0003050516,0.0003057942,0.0003062768,0.0003069970,0.0003074391,0.0003077871,0.0003085938,0.0003089564,0.0003094089,0.0003101669,0.0003106006,0.0003115013,0.0003111604,0.0003122179,0.0003125041,0.0003131435,0.0003133690,0.0003139146,0.0003143595,0.0003144354,0.0003149823,0.0003152693,0.0003158086,0.0003162236,0.0003168605,0.0003169168,0.0003175395,0.0003174333,0.0003181464,0.0003184157,0.0003186732,0.0003190865,0.0003193771,0.0003197309,0.0003200165,0.0003204788,0.0003205230,0.0003207673,0.0003209375,0.0003215684,0.0003217358,0.0003220898,0.0003220916,0.0003226279,0.0003229125,0.0003231430,0.0003232488,0.0003233265,0.0003241111,0.0003240716,0.0003241425,0.0003245543,0.0003244719,0.0003248841,0.0003249815,0.0003252740,0.0003254564,0.0003259110,0.0003262010,0.0003260492,0.0002325333,0.0002359226,0.0002396306,0.0002421421,0.0002448734,0.0002476248,0.0002496481,0.0002523901,0.0002543007,0.0002563683,0.0002582004,0.0002600339,0.0002618935,0.0002636008,0.0002652068,0.0002668895,0.0002678333,0.0002694503,0.0002705999,0.0002721286,0.0002729610,0.0002742025,0.0002751762,0.0002760113,0.0002772296,0.0002783236,0.0002791693,0.0002799259,0.0002808745,0.0002818761,0.0002829629,0.0002836405,0.0002845407,0.0002849894,0.0002856020,0.0002864437,0.0002869891,0.0002877753,0.0002887846,0.0002891153,0.0002896744,0.0002903628,0.0002908276,0.0002914843,0.0002919464,0.0002922862,0.0002930718,0.0002934166,0.0002938560,0.0002945650,0.0002949172,0.0002958183,0.0002954841,0.0002964995,0.0002967480,0.0002974303,0.0002975925,0.0002981653,0.0002985682,0.0002986380,0.0002992064,0.0002994257,0.0002999749,0.0003003040,0.0003009736,0.0003009985,0.0003015538,0.0003014289,0.0003021663,0.0003024243,0.0003026554,0.0003031909,0.0003033823,0.0003037124,0.0003039537,0.0003044032,0.0003044840,0.0003046921,0.0003047762,0.0003054569,0.0003055760,0.0003058149,0.0003059812,0.0003064353,0.0003066967,0.0003068977,0.0003069893,0.0003070764,0.0003078779,0.0003077742,0.0003079992,0.0003082900,0.0003082223,0.0003085739,0.0003087225,0.0003089975,0.0003090736,0.0003096007,0.0003097962,0.0003096806,0.0002207356,0.0002239307,0.0002274253,0.0002298392,0.0002324824,0.0002350472,0.0002369020,0.0002395358,0.0002413840,0.0002433925,0.0002451083,0.0002468414,0.0002485963,0.0002502983,0.0002517600,0.0002534155,0.0002542859,0.0002558211,0.0002568932,0.0002583372,0.0002591350,0.0002603307,0.0002611771,0.0002620745,0.0002632662,0.0002642336,0.0002650685,0.0002658145,0.0002666871,0.0002676468,0.0002686834,0.0002693512,0.0002701912,0.0002706371,0.0002712203,0.0002719811,0.0002725295,0.0002732214,0.0002742222,0.0002745454,0.0002750519,0.0002756953,0.0002761109,0.0002767740,0.0002772151,0.0002775972,0.0002783529,0.0002786696,0.0002791126,0.0002797495,0.0002801570,0.0002808972,0.0002806351,0.0002816196,0.0002818371,0.0002824239,0.0002826888,0.0002831400,0.0002835817,0.0002836076,0.0002842307,0.0002843286,0.0002848946,0.0002852009,0.0002859084,0.0002858056,0.0002864188,0.0002862876,0.0002869409,0.0002872377,0.0002874416,0.0002879733,0.0002881572,0.0002884273,0.0002886502,0.0002891600,0.0002892186,0.0002894116,0.0002894397,0.0002901586,0.0002901694,0.0002905059,0.0002906245,0.0002910251,0.0002913222,0.0002915671,0.0002916383,0.0002916200,0.0002924695,0.0002923333,0.0002925910,0.0002928229,0.0002927529,0.0002931447,0.0002933016,0.0002935374,0.0002935579,0.0002940798,0.0002942503,0.0002941312,0.0002095436,0.0002125268,0.0002158927,0.0002181454,0.0002207086,0.0002231023,0.0002249202,0.0002273786,0.0002291636,0.0002310205,0.0002326992,0.0002343529,0.0002360025,0.0002376163,0.0002390141,0.0002405705,0.0002414206,0.0002428887,0.0002438756,0.0002452485,0.0002460510,0.0002471738,0.0002479694,0.0002488833,0.0002499368,0.0002507911,0.0002516446,0.0002523782,0.0002531874,0.0002541914,0.0002551296,0.0002558208,0.0002565329,0.0002569785,0.0002574617,0.0002583413,0.0002588627,0.0002594580,0.0002604060,0.0002607261,0.0002612421,0.0002617613,0.0002622337,0.0002628683,0.0002632271,0.0002636187,0.0002643260,0.0002646505,0.0002651096,0.0002656995,0.0002660791,0.0002667486,0.0002664783,0.0002674555,0.0002676217,0.0002682701,0.0002684799,0.0002689355,0.0002693068,0.0002693721,0.0002698938,0.0002700408,0.0002705917,0.0002708486,0.0002715827,0.0002714024,0.0002719567,0.0002719094,0.0002724996,0.0002728589,0.0002729879,0.0002735482,0.0002736212,0.0002740021,0.0002741370,0.0002746674,0.0002747898,0.0002749244,0.0002749109,0.0002755606,0.0002755473,0.0002759586,0.0002761026,0.0002763916,0.0002767164,0.0002770112,0.0002769713,0.0002769721,0.0002777794,0.0002776726,0.0002778713,0.0002781537,0.0002780864,0.0002784269,0.0002786225,0.0002788488,0.0002788054,0.0002793292,0.0002794976,0.0002793476,0.0001988792,0.0002016680,0.0002049152,0.0002070651,0.0002095074,0.0002118447,0.0002135082,0.0002158643,0.0002174775,0.0002193353,0.0002209212,0.0002224914,0.0002241033,0.0002256108,0.0002269032,0.0002284339,0.0002291331,0.0002306092,0.0002315801,0.0002328318,0.0002336010,0.0002346588,0.0002354279,0.0002362980,0.0002372791,0.0002381470,0.0002389016,0.0002395985,0.0002404930,0.0002413828,0.0002422075,0.0002429280,0.0002436113,0.0002440991,0.0002444800,0.0002453316,0.0002458437,0.0002464356,0.0002473290,0.0002475493,0.0002480966,0.0002486331,0.0002489751,0.0002497383,0.0002499328,0.0002503177,0.0002509759,0.0002513444,0.0002518157,0.0002523947,0.0002527010,0.0002533746,0.0002530990,0.0002539640,0.0002542039,0.0002547516,0.0002549241,0.0002554697,0.0002557894,0.0002558243,0.0002563212,0.0002564922,0.0002569688,0.0002572213,0.0002579283,0.0002577401,0.0002583263,0.0002582026,0.0002588751,0.0002591170,0.0002593174,0.0002598023,0.0002598692,0.0002602705,0.0002603274,0.0002608704,0.0002610245,0.0002611474,0.0002610924,0.0002617447,0.0002616621,0.0002620963,0.0002622675,0.0002625011,0.0002628476,0.0002631259,0.0002630542,0.0002630370,0.0002638059,0.0002637597,0.0002639140,0.0002641901,0.0002641396,0.0002645071,0.0002646687,0.0002648247,0.0002648461,0.0002653403,0.0002654647,0.0002653255,0.0001887439,0.0001914175,0.0001944905,0.0001965792,0.0001988713,0.0002010405,0.0002026678,0.0002048738,0.0002064775,0.0002081729,0.0002097728,0.0002111779,0.0002127326,0.0002142017,0.0002154129,0.0002168769,0.0002175230,0.0002190085,0.0002198630,0.0002210034,0.0002218097,0.0002227707,0.0002235454,0.0002243451,0.0002253141,0.0002260922,0.0002267737,0.0002275452,0.0002283371,0.0002292192,0.0002299818,0.0002307024,0.0002313335,0.0002317700,0.0002321968,0.0002330138,0.0002334578,0.0002340239,0.0002348771,0.0002351214,0.0002356590,0.0002360607,0.0002364039,0.0002371740,0.0002374159,0.0002376761,0.0002383336,0.0002387084,0.0002391521,0.0002397387,0.0002399788,0.0002406637,0.0002403085,0.0002412076,0.0002414039,0.0002419898,0.0002421871,0.0002426462,0.0002429828,0.0002430037,0.0002433961,0.0002435539,0.0002440563,0.0002442318,0.0002449857,0.0002448050,0.0002453197,0.0002452805,0.0002458474,0.0002461184,0.0002463370,0.0002468023,0.0002468236,0.0002471966,0.0002472196,0.0002478259,0.0002478791,0.0002480309,0.0002480044,0.0002486358,0.0002485968,0.0002489385,0.0002490600,0.0002493062,0.0002496781,0.0002499585,0.0002498762,0.0002498086,0.0002505879,0.0002505023,0.0002506734,0.0002509759,0.0002509313,0.0002512347,0.0002513763,0.0002515559,0.0002516029,0.0002519862,0.0002522061,0.0002520394,0.0001791824,0.0001816987,0.0001846215,0.0001865799,0.0001887520,0.0001908491,0.0001923658,0.0001944885,0.0001960162,0.0001975987,0.0001991374,0.0002004566,0.0002019638,0.0002033439,0.0002044805,0.0002058924,0.0002065629,0.0002079033,0.0002086411,0.0002098227,0.0002105726,0.0002115331,0.0002122514,0.0002129872,0.0002139190,0.0002146860,0.0002152951,0.0002160749,0.0002168548,0.0002176509,0.0002183757,0.0002191129,0.0002196977,0.0002201259,0.0002205011,0.0002212562,0.0002216718,0.0002222451,0.0002230339,0.0002232855,0.0002237581,0.0002242201,0.0002244953,0.0002252394,0.0002254366,0.0002256599,0.0002263493,0.0002266525,0.0002271136,0.0002276802,0.0002278926,0.0002285188,0.0002282227,0.0002290394,0.0002292337,0.0002298667,0.0002300147,0.0002304871,0.0002307690,0.0002307946,0.0002311665,0.0002313839,0.0002317344,0.0002319765,0.0002325903,0.0002325048,0.0002330431,0.0002329610,0.0002334581,0.0002337708,0.0002340045,0.0002344003,0.0002345258,0.0002347475,0.0002347396,0.0002353758,0.0002354494,0.0002355730,0.0002355574,0.0002360783,0.0002361484,0.0002364824,0.0002366107,0.0002367633,0.0002371915,0.0002374556,0.0002373325,0.0002373102,0.0002380051,0.0002379444,0.0002381039,0.0002383995,0.0002383455,0.0002386376,0.0002388173,0.0002389738,0.0002390040,0.0002393513,0.0002395968,0.0002394592,0.0001700043,0.0001724736,0.0001752853,0.0001771198,0.0001792179,0.0001811322,0.0001826159,0.0001846820,0.0001860702,0.0001875970,0.0001890719,0.0001903163,0.0001918036,0.0001930306,0.0001941239,0.0001954416,0.0001960391,0.0001973409,0.0001981249,0.0001992196,0.0001998989,0.0002008274,0.0002015116,0.0002022700,0.0002031147,0.0002038978,0.0002044127,0.0002051188,0.0002059229,0.0002066718,0.0002073541,0.0002080502,0.0002086264,0.0002090487,0.0002093821,0.0002101067,0.0002104669,0.0002110312,0.0002117529,0.0002120211,0.0002124830,0.0002128904,0.0002132408,0.0002139413,0.0002141219,0.0002143497,0.0002150407,0.0002152095,0.0002157019,0.0002162287,0.0002164687,0.0002170587,0.0002167582,0.0002175283,0.0002177037,0.0002183743,0.0002185079,0.0002189270,0.0002191584,0.0002191796,0.0002195754,0.0002197296,0.0002200736,0.0002203277,0.0002209170,0.0002208647,0.0002213759,0.0002212467,0.0002217866,0.0002221261,0.0002222179,0.0002226096,0.0002227202,0.0002229978,0.0002229900,0.0002236143,0.0002236500,0.0002237440,0.0002237977,0.0002242547,0.0002243122,0.0002246415,0.0002247542,0.0002248919,0.0002253184,0.0002255550,0.0002254161,0.0002254570,0.0002260462,0.0002259991,0.0002262326,0.0002265135,0.0002264320,0.0002266823,0.0002268715,0.0002269699,0.0002270886,0.0002273362,0.0002275816,0.0002274795,0.0001613312,0.0001637124,0.0001663335,0.0001681506,0.0001701258,0.0001719154,0.0001733444,0.0001753271,0.0001766371,0.0001781239,0.0001794879,0.0001807056,0.0001820825,0.0001832130,0.0001843061,0.0001855431,0.0001861233,0.0001873876,0.0001881419,0.0001891266,0.0001897767,0.0001906669,0.0001912794,0.0001920758,0.0001928644,0.0001936191,0.0001941536,0.0001947933,0.0001955121,0.0001962711,0.0001968947,0.0001975445,0.0001981087,0.0001985018,0.0001988033,0.0001994912,0.0001998766,0.0002004388,0.0002011139,0.0002013130,0.0002018036,0.0002021752,0.0002025186,0.0002032004,0.0002033632,0.0002035890,0.0002042572,0.0002043437,0.0002048617,0.0002053451,0.0002056099,0.0002061762,0.0002059215,0.0002066251,0.0002067385,0.0002073865,0.0002075258,0.0002079291,0.0002081761,0.0002081261,0.0002085513,0.0002087026,0.0002089895,0.0002092733,0.0002097611,0.0002097383,0.0002102617,0.0002101328,0.0002107018,0.0002109686,0.0002110278,0.0002114170,0.0002115156,0.0002118067,0.0002118220,0.0002123730,0.0002124337,0.0002124566,0.0002125729,0.0002130640,0.0002130343,0.0002134386,0.0002135064,0.0002136814,0.0002140428,0.0002142358,0.0002141874,0.0002142366,0.0002147672,0.0002146612,0.0002148599,0.0002151617,0.0002150604,0.0002153222,0.0002155017,0.0002155983,0.0002157463,0.0002159477,0.0002162176,0.0002161319,0.0001531696,0.0001553750,0.0001578937,0.0001596577,0.0001615477,0.0001632307,0.0001645959,0.0001664695,0.0001677513,0.0001690965,0.0001704224,0.0001715055,0.0001728479,0.0001739631,0.0001749415,0.0001761570,0.0001767831,0.0001779332,0.0001786206,0.0001795340,0.0001801733,0.0001809884,0.0001816474,0.0001823307,0.0001831153,0.0001838355,0.0001843998,0.0001849737,0.0001856516,0.0001863895,0.0001870615,0.0001875608,0.0001881371,0.0001885034,0.0001888570,0.0001894148,0.0001897786,0.0001903172,0.0001910151,0.0001911930,0.0001916360,0.0001920295,0.0001923679,0.0001929826,0.0001931432,0.0001933225,0.0001939989,0.0001941031,0.0001945358,0.0001950407,0.0001952892,0.0001957548,0.0001955417,0.0001962165,0.0001964217,0.0001969915,0.0001971425,0.0001975004,0.0001977260,0.0001976248,0.0001980682,0.0001982181,0.0001985166,0.0001987714,0.0001991622,0.0001992236,0.0001996940,0.0001996310,0.0002001844,0.0002003961,0.0002004267,0.0002007961,0.0002008840,0.0002012216,0.0002011899,0.0002016955,0.0002017901,0.0002017803,0.0002019313,0.0002024443,0.0002023880,0.0002026996,0.0002027474,0.0002030431,0.0002032417,0.0002035023,0.0002034805,0.0002034748,0.0002040197,0.0002038771,0.0002040476,0.0002044386,0.0002043421,0.0002045052,0.0002047211,0.0002047824,0.0002049214,0.0002051106,0.0002053725,0.0002052959,0.0001454162,0.0001475076,0.0001498811,0.0001515458,0.0001533281,0.0001549268,0.0001562489,0.0001580332,0.0001592706,0.0001605122,0.0001618048,0.0001628351,0.0001641104,0.0001652008,0.0001661387,0.0001672441,0.0001678463,0.0001689696,0.0001695721,0.0001704910,0.0001711091,0.0001719021,0.0001725081,0.0001731282,0.0001738978,0.0001745289,0.0001750585,0.0001756320,0.0001763165,0.0001770261,0.0001776351,0.0001781082,0.0001787218,0.0001790227,0.0001793576,0.0001798467,0.0001802751,0.0001806699,0.0001814197,0.0001815427,0.0001820199,0.0001824176,0.0001826688,0.0001832201,0.0001834196,0.0001835849,0.0001842504,0.0001843651,0.0001847750,0.0001852725,0.0001854643,0.0001859201,0.0001857333,0.0001863618,0.0001865774,0.0001870131,0.0001872064,0.0001876111,0.0001878085,0.0001877351,0.0001881156,0.0001882843,0.0001885267,0.0001888115,0.0001892013,0.0001892475,0.0001897744,0.0001896077,0.0001901632,0.0001903876,0.0001903842,0.0001907640,0.0001908195,0.0001911383,0.0001910804,0.0001915384,0.0001916269,0.0001916430,0.0001918171,0.0001923200,0.0001922583,0.0001925230,0.0001925387,0.0001928486,0.0001930552,0.0001932914,0.0001933201,0.0001932823,0.0001938014,0.0001936955,0.0001938422,0.0001942013,0.0001941060,0.0001942855,0.0001944534,0.0001945568,0.0001946058,0.0001948543,0.0001950284,0.0001949522,0.0001380212,0.0001400491,0.0001422358,0.0001438816,0.0001455251,0.0001470985,0.0001483115,0.0001500131,0.0001512356,0.0001523717,0.0001536253,0.0001546081,0.0001558164,0.0001568540,0.0001577784,0.0001588077,0.0001593680,0.0001604381,0.0001610397,0.0001618433,0.0001624308,0.0001632330,0.0001637779,0.0001644451,0.0001651076,0.0001657219,0.0001661971,0.0001667476,0.0001674941,0.0001680933,0.0001686887,0.0001691380,0.0001697533,0.0001700022,0.0001703149,0.0001707979,0.0001711954,0.0001715822,0.0001722710,0.0001723292,0.0001728162,0.0001732211,0.0001734763,0.0001740196,0.0001742070,0.0001743516,0.0001750218,0.0001750752,0.0001754909,0.0001759417,0.0001761611,0.0001765810,0.0001764271,0.0001770089,0.0001772208,0.0001776295,0.0001778443,0.0001782521,0.0001783771,0.0001783232,0.0001786934,0.0001788482,0.0001790129,0.0001792983,0.0001797171,0.0001797697,0.0001802325,0.0001800935,0.0001805503,0.0001808492,0.0001808401,0.0001811609,0.0001812500,0.0001815086,0.0001815636,0.0001819758,0.0001819645,0.0001820466,0.0001822288,0.0001827066,0.0001826354,0.0001828800,0.0001828867,0.0001831947,0.0001834218,0.0001836403,0.0001836772,0.0001835918,0.0001840942,0.0001839906,0.0001841360,0.0001844435,0.0001843562,0.0001845826,0.0001847387,0.0001848313,0.0001848622,0.0001850970,0.0001852744,0.0001852028,0.0001310296,0.0001329516,0.0001350597,0.0001365778,0.0001382226,0.0001396093,0.0001408042,0.0001424042,0.0001435740,0.0001446524,0.0001458479,0.0001468125,0.0001479603,0.0001489002,0.0001497715,0.0001508352,0.0001513116,0.0001523293,0.0001529371,0.0001536805,0.0001542464,0.0001549797,0.0001555205,0.0001561518,0.0001568278,0.0001573673,0.0001578245,0.0001582853,0.0001590860,0.0001596314,0.0001602109,0.0001605863,0.0001612190,0.0001614803,0.0001617300,0.0001621720,0.0001625674,0.0001629537,0.0001635883,0.0001636985,0.0001641173,0.0001644578,0.0001647324,0.0001652367,0.0001654429,0.0001656258,0.0001662503,0.0001662817,0.0001666764,0.0001671168,0.0001672467,0.0001676936,0.0001675715,0.0001680808,0.0001683033,0.0001687321,0.0001689290,0.0001692995,0.0001693800,0.0001693618,0.0001697071,0.0001698996,0.0001700256,0.0001702986,0.0001707285,0.0001707284,0.0001711953,0.0001710278,0.0001715103,0.0001718157,0.0001717636,0.0001720849,0.0001721445,0.0001724232,0.0001724109,0.0001728588,0.0001728198,0.0001729407,0.0001731398,0.0001735518,0.0001734482,0.0001737256,0.0001737198,0.0001740188,0.0001742196,0.0001744181,0.0001745060,0.0001743672,0.0001748598,0.0001747864,0.0001748955,0.0001752249,0.0001751472,0.0001753660,0.0001755047,0.0001755756,0.0001756640,0.0001758275,0.0001759876,0.0001759248,0.0001243709,0.0001262191,0.0001281986,0.0001297215,0.0001312402,0.0001325340,0.0001336462,0.0001351999,0.0001362853,0.0001373699,0.0001384306,0.0001393925,0.0001404305,0.0001413781,0.0001422203,0.0001432164,0.0001436445,0.0001446340,0.0001452279,0.0001459211,0.0001464808,0.0001471323,0.0001477046,0.0001482537,0.0001489285,0.0001494554,0.0001498212,0.0001502590,0.0001510923,0.0001515233,0.0001521585,0.0001524794,0.0001531080,0.0001533451,0.0001535980,0.0001540208,0.0001543890,0.0001547270,0.0001553571,0.0001554872,0.0001558857,0.0001561887,0.0001564551,0.0001569897,0.0001571370,0.0001572839,0.0001579097,0.0001579019,0.0001582927,0.0001587063,0.0001588404,0.0001592993,0.0001591327,0.0001596965,0.0001598705,0.0001602286,0.0001605012,0.0001608177,0.0001608846,0.0001608276,0.0001611893,0.0001613792,0.0001614939,0.0001617192,0.0001622256,0.0001621746,0.0001625742,0.0001623947,0.0001628694,0.0001632194,0.0001632067,0.0001634504,0.0001635270,0.0001638420,0.0001637852,0.0001641861,0.0001641543,0.0001643013,0.0001644592,0.0001648556,0.0001647777,0.0001650190,0.0001649958,0.0001653618,0.0001654951,0.0001656678,0.0001656986,0.0001656674,0.0001660480,0.0001660161,0.0001661245,0.0001664259,0.0001664044,0.0001665864,0.0001667207,0.0001667938,0.0001669103,0.0001670155,0.0001671752,0.0001671141,0.0001180423,0.0001198135,0.0001216681,0.0001231565,0.0001245560,0.0001258087,0.0001268640,0.0001283674,0.0001293633,0.0001304145,0.0001314203,0.0001323821,0.0001333574,0.0001342343,0.0001350394,0.0001359837,0.0001363981,0.0001373557,0.0001378608,0.0001386112,0.0001391040,0.0001396869,0.0001402698,0.0001407711,0.0001414450,0.0001419060,0.0001422950,0.0001426675,0.0001434531,0.0001438940,0.0001445244,0.0001448380,0.0001453710,0.0001456313,0.0001459034,0.0001462421,0.0001465589,0.0001469718,0.0001475752,0.0001476151,0.0001480461,0.0001483945,0.0001485759,0.0001491243,0.0001492930,0.0001493827,0.0001500003,0.0001499733,0.0001503340,0.0001507933,0.0001509262,0.0001513349,0.0001511371,0.0001516720,0.0001518764,0.0001522006,0.0001524740,0.0001527646,0.0001527999,0.0001527242,0.0001530973,0.0001533172,0.0001533920,0.0001535811,0.0001540967,0.0001540585,0.0001544082,0.0001542474,0.0001546724,0.0001550282,0.0001550538,0.0001552517,0.0001553640,0.0001556012,0.0001555684,0.0001559855,0.0001559116,0.0001560586,0.0001562335,0.0001566164,0.0001564744,0.0001567505,0.0001567903,0.0001570928,0.0001572474,0.0001573368,0.0001573956,0.0001573737,0.0001577518,0.0001577428,0.0001578329,0.0001581544,0.0001580775,0.0001582405,0.0001584188,0.0001584751,0.0001585338,0.0001587037,0.0001588218,0.0001587881,0.0001120614,0.0001137629,0.0001154847,0.0001169145,0.0001182511,0.0001193843,0.0001204002,0.0001219334,0.0001227958,0.0001238053,0.0001247785,0.0001256740,0.0001266331,0.0001274469,0.0001282143,0.0001291488,0.0001294975,0.0001304250,0.0001308928,0.0001316247,0.0001321107,0.0001326737,0.0001332298,0.0001336357,0.0001342770,0.0001348178,0.0001351008,0.0001354652,0.0001362438,0.0001366689,0.0001372850,0.0001375615,0.0001380143,0.0001382936,0.0001385796,0.0001388975,0.0001391673,0.0001396108,0.0001401291,0.0001402255,0.0001406003,0.0001409471,0.0001410985,0.0001416117,0.0001417854,0.0001418750,0.0001424828,0.0001424556,0.0001427586,0.0001432586,0.0001433507,0.0001437322,0.0001435247,0.0001440799,0.0001442888,0.0001445726,0.0001448116,0.0001451192,0.0001451447,0.0001450245,0.0001454136,0.0001456480,0.0001457076,0.0001458503,0.0001463746,0.0001463570,0.0001466555,0.0001464857,0.0001469166,0.0001472639,0.0001472896,0.0001474630,0.0001475593,0.0001477927,0.0001477909,0.0001482319,0.0001480787,0.0001482353,0.0001483559,0.0001487804,0.0001486607,0.0001489364,0.0001490090,0.0001492765,0.0001493526,0.0001494119,0.0001495487,0.0001495011,0.0001498703,0.0001498271,0.0001500087,0.0001502344,0.0001501458,0.0001502656,0.0001505054,0.0001505561,0.0001505874,0.0001507334,0.0001508667,0.0001508257,0.0001063511,0.0001079578,0.0001096279,0.0001109613,0.0001122387,0.0001133474,0.0001143369,0.0001157711,0.0001165632,0.0001175260,0.0001184847,0.0001192965,0.0001202347,0.0001210070,0.0001217588,0.0001226359,0.0001229938,0.0001238317,0.0001242945,0.0001249977,0.0001254497,0.0001259894,0.0001264913,0.0001268670,0.0001274872,0.0001280341,0.0001283262,0.0001286680,0.0001294112,0.0001297863,0.0001303057,0.0001306040,0.0001310754,0.0001313179,0.0001315811,0.0001319072,0.0001321702,0.0001326208,0.0001330966,0.0001331763,0.0001335163,0.0001337987,0.0001340056,0.0001344267,0.0001346492,0.0001347917,0.0001353277,0.0001352974,0.0001355823,0.0001360650,0.0001361588,0.0001365525,0.0001362890,0.0001368644,0.0001370796,0.0001372636,0.0001375294,0.0001378382,0.0001378376,0.0001377645,0.0001380518,0.0001383403,0.0001383747,0.0001385358,0.0001390596,0.0001389797,0.0001392574,0.0001391325,0.0001395692,0.0001398596,0.0001398943,0.0001400794,0.0001401499,0.0001403898,0.0001404205,0.0001407796,0.0001406456,0.0001408095,0.0001408850,0.0001413093,0.0001412587,0.0001414745,0.0001415616,0.0001418362,0.0001418632,0.0001419461,0.0001420678,0.0001419995,0.0001423527,0.0001423442,0.0001425415,0.0001427562,0.0001426545,0.0001427644,0.0001429450,0.0001429967,0.0001430594,0.0001432025,0.0001433190,0.0001432652,0.0001009818,0.0001025113,0.0001040660,0.0001053303,0.0001065775,0.0001076137,0.0001085350,0.0001099278,0.0001106906,0.0001115679,0.0001124966,0.0001132506,0.0001141156,0.0001148516,0.0001156202,0.0001164161,0.0001167582,0.0001175722,0.0001180226,0.0001187048,0.0001191217,0.0001196297,0.0001201010,0.0001204717,0.0001210511,0.0001215928,0.0001218141,0.0001222389,0.0001228785,0.0001232446,0.0001237271,0.0001240477,0.0001244590,0.0001247563,0.0001249578,0.0001252776,0.0001255366,0.0001259407,0.0001264261,0.0001264910,0.0001268066,0.0001270836,0.0001272519,0.0001276297,0.0001279026,0.0001280541,0.0001285274,0.0001285081,0.0001287698,0.0001292271,0.0001292952,0.0001297163,0.0001294052,0.0001299971,0.0001302133,0.0001303886,0.0001306520,0.0001309626,0.0001309253,0.0001308282,0.0001310999,0.0001314140,0.0001314343,0.0001316106,0.0001320598,0.0001320431,0.0001322956,0.0001322114,0.0001326099,0.0001328276,0.0001328839,0.0001330192,0.0001331482,0.0001333665,0.0001334235,0.0001337319,0.0001336276,0.0001337997,0.0001338530,0.0001342419,0.0001342067,0.0001343870,0.0001344317,0.0001347490,0.0001347751,0.0001348285,0.0001349879,0.0001349256,0.0001352334,0.0001352082,0.0001354226,0.0001355980,0.0001355336,0.0001356025,0.0001357812,0.0001358538,0.0001359008,0.0001360538,0.0001361383,0.0001361370,0.0000958802,0.0000972934,0.0000988489,0.0001000055,0.0001011563,0.0001021688,0.0001030279,0.0001043304,0.0001051021,0.0001059330,0.0001067971,0.0001074837,0.0001083596,0.0001090917,0.0001097946,0.0001105197,0.0001108686,0.0001116448,0.0001120866,0.0001127083,0.0001131340,0.0001135969,0.0001140379,0.0001144049,0.0001149727,0.0001154618,0.0001156678,0.0001161256,0.0001167069,0.0001170709,0.0001174896,0.0001177815,0.0001181643,0.0001184789,0.0001186820,0.0001189911,0.0001192188,0.0001196153,0.0001200507,0.0001201575,0.0001204413,0.0001206698,0.0001208892,0.0001212395,0.0001214937,0.0001215998,0.0001220628,0.0001220869,0.0001223009,0.0001227745,0.0001227971,0.0001232123,0.0001229782,0.0001234437,0.0001236800,0.0001238442,0.0001240942,0.0001244144,0.0001243963,0.0001242842,0.0001244960,0.0001248163,0.0001248571,0.0001249947,0.0001254074,0.0001254227,0.0001256193,0.0001255918,0.0001259460,0.0001261397,0.0001262189,0.0001263854,0.0001265068,0.0001267076,0.0001267508,0.0001270782,0.0001269276,0.0001270496,0.0001271771,0.0001275425,0.0001274878,0.0001276161,0.0001277011,0.0001280256,0.0001280290,0.0001280770,0.0001282475,0.0001282135,0.0001284573,0.0001284809,0.0001286575,0.0001288373,0.0001288031,0.0001288000,0.0001289234,0.0001290380,0.0001290947,0.0001292697,0.0001293253,0.0001293424,0.0000910474,0.0000924018,0.0000938220,0.0000949278,0.0000960215,0.0000970019,0.0000977843,0.0000990338,0.0000997739,0.0001005784,0.0001013656,0.0001020645,0.0001028664,0.0001035362,0.0001042390,0.0001049568,0.0001052473,0.0001060504,0.0001064618,0.0001070442,0.0001074362,0.0001078832,0.0001083076,0.0001086714,0.0001091536,0.0001096617,0.0001098635,0.0001102854,0.0001108076,0.0001111415,0.0001115654,0.0001118379,0.0001122291,0.0001124992,0.0001127092,0.0001130098,0.0001132625,0.0001136103,0.0001140318,0.0001140860,0.0001143861,0.0001146302,0.0001148334,0.0001151464,0.0001153839,0.0001154904,0.0001158898,0.0001159833,0.0001161568,0.0001165785,0.0001166199,0.0001170206,0.0001168359,0.0001173051,0.0001174823,0.0001176318,0.0001179138,0.0001182179,0.0001181597,0.0001180409,0.0001182437,0.0001185626,0.0001185752,0.0001187370,0.0001191068,0.0001191418,0.0001193020,0.0001193131,0.0001196177,0.0001198224,0.0001199137,0.0001200013,0.0001201795,0.0001203557,0.0001204294,0.0001207193,0.0001206123,0.0001207105,0.0001207934,0.0001211307,0.0001211025,0.0001212215,0.0001213547,0.0001215940,0.0001216581,0.0001216569,0.0001218371,0.0001217934,0.0001220261,0.0001220972,0.0001223091,0.0001223997,0.0001223868,0.0001223597,0.0001224876,0.0001225682,0.0001225940,0.0001228333,0.0001228808,0.0001228809,0.0000864044,0.0000877219,0.0000890889,0.0000901318,0.0000911776,0.0000920803,0.0000928172,0.0000940417,0.0000947249,0.0000954920,0.0000962533,0.0000969246,0.0000976687,0.0000983190,0.0000989680,0.0000996596,0.0000999209,0.0001007434,0.0001011010,0.0001016098,0.0001020190,0.0001024530,0.0001028642,0.0001031929,0.0001036563,0.0001041592,0.0001043300,0.0001047445,0.0001052504,0.0001055643,0.0001059571,0.0001062078,0.0001066059,0.0001068548,0.0001070763,0.0001073270,0.0001075733,0.0001078962,0.0001082902,0.0001083472,0.0001086192,0.0001088684,0.0001090258,0.0001093744,0.0001095689,0.0001096641,0.0001100502,0.0001101735,0.0001103326,0.0001107227,0.0001107603,0.0001111420,0.0001109762,0.0001114488,0.0001115776,0.0001117183,0.0001120342,0.0001122871,0.0001122079,0.0001120838,0.0001123166,0.0001126107,0.0001126193,0.0001127871,0.0001131334,0.0001131738,0.0001133413,0.0001133170,0.0001136122,0.0001138738,0.0001139173,0.0001139858,0.0001141705,0.0001143383,0.0001144101,0.0001147012,0.0001145168,0.0001146810,0.0001147473,0.0001150657,0.0001150422,0.0001151279,0.0001153249,0.0001154684,0.0001155796,0.0001155365,0.0001157528,0.0001156983,0.0001159097,0.0001159861,0.0001161658,0.0001162889,0.0001162554,0.0001162402,0.0001163487,0.0001165061,0.0001164135,0.0001167064,0.0001167437,0.0001167340,0.0000820235,0.0000832634,0.0000845906,0.0000855870,0.0000866054,0.0000873992,0.0000881503,0.0000892943,0.0000899097,0.0000906666,0.0000914033,0.0000920371,0.0000927301,0.0000933264,0.0000939706,0.0000946536,0.0000948801,0.0000956868,0.0000960111,0.0000964908,0.0000968912,0.0000973246,0.0000976663,0.0000980120,0.0000983916,0.0000989016,0.0000990696,0.0000994986,0.0000999645,0.0001002303,0.0001006419,0.0001008457,0.0001012293,0.0001014723,0.0001016292,0.0001019448,0.0001021278,0.0001024658,0.0001028543,0.0001029058,0.0001031662,0.0001034512,0.0001035565,0.0001038973,0.0001041028,0.0001041989,0.0001045272,0.0001046707,0.0001048015,0.0001051665,0.0001051819,0.0001055384,0.0001054262,0.0001058668,0.0001060078,0.0001061279,0.0001064309,0.0001066443,0.0001065659,0.0001064710,0.0001067048,0.0001069639,0.0001069595,0.0001071549,0.0001074720,0.0001074792,0.0001076460,0.0001076518,0.0001079199,0.0001081691,0.0001082105,0.0001082712,0.0001084184,0.0001086198,0.0001086981,0.0001089621,0.0001088123,0.0001089415,0.0001090232,0.0001093279,0.0001092930,0.0001093597,0.0001095707,0.0001097327,0.0001097768,0.0001097322,0.0001099525,0.0001099331,0.0001101298,0.0001102457,0.0001103547,0.0001104806,0.0001104255,0.0001104315,0.0001105349,0.0001106636,0.0001105823,0.0001108477,0.0001109167,0.0001108783,0.0000778660,0.0000790586,0.0000803018,0.0000812740,0.0000822464,0.0000829739,0.0000837030,0.0000847966,0.0000854265,0.0000860769,0.0000867771,0.0000873851,0.0000880384,0.0000886507,0.0000892551,0.0000898654,0.0000901238,0.0000908549,0.0000911707,0.0000916638,0.0000919866,0.0000923964,0.0000927529,0.0000930959,0.0000934431,0.0000938986,0.0000940675,0.0000944984,0.0000949540,0.0000951332,0.0000955861,0.0000957468,0.0000961287,0.0000963703,0.0000965079,0.0000968265,0.0000969765,0.0000973320,0.0000976741,0.0000977476,0.0000979881,0.0000982110,0.0000983644,0.0000986905,0.0000988875,0.0000989773,0.0000992731,0.0000994592,0.0000995887,0.0000999119,0.0000998939,0.0001002577,0.0001001575,0.0001005417,0.0001006672,0.0001008187,0.0001011273,0.0001012908,0.0001012283,0.0001010983,0.0001013470,0.0001015999,0.0001015774,0.0001017987,0.0001020954,0.0001020624,0.0001022838,0.0001022704,0.0001024816,0.0001027558,0.0001028282,0.0001028738,0.0001029592,0.0001031761,0.0001032448,0.0001034898,0.0001033845,0.0001034494,0.0001035577,0.0001038464,0.0001038009,0.0001038773,0.0001040785,0.0001042817,0.0001042919,0.0001042283,0.0001044577,0.0001044827,0.0001045694,0.0001046890,0.0001048238,0.0001049268,0.0001048897,0.0001049324,0.0001050139,0.0001051083,0.0001050531,0.0001053180,0.0001053293,0.0001053490,0.0000739420,0.0000750660,0.0000762772,0.0000771201,0.0000780963,0.0000787919,0.0000794658,0.0000805168,0.0000811004,0.0000817038,0.0000824242,0.0000829580,0.0000835621,0.0000841802,0.0000847267,0.0000853168,0.0000856089,0.0000862604,0.0000865602,0.0000870356,0.0000873483,0.0000877120,0.0000881015,0.0000883952,0.0000887518,0.0000891559,0.0000893667,0.0000897456,0.0000902020,0.0000903748,0.0000907675,0.0000909120,0.0000913014,0.0000915412,0.0000916474,0.0000919687,0.0000920859,0.0000924218,0.0000927768,0.0000928477,0.0000930848,0.0000932766,0.0000934513,0.0000937437,0.0000938865,0.0000939822,0.0000943084,0.0000944829,0.0000945923,0.0000948933,0.0000948949,0.0000952537,0.0000951404,0.0000955271,0.0000956564,0.0000957662,0.0000960396,0.0000961916,0.0000961562,0.0000960155,0.0000962705,0.0000965062,0.0000965070,0.0000967204,0.0000969856,0.0000969005,0.0000971228,0.0000971452,0.0000973139,0.0000976224,0.0000976324,0.0000977240,0.0000978103,0.0000980121,0.0000980792,0.0000982700,0.0000982136,0.0000982495,0.0000983367,0.0000986435,0.0000985978,0.0000987216,0.0000988691,0.0000990827,0.0000991032,0.0000990338,0.0000992345,0.0000992894,0.0000993234,0.0000994253,0.0000995469,0.0000997116,0.0000996274,0.0000996855,0.0000997762,0.0000998451,0.0000998055,0.0001000326,0.0001000707,0.0001001038,0.0000701601,0.0000712507,0.0000724087,0.0000732187,0.0000741345,0.0000748387,0.0000754285,0.0000764750,0.0000769852,0.0000776031,0.0000782568,0.0000787639,0.0000793656,0.0000799236,0.0000804571,0.0000809795,0.0000813260,0.0000819027,0.0000822084,0.0000826291,0.0000829479,0.0000832832,0.0000836905,0.0000839707,0.0000843233,0.0000846623,0.0000848819,0.0000852591,0.0000856563,0.0000857945,0.0000861931,0.0000863464,0.0000867436,0.0000869484,0.0000870172,0.0000873709,0.0000874552,0.0000877665,0.0000881117,0.0000881759,0.0000883985,0.0000885836,0.0000887387,0.0000890521,0.0000891828,0.0000892568,0.0000895616,0.0000897442,0.0000898237,0.0000901442,0.0000901626,0.0000904551,0.0000903586,0.0000906942,0.0000908632,0.0000909532,0.0000912084,0.0000913555,0.0000913686,0.0000912036,0.0000914360,0.0000916606,0.0000916770,0.0000918767,0.0000921355,0.0000920438,0.0000922111,0.0000923100,0.0000924182,0.0000926850,0.0000927662,0.0000928339,0.0000929282,0.0000931169,0.0000931682,0.0000933750,0.0000933240,0.0000933603,0.0000934153,0.0000937192,0.0000936494,0.0000937889,0.0000939514,0.0000941281,0.0000941465,0.0000940572,0.0000942907,0.0000943559,0.0000943346,0.0000944803,0.0000945523,0.0000947329,0.0000946199,0.0000946896,0.0000947977,0.0000949029,0.0000948126,0.0000950073,0.0000950492,0.0000951043,0.0000666372,0.0000676786,0.0000687640,0.0000695301,0.0000704184,0.0000710710,0.0000716057,0.0000726034,0.0000730862,0.0000736636,0.0000743076,0.0000747983,0.0000753775,0.0000759009,0.0000764198,0.0000769163,0.0000771970,0.0000777917,0.0000780643,0.0000784708,0.0000787459,0.0000790763,0.0000795056,0.0000797239,0.0000801271,0.0000804239,0.0000806225,0.0000809476,0.0000813639,0.0000814533,0.0000818914,0.0000820109,0.0000823964,0.0000825714,0.0000826564,0.0000829760,0.0000830954,0.0000833661,0.0000836869,0.0000837557,0.0000839856,0.0000841440,0.0000842579,0.0000845790,0.0000846861,0.0000847395,0.0000850557,0.0000852331,0.0000853227,0.0000856229,0.0000856267,0.0000859123,0.0000858153,0.0000861453,0.0000863038,0.0000864192,0.0000866647,0.0000868379,0.0000868184,0.0000866344,0.0000868634,0.0000870753,0.0000870897,0.0000872806,0.0000875395,0.0000874133,0.0000875977,0.0000877140,0.0000877891,0.0000880628,0.0000881297,0.0000881936,0.0000882940,0.0000884551,0.0000885327,0.0000887219,0.0000886474,0.0000886818,0.0000887478,0.0000890504,0.0000889866,0.0000891111,0.0000892872,0.0000894334,0.0000894543,0.0000893215,0.0000895488,0.0000896053,0.0000896182,0.0000897480,0.0000898261,0.0000899839,0.0000899083,0.0000899805,0.0000900493,0.0000901584,0.0000900385,0.0000903183,0.0000903108,0.0000903458,0.0000632742,0.0000642563,0.0000652966,0.0000659761,0.0000668404,0.0000674978,0.0000679450,0.0000689225,0.0000693650,0.0000699260,0.0000705821,0.0000710373,0.0000715873,0.0000720805,0.0000725756,0.0000730250,0.0000732902,0.0000738823,0.0000741252,0.0000745453,0.0000747671,0.0000751159,0.0000755159,0.0000756892,0.0000760999,0.0000763940,0.0000765444,0.0000768679,0.0000772606,0.0000773429,0.0000777925,0.0000778952,0.0000782767,0.0000784147,0.0000784963,0.0000787718,0.0000788797,0.0000791870,0.0000795021,0.0000795608,0.0000797523,0.0000798959,0.0000799887,0.0000803441,0.0000804461,0.0000804801,0.0000807975,0.0000809211,0.0000810575,0.0000813254,0.0000813648,0.0000816178,0.0000815019,0.0000817950,0.0000819525,0.0000821048,0.0000823348,0.0000825322,0.0000824872,0.0000822918,0.0000825469,0.0000827353,0.0000827423,0.0000828956,0.0000831500,0.0000830485,0.0000832018,0.0000833373,0.0000834133,0.0000836423,0.0000837133,0.0000837445,0.0000839226,0.0000840338,0.0000841434,0.0000842938,0.0000841909,0.0000842411,0.0000843354,0.0000846003,0.0000845637,0.0000846946,0.0000847859,0.0000849773,0.0000849869,0.0000848593,0.0000850855,0.0000851428,0.0000851504,0.0000852636,0.0000853707,0.0000854954,0.0000854287,0.0000854799,0.0000855358,0.0000857202,0.0000855329,0.0000858363,0.0000857962,0.0000858249,0.0000600607,0.0000609908,0.0000619724,0.0000626144,0.0000634275,0.0000640895,0.0000644976,0.0000654254,0.0000658815,0.0000664341,0.0000670496,0.0000674335,0.0000679694,0.0000684594,0.0000689500,0.0000693351,0.0000696047,0.0000701174,0.0000703821,0.0000707887,0.0000709894,0.0000713585,0.0000716851,0.0000718460,0.0000722454,0.0000725387,0.0000727611,0.0000730166,0.0000733981,0.0000734410,0.0000738827,0.0000739833,0.0000743251,0.0000744592,0.0000745478,0.0000747848,0.0000749034,0.0000752105,0.0000755048,0.0000755871,0.0000757175,0.0000758992,0.0000759343,0.0000763393,0.0000764209,0.0000764171,0.0000767637,0.0000768437,0.0000770046,0.0000772530,0.0000773063,0.0000775058,0.0000774003,0.0000777070,0.0000778742,0.0000779995,0.0000781974,0.0000783834,0.0000783330,0.0000781769,0.0000784223,0.0000785674,0.0000786172,0.0000787771,0.0000789997,0.0000789032,0.0000790376,0.0000791753,0.0000792457,0.0000794338,0.0000795543,0.0000795552,0.0000797176,0.0000798280,0.0000799569,0.0000800959,0.0000799889,0.0000800384,0.0000801130,0.0000803743,0.0000803696,0.0000804615,0.0000805347,0.0000807412,0.0000807088,0.0000806368,0.0000808115,0.0000808752,0.0000809155,0.0000810172,0.0000811249,0.0000812014,0.0000811925,0.0000811930,0.0000812752,0.0000814313,0.0000812511,0.0000815408,0.0000815011,0.0000815575,0.0000570174,0.0000578770,0.0000588726,0.0000594369,0.0000602575,0.0000608452,0.0000612427,0.0000621329,0.0000625398,0.0000630603,0.0000636744,0.0000640532,0.0000645554,0.0000650205,0.0000654772,0.0000658423,0.0000661134,0.0000665831,0.0000668262,0.0000671981,0.0000673804,0.0000677737,0.0000681013,0.0000682594,0.0000685860,0.0000689235,0.0000691071,0.0000693602,0.0000697020,0.0000697196,0.0000701837,0.0000702600,0.0000705839,0.0000707154,0.0000707769,0.0000710385,0.0000711331,0.0000714241,0.0000716935,0.0000717903,0.0000719284,0.0000721186,0.0000721524,0.0000725266,0.0000726142,0.0000725684,0.0000729441,0.0000730025,0.0000731481,0.0000733897,0.0000734269,0.0000736369,0.0000735096,0.0000738051,0.0000740022,0.0000741144,0.0000742656,0.0000744318,0.0000743777,0.0000742742,0.0000744896,0.0000746214,0.0000746811,0.0000748600,0.0000750302,0.0000749422,0.0000750646,0.0000752294,0.0000752917,0.0000754818,0.0000755714,0.0000756055,0.0000757458,0.0000758201,0.0000759443,0.0000761187,0.0000759915,0.0000760122,0.0000760956,0.0000763535,0.0000763410,0.0000764764,0.0000764740,0.0000766710,0.0000767018,0.0000765829,0.0000767632,0.0000768125,0.0000768906,0.0000769614,0.0000771328,0.0000771219,0.0000771083,0.0000771395,0.0000772412,0.0000773267,0.0000772105,0.0000774459,0.0000774204,0.0000774607,0.0000541450,0.0000549575,0.0000558672,0.0000564233,0.0000571959,0.0000577601,0.0000581426,0.0000589840,0.0000593712,0.0000598839,0.0000604229,0.0000608222,0.0000613178,0.0000617414,0.0000621988,0.0000625197,0.0000627759,0.0000632448,0.0000634679,0.0000638080,0.0000639726,0.0000643750,0.0000646632,0.0000648621,0.0000651540,0.0000654668,0.0000656449,0.0000658724,0.0000661928,0.0000661991,0.0000666691,0.0000667255,0.0000670382,0.0000671700,0.0000672170,0.0000674814,0.0000675506,0.0000678691,0.0000680753,0.0000681815,0.0000683348,0.0000685167,0.0000685361,0.0000688610,0.0000689630,0.0000689255,0.0000692899,0.0000693546,0.0000694800,0.0000696912,0.0000697346,0.0000699364,0.0000698199,0.0000701075,0.0000702672,0.0000704185,0.0000705299,0.0000706977,0.0000706665,0.0000705772,0.0000707673,0.0000708984,0.0000709518,0.0000710796,0.0000712591,0.0000711991,0.0000713285,0.0000714666,0.0000715448,0.0000716904,0.0000717844,0.0000717868,0.0000719693,0.0000720257,0.0000721532,0.0000723076,0.0000722108,0.0000722331,0.0000722855,0.0000725765,0.0000725278,0.0000726309,0.0000726503,0.0000728525,0.0000728706,0.0000727454,0.0000729630,0.0000729662,0.0000730355,0.0000731231,0.0000733174,0.0000732568,0.0000732582,0.0000732905,0.0000733906,0.0000734680,0.0000733878,0.0000735877,0.0000735810,0.0000735795,0.0000514104,0.0000521762,0.0000530304,0.0000535710,0.0000543173,0.0000548387,0.0000552056,0.0000559923,0.0000563771,0.0000568637,0.0000573488,0.0000577386,0.0000582485,0.0000586043,0.0000590655,0.0000593812,0.0000595994,0.0000600480,0.0000602924,0.0000605945,0.0000607890,0.0000611321,0.0000614107,0.0000615861,0.0000618780,0.0000621807,0.0000623536,0.0000625582,0.0000628492,0.0000628832,0.0000633503,0.0000633902,0.0000636973,0.0000637984,0.0000638612,0.0000640972,0.0000641680,0.0000644348,0.0000646124,0.0000647434,0.0000649078,0.0000650632,0.0000651133,0.0000654161,0.0000654696,0.0000654492,0.0000658266,0.0000658801,0.0000659874,0.0000662029,0.0000662406,0.0000664234,0.0000663181,0.0000666106,0.0000667488,0.0000668767,0.0000669611,0.0000671803,0.0000671113,0.0000670699,0.0000672139,0.0000673688,0.0000673892,0.0000674937,0.0000676763,0.0000676209,0.0000677702,0.0000679012,0.0000679615,0.0000681120,0.0000682072,0.0000681544,0.0000683892,0.0000684166,0.0000685293,0.0000687188,0.0000685959,0.0000686059,0.0000687099,0.0000689443,0.0000688986,0.0000689773,0.0000690332,0.0000692096,0.0000692220,0.0000690745,0.0000693044,0.0000693182,0.0000693886,0.0000694610,0.0000696706,0.0000695905,0.0000696273,0.0000696396,0.0000697493,0.0000697837,0.0000697323,0.0000699018,0.0000699285,0.0000698965,0.0000488151,0.0000495346,0.0000503571,0.0000508650,0.0000515656,0.0000520563,0.0000524302,0.0000531472,0.0000535472,0.0000539650,0.0000544477,0.0000548435,0.0000553034,0.0000556349,0.0000560661,0.0000563740,0.0000566015,0.0000570207,0.0000572509,0.0000575465,0.0000577472,0.0000580372,0.0000583231,0.0000584795,0.0000587624,0.0000590843,0.0000592251,0.0000594308,0.0000597018,0.0000597239,0.0000601577,0.0000601910,0.0000604964,0.0000606162,0.0000606323,0.0000608809,0.0000609592,0.0000612093,0.0000613872,0.0000614934,0.0000616272,0.0000618256,0.0000618458,0.0000621182,0.0000621777,0.0000621544,0.0000624808,0.0000625637,0.0000627012,0.0000628993,0.0000629183,0.0000631047,0.0000630070,0.0000633350,0.0000634237,0.0000635306,0.0000636155,0.0000638377,0.0000637553,0.0000637254,0.0000638759,0.0000640121,0.0000640011,0.0000640968,0.0000643042,0.0000642260,0.0000643907,0.0000645067,0.0000645912,0.0000646734,0.0000647954,0.0000647613,0.0000649861,0.0000650084,0.0000651011,0.0000652858,0.0000651503,0.0000651732,0.0000652333,0.0000654672,0.0000654688,0.0000655703,0.0000655872,0.0000657439,0.0000657512,0.0000656356,0.0000658699,0.0000658728,0.0000659095,0.0000659822,0.0000662058,0.0000660978,0.0000661575,0.0000661627,0.0000663301,0.0000662687,0.0000662320,0.0000663993,0.0000664724,0.0000664031,0.0000463264,0.0000470407,0.0000478206,0.0000482906,0.0000489570,0.0000494240,0.0000497737,0.0000504684,0.0000508423,0.0000512735,0.0000517337,0.0000520990,0.0000525437,0.0000528282,0.0000532475,0.0000534995,0.0000537467,0.0000541499,0.0000543758,0.0000546377,0.0000548306,0.0000551304,0.0000553862,0.0000555458,0.0000558175,0.0000561232,0.0000562886,0.0000564540,0.0000566912,0.0000567536,0.0000571295,0.0000571706,0.0000574675,0.0000575630,0.0000575801,0.0000578322,0.0000579016,0.0000581608,0.0000583348,0.0000584057,0.0000585222,0.0000586929,0.0000587489,0.0000590042,0.0000590600,0.0000590552,0.0000593392,0.0000594326,0.0000595477,0.0000597505,0.0000597707,0.0000599303,0.0000598396,0.0000601757,0.0000602183,0.0000603357,0.0000604300,0.0000606556,0.0000605526,0.0000605252,0.0000606877,0.0000608258,0.0000608149,0.0000608885,0.0000610873,0.0000609806,0.0000611893,0.0000612216,0.0000613469,0.0000614210,0.0000615415,0.0000615094,0.0000617546,0.0000617624,0.0000618322,0.0000620199,0.0000618919,0.0000619299,0.0000620063,0.0000621952,0.0000622152,0.0000622773,0.0000623101,0.0000624616,0.0000624793,0.0000623727,0.0000625985,0.0000626133,0.0000626296,0.0000626796,0.0000628905,0.0000628011,0.0000628672,0.0000628777,0.0000630169,0.0000629494,0.0000629083,0.0000630945,0.0000631587,0.0000630724,0.0000439533,0.0000446506,0.0000454020,0.0000458529,0.0000464896,0.0000469378,0.0000472613,0.0000479312,0.0000482765,0.0000486938,0.0000491196,0.0000494683,0.0000499139,0.0000501726,0.0000505559,0.0000508279,0.0000510381,0.0000514233,0.0000516454,0.0000518922,0.0000520725,0.0000523385,0.0000525985,0.0000527451,0.0000530078,0.0000532894,0.0000534531,0.0000536101,0.0000538239,0.0000538998,0.0000542834,0.0000543112,0.0000545711,0.0000546814,0.0000546713,0.0000549372,0.0000549750,0.0000552621,0.0000554005,0.0000554677,0.0000555697,0.0000557759,0.0000558146,0.0000560601,0.0000561057,0.0000561146,0.0000563791,0.0000564664,0.0000565581,0.0000567652,0.0000567813,0.0000569252,0.0000568451,0.0000571942,0.0000571983,0.0000573134,0.0000573840,0.0000575842,0.0000575239,0.0000574884,0.0000576611,0.0000578057,0.0000577454,0.0000578371,0.0000580224,0.0000579187,0.0000581395,0.0000581412,0.0000582526,0.0000583435,0.0000584804,0.0000584275,0.0000586346,0.0000586710,0.0000587666,0.0000588983,0.0000587881,0.0000588427,0.0000589113,0.0000591341,0.0000591097,0.0000591560,0.0000592123,0.0000593203,0.0000593645,0.0000592489,0.0000594755,0.0000594892,0.0000594795,0.0000595508,0.0000597648,0.0000596420,0.0000597314,0.0000597192,0.0000598628,0.0000598150,0.0000598031,0.0000599488,0.0000600189,0.0000598994,0.0000417201,0.0000424051,0.0000431160,0.0000435472,0.0000441374,0.0000445725,0.0000448549,0.0000455053,0.0000458422,0.0000462178,0.0000466340,0.0000469905,0.0000473828,0.0000476284,0.0000480085,0.0000482677,0.0000484920,0.0000488165,0.0000490404,0.0000492723,0.0000494571,0.0000497049,0.0000499663,0.0000500641,0.0000503548,0.0000506111,0.0000507643,0.0000509144,0.0000511340,0.0000511729,0.0000515666,0.0000515721,0.0000518542,0.0000519517,0.0000519295,0.0000521547,0.0000522262,0.0000524552,0.0000526195,0.0000527034,0.0000527500,0.0000529792,0.0000529893,0.0000532463,0.0000533031,0.0000533137,0.0000535533,0.0000535964,0.0000537478,0.0000539348,0.0000539253,0.0000540795,0.0000540176,0.0000543130,0.0000543286,0.0000544392,0.0000544872,0.0000546898,0.0000546274,0.0000546037,0.0000547631,0.0000549317,0.0000548613,0.0000548967,0.0000551285,0.0000550104,0.0000552131,0.0000552059,0.0000553585,0.0000554411,0.0000555514,0.0000554834,0.0000556816,0.0000557267,0.0000558123,0.0000559591,0.0000558452,0.0000558978,0.0000559667,0.0000561715,0.0000561738,0.0000561839,0.0000562329,0.0000563594,0.0000563851,0.0000562528,0.0000565255,0.0000565610,0.0000564727,0.0000565828,0.0000568065,0.0000566514,0.0000567382,0.0000567302,0.0000568740,0.0000568269,0.0000567930,0.0000569446,0.0000570206,0.0000569055,0.0000395961,0.0000402458,0.0000409552,0.0000413467,0.0000419054,0.0000423520,0.0000425971,0.0000432103,0.0000435210,0.0000439075,0.0000442622,0.0000446131,0.0000449861,0.0000452014,0.0000455806,0.0000458170,0.0000460885,0.0000463617,0.0000465806,0.0000467832,0.0000469852,0.0000471842,0.0000474486,0.0000475536,0.0000478385,0.0000480589,0.0000482048,0.0000483896,0.0000485603,0.0000486161,0.0000490032,0.0000489687,0.0000492616,0.0000493437,0.0000492975,0.0000495518,0.0000496065,0.0000498366,0.0000499967,0.0000500513,0.0000501160,0.0000503317,0.0000503302,0.0000505936,0.0000506177,0.0000506637,0.0000508805,0.0000509185,0.0000510548,0.0000512304,0.0000512047,0.0000513465,0.0000512996,0.0000515813,0.0000516096,0.0000517184,0.0000517617,0.0000519418,0.0000518842,0.0000518645,0.0000519958,0.0000521705,0.0000521183,0.0000521704,0.0000523910,0.0000522789,0.0000524563,0.0000524232,0.0000525903,0.0000526597,0.0000528092,0.0000527046,0.0000528946,0.0000529214,0.0000530379,0.0000531625,0.0000530380,0.0000531020,0.0000531484,0.0000533373,0.0000533619,0.0000533750,0.0000533996,0.0000535849,0.0000535650,0.0000534418,0.0000537195,0.0000537462,0.0000536356,0.0000537531,0.0000539516,0.0000538436,0.0000539185,0.0000539026,0.0000540406,0.0000539888,0.0000539264,0.0000540967,0.0000541605,0.0000541004,0.0000375887,0.0000381932,0.0000388768,0.0000392408,0.0000397660,0.0000402163,0.0000404591,0.0000410012,0.0000413358,0.0000416960,0.0000420184,0.0000423732,0.0000427116,0.0000429389,0.0000432883,0.0000435054,0.0000437680,0.0000440430,0.0000442662,0.0000444375,0.0000446197,0.0000448069,0.0000450376,0.0000451539,0.0000454664,0.0000456257,0.0000458060,0.0000459724,0.0000461118,0.0000461412,0.0000465187,0.0000464958,0.0000467907,0.0000468988,0.0000468202,0.0000470939,0.0000471514,0.0000473299,0.0000475122,0.0000475475,0.0000476308,0.0000478069,0.0000477879,0.0000480677,0.0000480630,0.0000481267,0.0000483269,0.0000483790,0.0000485011,0.0000486551,0.0000486428,0.0000487777,0.0000487152,0.0000489961,0.0000490263,0.0000491091,0.0000491629,0.0000493584,0.0000492931,0.0000492902,0.0000493971,0.0000495636,0.0000495103,0.0000495713,0.0000497805,0.0000496815,0.0000498449,0.0000498131,0.0000499572,0.0000500403,0.0000501776,0.0000500734,0.0000502554,0.0000502802,0.0000504166,0.0000504788,0.0000503659,0.0000504349,0.0000505213,0.0000506803,0.0000507416,0.0000507133,0.0000507309,0.0000509009,0.0000508779,0.0000507866,0.0000510331,0.0000510397,0.0000509526,0.0000510573,0.0000512388,0.0000511360,0.0000512197,0.0000512374,0.0000513489,0.0000512799,0.0000512218,0.0000513916,0.0000514682,0.0000513799,0.0000356874,0.0000362901,0.0000368971,0.0000372487,0.0000377630,0.0000381845,0.0000384143,0.0000389238,0.0000392438,0.0000395993,0.0000399025,0.0000402290,0.0000405529,0.0000408008,0.0000410885,0.0000412940,0.0000415499,0.0000418368,0.0000420422,0.0000421799,0.0000423932,0.0000425196,0.0000427541,0.0000428914,0.0000431880,0.0000433297,0.0000434974,0.0000436656,0.0000437932,0.0000438114,0.0000441800,0.0000441614,0.0000444702,0.0000445358,0.0000444748,0.0000447180,0.0000447821,0.0000449607,0.0000451110,0.0000451681,0.0000452355,0.0000454038,0.0000454038,0.0000456297,0.0000456467,0.0000457127,0.0000458984,0.0000459762,0.0000460751,0.0000462295,0.0000462033,0.0000463403,0.0000462952,0.0000465382,0.0000465620,0.0000466549,0.0000467275,0.0000468929,0.0000468085,0.0000468184,0.0000468921,0.0000470947,0.0000470319,0.0000470668,0.0000472549,0.0000472077,0.0000473811,0.0000473166,0.0000474466,0.0000475446,0.0000476438,0.0000475629,0.0000477470,0.0000477669,0.0000478951,0.0000479516,0.0000478428,0.0000479080,0.0000479763,0.0000481389,0.0000481813,0.0000481645,0.0000481968,0.0000483566,0.0000483386,0.0000482248,0.0000484868,0.0000484876,0.0000484101,0.0000485215,0.0000486914,0.0000485824,0.0000486485,0.0000486936,0.0000487982,0.0000487282,0.0000486581,0.0000488432,0.0000488989,0.0000488282,0.0000338831,0.0000344541,0.0000350236,0.0000353616,0.0000358659,0.0000362512,0.0000364828,0.0000369657,0.0000372551,0.0000375929,0.0000378976,0.0000382031,0.0000385058,0.0000387647,0.0000390086,0.0000392067,0.0000394595,0.0000397496,0.0000399407,0.0000400748,0.0000402718,0.0000403734,0.0000406194,0.0000407339,0.0000410177,0.0000411798,0.0000413291,0.0000414529,0.0000415867,0.0000415962,0.0000419743,0.0000419533,0.0000422406,0.0000422870,0.0000422670,0.0000424824,0.0000425325,0.0000427065,0.0000428471,0.0000428945,0.0000429510,0.0000431445,0.0000430826,0.0000433526,0.0000433729,0.0000434123,0.0000435894,0.0000436942,0.0000437681,0.0000439323,0.0000438773,0.0000440479,0.0000439531,0.0000442047,0.0000442358,0.0000443380,0.0000444096,0.0000445284,0.0000444595,0.0000444868,0.0000445268,0.0000447387,0.0000446550,0.0000447109,0.0000448972,0.0000448442,0.0000449900,0.0000449333,0.0000450822,0.0000451399,0.0000452789,0.0000451933,0.0000453675,0.0000453607,0.0000455006,0.0000455561,0.0000454565,0.0000455393,0.0000455838,0.0000457183,0.0000457796,0.0000457447,0.0000458080,0.0000459360,0.0000459008,0.0000458250,0.0000460547,0.0000460834,0.0000460177,0.0000460892,0.0000462441,0.0000461601,0.0000462309,0.0000462664,0.0000463722,0.0000463162,0.0000462435,0.0000464136,0.0000464961,0.0000463904,0.0000321632,0.0000327258,0.0000332518,0.0000335667,0.0000340535,0.0000344119,0.0000346300,0.0000351084,0.0000353790,0.0000357024,0.0000359958,0.0000362857,0.0000365653,0.0000368319,0.0000370779,0.0000372440,0.0000374607,0.0000377501,0.0000379407,0.0000380589,0.0000382409,0.0000383542,0.0000385955,0.0000387035,0.0000389808,0.0000391031,0.0000392704,0.0000393841,0.0000394948,0.0000394851,0.0000398840,0.0000398112,0.0000401229,0.0000401890,0.0000401417,0.0000403486,0.0000404039,0.0000405512,0.0000407014,0.0000407516,0.0000407904,0.0000409881,0.0000409064,0.0000411704,0.0000412194,0.0000412736,0.0000414049,0.0000414906,0.0000415837,0.0000417769,0.0000416791,0.0000418573,0.0000417637,0.0000419651,0.0000420046,0.0000421345,0.0000421765,0.0000423107,0.0000422285,0.0000422703,0.0000423223,0.0000424865,0.0000424170,0.0000424749,0.0000426388,0.0000425996,0.0000427095,0.0000426894,0.0000428225,0.0000429044,0.0000430265,0.0000429298,0.0000431184,0.0000431053,0.0000432081,0.0000433051,0.0000432033,0.0000432700,0.0000433134,0.0000434589,0.0000434892,0.0000434379,0.0000435340,0.0000436722,0.0000436108,0.0000435535,0.0000437754,0.0000437682,0.0000437090,0.0000437846,0.0000439491,0.0000438508,0.0000439378,0.0000439261,0.0000440515,0.0000439818,0.0000439331,0.0000440988,0.0000441850,0.0000440661,0.0000305357,0.0000310743,0.0000315755,0.0000318818,0.0000323182,0.0000326632,0.0000328781,0.0000333390,0.0000335944,0.0000339186,0.0000341781,0.0000344658,0.0000347258,0.0000349830,0.0000351904,0.0000353630,0.0000355729,0.0000358432,0.0000360481,0.0000361653,0.0000363214,0.0000364489,0.0000366431,0.0000367749,0.0000370176,0.0000371408,0.0000373105,0.0000373966,0.0000375289,0.0000375155,0.0000378740,0.0000378168,0.0000381204,0.0000381818,0.0000381311,0.0000383058,0.0000383988,0.0000385269,0.0000386740,0.0000387168,0.0000387529,0.0000389172,0.0000388514,0.0000391274,0.0000391511,0.0000392058,0.0000393260,0.0000394143,0.0000395039,0.0000396754,0.0000396057,0.0000397834,0.0000396833,0.0000398910,0.0000399224,0.0000400246,0.0000400576,0.0000401685,0.0000400939,0.0000401560,0.0000402019,0.0000403568,0.0000402963,0.0000403529,0.0000405038,0.0000404682,0.0000405719,0.0000405629,0.0000406729,0.0000407530,0.0000408588,0.0000407767,0.0000409657,0.0000409642,0.0000410347,0.0000411598,0.0000410270,0.0000411119,0.0000411891,0.0000412987,0.0000413080,0.0000412600,0.0000413652,0.0000415019,0.0000414088,0.0000413992,0.0000415983,0.0000415920,0.0000415244,0.0000416074,0.0000417620,0.0000416694,0.0000417580,0.0000417271,0.0000418661,0.0000417673,0.0000417123,0.0000419027,0.0000420056,0.0000418660,0.0000289941,0.0000295113,0.0000299715,0.0000302685,0.0000306935,0.0000310121,0.0000312262,0.0000316629,0.0000318992,0.0000322263,0.0000324494,0.0000327270,0.0000329731,0.0000332195,0.0000334197,0.0000335928,0.0000337941,0.0000340224,0.0000342331,0.0000343520,0.0000345103,0.0000346164,0.0000348033,0.0000349050,0.0000351574,0.0000352763,0.0000354492,0.0000355349,0.0000356410,0.0000356364,0.0000359872,0.0000359312,0.0000362298,0.0000362592,0.0000362170,0.0000363805,0.0000364866,0.0000366139,0.0000367134,0.0000367817,0.0000368098,0.0000369822,0.0000369193,0.0000371612,0.0000372028,0.0000372522,0.0000373531,0.0000374436,0.0000375376,0.0000376831,0.0000376332,0.0000377667,0.0000376880,0.0000379053,0.0000379387,0.0000380092,0.0000380398,0.0000381422,0.0000380980,0.0000381434,0.0000382053,0.0000383363,0.0000382907,0.0000383551,0.0000384882,0.0000384672,0.0000385397,0.0000385375,0.0000386062,0.0000387235,0.0000388086,0.0000387556,0.0000389259,0.0000388931,0.0000389861,0.0000390993,0.0000389703,0.0000390519,0.0000391296,0.0000392538,0.0000392598,0.0000391829,0.0000392957,0.0000394508,0.0000393197,0.0000393209,0.0000395359,0.0000395099,0.0000394689,0.0000395131,0.0000396846,0.0000396179,0.0000396822,0.0000396273,0.0000397592,0.0000397029,0.0000396356,0.0000397992,0.0000399004,0.0000397676,0.0000275294,0.0000280220,0.0000284612,0.0000287258,0.0000291324,0.0000294431,0.0000296328,0.0000300648,0.0000303129,0.0000305999,0.0000308286,0.0000310599,0.0000313184,0.0000315558,0.0000317609,0.0000318882,0.0000320762,0.0000323159,0.0000325487,0.0000326201,0.0000327727,0.0000328896,0.0000330541,0.0000331451,0.0000334162,0.0000335218,0.0000336766,0.0000337409,0.0000338584,0.0000338498,0.0000341896,0.0000341252,0.0000344155,0.0000344587,0.0000344072,0.0000345460,0.0000346199,0.0000347792,0.0000348339,0.0000349473,0.0000349720,0.0000351273,0.0000350791,0.0000352987,0.0000353241,0.0000353864,0.0000354819,0.0000355725,0.0000356782,0.0000357835,0.0000357325,0.0000358703,0.0000358113,0.0000359959,0.0000360127,0.0000361192,0.0000361406,0.0000362465,0.0000361903,0.0000362383,0.0000363051,0.0000364424,0.0000363941,0.0000364384,0.0000365599,0.0000365655,0.0000366022,0.0000366119,0.0000366831,0.0000368074,0.0000368659,0.0000368202,0.0000369861,0.0000369586,0.0000370465,0.0000371563,0.0000370329,0.0000371120,0.0000371659,0.0000372692,0.0000373166,0.0000372333,0.0000373171,0.0000374659,0.0000373734,0.0000373665,0.0000375496,0.0000375130,0.0000375011,0.0000375122,0.0000377120,0.0000376501,0.0000377049,0.0000376295,0.0000377763,0.0000377248,0.0000376497,0.0000378303,0.0000379225,0.0000377725,0.0000261364,0.0000266032,0.0000270361,0.0000272874,0.0000276490,0.0000279634,0.0000281568,0.0000285599,0.0000287951,0.0000290586,0.0000292725,0.0000294836,0.0000297425,0.0000299474,0.0000301677,0.0000302569,0.0000304593,0.0000306915,0.0000309068,0.0000309713,0.0000311144,0.0000312294,0.0000313689,0.0000314624,0.0000317400,0.0000318586,0.0000320000,0.0000320601,0.0000321972,0.0000321411,0.0000324793,0.0000324083,0.0000326765,0.0000327455,0.0000326741,0.0000328087,0.0000328865,0.0000330179,0.0000330936,0.0000331993,0.0000332359,0.0000333529,0.0000333137,0.0000335163,0.0000335783,0.0000336044,0.0000336864,0.0000338004,0.0000338984,0.0000339917,0.0000339517,0.0000340718,0.0000340394,0.0000342013,0.0000342069,0.0000343098,0.0000343336,0.0000344296,0.0000343763,0.0000344253,0.0000344797,0.0000345986,0.0000345769,0.0000346216,0.0000347029,0.0000347301,0.0000347873,0.0000347976,0.0000348465,0.0000349806,0.0000350335,0.0000349825,0.0000351400,0.0000351049,0.0000351774,0.0000353162,0.0000351789,0.0000352576,0.0000353092,0.0000354290,0.0000354659,0.0000353830,0.0000354612,0.0000356025,0.0000355149,0.0000355230,0.0000356942,0.0000356585,0.0000356657,0.0000356288,0.0000358269,0.0000357684,0.0000358200,0.0000357626,0.0000359164,0.0000358605,0.0000357538,0.0000359573,0.0000360516,0.0000358891,0.0000248296,0.0000252596,0.0000256683,0.0000259013,0.0000262573,0.0000265567,0.0000267097,0.0000271085,0.0000273507,0.0000276132,0.0000277997,0.0000279786,0.0000282314,0.0000284609,0.0000286435,0.0000287295,0.0000289557,0.0000291500,0.0000293482,0.0000294183,0.0000295516,0.0000296637,0.0000297781,0.0000298777,0.0000301329,0.0000302779,0.0000303935,0.0000304476,0.0000305630,0.0000305285,0.0000308528,0.0000307516,0.0000310412,0.0000311114,0.0000310400,0.0000311735,0.0000312496,0.0000313684,0.0000314570,0.0000315489,0.0000315608,0.0000316860,0.0000316499,0.0000318407,0.0000318971,0.0000319167,0.0000319946,0.0000320984,0.0000322038,0.0000323207,0.0000322693,0.0000323765,0.0000323379,0.0000325195,0.0000324856,0.0000326036,0.0000326230,0.0000327166,0.0000326688,0.0000327226,0.0000327561,0.0000328438,0.0000328369,0.0000329231,0.0000329759,0.0000330125,0.0000330652,0.0000330778,0.0000331290,0.0000332261,0.0000332763,0.0000332251,0.0000333884,0.0000333499,0.0000334321,0.0000335552,0.0000334103,0.0000335255,0.0000335596,0.0000336828,0.0000336852,0.0000336139,0.0000336941,0.0000338364,0.0000337345,0.0000337604,0.0000339383,0.0000338894,0.0000338841,0.0000338546,0.0000340366,0.0000339783,0.0000340339,0.0000339884,0.0000341104,0.0000340653,0.0000339688,0.0000341874,0.0000342564,0.0000340950,0.0000235759,0.0000239873,0.0000243757,0.0000245962,0.0000249262,0.0000251910,0.0000253556,0.0000257385,0.0000259713,0.0000262199,0.0000263986,0.0000265685,0.0000267905,0.0000270338,0.0000272095,0.0000272741,0.0000274899,0.0000276793,0.0000278583,0.0000279331,0.0000280668,0.0000281647,0.0000282758,0.0000283598,0.0000286262,0.0000287645,0.0000288601,0.0000289088,0.0000290393,0.0000290006,0.0000293104,0.0000292146,0.0000294963,0.0000295597,0.0000294640,0.0000296102,0.0000296730,0.0000297949,0.0000298556,0.0000299927,0.0000299819,0.0000301018,0.0000300494,0.0000302510,0.0000303019,0.0000303088,0.0000304071,0.0000304796,0.0000306079,0.0000307095,0.0000306626,0.0000307719,0.0000307213,0.0000308812,0.0000308720,0.0000309496,0.0000309874,0.0000310878,0.0000310494,0.0000311143,0.0000311171,0.0000312179,0.0000312119,0.0000312820,0.0000313226,0.0000313646,0.0000314102,0.0000314346,0.0000314713,0.0000315636,0.0000316057,0.0000315499,0.0000317135,0.0000316998,0.0000317699,0.0000318673,0.0000317605,0.0000318553,0.0000318832,0.0000319851,0.0000320218,0.0000319452,0.0000320262,0.0000321553,0.0000320340,0.0000320581,0.0000322276,0.0000322129,0.0000321804,0.0000321703,0.0000323479,0.0000322443,0.0000323515,0.0000323175,0.0000323914,0.0000323669,0.0000322526,0.0000324963,0.0000325646,0.0000323942,0.0000223851,0.0000227691,0.0000231445,0.0000233488,0.0000236710,0.0000239275,0.0000240694,0.0000244334,0.0000246682,0.0000248891,0.0000250479,0.0000252279,0.0000254391,0.0000256686,0.0000258639,0.0000259359,0.0000261208,0.0000262578,0.0000264266,0.0000265344,0.0000266390,0.0000267426,0.0000268542,0.0000269379,0.0000271983,0.0000273272,0.0000274040,0.0000274669,0.0000275871,0.0000275596,0.0000278338,0.0000277670,0.0000280198,0.0000280906,0.0000279751,0.0000281226,0.0000281699,0.0000283001,0.0000283660,0.0000284842,0.0000284850,0.0000285914,0.0000285529,0.0000287324,0.0000287871,0.0000287824,0.0000288920,0.0000289467,0.0000290889,0.0000291738,0.0000291040,0.0000292296,0.0000292183,0.0000293253,0.0000293202,0.0000294035,0.0000294395,0.0000295351,0.0000294982,0.0000295537,0.0000295814,0.0000296535,0.0000296680,0.0000297254,0.0000297787,0.0000298184,0.0000298224,0.0000298796,0.0000299043,0.0000300031,0.0000300366,0.0000299700,0.0000301051,0.0000301279,0.0000301779,0.0000302833,0.0000302025,0.0000302921,0.0000302815,0.0000303905,0.0000304327,0.0000303414,0.0000304468,0.0000305474,0.0000304439,0.0000304512,0.0000306299,0.0000306119,0.0000305634,0.0000305697,0.0000307401,0.0000306411,0.0000307499,0.0000307188,0.0000307499,0.0000307626,0.0000306484,0.0000308641,0.0000309201,0.0000307807,0.0000212640,0.0000216147,0.0000219788,0.0000221716,0.0000225034,0.0000227212,0.0000228746,0.0000232002,0.0000234520,0.0000236358,0.0000237833,0.0000239626,0.0000241563,0.0000243689,0.0000245770,0.0000246350,0.0000248252,0.0000249327,0.0000250941,0.0000251978,0.0000253101,0.0000254025,0.0000255211,0.0000255914,0.0000258379,0.0000259467,0.0000260294,0.0000260873,0.0000262021,0.0000261635,0.0000264395,0.0000263637,0.0000266033,0.0000266825,0.0000265700,0.0000267058,0.0000267602,0.0000268731,0.0000269429,0.0000270663,0.0000270544,0.0000271379,0.0000271256,0.0000272914,0.0000273541,0.0000273547,0.0000274477,0.0000274984,0.0000276275,0.0000277075,0.0000276424,0.0000277978,0.0000277729,0.0000278657,0.0000278593,0.0000279516,0.0000279710,0.0000280741,0.0000279991,0.0000280771,0.0000281063,0.0000281606,0.0000281987,0.0000282494,0.0000282950,0.0000283464,0.0000283251,0.0000283942,0.0000284174,0.0000284952,0.0000285441,0.0000284757,0.0000286198,0.0000286256,0.0000286847,0.0000287588,0.0000286833,0.0000287792,0.0000287748,0.0000288616,0.0000289011,0.0000288359,0.0000289320,0.0000290234,0.0000289290,0.0000289394,0.0000291007,0.0000290798,0.0000290373,0.0000290283,0.0000291958,0.0000290976,0.0000292265,0.0000291674,0.0000292366,0.0000292300,0.0000291245,0.0000293154,0.0000293764,0.0000292738,0.0000201917,0.0000205231,0.0000208618,0.0000210678,0.0000213654,0.0000215666,0.0000217343,0.0000220379,0.0000222693,0.0000224579,0.0000225856,0.0000227736,0.0000229485,0.0000231497,0.0000233323,0.0000233819,0.0000235732,0.0000236726,0.0000238441,0.0000239399,0.0000240240,0.0000241272,0.0000242313,0.0000243174,0.0000245363,0.0000246202,0.0000247323,0.0000247711,0.0000248684,0.0000248550,0.0000251157,0.0000250399,0.0000252747,0.0000253574,0.0000252428,0.0000253762,0.0000254401,0.0000255267,0.0000255958,0.0000257312,0.0000256789,0.0000257433,0.0000257573,0.0000259479,0.0000259824,0.0000259676,0.0000260626,0.0000261308,0.0000262448,0.0000263229,0.0000262479,0.0000264054,0.0000263885,0.0000264485,0.0000264587,0.0000265549,0.0000265683,0.0000266786,0.0000266122,0.0000266683,0.0000267005,0.0000267475,0.0000267846,0.0000268344,0.0000268702,0.0000269404,0.0000269060,0.0000269657,0.0000270169,0.0000270679,0.0000271077,0.0000270547,0.0000271791,0.0000272131,0.0000272313,0.0000273266,0.0000272572,0.0000273496,0.0000273324,0.0000274185,0.0000274838,0.0000273925,0.0000274859,0.0000276000,0.0000274911,0.0000274950,0.0000276419,0.0000276299,0.0000275868,0.0000275750,0.0000277325,0.0000276491,0.0000277773,0.0000277155,0.0000277771,0.0000277714,0.0000276548,0.0000278416,0.0000279034,0.0000278025,0.0000191597,0.0000194890,0.0000198034,0.0000199930,0.0000202909,0.0000204752,0.0000206386,0.0000209315,0.0000211406,0.0000213300,0.0000214300,0.0000216277,0.0000218014,0.0000219865,0.0000221831,0.0000222133,0.0000224010,0.0000224837,0.0000226675,0.0000227546,0.0000228105,0.0000229195,0.0000229996,0.0000230944,0.0000233082,0.0000234049,0.0000235074,0.0000235190,0.0000236379,0.0000236197,0.0000238355,0.0000237829,0.0000240107,0.0000240990,0.0000239811,0.0000240942,0.0000241646,0.0000242424,0.0000243231,0.0000244490,0.0000243857,0.0000244534,0.0000244764,0.0000246539,0.0000246909,0.0000246643,0.0000247608,0.0000248199,0.0000249374,0.0000249958,0.0000249419,0.0000250940,0.0000250766,0.0000251389,0.0000251259,0.0000252226,0.0000252495,0.0000253504,0.0000252809,0.0000253230,0.0000253627,0.0000254006,0.0000254424,0.0000254941,0.0000255326,0.0000256155,0.0000255684,0.0000256258,0.0000256682,0.0000257006,0.0000257505,0.0000257169,0.0000258205,0.0000258496,0.0000258749,0.0000259728,0.0000259031,0.0000259743,0.0000259800,0.0000260457,0.0000261169,0.0000260186,0.0000261027,0.0000262300,0.0000261238,0.0000261378,0.0000262648,0.0000262491,0.0000261897,0.0000262306,0.0000263511,0.0000262728,0.0000263913,0.0000263317,0.0000264133,0.0000263822,0.0000262843,0.0000264739,0.0000265079,0.0000264295,0.0000181938,0.0000185027,0.0000188105,0.0000189777,0.0000192561,0.0000194505,0.0000196035,0.0000198569,0.0000200656,0.0000202711,0.0000203636,0.0000205412,0.0000207105,0.0000208939,0.0000210512,0.0000210974,0.0000212647,0.0000213491,0.0000215370,0.0000216076,0.0000216729,0.0000217682,0.0000218500,0.0000219481,0.0000221504,0.0000222212,0.0000223465,0.0000223357,0.0000224612,0.0000224364,0.0000226386,0.0000225945,0.0000228049,0.0000228763,0.0000227851,0.0000228909,0.0000229660,0.0000230147,0.0000231009,0.0000232269,0.0000231897,0.0000232201,0.0000232289,0.0000234303,0.0000234501,0.0000234257,0.0000235277,0.0000235861,0.0000236810,0.0000237344,0.0000237071,0.0000238356,0.0000238355,0.0000238949,0.0000238665,0.0000239649,0.0000239954,0.0000240800,0.0000240138,0.0000240574,0.0000240963,0.0000241375,0.0000241640,0.0000242272,0.0000242642,0.0000243377,0.0000242824,0.0000243372,0.0000243800,0.0000244386,0.0000244487,0.0000244267,0.0000245322,0.0000245571,0.0000245919,0.0000246747,0.0000246070,0.0000246982,0.0000247087,0.0000247275,0.0000247938,0.0000247230,0.0000247976,0.0000249242,0.0000248432,0.0000248181,0.0000249515,0.0000249166,0.0000248953,0.0000249226,0.0000250383,0.0000249458,0.0000250900,0.0000250170,0.0000250914,0.0000250755,0.0000249775,0.0000251498,0.0000251986,0.0000251255,0.0000172771,0.0000175680,0.0000178637,0.0000180137,0.0000182850,0.0000184781,0.0000186106,0.0000188569,0.0000190560,0.0000192388,0.0000193462,0.0000195080,0.0000196772,0.0000198555,0.0000199928,0.0000200350,0.0000201765,0.0000202661,0.0000204712,0.0000205182,0.0000205842,0.0000206740,0.0000207473,0.0000208531,0.0000210378,0.0000210774,0.0000212357,0.0000212047,0.0000213454,0.0000213061,0.0000215087,0.0000214548,0.0000216709,0.0000217370,0.0000216431,0.0000217295,0.0000218158,0.0000218523,0.0000219381,0.0000220614,0.0000220381,0.0000220676,0.0000220652,0.0000222476,0.0000222844,0.0000222446,0.0000223394,0.0000224089,0.0000224908,0.0000225799,0.0000225214,0.0000226655,0.0000226410,0.0000227032,0.0000226787,0.0000227503,0.0000227832,0.0000228858,0.0000228149,0.0000228641,0.0000228959,0.0000229347,0.0000229616,0.0000230114,0.0000230443,0.0000231199,0.0000230904,0.0000231129,0.0000231671,0.0000232184,0.0000232457,0.0000231997,0.0000233043,0.0000233542,0.0000233561,0.0000234341,0.0000233783,0.0000234637,0.0000234861,0.0000234992,0.0000235592,0.0000234975,0.0000235419,0.0000236920,0.0000235975,0.0000235868,0.0000236995,0.0000236890,0.0000236331,0.0000236803,0.0000237724,0.0000237114,0.0000238472,0.0000237548,0.0000238425,0.0000238175,0.0000237393,0.0000238981,0.0000239497,0.0000238695,0.0000163962,0.0000166914,0.0000169500,0.0000171195,0.0000173546,0.0000175348,0.0000176842,0.0000179101,0.0000181040,0.0000182634,0.0000183737,0.0000185198,0.0000186860,0.0000188539,0.0000190017,0.0000190081,0.0000191448,0.0000192491,0.0000194507,0.0000194997,0.0000195397,0.0000196382,0.0000196950,0.0000198172,0.0000199830,0.0000200256,0.0000201626,0.0000201418,0.0000202801,0.0000202344,0.0000204403,0.0000203796,0.0000205794,0.0000206542,0.0000205692,0.0000206579,0.0000207253,0.0000207588,0.0000208522,0.0000209727,0.0000209269,0.0000209644,0.0000209629,0.0000211271,0.0000211724,0.0000211250,0.0000212216,0.0000212884,0.0000213489,0.0000214687,0.0000213798,0.0000215444,0.0000215235,0.0000215875,0.0000215514,0.0000216131,0.0000216495,0.0000217348,0.0000216692,0.0000217151,0.0000217506,0.0000217906,0.0000218131,0.0000218631,0.0000218975,0.0000219753,0.0000219431,0.0000219405,0.0000220116,0.0000220607,0.0000220932,0.0000220524,0.0000221231,0.0000221815,0.0000221826,0.0000222717,0.0000222141,0.0000222896,0.0000223298,0.0000223167,0.0000223784,0.0000223439,0.0000223928,0.0000225038,0.0000224239,0.0000224125,0.0000225349,0.0000225055,0.0000224610,0.0000225093,0.0000225833,0.0000225198,0.0000226369,0.0000225697,0.0000226699,0.0000226328,0.0000225627,0.0000227328,0.0000227616,0.0000226633,0.0000155779,0.0000158466,0.0000160894,0.0000162603,0.0000164588,0.0000166561,0.0000168120,0.0000170068,0.0000172063,0.0000173563,0.0000174478,0.0000175949,0.0000177684,0.0000179135,0.0000180550,0.0000180622,0.0000181959,0.0000182768,0.0000184736,0.0000185348,0.0000185405,0.0000186548,0.0000186875,0.0000188296,0.0000189747,0.0000190208,0.0000191579,0.0000191377,0.0000192637,0.0000192223,0.0000194145,0.0000193667,0.0000195366,0.0000196151,0.0000195527,0.0000196191,0.0000196846,0.0000197100,0.0000198030,0.0000199190,0.0000198627,0.0000198963,0.0000199156,0.0000200775,0.0000200990,0.0000200552,0.0000201441,0.0000202263,0.0000202653,0.0000204175,0.0000203116,0.0000204559,0.0000204433,0.0000205129,0.0000204742,0.0000205258,0.0000205689,0.0000206444,0.0000205799,0.0000206245,0.0000206773,0.0000207103,0.0000207381,0.0000207706,0.0000208032,0.0000208986,0.0000208467,0.0000208463,0.0000209229,0.0000209463,0.0000209964,0.0000209445,0.0000210276,0.0000210692,0.0000210549,0.0000211568,0.0000211090,0.0000211942,0.0000212255,0.0000212020,0.0000212417,0.0000212171,0.0000212747,0.0000213829,0.0000213205,0.0000212994,0.0000214114,0.0000213718,0.0000213442,0.0000213963,0.0000214414,0.0000214053,0.0000215236,0.0000214590,0.0000215215,0.0000215107,0.0000214368,0.0000216068,0.0000216460,0.0000215416,0.0000147969,0.0000150528,0.0000152936,0.0000154488,0.0000156343,0.0000158103,0.0000159624,0.0000161505,0.0000163450,0.0000164999,0.0000165499,0.0000167066,0.0000168644,0.0000170143,0.0000171444,0.0000171643,0.0000172879,0.0000173582,0.0000175390,0.0000175964,0.0000176124,0.0000177301,0.0000177508,0.0000178778,0.0000180371,0.0000180629,0.0000181886,0.0000181745,0.0000183003,0.0000182653,0.0000184286,0.0000183958,0.0000185458,0.0000186235,0.0000185861,0.0000186310,0.0000186986,0.0000187043,0.0000188176,0.0000189356,0.0000188589,0.0000189227,0.0000189212,0.0000190702,0.0000190800,0.0000190537,0.0000191376,0.0000191941,0.0000192319,0.0000193963,0.0000193125,0.0000194309,0.0000194240,0.0000194828,0.0000194497,0.0000194891,0.0000195362,0.0000196193,0.0000195638,0.0000195925,0.0000196274,0.0000196813,0.0000196934,0.0000197454,0.0000197751,0.0000198548,0.0000198066,0.0000198141,0.0000198773,0.0000199156,0.0000199529,0.0000199085,0.0000199792,0.0000200251,0.0000200030,0.0000200942,0.0000200616,0.0000201338,0.0000201676,0.0000201272,0.0000201752,0.0000201672,0.0000202071,0.0000203223,0.0000202419,0.0000202330,0.0000203459,0.0000203251,0.0000202835,0.0000203281,0.0000203704,0.0000203411,0.0000204591,0.0000204008,0.0000204548,0.0000204544,0.0000203490,0.0000205395,0.0000205488,0.0000204698,0.0000140623,0.0000142937,0.0000145022,0.0000146719,0.0000148584,0.0000150165,0.0000151586,0.0000153389,0.0000155216,0.0000156767,0.0000157187,0.0000158493,0.0000160121,0.0000161683,0.0000162836,0.0000162975,0.0000164254,0.0000164901,0.0000166587,0.0000167148,0.0000167247,0.0000168406,0.0000168512,0.0000169693,0.0000171414,0.0000171491,0.0000172774,0.0000172625,0.0000173967,0.0000173461,0.0000175001,0.0000174842,0.0000176121,0.0000176885,0.0000176460,0.0000176958,0.0000177664,0.0000177812,0.0000178703,0.0000179904,0.0000179035,0.0000179853,0.0000179725,0.0000181174,0.0000181542,0.0000181107,0.0000181751,0.0000182265,0.0000182691,0.0000184294,0.0000183397,0.0000184607,0.0000184619,0.0000185185,0.0000184623,0.0000185173,0.0000185585,0.0000186466,0.0000186017,0.0000186257,0.0000186311,0.0000187007,0.0000186994,0.0000187545,0.0000187857,0.0000188768,0.0000188108,0.0000188214,0.0000188699,0.0000189321,0.0000189788,0.0000189105,0.0000189952,0.0000190257,0.0000190040,0.0000190924,0.0000190601,0.0000191409,0.0000191709,0.0000191058,0.0000191736,0.0000191597,0.0000191892,0.0000193036,0.0000192256,0.0000192358,0.0000193198,0.0000193265,0.0000192607,0.0000193166,0.0000193407,0.0000193228,0.0000194424,0.0000193826,0.0000194325,0.0000194272,0.0000193408,0.0000195238,0.0000195272,0.0000194700,0.0000133277,0.0000135710,0.0000137711,0.0000139360,0.0000141165,0.0000142659,0.0000144001,0.0000145563,0.0000147477,0.0000148942,0.0000149300,0.0000150485,0.0000151979,0.0000153823,0.0000154600,0.0000154551,0.0000156166,0.0000156536,0.0000158184,0.0000158871,0.0000158815,0.0000160116,0.0000159985,0.0000161195,0.0000162593,0.0000162995,0.0000164020,0.0000163986,0.0000165292,0.0000164842,0.0000166272,0.0000166033,0.0000167325,0.0000167970,0.0000167784,0.0000168125,0.0000168811,0.0000168741,0.0000169754,0.0000170883,0.0000170152,0.0000170794,0.0000170711,0.0000172128,0.0000172585,0.0000172077,0.0000172683,0.0000173200,0.0000173536,0.0000175083,0.0000174394,0.0000175369,0.0000175332,0.0000176067,0.0000175417,0.0000175940,0.0000176478,0.0000177093,0.0000176848,0.0000176912,0.0000176933,0.0000177562,0.0000177618,0.0000178227,0.0000178528,0.0000179345,0.0000178771,0.0000178739,0.0000179226,0.0000179917,0.0000180258,0.0000179718,0.0000180485,0.0000180808,0.0000180565,0.0000181530,0.0000180903,0.0000181967,0.0000181952,0.0000181513,0.0000182296,0.0000181964,0.0000182288,0.0000183360,0.0000182641,0.0000182760,0.0000183673,0.0000183444,0.0000182853,0.0000183505,0.0000183705,0.0000183454,0.0000184603,0.0000184194,0.0000184702,0.0000184516,0.0000183768,0.0000185565,0.0000185497,0.0000184939,0.0000126550,0.0000128843,0.0000130874,0.0000132335,0.0000134197,0.0000135482,0.0000136753,0.0000138183,0.0000140117,0.0000141390,0.0000141741,0.0000142939,0.0000144410,0.0000145963,0.0000146787,0.0000146781,0.0000148314,0.0000148599,0.0000150229,0.0000150844,0.0000150807,0.0000152057,0.0000151938,0.0000153065,0.0000154546,0.0000154849,0.0000155848,0.0000155658,0.0000156918,0.0000156556,0.0000157982,0.0000157743,0.0000158985,0.0000159472,0.0000159262,0.0000159769,0.0000160436,0.0000160238,0.0000161338,0.0000162390,0.0000161639,0.0000162366,0.0000162234,0.0000163494,0.0000163977,0.0000163547,0.0000164230,0.0000164504,0.0000164903,0.0000166372,0.0000165613,0.0000166388,0.0000166546,0.0000167432,0.0000166587,0.0000167140,0.0000167701,0.0000168122,0.0000168076,0.0000168069,0.0000168105,0.0000168715,0.0000168625,0.0000169427,0.0000169663,0.0000170422,0.0000170011,0.0000169826,0.0000170326,0.0000170992,0.0000171115,0.0000170593,0.0000171521,0.0000171865,0.0000171384,0.0000172463,0.0000171876,0.0000172923,0.0000172893,0.0000172352,0.0000173209,0.0000172858,0.0000173179,0.0000174225,0.0000173681,0.0000173596,0.0000174659,0.0000174300,0.0000173940,0.0000174410,0.0000174499,0.0000174315,0.0000175321,0.0000175130,0.0000175248,0.0000175148,0.0000174539,0.0000176388,0.0000176231,0.0000175706,0.0000120125,0.0000122363,0.0000124138,0.0000125795,0.0000127491,0.0000128698,0.0000129750,0.0000131284,0.0000133126,0.0000134332,0.0000134585,0.0000135681,0.0000137009,0.0000138719,0.0000139404,0.0000139330,0.0000140802,0.0000141042,0.0000142758,0.0000143177,0.0000143308,0.0000144484,0.0000144299,0.0000145465,0.0000146834,0.0000147004,0.0000148139,0.0000147943,0.0000149062,0.0000148555,0.0000150056,0.0000149907,0.0000151001,0.0000151448,0.0000151139,0.0000151763,0.0000152232,0.0000152218,0.0000153315,0.0000154413,0.0000153622,0.0000154207,0.0000154143,0.0000155463,0.0000155821,0.0000155273,0.0000156091,0.0000156284,0.0000156769,0.0000158061,0.0000157372,0.0000158215,0.0000158077,0.0000159066,0.0000158321,0.0000158971,0.0000159372,0.0000159824,0.0000159584,0.0000159590,0.0000159715,0.0000160322,0.0000160122,0.0000161029,0.0000161469,0.0000161911,0.0000161650,0.0000161432,0.0000161967,0.0000162252,0.0000162513,0.0000162078,0.0000162884,0.0000163380,0.0000162824,0.0000163976,0.0000163197,0.0000164281,0.0000164313,0.0000163714,0.0000164663,0.0000164189,0.0000164673,0.0000165559,0.0000164984,0.0000165015,0.0000165882,0.0000165632,0.0000165169,0.0000165645,0.0000165822,0.0000165528,0.0000166623,0.0000166086,0.0000166670,0.0000166420,0.0000165836,0.0000167602,0.0000167417,0.0000167016,0.0000114136,0.0000116153,0.0000117885,0.0000119508,0.0000121129,0.0000122053,0.0000123060,0.0000124754,0.0000126498,0.0000127686,0.0000127884,0.0000128815,0.0000130179,0.0000131805,0.0000132496,0.0000132385,0.0000133763,0.0000134024,0.0000135549,0.0000136054,0.0000136287,0.0000137319,0.0000136957,0.0000138309,0.0000139472,0.0000139636,0.0000140688,0.0000140581,0.0000141686,0.0000141027,0.0000142556,0.0000142351,0.0000143468,0.0000143851,0.0000143521,0.0000144271,0.0000144679,0.0000144554,0.0000145774,0.0000146630,0.0000145986,0.0000146418,0.0000146548,0.0000147625,0.0000148061,0.0000147669,0.0000148239,0.0000148475,0.0000148978,0.0000150051,0.0000149464,0.0000150345,0.0000150247,0.0000151157,0.0000150393,0.0000151051,0.0000151411,0.0000151822,0.0000151712,0.0000151548,0.0000151739,0.0000152230,0.0000152170,0.0000152982,0.0000153527,0.0000153850,0.0000153643,0.0000153340,0.0000154021,0.0000154089,0.0000154388,0.0000154047,0.0000154861,0.0000155207,0.0000154641,0.0000155768,0.0000154974,0.0000155963,0.0000156079,0.0000155575,0.0000156485,0.0000155857,0.0000156340,0.0000157285,0.0000156660,0.0000156718,0.0000157647,0.0000157360,0.0000157097,0.0000157393,0.0000157659,0.0000157391,0.0000158352,0.0000157890,0.0000158214,0.0000158078,0.0000157427,0.0000159284,0.0000159070,0.0000158772,0.0000108381,0.0000110228,0.0000111858,0.0000113520,0.0000114978,0.0000115933,0.0000116929,0.0000118430,0.0000120047,0.0000121151,0.0000121390,0.0000122346,0.0000123695,0.0000125104,0.0000125971,0.0000125792,0.0000126926,0.0000127476,0.0000128657,0.0000129111,0.0000129286,0.0000130424,0.0000130098,0.0000131260,0.0000132545,0.0000132686,0.0000133546,0.0000133464,0.0000134617,0.0000133835,0.0000135448,0.0000135248,0.0000136165,0.0000136726,0.0000136306,0.0000137006,0.0000137346,0.0000137390,0.0000138345,0.0000139508,0.0000138598,0.0000139114,0.0000139096,0.0000140183,0.0000140721,0.0000140215,0.0000140652,0.0000140882,0.0000141644,0.0000142595,0.0000141914,0.0000142719,0.0000142867,0.0000143650,0.0000142878,0.0000143512,0.0000143916,0.0000144292,0.0000144164,0.0000143943,0.0000144258,0.0000144574,0.0000144578,0.0000145457,0.0000145940,0.0000146229,0.0000146027,0.0000145548,0.0000146324,0.0000146314,0.0000146565,0.0000146389,0.0000147085,0.0000147622,0.0000146991,0.0000147916,0.0000147270,0.0000148055,0.0000148283,0.0000147706,0.0000148815,0.0000148169,0.0000148655,0.0000149493,0.0000148864,0.0000148971,0.0000149901,0.0000149562,0.0000149376,0.0000149459,0.0000149766,0.0000149557,0.0000150497,0.0000149992,0.0000150231,0.0000150301,0.0000149398,0.0000151484,0.0000151144,0.0000150865,0.0000103059,0.0000104635,0.0000106138,0.0000107890,0.0000109349,0.0000110148,0.0000111092,0.0000112508,0.0000114131,0.0000115065,0.0000115182,0.0000116312,0.0000117538,0.0000118751,0.0000119700,0.0000119506,0.0000120512,0.0000121027,0.0000122254,0.0000122602,0.0000122689,0.0000123948,0.0000123553,0.0000124719,0.0000126021,0.0000126064,0.0000126681,0.0000126847,0.0000127847,0.0000127062,0.0000128622,0.0000128469,0.0000129324,0.0000129989,0.0000129383,0.0000130128,0.0000130313,0.0000130432,0.0000131389,0.0000132508,0.0000131756,0.0000132317,0.0000131965,0.0000133113,0.0000133703,0.0000133222,0.0000133741,0.0000133687,0.0000134655,0.0000135428,0.0000134930,0.0000135454,0.0000135751,0.0000136519,0.0000135795,0.0000136334,0.0000136788,0.0000137079,0.0000137067,0.0000136771,0.0000137120,0.0000137308,0.0000137360,0.0000138256,0.0000138642,0.0000138943,0.0000138712,0.0000138297,0.0000139173,0.0000139044,0.0000139204,0.0000139125,0.0000139738,0.0000140324,0.0000139575,0.0000140594,0.0000139913,0.0000140711,0.0000140920,0.0000140305,0.0000141215,0.0000140802,0.0000141318,0.0000142052,0.0000141426,0.0000141431,0.0000142554,0.0000142019,0.0000142025,0.0000141872,0.0000142301,0.0000142100,0.0000142944,0.0000142508,0.0000142773,0.0000142699,0.0000141777,0.0000144009,0.0000143587,0.0000143410,0.0000097847,0.0000099364,0.0000100796,0.0000102466,0.0000103761,0.0000104619,0.0000105503,0.0000106921,0.0000108392,0.0000109204,0.0000109331,0.0000110473,0.0000111686,0.0000112783,0.0000113648,0.0000113490,0.0000114632,0.0000115100,0.0000116271,0.0000116446,0.0000116352,0.0000117876,0.0000117337,0.0000118487,0.0000119742,0.0000119752,0.0000120313,0.0000120508,0.0000121432,0.0000120711,0.0000122322,0.0000122108,0.0000122791,0.0000123600,0.0000122965,0.0000123586,0.0000123793,0.0000124001,0.0000124657,0.0000125960,0.0000125254,0.0000125824,0.0000125331,0.0000126583,0.0000127099,0.0000126649,0.0000127176,0.0000126967,0.0000127913,0.0000128687,0.0000128033,0.0000128680,0.0000128885,0.0000129682,0.0000129160,0.0000129558,0.0000130033,0.0000130147,0.0000130131,0.0000129980,0.0000130313,0.0000130384,0.0000130434,0.0000131281,0.0000131772,0.0000132006,0.0000131751,0.0000131397,0.0000132373,0.0000132066,0.0000132243,0.0000132132,0.0000132790,0.0000133184,0.0000132669,0.0000133671,0.0000132899,0.0000133672,0.0000133837,0.0000133377,0.0000134235,0.0000133721,0.0000134183,0.0000134916,0.0000134283,0.0000134354,0.0000135412,0.0000134834,0.0000134973,0.0000134824,0.0000135144,0.0000135103,0.0000135604,0.0000135288,0.0000135817,0.0000135481,0.0000134875,0.0000136772,0.0000136534,0.0000136247,0.0000092999,0.0000094378,0.0000095856,0.0000097287,0.0000098548,0.0000099417,0.0000100133,0.0000101476,0.0000102994,0.0000103610,0.0000103935,0.0000104930,0.0000106021,0.0000107188,0.0000107836,0.0000107737,0.0000108875,0.0000109279,0.0000110505,0.0000110580,0.0000110504,0.0000111927,0.0000111377,0.0000112481,0.0000113828,0.0000113739,0.0000114254,0.0000114564,0.0000115345,0.0000114792,0.0000116123,0.0000115986,0.0000116546,0.0000117497,0.0000116887,0.0000117420,0.0000117518,0.0000117687,0.0000118408,0.0000119621,0.0000118984,0.0000119557,0.0000119013,0.0000120357,0.0000120648,0.0000120375,0.0000120796,0.0000120562,0.0000121545,0.0000122252,0.0000121822,0.0000122192,0.0000122362,0.0000123170,0.0000122790,0.0000123088,0.0000123578,0.0000123766,0.0000123630,0.0000123514,0.0000123836,0.0000124011,0.0000124015,0.0000124925,0.0000125241,0.0000125574,0.0000125259,0.0000124694,0.0000125730,0.0000125626,0.0000125649,0.0000125486,0.0000126344,0.0000126570,0.0000126165,0.0000127095,0.0000126313,0.0000126953,0.0000127093,0.0000126727,0.0000127555,0.0000127125,0.0000127454,0.0000128231,0.0000127759,0.0000127706,0.0000128554,0.0000128206,0.0000128196,0.0000128182,0.0000128636,0.0000128356,0.0000128752,0.0000128503,0.0000129033,0.0000128833,0.0000128172,0.0000129879,0.0000129636,0.0000129438,0.0000088288,0.0000089664,0.0000090999,0.0000092331,0.0000093680,0.0000094411,0.0000095111,0.0000096378,0.0000097826,0.0000098385,0.0000098746,0.0000099674,0.0000100802,0.0000101871,0.0000102492,0.0000102241,0.0000103536,0.0000103877,0.0000105051,0.0000105025,0.0000104882,0.0000106366,0.0000105827,0.0000106772,0.0000108184,0.0000108046,0.0000108648,0.0000108857,0.0000109629,0.0000109008,0.0000110426,0.0000110190,0.0000110539,0.0000111597,0.0000111073,0.0000111559,0.0000111692,0.0000111757,0.0000112540,0.0000113743,0.0000113132,0.0000113595,0.0000113118,0.0000114297,0.0000114627,0.0000114274,0.0000114801,0.0000114569,0.0000115378,0.0000116164,0.0000115631,0.0000116125,0.0000116190,0.0000116953,0.0000116571,0.0000116873,0.0000117292,0.0000117651,0.0000117499,0.0000117370,0.0000117641,0.0000117776,0.0000117868,0.0000118770,0.0000119034,0.0000119250,0.0000118985,0.0000118437,0.0000119548,0.0000119401,0.0000119504,0.0000119150,0.0000120104,0.0000120197,0.0000119918,0.0000120780,0.0000120032,0.0000120577,0.0000120674,0.0000120497,0.0000121167,0.0000120823,0.0000121058,0.0000121777,0.0000121409,0.0000121290,0.0000122267,0.0000121846,0.0000121867,0.0000121846,0.0000122024,0.0000122004,0.0000122278,0.0000122166,0.0000122624,0.0000122430,0.0000121867,0.0000123546,0.0000123064,0.0000123027,0.0000083764,0.0000085110,0.0000086475,0.0000087670,0.0000089007,0.0000089683,0.0000090319,0.0000091483,0.0000092873,0.0000093556,0.0000093745,0.0000094658,0.0000095737,0.0000096796,0.0000097379,0.0000097107,0.0000098199,0.0000098635,0.0000099792,0.0000099696,0.0000099643,0.0000101023,0.0000100517,0.0000101430,0.0000102793,0.0000102569,0.0000103140,0.0000103401,0.0000104169,0.0000103703,0.0000104915,0.0000104644,0.0000105039,0.0000105939,0.0000105526,0.0000105905,0.0000106168,0.0000106192,0.0000106754,0.0000108030,0.0000107402,0.0000107845,0.0000107470,0.0000108658,0.0000108947,0.0000108488,0.0000109113,0.0000108834,0.0000109612,0.0000110337,0.0000109851,0.0000110326,0.0000110475,0.0000111029,0.0000110790,0.0000110979,0.0000111468,0.0000111676,0.0000111708,0.0000111447,0.0000111666,0.0000111887,0.0000111942,0.0000112893,0.0000113173,0.0000113201,0.0000112991,0.0000112414,0.0000113610,0.0000113393,0.0000113377,0.0000113264,0.0000114047,0.0000114283,0.0000113910,0.0000114711,0.0000114064,0.0000114574,0.0000114507,0.0000114405,0.0000114992,0.0000114753,0.0000115175,0.0000115750,0.0000115412,0.0000115201,0.0000116140,0.0000115731,0.0000115771,0.0000115685,0.0000115970,0.0000115958,0.0000116194,0.0000116133,0.0000116569,0.0000116350,0.0000115653,0.0000117399,0.0000116983,0.0000116903,0.0000079535,0.0000080750,0.0000082150,0.0000083170,0.0000084556,0.0000085139,0.0000085865,0.0000086891,0.0000088241,0.0000088854,0.0000089135,0.0000089838,0.0000090937,0.0000091928,0.0000092453,0.0000092344,0.0000093157,0.0000093559,0.0000094748,0.0000094678,0.0000094583,0.0000096049,0.0000095565,0.0000096323,0.0000097622,0.0000097340,0.0000097982,0.0000098275,0.0000098867,0.0000098562,0.0000099672,0.0000099397,0.0000099758,0.0000100622,0.0000100289,0.0000100606,0.0000100800,0.0000100811,0.0000101489,0.0000102655,0.0000101979,0.0000102305,0.0000102134,0.0000103256,0.0000103561,0.0000102976,0.0000103707,0.0000103443,0.0000104163,0.0000104938,0.0000104426,0.0000104874,0.0000104929,0.0000105410,0.0000105235,0.0000105402,0.0000105874,0.0000106124,0.0000106172,0.0000105821,0.0000106202,0.0000106306,0.0000106452,0.0000107368,0.0000107432,0.0000107673,0.0000107320,0.0000106754,0.0000107948,0.0000107863,0.0000107747,0.0000107633,0.0000108348,0.0000108605,0.0000108232,0.0000108978,0.0000108257,0.0000108917,0.0000108822,0.0000108698,0.0000109235,0.0000109130,0.0000109565,0.0000109860,0.0000109604,0.0000109371,0.0000110388,0.0000110074,0.0000110123,0.0000109990,0.0000110133,0.0000110102,0.0000110374,0.0000110301,0.0000110882,0.0000110447,0.0000110015,0.0000111399,0.0000111269,0.0000111137,0.0000075498,0.0000076753,0.0000078026,0.0000079024,0.0000080328,0.0000080822,0.0000081502,0.0000082560,0.0000083720,0.0000084438,0.0000084614,0.0000085223,0.0000086353,0.0000087354,0.0000087855,0.0000087597,0.0000088520,0.0000088804,0.0000090168,0.0000089881,0.0000089761,0.0000091155,0.0000090721,0.0000091444,0.0000092812,0.0000092482,0.0000093141,0.0000093326,0.0000093898,0.0000093617,0.0000094736,0.0000094456,0.0000094863,0.0000095601,0.0000095311,0.0000095560,0.0000095742,0.0000095789,0.0000096384,0.0000097586,0.0000096888,0.0000097231,0.0000097110,0.0000098154,0.0000098344,0.0000097794,0.0000098493,0.0000098316,0.0000099004,0.0000099630,0.0000099150,0.0000099587,0.0000099746,0.0000100156,0.0000099998,0.0000100150,0.0000100648,0.0000100899,0.0000100977,0.0000100614,0.0000100864,0.0000100951,0.0000101102,0.0000101888,0.0000102042,0.0000102381,0.0000101862,0.0000101496,0.0000102531,0.0000102480,0.0000102322,0.0000102242,0.0000102943,0.0000103162,0.0000102957,0.0000103587,0.0000102806,0.0000103430,0.0000103487,0.0000103298,0.0000103779,0.0000103661,0.0000104085,0.0000104429,0.0000104136,0.0000103927,0.0000104941,0.0000104626,0.0000104631,0.0000104459,0.0000104594,0.0000104641,0.0000104905,0.0000104874,0.0000105307,0.0000104958,0.0000104504,0.0000105770,0.0000105762,0.0000105662,0.0000071765,0.0000072954,0.0000074186,0.0000075055,0.0000076346,0.0000076770,0.0000077372,0.0000078542,0.0000079416,0.0000080219,0.0000080272,0.0000080853,0.0000082124,0.0000082938,0.0000083416,0.0000083242,0.0000083963,0.0000084434,0.0000085717,0.0000085340,0.0000085275,0.0000086572,0.0000086085,0.0000086835,0.0000088276,0.0000087818,0.0000088475,0.0000088604,0.0000089299,0.0000088944,0.0000090045,0.0000089567,0.0000090132,0.0000090983,0.0000090569,0.0000090897,0.0000090980,0.0000090894,0.0000091463,0.0000092741,0.0000092076,0.0000092364,0.0000092093,0.0000093253,0.0000093559,0.0000092935,0.0000093533,0.0000093333,0.0000094073,0.0000094625,0.0000094255,0.0000094693,0.0000094672,0.0000095094,0.0000094974,0.0000095178,0.0000095627,0.0000095852,0.0000095921,0.0000095719,0.0000095887,0.0000095877,0.0000096022,0.0000096754,0.0000096988,0.0000097209,0.0000096732,0.0000096494,0.0000097293,0.0000097438,0.0000097300,0.0000097182,0.0000097819,0.0000098052,0.0000097777,0.0000098364,0.0000097656,0.0000098280,0.0000098370,0.0000098081,0.0000098625,0.0000098480,0.0000099015,0.0000099062,0.0000099030,0.0000098735,0.0000099698,0.0000099422,0.0000099338,0.0000099199,0.0000099525,0.0000099515,0.0000099706,0.0000099668,0.0000100015,0.0000099733,0.0000099498,0.0000100527,0.0000100552,0.0000100419,0.0000068113,0.0000069203,0.0000070412,0.0000071400,0.0000072541,0.0000072946,0.0000073562,0.0000074625,0.0000075475,0.0000076244,0.0000076342,0.0000076766,0.0000077947,0.0000078815,0.0000079189,0.0000079152,0.0000079783,0.0000080220,0.0000081503,0.0000080990,0.0000081054,0.0000082237,0.0000081733,0.0000082500,0.0000083919,0.0000083349,0.0000084000,0.0000084231,0.0000084966,0.0000084498,0.0000085470,0.0000085049,0.0000085696,0.0000086397,0.0000086056,0.0000086380,0.0000086390,0.0000086276,0.0000086952,0.0000088069,0.0000087422,0.0000087781,0.0000087463,0.0000088528,0.0000088856,0.0000088298,0.0000088977,0.0000088720,0.0000089374,0.0000089882,0.0000089522,0.0000089995,0.0000089939,0.0000090397,0.0000090254,0.0000090363,0.0000090916,0.0000091154,0.0000091096,0.0000090963,0.0000091015,0.0000091056,0.0000091270,0.0000092095,0.0000092232,0.0000092282,0.0000091838,0.0000091819,0.0000092476,0.0000092553,0.0000092597,0.0000092413,0.0000092879,0.0000093168,0.0000092885,0.0000093463,0.0000092810,0.0000093422,0.0000093445,0.0000093150,0.0000093712,0.0000093593,0.0000094016,0.0000094136,0.0000094082,0.0000093724,0.0000094750,0.0000094361,0.0000094249,0.0000094318,0.0000094570,0.0000094572,0.0000094741,0.0000094686,0.0000095147,0.0000094701,0.0000094530,0.0000095699,0.0000095615,0.0000095359,0.0000064628,0.0000065685,0.0000066838,0.0000067914,0.0000068815,0.0000069252,0.0000069917,0.0000070874,0.0000071714,0.0000072391,0.0000072449,0.0000072975,0.0000074022,0.0000074763,0.0000075227,0.0000075154,0.0000075662,0.0000076264,0.0000077403,0.0000076849,0.0000076928,0.0000078116,0.0000077579,0.0000078395,0.0000079746,0.0000079165,0.0000079652,0.0000080094,0.0000080740,0.0000080177,0.0000081180,0.0000080682,0.0000081431,0.0000082160,0.0000081653,0.0000082085,0.0000082092,0.0000082045,0.0000082634,0.0000083532,0.0000083070,0.0000083461,0.0000083054,0.0000084080,0.0000084506,0.0000083901,0.0000084577,0.0000084257,0.0000084951,0.0000085545,0.0000084980,0.0000085549,0.0000085351,0.0000086005,0.0000085655,0.0000085843,0.0000086385,0.0000086623,0.0000086601,0.0000086462,0.0000086473,0.0000086526,0.0000086733,0.0000087547,0.0000087646,0.0000087681,0.0000087306,0.0000087235,0.0000087943,0.0000087815,0.0000088072,0.0000087711,0.0000088088,0.0000088495,0.0000088246,0.0000088846,0.0000088162,0.0000088827,0.0000088764,0.0000088585,0.0000089054,0.0000088937,0.0000089305,0.0000089418,0.0000089395,0.0000088936,0.0000089943,0.0000089707,0.0000089482,0.0000089570,0.0000089851,0.0000089884,0.0000090043,0.0000090081,0.0000090444,0.0000089927,0.0000089859,0.0000090971,0.0000090908,0.0000090621,0.0000061331,0.0000062433,0.0000063534,0.0000064505,0.0000065405,0.0000065762,0.0000066482,0.0000067220,0.0000067997,0.0000068751,0.0000068700,0.0000069335,0.0000070275,0.0000071143,0.0000071546,0.0000071462,0.0000072006,0.0000072377,0.0000073537,0.0000073062,0.0000072981,0.0000074194,0.0000073675,0.0000074509,0.0000075725,0.0000075211,0.0000075688,0.0000076008,0.0000076650,0.0000076152,0.0000077039,0.0000076636,0.0000077371,0.0000078121,0.0000077593,0.0000077949,0.0000077998,0.0000077871,0.0000078497,0.0000079165,0.0000078968,0.0000079428,0.0000078935,0.0000079902,0.0000080319,0.0000079672,0.0000080399,0.0000079968,0.0000080692,0.0000081388,0.0000080823,0.0000081335,0.0000081085,0.0000081647,0.0000081397,0.0000081603,0.0000082071,0.0000082337,0.0000082131,0.0000082185,0.0000082142,0.0000082161,0.0000082245,0.0000083328,0.0000083307,0.0000083269,0.0000082971,0.0000082978,0.0000083560,0.0000083374,0.0000083657,0.0000083302,0.0000083618,0.0000084071,0.0000083781,0.0000084362,0.0000083755,0.0000084271,0.0000084443,0.0000084122,0.0000084613,0.0000084564,0.0000084890,0.0000085049,0.0000084935,0.0000084499,0.0000085490,0.0000085255,0.0000085082,0.0000085172,0.0000085467,0.0000085287,0.0000085638,0.0000085553,0.0000085884,0.0000085473,0.0000085390,0.0000086447,0.0000086451,0.0000086001,0.0000058260,0.0000059216,0.0000060283,0.0000061277,0.0000062140,0.0000062464,0.0000063141,0.0000063902,0.0000064639,0.0000065299,0.0000065220,0.0000065996,0.0000066689,0.0000067574,0.0000067914,0.0000067814,0.0000068389,0.0000068717,0.0000069845,0.0000069282,0.0000069293,0.0000070480,0.0000070032,0.0000070761,0.0000072007,0.0000071510,0.0000071927,0.0000072127,0.0000072812,0.0000072375,0.0000073231,0.0000072775,0.0000073533,0.0000074207,0.0000073657,0.0000073953,0.0000074140,0.0000073883,0.0000074593,0.0000075319,0.0000074947,0.0000075455,0.0000074960,0.0000075915,0.0000076294,0.0000075693,0.0000076472,0.0000075968,0.0000076658,0.0000077319,0.0000076870,0.0000077238,0.0000077020,0.0000077512,0.0000077232,0.0000077623,0.0000077927,0.0000078317,0.0000078044,0.0000078070,0.0000078075,0.0000078029,0.0000078163,0.0000079223,0.0000079085,0.0000079120,0.0000078786,0.0000078958,0.0000079506,0.0000079227,0.0000079463,0.0000079245,0.0000079546,0.0000079798,0.0000079602,0.0000080050,0.0000079581,0.0000080016,0.0000080211,0.0000079916,0.0000080350,0.0000080247,0.0000080684,0.0000080862,0.0000080743,0.0000080253,0.0000081272,0.0000081017,0.0000080785,0.0000080838,0.0000081209,0.0000081052,0.0000081398,0.0000081246,0.0000081574,0.0000081178,0.0000081116,0.0000082011,0.0000082209,0.0000081714,0.0000055327,0.0000056292,0.0000057269,0.0000058233,0.0000059078,0.0000059337,0.0000059968,0.0000060684,0.0000061353,0.0000062101,0.0000062065,0.0000062726,0.0000063238,0.0000064252,0.0000064429,0.0000064418,0.0000064875,0.0000065257,0.0000066341,0.0000065732,0.0000065873,0.0000066787,0.0000066568,0.0000067304,0.0000068549,0.0000067998,0.0000068361,0.0000068540,0.0000069054,0.0000068737,0.0000069611,0.0000069162,0.0000069901,0.0000070474,0.0000069981,0.0000070446,0.0000070371,0.0000070070,0.0000070859,0.0000071480,0.0000071273,0.0000071629,0.0000071223,0.0000072095,0.0000072409,0.0000071966,0.0000072576,0.0000072276,0.0000072910,0.0000073517,0.0000072940,0.0000073263,0.0000073205,0.0000073669,0.0000073434,0.0000073726,0.0000074115,0.0000074365,0.0000074125,0.0000074093,0.0000074178,0.0000074063,0.0000074259,0.0000075240,0.0000075146,0.0000075147,0.0000074803,0.0000075092,0.0000075544,0.0000075239,0.0000075481,0.0000075287,0.0000075488,0.0000075854,0.0000075657,0.0000076083,0.0000075656,0.0000076023,0.0000076178,0.0000075892,0.0000076369,0.0000076289,0.0000076680,0.0000076820,0.0000076781,0.0000076166,0.0000077250,0.0000077004,0.0000076757,0.0000076806,0.0000077211,0.0000077031,0.0000077374,0.0000077227,0.0000077472,0.0000077231,0.0000076980,0.0000077922,0.0000078108,0.0000077625,0.0000052650,0.0000053445,0.0000054364,0.0000055366,0.0000056039,0.0000056272,0.0000056918,0.0000057615,0.0000058282,0.0000059042,0.0000058890,0.0000059567,0.0000060011,0.0000061067,0.0000061180,0.0000061140,0.0000061747,0.0000061974,0.0000063081,0.0000062462,0.0000062589,0.0000063423,0.0000063242,0.0000063974,0.0000065055,0.0000064458,0.0000064895,0.0000065176,0.0000065654,0.0000065347,0.0000066134,0.0000065707,0.0000066376,0.0000067022,0.0000066594,0.0000066827,0.0000066885,0.0000066554,0.0000067386,0.0000067863,0.0000067674,0.0000068003,0.0000067763,0.0000068433,0.0000068826,0.0000068511,0.0000068935,0.0000068702,0.0000069322,0.0000069853,0.0000069210,0.0000069568,0.0000069591,0.0000070061,0.0000069745,0.0000070025,0.0000070412,0.0000070687,0.0000070435,0.0000070335,0.0000070467,0.0000070248,0.0000070659,0.0000071518,0.0000071460,0.0000071463,0.0000071152,0.0000071287,0.0000071738,0.0000071455,0.0000071651,0.0000071508,0.0000071778,0.0000072124,0.0000071802,0.0000072351,0.0000071878,0.0000072180,0.0000072319,0.0000072075,0.0000072618,0.0000072552,0.0000072988,0.0000072998,0.0000072913,0.0000072403,0.0000073420,0.0000073202,0.0000073107,0.0000073056,0.0000073306,0.0000073184,0.0000073490,0.0000073437,0.0000073630,0.0000073399,0.0000073151,0.0000074139,0.0000074296,0.0000073764,0.0000050007,0.0000050785,0.0000051685,0.0000052541,0.0000053184,0.0000053371,0.0000054065,0.0000054719,0.0000055321,0.0000056071,0.0000055835,0.0000056571,0.0000057083,0.0000058008,0.0000058152,0.0000058111,0.0000058649,0.0000058866,0.0000059917,0.0000059269,0.0000059460,0.0000060342,0.0000060176,0.0000060868,0.0000061784,0.0000061285,0.0000061668,0.0000061897,0.0000062429,0.0000062020,0.0000062772,0.0000062435,0.0000063112,0.0000063763,0.0000063425,0.0000063524,0.0000063563,0.0000063145,0.0000064008,0.0000064418,0.0000064343,0.0000064634,0.0000064361,0.0000065068,0.0000065453,0.0000065139,0.0000065395,0.0000065257,0.0000065869,0.0000066335,0.0000065837,0.0000066000,0.0000066164,0.0000066564,0.0000066278,0.0000066560,0.0000066905,0.0000067076,0.0000066906,0.0000066788,0.0000066899,0.0000066823,0.0000067185,0.0000067984,0.0000067884,0.0000067910,0.0000067637,0.0000067820,0.0000068157,0.0000067992,0.0000068081,0.0000067933,0.0000068290,0.0000068548,0.0000068178,0.0000068787,0.0000068355,0.0000068654,0.0000068816,0.0000068616,0.0000068989,0.0000068924,0.0000069384,0.0000069354,0.0000069271,0.0000068885,0.0000069728,0.0000069644,0.0000069367,0.0000069374,0.0000069760,0.0000069579,0.0000069866,0.0000069751,0.0000069996,0.0000069711,0.0000069424,0.0000070491,0.0000070678,0.0000070031,0.0000047507,0.0000048199,0.0000049078,0.0000049946,0.0000050520,0.0000050701,0.0000051406,0.0000051967,0.0000052605,0.0000053258,0.0000053068,0.0000053768,0.0000054272,0.0000055000,0.0000055242,0.0000055178,0.0000055716,0.0000055884,0.0000056902,0.0000056345,0.0000056398,0.0000057251,0.0000057177,0.0000057839,0.0000058703,0.0000058207,0.0000058644,0.0000058823,0.0000059374,0.0000058936,0.0000059585,0.0000059270,0.0000059997,0.0000060641,0.0000060284,0.0000060339,0.0000060366,0.0000059991,0.0000060841,0.0000061208,0.0000061059,0.0000061367,0.0000061146,0.0000061764,0.0000062177,0.0000061841,0.0000062115,0.0000061940,0.0000062576,0.0000063048,0.0000062519,0.0000062751,0.0000062906,0.0000063164,0.0000062997,0.0000063234,0.0000063523,0.0000063630,0.0000063661,0.0000063402,0.0000063568,0.0000063527,0.0000063782,0.0000064565,0.0000064526,0.0000064523,0.0000064302,0.0000064482,0.0000064789,0.0000064625,0.0000064680,0.0000064613,0.0000064841,0.0000065163,0.0000064879,0.0000065307,0.0000065004,0.0000065291,0.0000065313,0.0000065247,0.0000065508,0.0000065533,0.0000065880,0.0000065829,0.0000065835,0.0000065415,0.0000066237,0.0000066136,0.0000065922,0.0000065922,0.0000066261,0.0000066101,0.0000066247,0.0000066244,0.0000066337,0.0000066224,0.0000065956,0.0000066974,0.0000067009,0.0000066556,0.0000045156,0.0000045704,0.0000046570,0.0000047541,0.0000047980,0.0000048229,0.0000048843,0.0000049444,0.0000049937,0.0000050559,0.0000050346,0.0000051094,0.0000051577,0.0000052253,0.0000052500,0.0000052335,0.0000052990,0.0000053167,0.0000054023,0.0000053447,0.0000053561,0.0000054334,0.0000054279,0.0000054999,0.0000055680,0.0000055338,0.0000055701,0.0000055880,0.0000056368,0.0000055920,0.0000056626,0.0000056303,0.0000057060,0.0000057587,0.0000057267,0.0000057387,0.0000057265,0.0000057003,0.0000057839,0.0000058218,0.0000058062,0.0000058280,0.0000058108,0.0000058706,0.0000059069,0.0000058705,0.0000059052,0.0000058924,0.0000059476,0.0000059912,0.0000059446,0.0000059597,0.0000059784,0.0000059974,0.0000059900,0.0000060088,0.0000060347,0.0000060514,0.0000060541,0.0000060311,0.0000060314,0.0000060234,0.0000060556,0.0000061366,0.0000061390,0.0000061369,0.0000061143,0.0000061261,0.0000061613,0.0000061490,0.0000061482,0.0000061281,0.0000061661,0.0000061924,0.0000061666,0.0000062008,0.0000061764,0.0000062046,0.0000062030,0.0000062067,0.0000062186,0.0000062345,0.0000062511,0.0000062513,0.0000062577,0.0000062180,0.0000062976,0.0000062833,0.0000062655,0.0000062666,0.0000062986,0.0000062812,0.0000062996,0.0000062964,0.0000063039,0.0000062809,0.0000062695,0.0000063651,0.0000063684,0.0000063294,0.0000042856,0.0000043433,0.0000044276,0.0000045191,0.0000045539,0.0000045818,0.0000046425,0.0000046908,0.0000047459,0.0000048034,0.0000047795,0.0000048576,0.0000049007,0.0000049616,0.0000049838,0.0000049710,0.0000050337,0.0000050536,0.0000051382,0.0000050806,0.0000050887,0.0000051598,0.0000051609,0.0000052206,0.0000052872,0.0000052559,0.0000052968,0.0000053098,0.0000053578,0.0000053135,0.0000053799,0.0000053503,0.0000054230,0.0000054726,0.0000054351,0.0000054551,0.0000054392,0.0000054131,0.0000054919,0.0000055336,0.0000055194,0.0000055389,0.0000055127,0.0000055744,0.0000056127,0.0000055824,0.0000056186,0.0000056082,0.0000056478,0.0000056870,0.0000056456,0.0000056599,0.0000056805,0.0000057014,0.0000056916,0.0000057125,0.0000057316,0.0000057582,0.0000057521,0.0000057349,0.0000057353,0.0000057329,0.0000057512,0.0000058361,0.0000058314,0.0000058329,0.0000058076,0.0000058189,0.0000058522,0.0000058457,0.0000058302,0.0000058218,0.0000058484,0.0000058866,0.0000058605,0.0000058873,0.0000058727,0.0000058972,0.0000058955,0.0000058952,0.0000059045,0.0000059240,0.0000059320,0.0000059405,0.0000059484,0.0000059079,0.0000059799,0.0000059694,0.0000059505,0.0000059504,0.0000059864,0.0000059770,0.0000059934,0.0000059910,0.0000059943,0.0000059660,0.0000059695,0.0000060468,0.0000060519,0.0000060022,0.0000040684,0.0000041279,0.0000042064,0.0000042929,0.0000043223,0.0000043485,0.0000044041,0.0000044441,0.0000045080,0.0000045624,0.0000045366,0.0000046126,0.0000046586,0.0000047088,0.0000047343,0.0000047235,0.0000047800,0.0000048050,0.0000048881,0.0000048297,0.0000048347,0.0000048976,0.0000049004,0.0000049564,0.0000050293,0.0000049898,0.0000050391,0.0000050421,0.0000050925,0.0000050575,0.0000051156,0.0000050798,0.0000051440,0.0000051935,0.0000051730,0.0000051855,0.0000051686,0.0000051328,0.0000052160,0.0000052494,0.0000052315,0.0000052547,0.0000052410,0.0000052963,0.0000053302,0.0000052947,0.0000053397,0.0000053262,0.0000053654,0.0000054053,0.0000053695,0.0000053793,0.0000054005,0.0000054280,0.0000054152,0.0000054262,0.0000054381,0.0000054812,0.0000054654,0.0000054540,0.0000054617,0.0000054525,0.0000054672,0.0000055456,0.0000055442,0.0000055410,0.0000055131,0.0000055296,0.0000055558,0.0000055482,0.0000055426,0.0000055340,0.0000055548,0.0000055937,0.0000055739,0.0000055876,0.0000055889,0.0000056039,0.0000056031,0.0000055955,0.0000056129,0.0000056308,0.0000056383,0.0000056517,0.0000056471,0.0000056247,0.0000056812,0.0000056790,0.0000056535,0.0000056378,0.0000056930,0.0000056810,0.0000056881,0.0000056961,0.0000056994,0.0000056785,0.0000056798,0.0000057374,0.0000057508,0.0000057165,0.0000038595,0.0000039082,0.0000039993,0.0000040736,0.0000041074,0.0000041239,0.0000041850,0.0000042181,0.0000042828,0.0000043376,0.0000043080,0.0000043868,0.0000044289,0.0000044693,0.0000045026,0.0000044853,0.0000045369,0.0000045644,0.0000046474,0.0000045785,0.0000045862,0.0000046564,0.0000046542,0.0000047060,0.0000047764,0.0000047280,0.0000047853,0.0000047894,0.0000048334,0.0000048038,0.0000048614,0.0000048244,0.0000048806,0.0000049322,0.0000049108,0.0000049324,0.0000049084,0.0000048774,0.0000049515,0.0000049896,0.0000049813,0.0000049889,0.0000049749,0.0000050328,0.0000050706,0.0000050262,0.0000050623,0.0000050602,0.0000050937,0.0000051339,0.0000050970,0.0000051145,0.0000051359,0.0000051535,0.0000051573,0.0000051569,0.0000051647,0.0000052114,0.0000052011,0.0000051746,0.0000051988,0.0000051848,0.0000051931,0.0000052751,0.0000052666,0.0000052626,0.0000052327,0.0000052496,0.0000052872,0.0000052788,0.0000052614,0.0000052561,0.0000052780,0.0000053194,0.0000053060,0.0000053108,0.0000053129,0.0000053261,0.0000053214,0.0000053231,0.0000053348,0.0000053479,0.0000053606,0.0000053676,0.0000053691,0.0000053448,0.0000054017,0.0000053980,0.0000053657,0.0000053520,0.0000054103,0.0000054036,0.0000054020,0.0000054205,0.0000054148,0.0000053963,0.0000053981,0.0000054472,0.0000054626,0.0000054363,0.0000036645,0.0000037087,0.0000038038,0.0000038707,0.0000039043,0.0000039125,0.0000039744,0.0000040112,0.0000040709,0.0000041223,0.0000040879,0.0000041613,0.0000042029,0.0000042439,0.0000042709,0.0000042624,0.0000043134,0.0000043397,0.0000044093,0.0000043387,0.0000043558,0.0000044214,0.0000044188,0.0000044719,0.0000045462,0.0000044831,0.0000045525,0.0000045551,0.0000045895,0.0000045597,0.0000046210,0.0000045822,0.0000046427,0.0000046771,0.0000046565,0.0000046880,0.0000046739,0.0000046277,0.0000047000,0.0000047393,0.0000047394,0.0000047376,0.0000047298,0.0000047811,0.0000048171,0.0000047718,0.0000048082,0.0000048156,0.0000048414,0.0000048773,0.0000048510,0.0000048602,0.0000048866,0.0000049015,0.0000048995,0.0000048997,0.0000049070,0.0000049522,0.0000049421,0.0000049180,0.0000049406,0.0000049235,0.0000049321,0.0000050161,0.0000050047,0.0000049906,0.0000049714,0.0000049881,0.0000050245,0.0000050205,0.0000049960,0.0000049886,0.0000050113,0.0000050612,0.0000050399,0.0000050489,0.0000050465,0.0000050625,0.0000050672,0.0000050569,0.0000050680,0.0000050809,0.0000050857,0.0000050983,0.0000051025,0.0000050787,0.0000051321,0.0000051286,0.0000050920,0.0000050828,0.0000051365,0.0000051329,0.0000051355,0.0000051472,0.0000051448,0.0000051212,0.0000051325,0.0000051757,0.0000051931,0.0000051812,0.0000034780,0.0000035141,0.0000036173,0.0000036792,0.0000037157,0.0000037179,0.0000037728,0.0000038052,0.0000038692,0.0000039141,0.0000038849,0.0000039479,0.0000039911,0.0000040248,0.0000040572,0.0000040422,0.0000041053,0.0000041298,0.0000041916,0.0000041220,0.0000041390,0.0000041953,0.0000041979,0.0000042546,0.0000043151,0.0000042614,0.0000043248,0.0000043359,0.0000043488,0.0000043342,0.0000043830,0.0000043579,0.0000044181,0.0000044462,0.0000044219,0.0000044623,0.0000044474,0.0000044071,0.0000044553,0.0000045001,0.0000045152,0.0000044984,0.0000044985,0.0000045475,0.0000045702,0.0000045300,0.0000045634,0.0000045792,0.0000045878,0.0000046336,0.0000046035,0.0000046153,0.0000046520,0.0000046619,0.0000046499,0.0000046552,0.0000046682,0.0000047054,0.0000046968,0.0000046766,0.0000046968,0.0000046796,0.0000046845,0.0000047593,0.0000047565,0.0000047489,0.0000047178,0.0000047458,0.0000047623,0.0000047765,0.0000047503,0.0000047333,0.0000047583,0.0000048063,0.0000047893,0.0000048006,0.0000048092,0.0000048053,0.0000048161,0.0000048047,0.0000048175,0.0000048250,0.0000048324,0.0000048538,0.0000048460,0.0000048323,0.0000048731,0.0000048726,0.0000048439,0.0000048199,0.0000048777,0.0000048916,0.0000048902,0.0000048883,0.0000048877,0.0000048726,0.0000048827,0.0000049141,0.0000049397,0.0000049273,0.0000033027,0.0000033338,0.0000034345,0.0000034928,0.0000035316,0.0000035346,0.0000035914,0.0000036072,0.0000036727,0.0000037146,0.0000036872,0.0000037435,0.0000037964,0.0000038189,0.0000038566,0.0000038395,0.0000039040,0.0000039221,0.0000039913,0.0000039146,0.0000039382,0.0000039894,0.0000039882,0.0000040421,0.0000040997,0.0000040386,0.0000040991,0.0000041210,0.0000041241,0.0000041169,0.0000041722,0.0000041341,0.0000041997,0.0000042258,0.0000042060,0.0000042355,0.0000042284,0.0000041846,0.0000042395,0.0000042771,0.0000042905,0.0000042766,0.0000042644,0.0000043276,0.0000043466,0.0000043073,0.0000043389,0.0000043501,0.0000043620,0.0000044096,0.0000043764,0.0000043868,0.0000044250,0.0000044231,0.0000044143,0.0000044237,0.0000044322,0.0000044805,0.0000044626,0.0000044467,0.0000044740,0.0000044525,0.0000044516,0.0000045256,0.0000045162,0.0000045198,0.0000044843,0.0000045165,0.0000045141,0.0000045354,0.0000045156,0.0000044938,0.0000045252,0.0000045665,0.0000045510,0.0000045589,0.0000045717,0.0000045689,0.0000045699,0.0000045739,0.0000045748,0.0000045865,0.0000045901,0.0000046121,0.0000046060,0.0000045949,0.0000046321,0.0000046302,0.0000046048,0.0000045851,0.0000046247,0.0000046439,0.0000046449,0.0000046418,0.0000046471,0.0000046271,0.0000046445,0.0000046692,0.0000046922,0.0000046868,0.0000031362,0.0000031670,0.0000032547,0.0000033200,0.0000033570,0.0000033584,0.0000034115,0.0000034246,0.0000034867,0.0000035293,0.0000035030,0.0000035526,0.0000036064,0.0000036260,0.0000036642,0.0000036362,0.0000037094,0.0000037239,0.0000037962,0.0000037170,0.0000037480,0.0000037950,0.0000037865,0.0000038378,0.0000038987,0.0000038372,0.0000038919,0.0000039173,0.0000039224,0.0000039169,0.0000039707,0.0000039229,0.0000039943,0.0000040185,0.0000040025,0.0000040215,0.0000040178,0.0000039751,0.0000040248,0.0000040646,0.0000040738,0.0000040719,0.0000040481,0.0000041145,0.0000041343,0.0000040969,0.0000041219,0.0000041223,0.0000041414,0.0000041939,0.0000041606,0.0000041577,0.0000042024,0.0000042077,0.0000041933,0.0000042053,0.0000042132,0.0000042597,0.0000042426,0.0000042253,0.0000042551,0.0000042239,0.0000042318,0.0000043100,0.0000042878,0.0000042914,0.0000042655,0.0000042909,0.0000042911,0.0000043027,0.0000042919,0.0000042697,0.0000042970,0.0000043364,0.0000043325,0.0000043297,0.0000043464,0.0000043530,0.0000043383,0.0000043438,0.0000043469,0.0000043547,0.0000043691,0.0000043853,0.0000043776,0.0000043645,0.0000043990,0.0000044088,0.0000043799,0.0000043508,0.0000043944,0.0000044099,0.0000044148,0.0000044161,0.0000044139,0.0000043974,0.0000044083,0.0000044357,0.0000044602,0.0000044586,0.0000029785,0.0000030066,0.0000030943,0.0000031533,0.0000031860,0.0000031930,0.0000032385,0.0000032537,0.0000033103,0.0000033494,0.0000033257,0.0000033758,0.0000034243,0.0000034364,0.0000034821,0.0000034550,0.0000035255,0.0000035393,0.0000036024,0.0000035297,0.0000035640,0.0000036075,0.0000035980,0.0000036406,0.0000036980,0.0000036418,0.0000036939,0.0000037260,0.0000037311,0.0000037295,0.0000037796,0.0000037246,0.0000037936,0.0000038103,0.0000038005,0.0000038218,0.0000038215,0.0000037765,0.0000038194,0.0000038595,0.0000038700,0.0000038676,0.0000038376,0.0000039084,0.0000039328,0.0000038868,0.0000039127,0.0000039072,0.0000039316,0.0000039909,0.0000039462,0.0000039448,0.0000039884,0.0000040061,0.0000039778,0.0000040031,0.0000040033,0.0000040520,0.0000040300,0.0000040155,0.0000040443,0.0000040166,0.0000040165,0.0000040972,0.0000040687,0.0000040794,0.0000040567,0.0000040821,0.0000040760,0.0000040928,0.0000040677,0.0000040532,0.0000040886,0.0000041138,0.0000041116,0.0000041121,0.0000041316,0.0000041361,0.0000041226,0.0000041196,0.0000041325,0.0000041355,0.0000041530,0.0000041737,0.0000041596,0.0000041540,0.0000041804,0.0000041938,0.0000041615,0.0000041395,0.0000041782,0.0000041838,0.0000042018,0.0000041944,0.0000041896,0.0000041753,0.0000041891,0.0000042187,0.0000042348,0.0000042329,0.0000028428,0.0000028573,0.0000029452,0.0000029955,0.0000030211,0.0000030302,0.0000030757,0.0000030921,0.0000031385,0.0000031800,0.0000031629,0.0000032097,0.0000032546,0.0000032661,0.0000033031,0.0000032826,0.0000033552,0.0000033656,0.0000034225,0.0000033521,0.0000033850,0.0000034189,0.0000034205,0.0000034706,0.0000035137,0.0000034557,0.0000035102,0.0000035433,0.0000035447,0.0000035444,0.0000035805,0.0000035391,0.0000036089,0.0000036125,0.0000036071,0.0000036273,0.0000036306,0.0000035887,0.0000036333,0.0000036677,0.0000036789,0.0000036725,0.0000036464,0.0000037149,0.0000037341,0.0000037006,0.0000037162,0.0000037110,0.0000037337,0.0000037879,0.0000037513,0.0000037578,0.0000037829,0.0000038049,0.0000037763,0.0000037956,0.0000038062,0.0000038532,0.0000038403,0.0000038167,0.0000038469,0.0000038179,0.0000038041,0.0000038919,0.0000038679,0.0000038871,0.0000038630,0.0000038812,0.0000038801,0.0000038905,0.0000038682,0.0000038518,0.0000038931,0.0000039073,0.0000039077,0.0000039153,0.0000039279,0.0000039324,0.0000039187,0.0000039120,0.0000039310,0.0000039278,0.0000039474,0.0000039672,0.0000039543,0.0000039444,0.0000039743,0.0000039816,0.0000039531,0.0000039410,0.0000039676,0.0000039714,0.0000039935,0.0000039867,0.0000039863,0.0000039717,0.0000039725,0.0000040106,0.0000040286,0.0000040276,0.0000027021,0.0000027115,0.0000027943,0.0000028462,0.0000028743,0.0000028789,0.0000029188,0.0000029358,0.0000029859,0.0000030208,0.0000030036,0.0000030492,0.0000030918,0.0000031054,0.0000031378,0.0000031239,0.0000031926,0.0000031942,0.0000032552,0.0000031826,0.0000032208,0.0000032434,0.0000032544,0.0000032991,0.0000033355,0.0000032809,0.0000033390,0.0000033653,0.0000033742,0.0000033685,0.0000034001,0.0000033593,0.0000034297,0.0000034342,0.0000034320,0.0000034439,0.0000034517,0.0000034039,0.0000034496,0.0000034815,0.0000034994,0.0000034918,0.0000034594,0.0000035344,0.0000035466,0.0000035162,0.0000035302,0.0000035287,0.0000035472,0.0000035936,0.0000035685,0.0000035687,0.0000035951,0.0000036095,0.0000035906,0.0000036035,0.0000036139,0.0000036571,0.0000036500,0.0000036234,0.0000036481,0.0000036208,0.0000036129,0.0000036917,0.0000036716,0.0000036973,0.0000036705,0.0000036885,0.0000036870,0.0000037001,0.0000036734,0.0000036678,0.0000036984,0.0000037088,0.0000037207,0.0000037176,0.0000037303,0.0000037373,0.0000037225,0.0000037193,0.0000037383,0.0000037302,0.0000037455,0.0000037687,0.0000037546,0.0000037479,0.0000037835,0.0000037752,0.0000037619,0.0000037393,0.0000037673,0.0000037749,0.0000037921,0.0000037960,0.0000037935,0.0000037753,0.0000037815,0.0000038169,0.0000038280,0.0000038268,0.0000025680,0.0000025748,0.0000026571,0.0000027044,0.0000027310,0.0000027418,0.0000027731,0.0000027886,0.0000028312,0.0000028702,0.0000028489,0.0000029041,0.0000029336,0.0000029464,0.0000029805,0.0000029749,0.0000030308,0.0000030318,0.0000030919,0.0000030319,0.0000030573,0.0000030765,0.0000030936,0.0000031376,0.0000031754,0.0000031181,0.0000031678,0.0000031947,0.0000032075,0.0000031954,0.0000032294,0.0000031837,0.0000032631,0.0000032644,0.0000032562,0.0000032766,0.0000032872,0.0000032390,0.0000032830,0.0000033037,0.0000033220,0.0000033227,0.0000032862,0.0000033571,0.0000033720,0.0000033399,0.0000033585,0.0000033576,0.0000033699,0.0000034167,0.0000033909,0.0000033969,0.0000034213,0.0000034248,0.0000034043,0.0000034333,0.0000034249,0.0000034739,0.0000034671,0.0000034390,0.0000034678,0.0000034423,0.0000034286,0.0000035052,0.0000034927,0.0000035145,0.0000034882,0.0000035066,0.0000035018,0.0000035119,0.0000034835,0.0000034835,0.0000035167,0.0000035247,0.0000035432,0.0000035288,0.0000035443,0.0000035507,0.0000035428,0.0000035345,0.0000035556,0.0000035425,0.0000035590,0.0000035722,0.0000035615,0.0000035621,0.0000035949,0.0000035904,0.0000035751,0.0000035562,0.0000035771,0.0000035873,0.0000036002,0.0000036035,0.0000036041,0.0000035924,0.0000035975,0.0000036258,0.0000036404,0.0000036363,0.0000024381,0.0000024458,0.0000025228,0.0000025707,0.0000026015,0.0000026046,0.0000026370,0.0000026477,0.0000026888,0.0000027261,0.0000027106,0.0000027507,0.0000027911,0.0000027975,0.0000028342,0.0000028284,0.0000028786,0.0000028816,0.0000029336,0.0000028782,0.0000029036,0.0000029233,0.0000029396,0.0000029834,0.0000030200,0.0000029658,0.0000030017,0.0000030418,0.0000030431,0.0000030347,0.0000030631,0.0000030245,0.0000030986,0.0000030937,0.0000030911,0.0000031048,0.0000031223,0.0000030754,0.0000031158,0.0000031418,0.0000031574,0.0000031648,0.0000031161,0.0000031945,0.0000032027,0.0000031720,0.0000031878,0.0000031902,0.0000031932,0.0000032427,0.0000032193,0.0000032245,0.0000032448,0.0000032536,0.0000032376,0.0000032582,0.0000032537,0.0000032943,0.0000032988,0.0000032722,0.0000032906,0.0000032693,0.0000032502,0.0000033313,0.0000033240,0.0000033417,0.0000033152,0.0000033317,0.0000033326,0.0000033423,0.0000033063,0.0000033081,0.0000033442,0.0000033535,0.0000033637,0.0000033517,0.0000033633,0.0000033675,0.0000033638,0.0000033532,0.0000033834,0.0000033641,0.0000033849,0.0000034027,0.0000033853,0.0000033897,0.0000034227,0.0000034083,0.0000033956,0.0000033806,0.0000033976,0.0000034086,0.0000034248,0.0000034303,0.0000034253,0.0000034205,0.0000034280,0.0000034428,0.0000034667,0.0000034543,0.0000023172,0.0000023246,0.0000024012,0.0000024392,0.0000024622,0.0000024760,0.0000025112,0.0000025121,0.0000025487,0.0000025839,0.0000025715,0.0000026096,0.0000026481,0.0000026553,0.0000026894,0.0000026827,0.0000027269,0.0000027412,0.0000027877,0.0000027313,0.0000027570,0.0000027787,0.0000027965,0.0000028353,0.0000028648,0.0000028192,0.0000028525,0.0000028929,0.0000028849,0.0000028846,0.0000029071,0.0000028738,0.0000029443,0.0000029398,0.0000029315,0.0000029609,0.0000029728,0.0000029294,0.0000029602,0.0000029837,0.0000029973,0.0000030122,0.0000029636,0.0000030281,0.0000030490,0.0000030159,0.0000030255,0.0000030345,0.0000030437,0.0000030861,0.0000030580,0.0000030573,0.0000030823,0.0000030904,0.0000030739,0.0000030953,0.0000030955,0.0000031272,0.0000031300,0.0000031076,0.0000031259,0.0000031025,0.0000030809,0.0000031680,0.0000031581,0.0000031707,0.0000031484,0.0000031664,0.0000031621,0.0000031758,0.0000031483,0.0000031469,0.0000031778,0.0000031869,0.0000032020,0.0000031856,0.0000031986,0.0000032022,0.0000031902,0.0000031845,0.0000032128,0.0000032007,0.0000032177,0.0000032358,0.0000032119,0.0000032179,0.0000032495,0.0000032411,0.0000032217,0.0000032092,0.0000032279,0.0000032392,0.0000032553,0.0000032613,0.0000032559,0.0000032525,0.0000032573,0.0000032724,0.0000032950,0.0000032859,0.0000022037,0.0000022120,0.0000022787,0.0000023175,0.0000023392,0.0000023486,0.0000023884,0.0000023789,0.0000024183,0.0000024503,0.0000024423,0.0000024805,0.0000025139,0.0000025190,0.0000025579,0.0000025485,0.0000025934,0.0000026026,0.0000026461,0.0000025922,0.0000026180,0.0000026395,0.0000026558,0.0000026941,0.0000027243,0.0000026802,0.0000027134,0.0000027444,0.0000027304,0.0000027416,0.0000027600,0.0000027280,0.0000027958,0.0000027913,0.0000027874,0.0000028215,0.0000028220,0.0000027811,0.0000028052,0.0000028308,0.0000028425,0.0000028596,0.0000028206,0.0000028743,0.0000028961,0.0000028700,0.0000028801,0.0000028832,0.0000028879,0.0000029356,0.0000028962,0.0000029038,0.0000029275,0.0000029320,0.0000029228,0.0000029429,0.0000029381,0.0000029769,0.0000029740,0.0000029584,0.0000029748,0.0000029556,0.0000029255,0.0000030122,0.0000030002,0.0000030143,0.0000029910,0.0000030079,0.0000030043,0.0000030258,0.0000029923,0.0000029936,0.0000030197,0.0000030295,0.0000030369,0.0000030255,0.0000030439,0.0000030433,0.0000030364,0.0000030297,0.0000030519,0.0000030404,0.0000030565,0.0000030727,0.0000030461,0.0000030589,0.0000030882,0.0000030746,0.0000030699,0.0000030533,0.0000030701,0.0000030749,0.0000030910,0.0000030898,0.0000030948,0.0000030933,0.0000030990,0.0000031109,0.0000031314,0.0000031209,0.0000020845,0.0000020999,0.0000021662,0.0000021970,0.0000022221,0.0000022311,0.0000022638,0.0000022593,0.0000022979,0.0000023310,0.0000023193,0.0000023557,0.0000023957,0.0000023892,0.0000024278,0.0000024281,0.0000024703,0.0000024753,0.0000025149,0.0000024582,0.0000024859,0.0000025075,0.0000025131,0.0000025592,0.0000025837,0.0000025429,0.0000025774,0.0000026086,0.0000025900,0.0000026031,0.0000026153,0.0000025948,0.0000026498,0.0000026504,0.0000026493,0.0000026801,0.0000026772,0.0000026475,0.0000026655,0.0000026902,0.0000027039,0.0000027156,0.0000026745,0.0000027358,0.0000027537,0.0000027234,0.0000027385,0.0000027420,0.0000027460,0.0000027925,0.0000027531,0.0000027556,0.0000027707,0.0000027810,0.0000027814,0.0000027968,0.0000027957,0.0000028264,0.0000028275,0.0000028167,0.0000028246,0.0000028034,0.0000027797,0.0000028646,0.0000028600,0.0000028614,0.0000028415,0.0000028624,0.0000028532,0.0000028705,0.0000028404,0.0000028501,0.0000028736,0.0000028720,0.0000028877,0.0000028791,0.0000028944,0.0000028864,0.0000028832,0.0000028736,0.0000029046,0.0000028936,0.0000029091,0.0000029195,0.0000029003,0.0000029025,0.0000029306,0.0000029264,0.0000029164,0.0000029012,0.0000029127,0.0000029254,0.0000029328,0.0000029365,0.0000029396,0.0000029434,0.0000029464,0.0000029616,0.0000029696,0.0000029626,0.0000019817,0.0000019903,0.0000020597,0.0000020803,0.0000021160,0.0000021225,0.0000021485,0.0000021455,0.0000021855,0.0000022120,0.0000022076,0.0000022420,0.0000022748,0.0000022674,0.0000023036,0.0000023106,0.0000023497,0.0000023518,0.0000023881,0.0000023362,0.0000023594,0.0000023800,0.0000023906,0.0000024312,0.0000024483,0.0000024120,0.0000024443,0.0000024755,0.0000024621,0.0000024671,0.0000024840,0.0000024679,0.0000025171,0.0000025194,0.0000025205,0.0000025443,0.0000025444,0.0000025132,0.0000025305,0.0000025595,0.0000025709,0.0000025765,0.0000025415,0.0000025974,0.0000026200,0.0000025866,0.0000026046,0.0000026064,0.0000026055,0.0000026539,0.0000026169,0.0000026248,0.0000026362,0.0000026456,0.0000026440,0.0000026542,0.0000026586,0.0000026843,0.0000026933,0.0000026768,0.0000026852,0.0000026610,0.0000026403,0.0000027218,0.0000027171,0.0000027125,0.0000026971,0.0000027282,0.0000027066,0.0000027317,0.0000027061,0.0000027108,0.0000027301,0.0000027310,0.0000027401,0.0000027439,0.0000027524,0.0000027444,0.0000027397,0.0000027342,0.0000027546,0.0000027499,0.0000027611,0.0000027801,0.0000027535,0.0000027533,0.0000027824,0.0000027800,0.0000027733,0.0000027566,0.0000027567,0.0000027817,0.0000027858,0.0000027886,0.0000027927,0.0000027986,0.0000028013,0.0000028143,0.0000028246,0.0000028165,0.0000018784,0.0000018915,0.0000019561,0.0000019714,0.0000020051,0.0000020141,0.0000020394,0.0000020374,0.0000020803,0.0000020980,0.0000020979,0.0000021320,0.0000021626,0.0000021580,0.0000021781,0.0000021962,0.0000022306,0.0000022366,0.0000022627,0.0000022182,0.0000022481,0.0000022599,0.0000022739,0.0000023106,0.0000023271,0.0000022958,0.0000023176,0.0000023495,0.0000023403,0.0000023451,0.0000023616,0.0000023496,0.0000023845,0.0000023911,0.0000023981,0.0000024067,0.0000024224,0.0000023894,0.0000024048,0.0000024299,0.0000024414,0.0000024473,0.0000024200,0.0000024665,0.0000024843,0.0000024575,0.0000024761,0.0000024740,0.0000024772,0.0000025210,0.0000024905,0.0000024958,0.0000025074,0.0000025113,0.0000025128,0.0000025217,0.0000025259,0.0000025516,0.0000025613,0.0000025404,0.0000025534,0.0000025295,0.0000025112,0.0000025874,0.0000025779,0.0000025758,0.0000025622,0.0000025949,0.0000025719,0.0000025954,0.0000025745,0.0000025777,0.0000025956,0.0000025939,0.0000026058,0.0000026063,0.0000026134,0.0000026121,0.0000026067,0.0000025996,0.0000026178,0.0000026101,0.0000026252,0.0000026419,0.0000026138,0.0000026165,0.0000026469,0.0000026357,0.0000026369,0.0000026191,0.0000026170,0.0000026454,0.0000026459,0.0000026531,0.0000026543,0.0000026599,0.0000026584,0.0000026751,0.0000026797,0.0000026715,0.0000017821,0.0000017950,0.0000018568,0.0000018747,0.0000019119,0.0000019128,0.0000019420,0.0000019392,0.0000019776,0.0000019856,0.0000019935,0.0000020248,0.0000020522,0.0000020485,0.0000020650,0.0000020911,0.0000021184,0.0000021262,0.0000021424,0.0000021094,0.0000021375,0.0000021493,0.0000021584,0.0000021914,0.0000022158,0.0000021790,0.0000022007,0.0000022382,0.0000022236,0.0000022297,0.0000022452,0.0000022310,0.0000022664,0.0000022695,0.0000022747,0.0000022909,0.0000022978,0.0000022689,0.0000022832,0.0000023056,0.0000023246,0.0000023250,0.0000022999,0.0000023483,0.0000023639,0.0000023418,0.0000023573,0.0000023546,0.0000023588,0.0000023925,0.0000023658,0.0000023711,0.0000023819,0.0000023862,0.0000023878,0.0000023967,0.0000024043,0.0000024259,0.0000024340,0.0000024163,0.0000024289,0.0000024016,0.0000023833,0.0000024627,0.0000024510,0.0000024483,0.0000024344,0.0000024666,0.0000024422,0.0000024641,0.0000024412,0.0000024518,0.0000024683,0.0000024655,0.0000024750,0.0000024774,0.0000024879,0.0000024800,0.0000024748,0.0000024682,0.0000024826,0.0000024813,0.0000024950,0.0000025124,0.0000024842,0.0000024794,0.0000025215,0.0000025043,0.0000025095,0.0000024813,0.0000024941,0.0000025094,0.0000025100,0.0000025212,0.0000025242,0.0000025289,0.0000025183,0.0000025427,0.0000025487,0.0000025392,0.0000016841,0.0000017101,0.0000017678,0.0000017843,0.0000018161,0.0000018171,0.0000018370,0.0000018462,0.0000018807,0.0000018864,0.0000019004,0.0000019247,0.0000019512,0.0000019431,0.0000019592,0.0000019826,0.0000020139,0.0000020242,0.0000020355,0.0000020097,0.0000020334,0.0000020378,0.0000020501,0.0000020832,0.0000021030,0.0000020655,0.0000020913,0.0000021252,0.0000021043,0.0000021174,0.0000021335,0.0000021170,0.0000021545,0.0000021598,0.0000021649,0.0000021790,0.0000021860,0.0000021529,0.0000021682,0.0000021879,0.0000022042,0.0000022094,0.0000021831,0.0000022330,0.0000022456,0.0000022246,0.0000022373,0.0000022465,0.0000022412,0.0000022727,0.0000022505,0.0000022507,0.0000022654,0.0000022701,0.0000022619,0.0000022817,0.0000022921,0.0000023019,0.0000023069,0.0000022926,0.0000023095,0.0000022767,0.0000022630,0.0000023395,0.0000023336,0.0000023251,0.0000023102,0.0000023388,0.0000023145,0.0000023356,0.0000023170,0.0000023271,0.0000023465,0.0000023450,0.0000023564,0.0000023530,0.0000023588,0.0000023567,0.0000023494,0.0000023495,0.0000023593,0.0000023640,0.0000023735,0.0000023927,0.0000023623,0.0000023646,0.0000023906,0.0000023730,0.0000023820,0.0000023583,0.0000023737,0.0000023790,0.0000023880,0.0000023994,0.0000024083,0.0000024014,0.0000023926,0.0000024162,0.0000024267,0.0000024212,0.0000016014,0.0000016201,0.0000016844,0.0000016988,0.0000017244,0.0000017275,0.0000017409,0.0000017517,0.0000017898,0.0000017864,0.0000018037,0.0000018286,0.0000018539,0.0000018434,0.0000018583,0.0000018868,0.0000019142,0.0000019298,0.0000019373,0.0000019089,0.0000019277,0.0000019346,0.0000019550,0.0000019760,0.0000019936,0.0000019625,0.0000019851,0.0000020192,0.0000019973,0.0000020129,0.0000020286,0.0000020086,0.0000020495,0.0000020543,0.0000020533,0.0000020728,0.0000020764,0.0000020496,0.0000020558,0.0000020811,0.0000021002,0.0000020966,0.0000020742,0.0000021242,0.0000021305,0.0000021178,0.0000021304,0.0000021348,0.0000021279,0.0000021594,0.0000021395,0.0000021389,0.0000021486,0.0000021560,0.0000021467,0.0000021660,0.0000021780,0.0000021858,0.0000021937,0.0000021776,0.0000021948,0.0000021633,0.0000021540,0.0000022298,0.0000022167,0.0000022059,0.0000021963,0.0000022165,0.0000022033,0.0000022121,0.0000022026,0.0000022118,0.0000022263,0.0000022303,0.0000022398,0.0000022376,0.0000022413,0.0000022402,0.0000022278,0.0000022299,0.0000022394,0.0000022421,0.0000022565,0.0000022731,0.0000022485,0.0000022461,0.0000022713,0.0000022552,0.0000022620,0.0000022387,0.0000022572,0.0000022624,0.0000022703,0.0000022856,0.0000022883,0.0000022788,0.0000022702,0.0000022991,0.0000023041,0.0000022954,0.0000015211,0.0000015341,0.0000015981,0.0000016109,0.0000016352,0.0000016452,0.0000016516,0.0000016614,0.0000017021,0.0000016974,0.0000017113,0.0000017400,0.0000017577,0.0000017519,0.0000017651,0.0000017925,0.0000018149,0.0000018323,0.0000018409,0.0000018179,0.0000018324,0.0000018390,0.0000018522,0.0000018828,0.0000018930,0.0000018679,0.0000018862,0.0000019232,0.0000019005,0.0000019106,0.0000019275,0.0000019088,0.0000019515,0.0000019522,0.0000019454,0.0000019754,0.0000019794,0.0000019491,0.0000019493,0.0000019770,0.0000019968,0.0000019915,0.0000019714,0.0000020154,0.0000020280,0.0000020143,0.0000020294,0.0000020277,0.0000020263,0.0000020572,0.0000020301,0.0000020365,0.0000020371,0.0000020530,0.0000020354,0.0000020596,0.0000020729,0.0000020750,0.0000020802,0.0000020624,0.0000020872,0.0000020551,0.0000020452,0.0000021184,0.0000021054,0.0000021031,0.0000020853,0.0000021063,0.0000020998,0.0000021018,0.0000020904,0.0000021012,0.0000021220,0.0000021177,0.0000021247,0.0000021286,0.0000021302,0.0000021224,0.0000021160,0.0000021211,0.0000021191,0.0000021337,0.0000021442,0.0000021651,0.0000021341,0.0000021373,0.0000021550,0.0000021463,0.0000021459,0.0000021265,0.0000021422,0.0000021491,0.0000021580,0.0000021677,0.0000021734,0.0000021594,0.0000021563,0.0000021944,0.0000021911,0.0000021886,0.0000014443,0.0000014616,0.0000015176,0.0000015270,0.0000015526,0.0000015640,0.0000015664,0.0000015785,0.0000016200,0.0000016118,0.0000016242,0.0000016547,0.0000016698,0.0000016683,0.0000016802,0.0000017063,0.0000017275,0.0000017447,0.0000017492,0.0000017299,0.0000017437,0.0000017395,0.0000017561,0.0000017831,0.0000017977,0.0000017771,0.0000017928,0.0000018311,0.0000018063,0.0000018078,0.0000018304,0.0000018149,0.0000018541,0.0000018560,0.0000018488,0.0000018781,0.0000018806,0.0000018522,0.0000018540,0.0000018753,0.0000018980,0.0000018841,0.0000018766,0.0000019160,0.0000019245,0.0000019159,0.0000019323,0.0000019261,0.0000019211,0.0000019507,0.0000019279,0.0000019358,0.0000019332,0.0000019481,0.0000019351,0.0000019657,0.0000019730,0.0000019689,0.0000019845,0.0000019597,0.0000019847,0.0000019519,0.0000019485,0.0000020135,0.0000019933,0.0000020014,0.0000019777,0.0000020035,0.0000019971,0.0000019949,0.0000019896,0.0000019949,0.0000020197,0.0000020078,0.0000020231,0.0000020234,0.0000020198,0.0000020215,0.0000020092,0.0000020171,0.0000020133,0.0000020301,0.0000020348,0.0000020568,0.0000020333,0.0000020308,0.0000020471,0.0000020359,0.0000020379,0.0000020223,0.0000020408,0.0000020493,0.0000020548,0.0000020639,0.0000020637,0.0000020531,0.0000020505,0.0000020802,0.0000020843,0.0000020800,0.0000013706,0.0000013905,0.0000014389,0.0000014468,0.0000014757,0.0000014875,0.0000014814,0.0000014991,0.0000015403,0.0000015300,0.0000015391,0.0000015702,0.0000015915,0.0000015817,0.0000015957,0.0000016130,0.0000016393,0.0000016568,0.0000016582,0.0000016486,0.0000016492,0.0000016521,0.0000016710,0.0000016943,0.0000017050,0.0000016853,0.0000017042,0.0000017392,0.0000017147,0.0000017147,0.0000017389,0.0000017243,0.0000017645,0.0000017608,0.0000017608,0.0000017853,0.0000017899,0.0000017592,0.0000017598,0.0000017805,0.0000018012,0.0000017943,0.0000017844,0.0000018225,0.0000018217,0.0000018204,0.0000018376,0.0000018310,0.0000018227,0.0000018606,0.0000018275,0.0000018445,0.0000018383,0.0000018568,0.0000018393,0.0000018694,0.0000018744,0.0000018715,0.0000018837,0.0000018610,0.0000018851,0.0000018478,0.0000018494,0.0000019154,0.0000018968,0.0000019072,0.0000018798,0.0000019066,0.0000018983,0.0000018933,0.0000018936,0.0000018918,0.0000019176,0.0000019096,0.0000019235,0.0000019268,0.0000019185,0.0000019186,0.0000019065,0.0000019162,0.0000019120,0.0000019298,0.0000019357,0.0000019544,0.0000019354,0.0000019264,0.0000019448,0.0000019356,0.0000019388,0.0000019205,0.0000019393,0.0000019472,0.0000019510,0.0000019582,0.0000019613,0.0000019529,0.0000019461,0.0000019758,0.0000019788,0.0000019718,0.0000013064,0.0000013226,0.0000013639,0.0000013746,0.0000014027,0.0000014137,0.0000014022,0.0000014225,0.0000014628,0.0000014501,0.0000014667,0.0000014907,0.0000015056,0.0000015073,0.0000015140,0.0000015324,0.0000015590,0.0000015736,0.0000015773,0.0000015674,0.0000015713,0.0000015687,0.0000015878,0.0000016061,0.0000016171,0.0000015995,0.0000016125,0.0000016558,0.0000016354,0.0000016307,0.0000016523,0.0000016397,0.0000016776,0.0000016692,0.0000016738,0.0000016967,0.0000016997,0.0000016743,0.0000016701,0.0000016898,0.0000017114,0.0000017062,0.0000016968,0.0000017320,0.0000017327,0.0000017301,0.0000017407,0.0000017420,0.0000017281,0.0000017667,0.0000017374,0.0000017528,0.0000017494,0.0000017603,0.0000017496,0.0000017771,0.0000017850,0.0000017790,0.0000017895,0.0000017684,0.0000017908,0.0000017510,0.0000017547,0.0000018266,0.0000018011,0.0000018137,0.0000017857,0.0000018135,0.0000018007,0.0000017997,0.0000018022,0.0000018015,0.0000018227,0.0000018159,0.0000018280,0.0000018331,0.0000018228,0.0000018166,0.0000018123,0.0000018238,0.0000018214,0.0000018362,0.0000018420,0.0000018565,0.0000018357,0.0000018272,0.0000018457,0.0000018419,0.0000018412,0.0000018276,0.0000018493,0.0000018498,0.0000018568,0.0000018635,0.0000018632,0.0000018541,0.0000018458,0.0000018781,0.0000018855,0.0000018762,0.0000012371,0.0000012595,0.0000012969,0.0000013060,0.0000013310,0.0000013436,0.0000013282,0.0000013532,0.0000013891,0.0000013751,0.0000013931,0.0000014193,0.0000014312,0.0000014326,0.0000014368,0.0000014562,0.0000014837,0.0000014971,0.0000014939,0.0000014872,0.0000014918,0.0000014910,0.0000015114,0.0000015309,0.0000015349,0.0000015188,0.0000015315,0.0000015753,0.0000015532,0.0000015501,0.0000015684,0.0000015578,0.0000015962,0.0000015841,0.0000015908,0.0000016147,0.0000016100,0.0000015892,0.0000015856,0.0000016032,0.0000016243,0.0000016195,0.0000016129,0.0000016463,0.0000016459,0.0000016422,0.0000016506,0.0000016549,0.0000016428,0.0000016762,0.0000016483,0.0000016705,0.0000016626,0.0000016733,0.0000016576,0.0000016910,0.0000016954,0.0000016874,0.0000017006,0.0000016824,0.0000017029,0.0000016705,0.0000016662,0.0000017335,0.0000017072,0.0000017225,0.0000016998,0.0000017258,0.0000017139,0.0000017142,0.0000017131,0.0000017078,0.0000017284,0.0000017311,0.0000017365,0.0000017414,0.0000017352,0.0000017250,0.0000017276,0.0000017290,0.0000017351,0.0000017423,0.0000017525,0.0000017619,0.0000017442,0.0000017402,0.0000017537,0.0000017506,0.0000017515,0.0000017361,0.0000017533,0.0000017575,0.0000017640,0.0000017619,0.0000017676,0.0000017622,0.0000017555,0.0000017872,0.0000017953,0.0000017873,0.0000011699,0.0000011958,0.0000012334,0.0000012406,0.0000012641,0.0000012759,0.0000012597,0.0000012854,0.0000013185,0.0000013077,0.0000013244,0.0000013499,0.0000013625,0.0000013594,0.0000013695,0.0000013804,0.0000014134,0.0000014287,0.0000014139,0.0000014123,0.0000014126,0.0000014173,0.0000014406,0.0000014557,0.0000014593,0.0000014410,0.0000014508,0.0000014999,0.0000014752,0.0000014760,0.0000014932,0.0000014795,0.0000015156,0.0000015034,0.0000015139,0.0000015358,0.0000015335,0.0000015078,0.0000015045,0.0000015231,0.0000015508,0.0000015377,0.0000015306,0.0000015636,0.0000015651,0.0000015601,0.0000015716,0.0000015709,0.0000015641,0.0000015924,0.0000015666,0.0000015868,0.0000015802,0.0000015861,0.0000015833,0.0000016089,0.0000016102,0.0000016073,0.0000016185,0.0000015947,0.0000016163,0.0000015889,0.0000015841,0.0000016485,0.0000016218,0.0000016396,0.0000016144,0.0000016447,0.0000016256,0.0000016268,0.0000016287,0.0000016225,0.0000016468,0.0000016417,0.0000016490,0.0000016585,0.0000016462,0.0000016371,0.0000016472,0.0000016427,0.0000016483,0.0000016548,0.0000016657,0.0000016760,0.0000016567,0.0000016540,0.0000016720,0.0000016666,0.0000016659,0.0000016478,0.0000016683,0.0000016725,0.0000016762,0.0000016802,0.0000016799,0.0000016742,0.0000016645,0.0000016973,0.0000017063,0.0000016974,0.0000011057,0.0000011297,0.0000011697,0.0000011796,0.0000011982,0.0000012114,0.0000011965,0.0000012193,0.0000012464,0.0000012400,0.0000012596,0.0000012862,0.0000012928,0.0000012895,0.0000012981,0.0000013127,0.0000013405,0.0000013575,0.0000013460,0.0000013439,0.0000013457,0.0000013487,0.0000013719,0.0000013880,0.0000013920,0.0000013661,0.0000013759,0.0000014239,0.0000014026,0.0000014025,0.0000014191,0.0000014082,0.0000014353,0.0000014254,0.0000014344,0.0000014596,0.0000014540,0.0000014321,0.0000014310,0.0000014456,0.0000014753,0.0000014605,0.0000014535,0.0000014879,0.0000014889,0.0000014802,0.0000014942,0.0000014930,0.0000014911,0.0000015152,0.0000014871,0.0000015110,0.0000014990,0.0000015054,0.0000015028,0.0000015286,0.0000015295,0.0000015289,0.0000015353,0.0000015169,0.0000015327,0.0000015037,0.0000015023,0.0000015691,0.0000015373,0.0000015609,0.0000015321,0.0000015621,0.0000015472,0.0000015436,0.0000015433,0.0000015494,0.0000015661,0.0000015605,0.0000015668,0.0000015762,0.0000015647,0.0000015515,0.0000015598,0.0000015569,0.0000015641,0.0000015767,0.0000015822,0.0000015952,0.0000015771,0.0000015752,0.0000015914,0.0000015863,0.0000015795,0.0000015635,0.0000015864,0.0000015900,0.0000015884,0.0000015968,0.0000015938,0.0000015944,0.0000015815,0.0000016178,0.0000016253,0.0000016125,0.0000010473,0.0000010721,0.0000011108,0.0000011219,0.0000011373,0.0000011525,0.0000011369,0.0000011591,0.0000011845,0.0000011854,0.0000011947,0.0000012231,0.0000012255,0.0000012248,0.0000012330,0.0000012477,0.0000012766,0.0000012870,0.0000012786,0.0000012797,0.0000012790,0.0000012805,0.0000013055,0.0000013193,0.0000013166,0.0000012996,0.0000013090,0.0000013516,0.0000013324,0.0000013374,0.0000013507,0.0000013360,0.0000013654,0.0000013540,0.0000013668,0.0000013837,0.0000013803,0.0000013614,0.0000013633,0.0000013714,0.0000014041,0.0000013867,0.0000013769,0.0000014161,0.0000014137,0.0000014101,0.0000014201,0.0000014168,0.0000014214,0.0000014372,0.0000014138,0.0000014354,0.0000014207,0.0000014313,0.0000014287,0.0000014589,0.0000014564,0.0000014511,0.0000014651,0.0000014446,0.0000014577,0.0000014292,0.0000014221,0.0000014909,0.0000014620,0.0000014832,0.0000014523,0.0000014880,0.0000014739,0.0000014722,0.0000014599,0.0000014756,0.0000014963,0.0000014792,0.0000014918,0.0000014925,0.0000014867,0.0000014712,0.0000014809,0.0000014784,0.0000014847,0.0000014976,0.0000015056,0.0000015153,0.0000014975,0.0000014955,0.0000015114,0.0000015082,0.0000015035,0.0000014875,0.0000015075,0.0000015083,0.0000015091,0.0000015193,0.0000015128,0.0000015132,0.0000015036,0.0000015366,0.0000015433,0.0000015317,0.0000009954,0.0000010162,0.0000010580,0.0000010664,0.0000010806,0.0000010932,0.0000010829,0.0000010991,0.0000011247,0.0000011257,0.0000011358,0.0000011607,0.0000011633,0.0000011654,0.0000011711,0.0000011894,0.0000012107,0.0000012224,0.0000012115,0.0000012202,0.0000012134,0.0000012170,0.0000012421,0.0000012540,0.0000012474,0.0000012344,0.0000012441,0.0000012845,0.0000012674,0.0000012660,0.0000012861,0.0000012682,0.0000012982,0.0000012864,0.0000013037,0.0000013132,0.0000013093,0.0000012931,0.0000012951,0.0000012994,0.0000013365,0.0000013185,0.0000013090,0.0000013379,0.0000013452,0.0000013428,0.0000013489,0.0000013480,0.0000013505,0.0000013696,0.0000013457,0.0000013634,0.0000013518,0.0000013596,0.0000013525,0.0000013840,0.0000013842,0.0000013812,0.0000013896,0.0000013766,0.0000013863,0.0000013597,0.0000013512,0.0000014159,0.0000013876,0.0000014127,0.0000013755,0.0000014112,0.0000014013,0.0000013943,0.0000013887,0.0000014031,0.0000014238,0.0000014047,0.0000014162,0.0000014143,0.0000014125,0.0000013977,0.0000014086,0.0000014064,0.0000014094,0.0000014240,0.0000014277,0.0000014335,0.0000014223,0.0000014235,0.0000014365,0.0000014365,0.0000014321,0.0000014110,0.0000014303,0.0000014337,0.0000014340,0.0000014442,0.0000014412,0.0000014363,0.0000014223,0.0000014594,0.0000014652,0.0000014562,0.0000009461,0.0000009669,0.0000010038,0.0000010130,0.0000010282,0.0000010377,0.0000010307,0.0000010463,0.0000010702,0.0000010682,0.0000010775,0.0000011015,0.0000011084,0.0000011097,0.0000011113,0.0000011260,0.0000011460,0.0000011619,0.0000011505,0.0000011577,0.0000011531,0.0000011632,0.0000011790,0.0000011879,0.0000011839,0.0000011737,0.0000011814,0.0000012201,0.0000012026,0.0000012010,0.0000012198,0.0000012059,0.0000012321,0.0000012252,0.0000012365,0.0000012483,0.0000012424,0.0000012268,0.0000012342,0.0000012337,0.0000012694,0.0000012548,0.0000012402,0.0000012715,0.0000012824,0.0000012752,0.0000012797,0.0000012833,0.0000012814,0.0000013050,0.0000012795,0.0000012956,0.0000012869,0.0000012946,0.0000012848,0.0000013181,0.0000013166,0.0000013122,0.0000013253,0.0000013049,0.0000013165,0.0000012938,0.0000012859,0.0000013452,0.0000013157,0.0000013414,0.0000013073,0.0000013425,0.0000013313,0.0000013282,0.0000013227,0.0000013300,0.0000013534,0.0000013384,0.0000013514,0.0000013455,0.0000013416,0.0000013324,0.0000013379,0.0000013321,0.0000013395,0.0000013562,0.0000013581,0.0000013605,0.0000013497,0.0000013528,0.0000013664,0.0000013670,0.0000013623,0.0000013377,0.0000013616,0.0000013627,0.0000013630,0.0000013727,0.0000013708,0.0000013649,0.0000013554,0.0000013856,0.0000013932,0.0000013854,0.0000008986,0.0000009167,0.0000009519,0.0000009626,0.0000009780,0.0000009867,0.0000009802,0.0000009923,0.0000010171,0.0000010156,0.0000010229,0.0000010428,0.0000010524,0.0000010572,0.0000010550,0.0000010670,0.0000010871,0.0000011058,0.0000010905,0.0000011001,0.0000010933,0.0000011084,0.0000011196,0.0000011281,0.0000011254,0.0000011168,0.0000011167,0.0000011625,0.0000011457,0.0000011429,0.0000011640,0.0000011478,0.0000011678,0.0000011657,0.0000011728,0.0000011881,0.0000011790,0.0000011655,0.0000011718,0.0000011715,0.0000012060,0.0000011940,0.0000011745,0.0000012089,0.0000012192,0.0000012148,0.0000012139,0.0000012150,0.0000012169,0.0000012305,0.0000012133,0.0000012319,0.0000012268,0.0000012292,0.0000012228,0.0000012490,0.0000012535,0.0000012495,0.0000012601,0.0000012429,0.0000012519,0.0000012262,0.0000012205,0.0000012781,0.0000012527,0.0000012722,0.0000012445,0.0000012752,0.0000012606,0.0000012583,0.0000012565,0.0000012650,0.0000012869,0.0000012722,0.0000012823,0.0000012816,0.0000012766,0.0000012632,0.0000012679,0.0000012682,0.0000012763,0.0000012902,0.0000012865,0.0000012943,0.0000012875,0.0000012887,0.0000013001,0.0000012962,0.0000012955,0.0000012672,0.0000012987,0.0000012940,0.0000012956,0.0000013040,0.0000013053,0.0000012923,0.0000012881,0.0000013198,0.0000013247,0.0000013164,0.0000008541,0.0000008705,0.0000009051,0.0000009138,0.0000009328,0.0000009388,0.0000009303,0.0000009407,0.0000009663,0.0000009665,0.0000009693,0.0000009922,0.0000010015,0.0000010016,0.0000010035,0.0000010158,0.0000010368,0.0000010542,0.0000010363,0.0000010474,0.0000010377,0.0000010482,0.0000010638,0.0000010726,0.0000010671,0.0000010606,0.0000010642,0.0000011002,0.0000010894,0.0000010865,0.0000011069,0.0000010895,0.0000011145,0.0000011012,0.0000011133,0.0000011270,0.0000011192,0.0000011085,0.0000011139,0.0000011147,0.0000011473,0.0000011334,0.0000011162,0.0000011525,0.0000011596,0.0000011521,0.0000011526,0.0000011513,0.0000011565,0.0000011707,0.0000011544,0.0000011700,0.0000011675,0.0000011654,0.0000011595,0.0000011859,0.0000011953,0.0000011843,0.0000012021,0.0000011821,0.0000011919,0.0000011648,0.0000011589,0.0000012148,0.0000011866,0.0000012094,0.0000011843,0.0000012106,0.0000012020,0.0000011922,0.0000011944,0.0000011987,0.0000012214,0.0000012059,0.0000012173,0.0000012209,0.0000012132,0.0000012072,0.0000012059,0.0000012015,0.0000012139,0.0000012245,0.0000012225,0.0000012298,0.0000012281,0.0000012278,0.0000012393,0.0000012349,0.0000012297,0.0000012063,0.0000012359,0.0000012269,0.0000012312,0.0000012433,0.0000012378,0.0000012274,0.0000012214,0.0000012526,0.0000012613,0.0000012470,0.0000008122,0.0000008285,0.0000008608,0.0000008700,0.0000008892,0.0000008916,0.0000008816,0.0000008893,0.0000009168,0.0000009199,0.0000009183,0.0000009473,0.0000009486,0.0000009493,0.0000009526,0.0000009648,0.0000009815,0.0000010031,0.0000009853,0.0000009955,0.0000009860,0.0000009978,0.0000010056,0.0000010174,0.0000010111,0.0000010047,0.0000010150,0.0000010450,0.0000010318,0.0000010355,0.0000010478,0.0000010359,0.0000010607,0.0000010473,0.0000010550,0.0000010706,0.0000010650,0.0000010515,0.0000010549,0.0000010625,0.0000010926,0.0000010752,0.0000010631,0.0000010932,0.0000011009,0.0000010945,0.0000010947,0.0000010921,0.0000010961,0.0000011093,0.0000010973,0.0000011154,0.0000011084,0.0000011051,0.0000010998,0.0000011296,0.0000011338,0.0000011246,0.0000011405,0.0000011281,0.0000011315,0.0000011099,0.0000010997,0.0000011577,0.0000011282,0.0000011473,0.0000011280,0.0000011494,0.0000011379,0.0000011332,0.0000011324,0.0000011402,0.0000011613,0.0000011463,0.0000011560,0.0000011600,0.0000011562,0.0000011498,0.0000011454,0.0000011429,0.0000011543,0.0000011632,0.0000011616,0.0000011681,0.0000011679,0.0000011677,0.0000011803,0.0000011694,0.0000011665,0.0000011467,0.0000011718,0.0000011629,0.0000011723,0.0000011813,0.0000011768,0.0000011655,0.0000011647,0.0000011926,0.0000011983,0.0000011829,0.0000007698,0.0000007863,0.0000008159,0.0000008264,0.0000008436,0.0000008489,0.0000008329,0.0000008437,0.0000008707,0.0000008748,0.0000008726,0.0000008997,0.0000009070,0.0000009002,0.0000009045,0.0000009165,0.0000009303,0.0000009532,0.0000009390,0.0000009442,0.0000009381,0.0000009474,0.0000009532,0.0000009686,0.0000009648,0.0000009534,0.0000009674,0.0000009915,0.0000009787,0.0000009830,0.0000009964,0.0000009859,0.0000010060,0.0000009910,0.0000009998,0.0000010189,0.0000010146,0.0000010009,0.0000010039,0.0000010132,0.0000010375,0.0000010229,0.0000010133,0.0000010417,0.0000010451,0.0000010357,0.0000010405,0.0000010394,0.0000010368,0.0000010524,0.0000010437,0.0000010621,0.0000010496,0.0000010524,0.0000010469,0.0000010805,0.0000010780,0.0000010686,0.0000010839,0.0000010734,0.0000010731,0.0000010496,0.0000010441,0.0000010988,0.0000010702,0.0000010916,0.0000010724,0.0000010955,0.0000010810,0.0000010734,0.0000010778,0.0000010820,0.0000011030,0.0000010890,0.0000010983,0.0000011013,0.0000011014,0.0000010933,0.0000010875,0.0000010858,0.0000010994,0.0000011043,0.0000011025,0.0000011126,0.0000011107,0.0000011090,0.0000011181,0.0000011125,0.0000011045,0.0000010904,0.0000011137,0.0000011011,0.0000011136,0.0000011215,0.0000011168,0.0000011061,0.0000011059,0.0000011327,0.0000011384,0.0000011224,0.0000007324,0.0000007452,0.0000007753,0.0000007847,0.0000008007,0.0000008056,0.0000007921,0.0000008022,0.0000008202,0.0000008306,0.0000008261,0.0000008533,0.0000008619,0.0000008561,0.0000008602,0.0000008726,0.0000008809,0.0000009051,0.0000008935,0.0000009001,0.0000008922,0.0000009022,0.0000009039,0.0000009218,0.0000009158,0.0000009027,0.0000009176,0.0000009396,0.0000009334,0.0000009344,0.0000009448,0.0000009336,0.0000009584,0.0000009426,0.0000009488,0.0000009699,0.0000009625,0.0000009523,0.0000009563,0.0000009658,0.0000009869,0.0000009722,0.0000009627,0.0000009886,0.0000009944,0.0000009839,0.0000009888,0.0000009884,0.0000009848,0.0000009982,0.0000009893,0.0000010118,0.0000009995,0.0000009969,0.0000009942,0.0000010302,0.0000010218,0.0000010161,0.0000010325,0.0000010169,0.0000010180,0.0000009967,0.0000009955,0.0000010421,0.0000010148,0.0000010354,0.0000010187,0.0000010367,0.0000010274,0.0000010202,0.0000010255,0.0000010292,0.0000010519,0.0000010332,0.0000010433,0.0000010485,0.0000010521,0.0000010423,0.0000010364,0.0000010335,0.0000010461,0.0000010490,0.0000010491,0.0000010573,0.0000010565,0.0000010562,0.0000010619,0.0000010532,0.0000010495,0.0000010390,0.0000010599,0.0000010459,0.0000010518,0.0000010615,0.0000010628,0.0000010487,0.0000010498,0.0000010777,0.0000010806,0.0000010717,0.0000006941,0.0000007101,0.0000007390,0.0000007480,0.0000007608,0.0000007636,0.0000007520,0.0000007623,0.0000007788,0.0000007864,0.0000007852,0.0000008115,0.0000008184,0.0000008139,0.0000008172,0.0000008308,0.0000008378,0.0000008556,0.0000008507,0.0000008553,0.0000008513,0.0000008578,0.0000008603,0.0000008759,0.0000008699,0.0000008565,0.0000008679,0.0000008946,0.0000008836,0.0000008898,0.0000008988,0.0000008886,0.0000009133,0.0000008920,0.0000009028,0.0000009208,0.0000009098,0.0000009059,0.0000009093,0.0000009154,0.0000009382,0.0000009252,0.0000009148,0.0000009380,0.0000009455,0.0000009328,0.0000009391,0.0000009402,0.0000009345,0.0000009468,0.0000009370,0.0000009607,0.0000009510,0.0000009433,0.0000009415,0.0000009779,0.0000009700,0.0000009640,0.0000009753,0.0000009674,0.0000009656,0.0000009492,0.0000009470,0.0000009886,0.0000009593,0.0000009851,0.0000009641,0.0000009822,0.0000009756,0.0000009646,0.0000009754,0.0000009740,0.0000009966,0.0000009831,0.0000009935,0.0000009913,0.0000010010,0.0000009857,0.0000009842,0.0000009799,0.0000009959,0.0000010004,0.0000009968,0.0000010059,0.0000010062,0.0000010063,0.0000010127,0.0000010031,0.0000009955,0.0000009886,0.0000010050,0.0000009909,0.0000010006,0.0000010093,0.0000010090,0.0000009937,0.0000009949,0.0000010215,0.0000010265,0.0000010175,0.0000006590,0.0000006728,0.0000007031,0.0000007118,0.0000007244,0.0000007270,0.0000007135,0.0000007266,0.0000007407,0.0000007468,0.0000007485,0.0000007704,0.0000007775,0.0000007702,0.0000007742,0.0000007920,0.0000007947,0.0000008104,0.0000008081,0.0000008138,0.0000008083,0.0000008177,0.0000008170,0.0000008349,0.0000008298,0.0000008140,0.0000008253,0.0000008509,0.0000008375,0.0000008460,0.0000008525,0.0000008429,0.0000008672,0.0000008470,0.0000008563,0.0000008769,0.0000008667,0.0000008609,0.0000008635,0.0000008691,0.0000008903,0.0000008796,0.0000008685,0.0000008883,0.0000008964,0.0000008866,0.0000008907,0.0000008908,0.0000008869,0.0000009030,0.0000008930,0.0000009130,0.0000009057,0.0000008957,0.0000008984,0.0000009265,0.0000009225,0.0000009132,0.0000009231,0.0000009211,0.0000009182,0.0000009032,0.0000008993,0.0000009404,0.0000009112,0.0000009374,0.0000009164,0.0000009334,0.0000009271,0.0000009145,0.0000009246,0.0000009239,0.0000009481,0.0000009346,0.0000009451,0.0000009416,0.0000009491,0.0000009397,0.0000009313,0.0000009314,0.0000009443,0.0000009507,0.0000009485,0.0000009555,0.0000009565,0.0000009587,0.0000009585,0.0000009495,0.0000009446,0.0000009383,0.0000009521,0.0000009416,0.0000009513,0.0000009575,0.0000009614,0.0000009478,0.0000009448,0.0000009687,0.0000009745,0.0000009631,0.0000006263,0.0000006399,0.0000006706,0.0000006769,0.0000006841,0.0000006923,0.0000006790,0.0000006893,0.0000007057,0.0000007107,0.0000007104,0.0000007328,0.0000007437,0.0000007343,0.0000007350,0.0000007541,0.0000007580,0.0000007719,0.0000007680,0.0000007736,0.0000007650,0.0000007797,0.0000007735,0.0000007954,0.0000007887,0.0000007729,0.0000007814,0.0000008112,0.0000007950,0.0000008048,0.0000008131,0.0000008046,0.0000008226,0.0000008052,0.0000008164,0.0000008288,0.0000008204,0.0000008150,0.0000008201,0.0000008287,0.0000008489,0.0000008368,0.0000008281,0.0000008398,0.0000008520,0.0000008403,0.0000008456,0.0000008446,0.0000008444,0.0000008594,0.0000008490,0.0000008698,0.0000008586,0.0000008500,0.0000008535,0.0000008810,0.0000008774,0.0000008683,0.0000008745,0.0000008771,0.0000008710,0.0000008602,0.0000008540,0.0000008926,0.0000008659,0.0000008945,0.0000008730,0.0000008873,0.0000008866,0.0000008714,0.0000008783,0.0000008795,0.0000009004,0.0000008863,0.0000008974,0.0000008957,0.0000008979,0.0000008963,0.0000008887,0.0000008847,0.0000008993,0.0000009008,0.0000009012,0.0000009075,0.0000009077,0.0000009108,0.0000009107,0.0000009054,0.0000008960,0.0000008934,0.0000009050,0.0000008912,0.0000009004,0.0000009080,0.0000009131,0.0000009004,0.0000008984,0.0000009234,0.0000009290,0.0000009161,0.0000005956,0.0000006090,0.0000006366,0.0000006428,0.0000006508,0.0000006577,0.0000006471,0.0000006563,0.0000006711,0.0000006760,0.0000006758,0.0000007001,0.0000007028,0.0000006919,0.0000006982,0.0000007160,0.0000007187,0.0000007359,0.0000007297,0.0000007385,0.0000007266,0.0000007403,0.0000007339,0.0000007564,0.0000007492,0.0000007376,0.0000007432,0.0000007698,0.0000007534,0.0000007641,0.0000007708,0.0000007635,0.0000007833,0.0000007663,0.0000007763,0.0000007903,0.0000007769,0.0000007719,0.0000007797,0.0000007878,0.0000008072,0.0000007948,0.0000007870,0.0000007977,0.0000008120,0.0000008008,0.0000008046,0.0000008017,0.0000008017,0.0000008145,0.0000008059,0.0000008305,0.0000008125,0.0000008076,0.0000008105,0.0000008388,0.0000008344,0.0000008239,0.0000008337,0.0000008352,0.0000008290,0.0000008176,0.0000008116,0.0000008467,0.0000008211,0.0000008494,0.0000008276,0.0000008450,0.0000008442,0.0000008327,0.0000008344,0.0000008399,0.0000008580,0.0000008422,0.0000008523,0.0000008533,0.0000008478,0.0000008541,0.0000008474,0.0000008382,0.0000008561,0.0000008557,0.0000008573,0.0000008659,0.0000008615,0.0000008694,0.0000008695,0.0000008619,0.0000008526,0.0000008512,0.0000008591,0.0000008500,0.0000008538,0.0000008637,0.0000008672,0.0000008531,0.0000008552,0.0000008804,0.0000008813,0.0000008677,0.0000005657,0.0000005784,0.0000006017,0.0000006134,0.0000006179,0.0000006221,0.0000006149,0.0000006243,0.0000006404,0.0000006428,0.0000006414,0.0000006647,0.0000006675,0.0000006580,0.0000006643,0.0000006811,0.0000006823,0.0000007001,0.0000006941,0.0000007018,0.0000006879,0.0000007059,0.0000006968,0.0000007187,0.0000007131,0.0000007007,0.0000007029,0.0000007307,0.0000007196,0.0000007291,0.0000007312,0.0000007223,0.0000007467,0.0000007274,0.0000007377,0.0000007507,0.0000007378,0.0000007331,0.0000007389,0.0000007465,0.0000007665,0.0000007567,0.0000007497,0.0000007568,0.0000007680,0.0000007578,0.0000007656,0.0000007631,0.0000007614,0.0000007752,0.0000007617,0.0000007913,0.0000007711,0.0000007668,0.0000007735,0.0000007920,0.0000007949,0.0000007809,0.0000007910,0.0000007937,0.0000007867,0.0000007807,0.0000007746,0.0000008073,0.0000007810,0.0000008095,0.0000007857,0.0000008031,0.0000008011,0.0000007936,0.0000007905,0.0000007986,0.0000008160,0.0000008026,0.0000008103,0.0000008095,0.0000008042,0.0000008132,0.0000008061,0.0000007941,0.0000008131,0.0000008148,0.0000008180,0.0000008222,0.0000008209,0.0000008230,0.0000008264,0.0000008225,0.0000008117,0.0000008104,0.0000008179,0.0000008078,0.0000008148,0.0000008261,0.0000008257,0.0000008068,0.0000008121,0.0000008397,0.0000008417,0.0000008236,0.0000005375,0.0000005485,0.0000005721,0.0000005813,0.0000005877,0.0000005892,0.0000005865,0.0000005946,0.0000006091,0.0000006115,0.0000006080,0.0000006289,0.0000006366,0.0000006223,0.0000006321,0.0000006462,0.0000006469,0.0000006643,0.0000006579,0.0000006688,0.0000006543,0.0000006687,0.0000006594,0.0000006854,0.0000006740,0.0000006643,0.0000006653,0.0000006943,0.0000006815,0.0000006932,0.0000006950,0.0000006854,0.0000007107,0.0000006915,0.0000007057,0.0000007138,0.0000006971,0.0000007003,0.0000007017,0.0000007085,0.0000007266,0.0000007183,0.0000007087,0.0000007213,0.0000007319,0.0000007231,0.0000007252,0.0000007219,0.0000007200,0.0000007360,0.0000007274,0.0000007528,0.0000007348,0.0000007327,0.0000007366,0.0000007513,0.0000007576,0.0000007389,0.0000007544,0.0000007530,0.0000007488,0.0000007419,0.0000007331,0.0000007711,0.0000007422,0.0000007688,0.0000007464,0.0000007630,0.0000007621,0.0000007524,0.0000007507,0.0000007569,0.0000007762,0.0000007656,0.0000007724,0.0000007665,0.0000007611,0.0000007723,0.0000007642,0.0000007563,0.0000007746,0.0000007759,0.0000007741,0.0000007800,0.0000007794,0.0000007809,0.0000007820,0.0000007807,0.0000007706,0.0000007689,0.0000007761,0.0000007680,0.0000007739,0.0000007848,0.0000007868,0.0000007682,0.0000007718,0.0000007965,0.0000007965,0.0000007825,0.0000005113,0.0000005207,0.0000005451,0.0000005522,0.0000005580,0.0000005587,0.0000005539,0.0000005651,0.0000005783,0.0000005822,0.0000005789,0.0000005975,0.0000006048,0.0000005926,0.0000005990,0.0000006111,0.0000006157,0.0000006278,0.0000006255,0.0000006368,0.0000006229,0.0000006333,0.0000006260,0.0000006519,0.0000006435,0.0000006313,0.0000006316,0.0000006616,0.0000006479,0.0000006611,0.0000006629,0.0000006520,0.0000006714,0.0000006564,0.0000006707,0.0000006757,0.0000006623,0.0000006654,0.0000006687,0.0000006739,0.0000006897,0.0000006834,0.0000006719,0.0000006863,0.0000006976,0.0000006880,0.0000006842,0.0000006844,0.0000006838,0.0000007002,0.0000006936,0.0000007134,0.0000006976,0.0000006959,0.0000007014,0.0000007139,0.0000007201,0.0000007023,0.0000007155,0.0000007167,0.0000007122,0.0000007030,0.0000006980,0.0000007335,0.0000007048,0.0000007302,0.0000007083,0.0000007251,0.0000007255,0.0000007143,0.0000007153,0.0000007213,0.0000007364,0.0000007248,0.0000007344,0.0000007322,0.0000007238,0.0000007335,0.0000007262,0.0000007178,0.0000007372,0.0000007388,0.0000007339,0.0000007430,0.0000007415,0.0000007454,0.0000007439,0.0000007411,0.0000007336,0.0000007329,0.0000007357,0.0000007298,0.0000007329,0.0000007450,0.0000007444,0.0000007295,0.0000007324,0.0000007604,0.0000007566,0.0000007420,0.0000004832,0.0000004946,0.0000005186,0.0000005265,0.0000005299,0.0000005316,0.0000005222,0.0000005355,0.0000005462,0.0000005507,0.0000005485,0.0000005669,0.0000005721,0.0000005653,0.0000005693,0.0000005817,0.0000005837,0.0000005932,0.0000005940,0.0000006025,0.0000005869,0.0000006020,0.0000005958,0.0000006188,0.0000006116,0.0000005996,0.0000005983,0.0000006280,0.0000006118,0.0000006237,0.0000006331,0.0000006199,0.0000006395,0.0000006236,0.0000006372,0.0000006437,0.0000006293,0.0000006332,0.0000006369,0.0000006395,0.0000006565,0.0000006475,0.0000006340,0.0000006507,0.0000006629,0.0000006523,0.0000006521,0.0000006473,0.0000006496,0.0000006651,0.0000006593,0.0000006797,0.0000006595,0.0000006605,0.0000006653,0.0000006743,0.0000006832,0.0000006675,0.0000006783,0.0000006797,0.0000006744,0.0000006677,0.0000006609,0.0000006979,0.0000006699,0.0000006921,0.0000006721,0.0000006900,0.0000006926,0.0000006810,0.0000006780,0.0000006824,0.0000007024,0.0000006897,0.0000006998,0.0000006962,0.0000006892,0.0000007018,0.0000006921,0.0000006839,0.0000007001,0.0000007035,0.0000006971,0.0000007078,0.0000007062,0.0000007095,0.0000007073,0.0000007037,0.0000006972,0.0000006977,0.0000006991,0.0000006956,0.0000006936,0.0000007119,0.0000007069,0.0000006928,0.0000006933,0.0000007231,0.0000007147,0.0000007054,0.0000004597,0.0000004712,0.0000004933,0.0000005000,0.0000005053,0.0000005038,0.0000004971,0.0000005080,0.0000005210,0.0000005235,0.0000005218,0.0000005394,0.0000005439,0.0000005351,0.0000005425,0.0000005567,0.0000005526,0.0000005655,0.0000005653,0.0000005706,0.0000005581,0.0000005723,0.0000005648,0.0000005899,0.0000005825,0.0000005702,0.0000005670,0.0000005958,0.0000005807,0.0000005915,0.0000006057,0.0000005905,0.0000006073,0.0000005896,0.0000006039,0.0000006103,0.0000005993,0.0000006038,0.0000006052,0.0000006080,0.0000006210,0.0000006171,0.0000006016,0.0000006196,0.0000006297,0.0000006183,0.0000006175,0.0000006133,0.0000006189,0.0000006278,0.0000006241,0.0000006476,0.0000006298,0.0000006272,0.0000006343,0.0000006378,0.0000006518,0.0000006335,0.0000006468,0.0000006457,0.0000006395,0.0000006324,0.0000006270,0.0000006651,0.0000006382,0.0000006596,0.0000006351,0.0000006521,0.0000006564,0.0000006470,0.0000006470,0.0000006476,0.0000006686,0.0000006563,0.0000006630,0.0000006627,0.0000006542,0.0000006657,0.0000006606,0.0000006460,0.0000006645,0.0000006688,0.0000006645,0.0000006698,0.0000006731,0.0000006735,0.0000006755,0.0000006666,0.0000006646,0.0000006640,0.0000006647,0.0000006582,0.0000006616,0.0000006764,0.0000006698,0.0000006605,0.0000006621,0.0000006891,0.0000006810,0.0000006729,0.0000004387,0.0000004484,0.0000004687,0.0000004763,0.0000004797,0.0000004794,0.0000004722,0.0000004832,0.0000004959,0.0000004976,0.0000004968,0.0000005129,0.0000005189,0.0000005082,0.0000005164,0.0000005325,0.0000005233,0.0000005377,0.0000005404,0.0000005417,0.0000005290,0.0000005449,0.0000005381,0.0000005641,0.0000005549,0.0000005440,0.0000005387,0.0000005655,0.0000005518,0.0000005633,0.0000005746,0.0000005630,0.0000005755,0.0000005609,0.0000005740,0.0000005818,0.0000005683,0.0000005714,0.0000005729,0.0000005767,0.0000005898,0.0000005869,0.0000005726,0.0000005886,0.0000005970,0.0000005840,0.0000005870,0.0000005859,0.0000005870,0.0000005986,0.0000005917,0.0000006143,0.0000005994,0.0000005954,0.0000006042,0.0000006044,0.0000006213,0.0000006036,0.0000006185,0.0000006147,0.0000006084,0.0000006001,0.0000005950,0.0000006304,0.0000006055,0.0000006251,0.0000006043,0.0000006208,0.0000006215,0.0000006122,0.0000006127,0.0000006176,0.0000006363,0.0000006251,0.0000006303,0.0000006318,0.0000006212,0.0000006327,0.0000006278,0.0000006128,0.0000006324,0.0000006349,0.0000006306,0.0000006367,0.0000006384,0.0000006412,0.0000006441,0.0000006364,0.0000006324,0.0000006327,0.0000006316,0.0000006244,0.0000006302,0.0000006452,0.0000006358,0.0000006265,0.0000006274,0.0000006554,0.0000006482,0.0000006410,0.0000004168,0.0000004276,0.0000004417,0.0000004536,0.0000004585,0.0000004546,0.0000004504,0.0000004598,0.0000004710,0.0000004732,0.0000004748,0.0000004863,0.0000004927,0.0000004808,0.0000004879,0.0000005063,0.0000004965,0.0000005126,0.0000005129,0.0000005168,0.0000005053,0.0000005185,0.0000005121,0.0000005352,0.0000005284,0.0000005145,0.0000005100,0.0000005370,0.0000005215,0.0000005326,0.0000005473,0.0000005341,0.0000005480,0.0000005333,0.0000005430,0.0000005565,0.0000005407,0.0000005443,0.0000005434,0.0000005492,0.0000005633,0.0000005558,0.0000005447,0.0000005607,0.0000005644,0.0000005549,0.0000005575,0.0000005589,0.0000005587,0.0000005694,0.0000005613,0.0000005816,0.0000005675,0.0000005669,0.0000005730,0.0000005749,0.0000005916,0.0000005740,0.0000005915,0.0000005858,0.0000005769,0.0000005718,0.0000005652,0.0000006026,0.0000005782,0.0000005966,0.0000005751,0.0000005908,0.0000005922,0.0000005778,0.0000005812,0.0000005878,0.0000006067,0.0000005939,0.0000005981,0.0000006000,0.0000005901,0.0000006005,0.0000005955,0.0000005850,0.0000006018,0.0000005989,0.0000006006,0.0000006068,0.0000006069,0.0000006082,0.0000006118,0.0000006048,0.0000005993,0.0000006037,0.0000005998,0.0000005929,0.0000006004,0.0000006121,0.0000006031,0.0000005968,0.0000005946,0.0000006216,0.0000006181,0.0000006085,0.0000003968,0.0000004082,0.0000004170,0.0000004302,0.0000004346,0.0000004317,0.0000004263,0.0000004369,0.0000004471,0.0000004501,0.0000004514,0.0000004621,0.0000004680,0.0000004571,0.0000004639,0.0000004820,0.0000004707,0.0000004852,0.0000004879,0.0000004899,0.0000004796,0.0000004915,0.0000004862,0.0000005065,0.0000005038,0.0000004868,0.0000004860,0.0000005089,0.0000004954,0.0000005067,0.0000005181,0.0000005098,0.0000005216,0.0000005095,0.0000005153,0.0000005299,0.0000005141,0.0000005189,0.0000005189,0.0000005207,0.0000005364,0.0000005284,0.0000005167,0.0000005303,0.0000005354,0.0000005302,0.0000005313,0.0000005299,0.0000005312,0.0000005400,0.0000005319,0.0000005528,0.0000005381,0.0000005389,0.0000005450,0.0000005472,0.0000005652,0.0000005450,0.0000005638,0.0000005587,0.0000005483,0.0000005434,0.0000005351,0.0000005759,0.0000005491,0.0000005685,0.0000005458,0.0000005594,0.0000005642,0.0000005499,0.0000005491,0.0000005591,0.0000005739,0.0000005643,0.0000005658,0.0000005734,0.0000005622,0.0000005694,0.0000005657,0.0000005589,0.0000005719,0.0000005689,0.0000005725,0.0000005752,0.0000005772,0.0000005802,0.0000005805,0.0000005741,0.0000005676,0.0000005736,0.0000005707,0.0000005647,0.0000005708,0.0000005809,0.0000005709,0.0000005669,0.0000005642,0.0000005907,0.0000005864,0.0000005775,0.0000003744,0.0000003859,0.0000003971,0.0000004093,0.0000004143,0.0000004100,0.0000004023,0.0000004152,0.0000004232,0.0000004275,0.0000004263,0.0000004405,0.0000004441,0.0000004336,0.0000004415,0.0000004589,0.0000004465,0.0000004630,0.0000004617,0.0000004685,0.0000004561,0.0000004683,0.0000004623,0.0000004818,0.0000004757,0.0000004638,0.0000004594,0.0000004844,0.0000004716,0.0000004811,0.0000004938,0.0000004840,0.0000004978,0.0000004852,0.0000004894,0.0000005015,0.0000004900,0.0000004946,0.0000004916,0.0000004947,0.0000005118,0.0000004986,0.0000004921,0.0000005050,0.0000005104,0.0000005048,0.0000005058,0.0000005020,0.0000005040,0.0000005125,0.0000005040,0.0000005255,0.0000005109,0.0000005121,0.0000005173,0.0000005194,0.0000005363,0.0000005187,0.0000005358,0.0000005336,0.0000005225,0.0000005150,0.0000005091,0.0000005483,0.0000005199,0.0000005395,0.0000005169,0.0000005319,0.0000005356,0.0000005218,0.0000005241,0.0000005341,0.0000005453,0.0000005367,0.0000005381,0.0000005448,0.0000005341,0.0000005427,0.0000005391,0.0000005279,0.0000005453,0.0000005394,0.0000005441,0.0000005468,0.0000005481,0.0000005500,0.0000005499,0.0000005451,0.0000005390,0.0000005421,0.0000005422,0.0000005379,0.0000005414,0.0000005527,0.0000005399,0.0000005369,0.0000005336,0.0000005620,0.0000005599,0.0000005491,0.0000003545,0.0000003677,0.0000003771,0.0000003896,0.0000003928,0.0000003879,0.0000003838,0.0000003950,0.0000004020,0.0000004041,0.0000004065,0.0000004184,0.0000004213,0.0000004149,0.0000004223,0.0000004348,0.0000004267,0.0000004381,0.0000004404,0.0000004469,0.0000004320,0.0000004484,0.0000004386,0.0000004590,0.0000004524,0.0000004403,0.0000004374,0.0000004605,0.0000004495,0.0000004575,0.0000004694,0.0000004612,0.0000004737,0.0000004631,0.0000004634,0.0000004756,0.0000004654,0.0000004687,0.0000004684,0.0000004710,0.0000004864,0.0000004723,0.0000004685,0.0000004809,0.0000004839,0.0000004799,0.0000004808,0.0000004755,0.0000004791,0.0000004874,0.0000004793,0.0000004997,0.0000004846,0.0000004882,0.0000004916,0.0000004932,0.0000005120,0.0000004935,0.0000005106,0.0000005051,0.0000004962,0.0000004876,0.0000004866,0.0000005218,0.0000004921,0.0000005139,0.0000004918,0.0000005061,0.0000005103,0.0000004950,0.0000004981,0.0000005080,0.0000005202,0.0000005116,0.0000005124,0.0000005139,0.0000005055,0.0000005165,0.0000005127,0.0000005037,0.0000005199,0.0000005102,0.0000005184,0.0000005182,0.0000005239,0.0000005219,0.0000005231,0.0000005182,0.0000005143,0.0000005133,0.0000005143,0.0000005086,0.0000005168,0.0000005287,0.0000005130,0.0000005092,0.0000005079,0.0000005322,0.0000005325,0.0000005208,0.0000003346,0.0000003487,0.0000003592,0.0000003708,0.0000003727,0.0000003688,0.0000003636,0.0000003754,0.0000003816,0.0000003863,0.0000003863,0.0000003982,0.0000004012,0.0000003949,0.0000004009,0.0000004126,0.0000004048,0.0000004169,0.0000004188,0.0000004255,0.0000004100,0.0000004253,0.0000004153,0.0000004349,0.0000004276,0.0000004168,0.0000004168,0.0000004358,0.0000004281,0.0000004364,0.0000004454,0.0000004362,0.0000004491,0.0000004374,0.0000004418,0.0000004533,0.0000004443,0.0000004448,0.0000004443,0.0000004480,0.0000004638,0.0000004497,0.0000004438,0.0000004544,0.0000004607,0.0000004565,0.0000004579,0.0000004533,0.0000004555,0.0000004644,0.0000004559,0.0000004785,0.0000004584,0.0000004613,0.0000004643,0.0000004707,0.0000004853,0.0000004687,0.0000004850,0.0000004782,0.0000004712,0.0000004660,0.0000004627,0.0000004921,0.0000004674,0.0000004908,0.0000004692,0.0000004818,0.0000004833,0.0000004698,0.0000004745,0.0000004838,0.0000004927,0.0000004874,0.0000004855,0.0000004882,0.0000004822,0.0000004935,0.0000004852,0.0000004793,0.0000004952,0.0000004814,0.0000004924,0.0000004896,0.0000004966,0.0000004958,0.0000004958,0.0000004934,0.0000004878,0.0000004874,0.0000004886,0.0000004841,0.0000004916,0.0000005045,0.0000004880,0.0000004823,0.0000004816,0.0000005064,0.0000005052,0.0000004943,0.0000003185,0.0000003306,0.0000003428,0.0000003515,0.0000003523,0.0000003522,0.0000003482,0.0000003592,0.0000003606,0.0000003666,0.0000003659,0.0000003796,0.0000003822,0.0000003735,0.0000003803,0.0000003936,0.0000003848,0.0000003969,0.0000003996,0.0000004047,0.0000003916,0.0000004048,0.0000003964,0.0000004141,0.0000004064,0.0000003955,0.0000003932,0.0000004123,0.0000004078,0.0000004140,0.0000004207,0.0000004152,0.0000004290,0.0000004154,0.0000004214,0.0000004313,0.0000004208,0.0000004242,0.0000004214,0.0000004248,0.0000004385,0.0000004287,0.0000004218,0.0000004338,0.0000004398,0.0000004326,0.0000004322,0.0000004289,0.0000004338,0.0000004413,0.0000004327,0.0000004548,0.0000004357,0.0000004406,0.0000004419,0.0000004471,0.0000004585,0.0000004468,0.0000004592,0.0000004530,0.0000004469,0.0000004433,0.0000004399,0.0000004654,0.0000004448,0.0000004656,0.0000004456,0.0000004552,0.0000004609,0.0000004482,0.0000004506,0.0000004624,0.0000004685,0.0000004635,0.0000004614,0.0000004627,0.0000004582,0.0000004674,0.0000004622,0.0000004560,0.0000004708,0.0000004554,0.0000004691,0.0000004641,0.0000004700,0.0000004697,0.0000004723,0.0000004698,0.0000004670,0.0000004617,0.0000004638,0.0000004628,0.0000004664,0.0000004772,0.0000004616,0.0000004611,0.0000004574,0.0000004812,0.0000004775,0.0000004702,0.0000003035,0.0000003126,0.0000003246,0.0000003341,0.0000003340,0.0000003344,0.0000003288,0.0000003438,0.0000003428,0.0000003494,0.0000003460,0.0000003589,0.0000003664,0.0000003531,0.0000003598,0.0000003731,0.0000003673,0.0000003765,0.0000003810,0.0000003848,0.0000003716,0.0000003870,0.0000003772,0.0000003939,0.0000003867,0.0000003754,0.0000003748,0.0000003935,0.0000003892,0.0000003922,0.0000004013,0.0000003938,0.0000004098,0.0000003953,0.0000003997,0.0000004082,0.0000003989,0.0000004040,0.0000003996,0.0000004003,0.0000004155,0.0000004082,0.0000004011,0.0000004124,0.0000004180,0.0000004110,0.0000004095,0.0000004095,0.0000004113,0.0000004168,0.0000004098,0.0000004338,0.0000004149,0.0000004193,0.0000004202,0.0000004234,0.0000004350,0.0000004253,0.0000004364,0.0000004308,0.0000004262,0.0000004194,0.0000004196,0.0000004438,0.0000004257,0.0000004434,0.0000004230,0.0000004328,0.0000004365,0.0000004257,0.0000004277,0.0000004379,0.0000004471,0.0000004410,0.0000004383,0.0000004396,0.0000004362,0.0000004466,0.0000004406,0.0000004348,0.0000004469,0.0000004326,0.0000004468,0.0000004411,0.0000004461,0.0000004443,0.0000004482,0.0000004450,0.0000004442,0.0000004368,0.0000004378,0.0000004383,0.0000004447,0.0000004551,0.0000004398,0.0000004385,0.0000004362,0.0000004567,0.0000004535,0.0000004456,0.0000002861,0.0000002951,0.0000003106,0.0000003170,0.0000003180,0.0000003178,0.0000003154,0.0000003256,0.0000003250,0.0000003320,0.0000003292,0.0000003384,0.0000003490,0.0000003363,0.0000003435,0.0000003566,0.0000003489,0.0000003561,0.0000003632,0.0000003675,0.0000003537,0.0000003685,0.0000003597,0.0000003724,0.0000003657,0.0000003558,0.0000003545,0.0000003756,0.0000003705,0.0000003725,0.0000003819,0.0000003738,0.0000003905,0.0000003753,0.0000003782,0.0000003856,0.0000003783,0.0000003851,0.0000003776,0.0000003803,0.0000003984,0.0000003896,0.0000003858,0.0000003907,0.0000003972,0.0000003920,0.0000003863,0.0000003874,0.0000003924,0.0000003955,0.0000003921,0.0000004118,0.0000003947,0.0000003991,0.0000003986,0.0000004019,0.0000004116,0.0000004050,0.0000004142,0.0000004080,0.0000004057,0.0000003988,0.0000004007,0.0000004216,0.0000004032,0.0000004192,0.0000004039,0.0000004119,0.0000004161,0.0000004022,0.0000004078,0.0000004134,0.0000004217,0.0000004196,0.0000004155,0.0000004185,0.0000004137,0.0000004244,0.0000004166,0.0000004133,0.0000004212,0.0000004115,0.0000004251,0.0000004189,0.0000004237,0.0000004214,0.0000004262,0.0000004217,0.0000004204,0.0000004170,0.0000004158,0.0000004163,0.0000004229,0.0000004329,0.0000004140,0.0000004179,0.0000004136,0.0000004331,0.0000004322,0.0000004238,0.0000002724,0.0000002778,0.0000002944,0.0000003016,0.0000003027,0.0000003013,0.0000003010,0.0000003092,0.0000003101,0.0000003161,0.0000003123,0.0000003201,0.0000003317,0.0000003184,0.0000003252,0.0000003373,0.0000003323,0.0000003367,0.0000003446,0.0000003484,0.0000003362,0.0000003492,0.0000003422,0.0000003560,0.0000003470,0.0000003369,0.0000003355,0.0000003565,0.0000003501,0.0000003535,0.0000003630,0.0000003573,0.0000003721,0.0000003578,0.0000003570,0.0000003661,0.0000003612,0.0000003665,0.0000003585,0.0000003617,0.0000003793,0.0000003714,0.0000003660,0.0000003712,0.0000003782,0.0000003724,0.0000003665,0.0000003647,0.0000003731,0.0000003785,0.0000003722,0.0000003912,0.0000003729,0.0000003780,0.0000003795,0.0000003839,0.0000003935,0.0000003861,0.0000003944,0.0000003881,0.0000003848,0.0000003812,0.0000003812,0.0000004027,0.0000003816,0.0000003978,0.0000003839,0.0000003926,0.0000003968,0.0000003843,0.0000003871,0.0000003922,0.0000003989,0.0000003962,0.0000003974,0.0000003988,0.0000003916,0.0000004016,0.0000003937,0.0000003918,0.0000004020,0.0000003913,0.0000004048,0.0000003995,0.0000004039,0.0000003983,0.0000004049,0.0000003998,0.0000004006,0.0000003976,0.0000003955,0.0000003951,0.0000004011,0.0000004118,0.0000003958,0.0000003987,0.0000003927,0.0000004132,0.0000004126,0.0000004043,0.0000002571,0.0000002626,0.0000002822,0.0000002874,0.0000002883,0.0000002871,0.0000002853,0.0000002939,0.0000002944,0.0000003013,0.0000002983,0.0000003025,0.0000003160,0.0000003043,0.0000003119,0.0000003220,0.0000003161,0.0000003192,0.0000003279,0.0000003315,0.0000003224,0.0000003320,0.0000003241,0.0000003377,0.0000003303,0.0000003213,0.0000003189,0.0000003384,0.0000003324,0.0000003347,0.0000003478,0.0000003394,0.0000003517,0.0000003404,0.0000003410,0.0000003489,0.0000003420,0.0000003488,0.0000003403,0.0000003437,0.0000003610,0.0000003526,0.0000003486,0.0000003527,0.0000003583,0.0000003551,0.0000003458,0.0000003455,0.0000003549,0.0000003587,0.0000003542,0.0000003728,0.0000003558,0.0000003611,0.0000003622,0.0000003645,0.0000003745,0.0000003667,0.0000003757,0.0000003693,0.0000003649,0.0000003622,0.0000003626,0.0000003831,0.0000003624,0.0000003789,0.0000003641,0.0000003727,0.0000003768,0.0000003649,0.0000003665,0.0000003714,0.0000003805,0.0000003764,0.0000003768,0.0000003778,0.0000003736,0.0000003810,0.0000003716,0.0000003732,0.0000003859,0.0000003746,0.0000003820,0.0000003824,0.0000003815,0.0000003805,0.0000003835,0.0000003792,0.0000003798,0.0000003792,0.0000003750,0.0000003751,0.0000003833,0.0000003913,0.0000003747,0.0000003807,0.0000003725,0.0000003939,0.0000003907,0.0000003855,0.0000002452,0.0000002515,0.0000002698,0.0000002726,0.0000002726,0.0000002740,0.0000002694,0.0000002802,0.0000002801,0.0000002882,0.0000002842,0.0000002878,0.0000002998,0.0000002907,0.0000002968,0.0000003066,0.0000003018,0.0000003053,0.0000003118,0.0000003151,0.0000003063,0.0000003146,0.0000003087,0.0000003198,0.0000003141,0.0000003039,0.0000003034,0.0000003206,0.0000003139,0.0000003175,0.0000003302,0.0000003221,0.0000003333,0.0000003260,0.0000003219,0.0000003286,0.0000003255,0.0000003301,0.0000003247,0.0000003257,0.0000003425,0.0000003359,0.0000003320,0.0000003348,0.0000003404,0.0000003386,0.0000003280,0.0000003304,0.0000003382,0.0000003410,0.0000003368,0.0000003546,0.0000003368,0.0000003450,0.0000003428,0.0000003471,0.0000003562,0.0000003491,0.0000003566,0.0000003503,0.0000003453,0.0000003426,0.0000003454,0.0000003636,0.0000003454,0.0000003610,0.0000003464,0.0000003543,0.0000003594,0.0000003462,0.0000003486,0.0000003530,0.0000003624,0.0000003589,0.0000003589,0.0000003575,0.0000003567,0.0000003612,0.0000003526,0.0000003556,0.0000003650,0.0000003558,0.0000003660,0.0000003621,0.0000003625,0.0000003592,0.0000003651,0.0000003613,0.0000003588,0.0000003598,0.0000003567,0.0000003561,0.0000003660,0.0000003701,0.0000003565,0.0000003605,0.0000003544,0.0000003750,0.0000003708,0.0000003655,0.0000002322,0.0000002386,0.0000002551,0.0000002598,0.0000002560,0.0000002590,0.0000002560,0.0000002666,0.0000002662,0.0000002732,0.0000002692,0.0000002743,0.0000002858,0.0000002763,0.0000002815,0.0000002903,0.0000002838,0.0000002899,0.0000002963,0.0000002988,0.0000002916,0.0000002994,0.0000002938,0.0000003046,0.0000002998,0.0000002891,0.0000002872,0.0000003041,0.0000003004,0.0000003024,0.0000003155,0.0000003054,0.0000003165,0.0000003091,0.0000003066,0.0000003136,0.0000003110,0.0000003147,0.0000003102,0.0000003093,0.0000003273,0.0000003199,0.0000003119,0.0000003177,0.0000003221,0.0000003218,0.0000003122,0.0000003119,0.0000003223,0.0000003244,0.0000003197,0.0000003356,0.0000003204,0.0000003262,0.0000003255,0.0000003297,0.0000003376,0.0000003310,0.0000003386,0.0000003315,0.0000003276,0.0000003280,0.0000003279,0.0000003465,0.0000003304,0.0000003421,0.0000003299,0.0000003361,0.0000003441,0.0000003283,0.0000003300,0.0000003360,0.0000003474,0.0000003431,0.0000003424,0.0000003383,0.0000003363,0.0000003406,0.0000003378,0.0000003381,0.0000003476,0.0000003364,0.0000003469,0.0000003445,0.0000003451,0.0000003426,0.0000003480,0.0000003423,0.0000003420,0.0000003416,0.0000003372,0.0000003382,0.0000003478,0.0000003491,0.0000003395,0.0000003431,0.0000003352,0.0000003591,0.0000003521,0.0000003456,0.0000002185,0.0000002263,0.0000002423,0.0000002468,0.0000002440,0.0000002441,0.0000002416,0.0000002528,0.0000002530,0.0000002582,0.0000002555,0.0000002601,0.0000002739,0.0000002600,0.0000002668,0.0000002764,0.0000002708,0.0000002735,0.0000002833,0.0000002861,0.0000002800,0.0000002847,0.0000002806,0.0000002881,0.0000002849,0.0000002742,0.0000002739,0.0000002890,0.0000002865,0.0000002868,0.0000002991,0.0000002897,0.0000003011,0.0000002943,0.0000002923,0.0000002961,0.0000002954,0.0000003008,0.0000002941,0.0000002944,0.0000003131,0.0000003047,0.0000002958,0.0000003020,0.0000003043,0.0000003053,0.0000002962,0.0000002976,0.0000003068,0.0000003085,0.0000003036,0.0000003168,0.0000003038,0.0000003104,0.0000003092,0.0000003122,0.0000003205,0.0000003162,0.0000003207,0.0000003149,0.0000003119,0.0000003134,0.0000003110,0.0000003294,0.0000003164,0.0000003259,0.0000003138,0.0000003175,0.0000003284,0.0000003108,0.0000003138,0.0000003188,0.0000003275,0.0000003264,0.0000003265,0.0000003199,0.0000003187,0.0000003224,0.0000003216,0.0000003221,0.0000003310,0.0000003160,0.0000003293,0.0000003293,0.0000003268,0.0000003266,0.0000003303,0.0000003240,0.0000003247,0.0000003250,0.0000003204,0.0000003203,0.0000003308,0.0000003311,0.0000003227,0.0000003271,0.0000003157,0.0000003407,0.0000003344,0.0000003289,0.0000002075,0.0000002145,0.0000002311,0.0000002348,0.0000002322,0.0000002331,0.0000002277,0.0000002395,0.0000002403,0.0000002468,0.0000002431,0.0000002471,0.0000002616,0.0000002462,0.0000002533,0.0000002620,0.0000002558,0.0000002586,0.0000002715,0.0000002723,0.0000002655,0.0000002685,0.0000002668,0.0000002740,0.0000002684,0.0000002591,0.0000002598,0.0000002742,0.0000002718,0.0000002731,0.0000002823,0.0000002764,0.0000002865,0.0000002822,0.0000002794,0.0000002826,0.0000002813,0.0000002876,0.0000002790,0.0000002801,0.0000002972,0.0000002896,0.0000002807,0.0000002878,0.0000002879,0.0000002915,0.0000002793,0.0000002833,0.0000002931,0.0000002914,0.0000002888,0.0000003023,0.0000002869,0.0000002942,0.0000002965,0.0000002961,0.0000003040,0.0000003006,0.0000003055,0.0000002985,0.0000002967,0.0000002973,0.0000002943,0.0000003137,0.0000003008,0.0000003110,0.0000002972,0.0000003037,0.0000003081,0.0000002966,0.0000002988,0.0000003028,0.0000003111,0.0000003107,0.0000003131,0.0000003038,0.0000003049,0.0000003054,0.0000003066,0.0000003050,0.0000003134,0.0000003030,0.0000003128,0.0000003140,0.0000003116,0.0000003103,0.0000003125,0.0000003077,0.0000003097,0.0000003093,0.0000003045,0.0000003022,0.0000003147,0.0000003147,0.0000003078,0.0000003102,0.0000002999,0.0000003246,0.0000003163,0.0000003120,0.0000001964,0.0000002017,0.0000002189,0.0000002245,0.0000002214,0.0000002195,0.0000002159,0.0000002276,0.0000002279,0.0000002337,0.0000002324,0.0000002340,0.0000002479,0.0000002343,0.0000002412,0.0000002482,0.0000002427,0.0000002464,0.0000002581,0.0000002591,0.0000002515,0.0000002542,0.0000002525,0.0000002594,0.0000002530,0.0000002474,0.0000002478,0.0000002607,0.0000002590,0.0000002589,0.0000002673,0.0000002630,0.0000002742,0.0000002668,0.0000002642,0.0000002681,0.0000002676,0.0000002743,0.0000002639,0.0000002658,0.0000002823,0.0000002747,0.0000002659,0.0000002723,0.0000002726,0.0000002797,0.0000002646,0.0000002682,0.0000002785,0.0000002767,0.0000002744,0.0000002877,0.0000002733,0.0000002795,0.0000002799,0.0000002796,0.0000002879,0.0000002853,0.0000002898,0.0000002846,0.0000002821,0.0000002815,0.0000002799,0.0000002978,0.0000002859,0.0000002938,0.0000002829,0.0000002900,0.0000002928,0.0000002819,0.0000002834,0.0000002845,0.0000002946,0.0000002949,0.0000002996,0.0000002880,0.0000002900,0.0000002890,0.0000002928,0.0000002891,0.0000002982,0.0000002887,0.0000002982,0.0000002988,0.0000002989,0.0000002930,0.0000002976,0.0000002936,0.0000002944,0.0000002917,0.0000002904,0.0000002872,0.0000003000,0.0000002981,0.0000002923,0.0000002942,0.0000002873,0.0000003076,0.0000002997,0.0000002970,0.0000001870,0.0000001920,0.0000002076,0.0000002134,0.0000002107,0.0000002069,0.0000002061,0.0000002177,0.0000002168,0.0000002225,0.0000002225,0.0000002223,0.0000002360,0.0000002233,0.0000002298,0.0000002349,0.0000002313,0.0000002344,0.0000002444,0.0000002468,0.0000002390,0.0000002436,0.0000002399,0.0000002448,0.0000002414,0.0000002339,0.0000002359,0.0000002481,0.0000002470,0.0000002464,0.0000002551,0.0000002490,0.0000002610,0.0000002524,0.0000002511,0.0000002548,0.0000002565,0.0000002606,0.0000002501,0.0000002533,0.0000002703,0.0000002628,0.0000002539,0.0000002591,0.0000002591,0.0000002679,0.0000002528,0.0000002552,0.0000002654,0.0000002618,0.0000002593,0.0000002750,0.0000002589,0.0000002666,0.0000002662,0.0000002670,0.0000002732,0.0000002715,0.0000002759,0.0000002692,0.0000002691,0.0000002667,0.0000002668,0.0000002838,0.0000002707,0.0000002799,0.0000002702,0.0000002755,0.0000002766,0.0000002673,0.0000002717,0.0000002713,0.0000002797,0.0000002800,0.0000002860,0.0000002738,0.0000002765,0.0000002765,0.0000002788,0.0000002758,0.0000002815,0.0000002742,0.0000002822,0.0000002838,0.0000002841,0.0000002790,0.0000002846,0.0000002802,0.0000002793,0.0000002759,0.0000002761,0.0000002718,0.0000002849,0.0000002847,0.0000002757,0.0000002790,0.0000002743,0.0000002934,0.0000002830,0.0000002829,0.0000001768,0.0000001826,0.0000001994,0.0000002027,0.0000001994,0.0000001953,0.0000001942,0.0000002066,0.0000002063,0.0000002113,0.0000002112,0.0000002127,0.0000002255,0.0000002127,0.0000002183,0.0000002244,0.0000002199,0.0000002229,0.0000002329,0.0000002365,0.0000002266,0.0000002315,0.0000002294,0.0000002331,0.0000002299,0.0000002246,0.0000002235,0.0000002365,0.0000002367,0.0000002344,0.0000002443,0.0000002377,0.0000002459,0.0000002405,0.0000002383,0.0000002413,0.0000002452,0.0000002485,0.0000002373,0.0000002420,0.0000002583,0.0000002499,0.0000002422,0.0000002457,0.0000002452,0.0000002548,0.0000002403,0.0000002424,0.0000002529,0.0000002498,0.0000002460,0.0000002629,0.0000002457,0.0000002534,0.0000002530,0.0000002516,0.0000002603,0.0000002568,0.0000002611,0.0000002583,0.0000002540,0.0000002539,0.0000002535,0.0000002691,0.0000002547,0.0000002668,0.0000002561,0.0000002639,0.0000002618,0.0000002518,0.0000002558,0.0000002581,0.0000002661,0.0000002659,0.0000002734,0.0000002595,0.0000002613,0.0000002626,0.0000002665,0.0000002650,0.0000002688,0.0000002614,0.0000002680,0.0000002710,0.0000002715,0.0000002669,0.0000002694,0.0000002655,0.0000002672,0.0000002616,0.0000002625,0.0000002588,0.0000002713,0.0000002703,0.0000002600,0.0000002652,0.0000002589,0.0000002788,0.0000002669,0.0000002682,0.0000001671,0.0000001746,0.0000001902,0.0000001930,0.0000001916,0.0000001863,0.0000001845,0.0000001984,0.0000001957,0.0000002017,0.0000002015,0.0000002018,0.0000002141,0.0000002047,0.0000002096,0.0000002132,0.0000002080,0.0000002123,0.0000002197,0.0000002231,0.0000002161,0.0000002213,0.0000002171,0.0000002201,0.0000002175,0.0000002127,0.0000002126,0.0000002249,0.0000002257,0.0000002229,0.0000002321,0.0000002264,0.0000002343,0.0000002290,0.0000002268,0.0000002294,0.0000002317,0.0000002353,0.0000002256,0.0000002297,0.0000002464,0.0000002383,0.0000002319,0.0000002328,0.0000002338,0.0000002413,0.0000002270,0.0000002308,0.0000002412,0.0000002372,0.0000002325,0.0000002480,0.0000002335,0.0000002411,0.0000002401,0.0000002377,0.0000002460,0.0000002443,0.0000002476,0.0000002468,0.0000002414,0.0000002420,0.0000002420,0.0000002565,0.0000002407,0.0000002542,0.0000002443,0.0000002489,0.0000002484,0.0000002395,0.0000002428,0.0000002459,0.0000002507,0.0000002521,0.0000002591,0.0000002455,0.0000002477,0.0000002493,0.0000002546,0.0000002510,0.0000002542,0.0000002502,0.0000002550,0.0000002560,0.0000002587,0.0000002538,0.0000002547,0.0000002511,0.0000002543,0.0000002481,0.0000002505,0.0000002460,0.0000002584,0.0000002562,0.0000002471,0.0000002515,0.0000002449,0.0000002660,0.0000002534,0.0000002558,0.0000001590,0.0000001665,0.0000001793,0.0000001837,0.0000001819,0.0000001770,0.0000001764,0.0000001875,0.0000001861,0.0000001910,0.0000001915,0.0000001921,0.0000002028,0.0000001947,0.0000001983,0.0000002026,0.0000001976,0.0000002015,0.0000002098,0.0000002108,0.0000002047,0.0000002113,0.0000002053,0.0000002099,0.0000002074,0.0000002011,0.0000002025,0.0000002131,0.0000002143,0.0000002105,0.0000002181,0.0000002162,0.0000002222,0.0000002180,0.0000002153,0.0000002176,0.0000002190,0.0000002242,0.0000002150,0.0000002183,0.0000002350,0.0000002264,0.0000002216,0.0000002212,0.0000002238,0.0000002300,0.0000002147,0.0000002180,0.0000002279,0.0000002244,0.0000002201,0.0000002360,0.0000002212,0.0000002299,0.0000002286,0.0000002253,0.0000002330,0.0000002316,0.0000002350,0.0000002336,0.0000002281,0.0000002318,0.0000002305,0.0000002438,0.0000002294,0.0000002414,0.0000002328,0.0000002367,0.0000002357,0.0000002247,0.0000002306,0.0000002327,0.0000002376,0.0000002390,0.0000002460,0.0000002338,0.0000002370,0.0000002369,0.0000002414,0.0000002385,0.0000002406,0.0000002380,0.0000002423,0.0000002425,0.0000002443,0.0000002416,0.0000002403,0.0000002384,0.0000002412,0.0000002367,0.0000002384,0.0000002352,0.0000002454,0.0000002452,0.0000002349,0.0000002406,0.0000002327,0.0000002513,0.0000002409,0.0000002424,0.0000001501,0.0000001582,0.0000001712,0.0000001728,0.0000001735,0.0000001683,0.0000001672,0.0000001780,0.0000001769,0.0000001816,0.0000001812,0.0000001819,0.0000001933,0.0000001851,0.0000001887,0.0000001927,0.0000001876,0.0000001938,0.0000001999,0.0000001993,0.0000001941,0.0000001993,0.0000001957,0.0000001995,0.0000001980,0.0000001909,0.0000001918,0.0000002036,0.0000002034,0.0000002008,0.0000002083,0.0000002050,0.0000002110,0.0000002061,0.0000002038,0.0000002061,0.0000002072,0.0000002131,0.0000002048,0.0000002069,0.0000002219,0.0000002148,0.0000002095,0.0000002123,0.0000002108,0.0000002185,0.0000002046,0.0000002071,0.0000002175,0.0000002139,0.0000002085,0.0000002235,0.0000002090,0.0000002172,0.0000002188,0.0000002142,0.0000002218,0.0000002189,0.0000002228,0.0000002235,0.0000002184,0.0000002209,0.0000002183,0.0000002303,0.0000002166,0.0000002311,0.0000002194,0.0000002255,0.0000002228,0.0000002142,0.0000002208,0.0000002215,0.0000002267,0.0000002274,0.0000002333,0.0000002213,0.0000002257,0.0000002256,0.0000002285,0.0000002236,0.0000002307,0.0000002256,0.0000002316,0.0000002308,0.0000002320,0.0000002286,0.0000002297,0.0000002259,0.0000002294,0.0000002251,0.0000002260,0.0000002243,0.0000002335,0.0000002323,0.0000002226,0.0000002296,0.0000002208,0.0000002388,0.0000002297,0.0000002320,0.0000001419,0.0000001497,0.0000001619,0.0000001655,0.0000001655,0.0000001605,0.0000001571,0.0000001684,0.0000001682,0.0000001723,0.0000001726,0.0000001720,0.0000001834,0.0000001757,0.0000001797,0.0000001826,0.0000001780,0.0000001841,0.0000001910,0.0000001884,0.0000001839,0.0000001886,0.0000001847,0.0000001908,0.0000001876,0.0000001812,0.0000001806,0.0000001951,0.0000001940,0.0000001903,0.0000001983,0.0000001954,0.0000001990,0.0000001951,0.0000001926,0.0000001960,0.0000001963,0.0000002013,0.0000001937,0.0000001963,0.0000002116,0.0000002044,0.0000001994,0.0000002012,0.0000001998,0.0000002067,0.0000001949,0.0000001989,0.0000002078,0.0000002028,0.0000001990,0.0000002126,0.0000001992,0.0000002063,0.0000002086,0.0000002035,0.0000002120,0.0000002073,0.0000002124,0.0000002130,0.0000002083,0.0000002109,0.0000002066,0.0000002182,0.0000002065,0.0000002204,0.0000002077,0.0000002147,0.0000002123,0.0000002041,0.0000002100,0.0000002119,0.0000002152,0.0000002173,0.0000002233,0.0000002109,0.0000002149,0.0000002134,0.0000002173,0.0000002119,0.0000002195,0.0000002146,0.0000002205,0.0000002194,0.0000002207,0.0000002186,0.0000002184,0.0000002147,0.0000002171,0.0000002140,0.0000002155,0.0000002125,0.0000002230,0.0000002195,0.0000002113,0.0000002162,0.0000002106,0.0000002288,0.0000002187,0.0000002210,0.0000001351,0.0000001435,0.0000001527,0.0000001572,0.0000001582,0.0000001524,0.0000001491,0.0000001590,0.0000001593,0.0000001634,0.0000001640,0.0000001634,0.0000001739,0.0000001660,0.0000001713,0.0000001716,0.0000001678,0.0000001758,0.0000001821,0.0000001796,0.0000001748,0.0000001813,0.0000001735,0.0000001825,0.0000001795,0.0000001723,0.0000001718,0.0000001848,0.0000001833,0.0000001781,0.0000001893,0.0000001860,0.0000001885,0.0000001844,0.0000001824,0.0000001847,0.0000001866,0.0000001907,0.0000001848,0.0000001857,0.0000002004,0.0000001938,0.0000001887,0.0000001915,0.0000001893,0.0000001966,0.0000001853,0.0000001891,0.0000001985,0.0000001911,0.0000001875,0.0000002042,0.0000001890,0.0000001967,0.0000001989,0.0000001921,0.0000002006,0.0000001984,0.0000002019,0.0000002039,0.0000001988,0.0000001995,0.0000001942,0.0000002088,0.0000001966,0.0000002086,0.0000001980,0.0000002036,0.0000002003,0.0000001942,0.0000001993,0.0000002006,0.0000002035,0.0000002068,0.0000002118,0.0000001992,0.0000002044,0.0000002015,0.0000002078,0.0000002014,0.0000002088,0.0000002042,0.0000002084,0.0000002077,0.0000002095,0.0000002075,0.0000002051,0.0000002016,0.0000002078,0.0000002053,0.0000002037,0.0000002011,0.0000002130,0.0000002066,0.0000002012,0.0000002061,0.0000001991,0.0000002178,0.0000002062,0.0000002095,0.0000001275,0.0000001375,0.0000001456,0.0000001488,0.0000001504,0.0000001457,0.0000001408,0.0000001526,0.0000001508,0.0000001556,0.0000001557,0.0000001557,0.0000001673,0.0000001571,0.0000001633,0.0000001620,0.0000001592,0.0000001686,0.0000001723,0.0000001701,0.0000001655,0.0000001709,0.0000001641,0.0000001736,0.0000001708,0.0000001636,0.0000001621,0.0000001757,0.0000001724,0.0000001679,0.0000001806,0.0000001761,0.0000001797,0.0000001738,0.0000001736,0.0000001771,0.0000001783,0.0000001811,0.0000001748,0.0000001766,0.0000001903,0.0000001831,0.0000001803,0.0000001802,0.0000001794,0.0000001866,0.0000001765,0.0000001784,0.0000001884,0.0000001823,0.0000001793,0.0000001936,0.0000001813,0.0000001877,0.0000001891,0.0000001830,0.0000001917,0.0000001895,0.0000001916,0.0000001935,0.0000001890,0.0000001899,0.0000001833,0.0000001979,0.0000001870,0.0000002014,0.0000001861,0.0000001940,0.0000001900,0.0000001847,0.0000001894,0.0000001891,0.0000001939,0.0000001954,0.0000002025,0.0000001901,0.0000001959,0.0000001913,0.0000001968,0.0000001908,0.0000001985,0.0000001950,0.0000001958,0.0000001970,0.0000001968,0.0000001978,0.0000001946,0.0000001903,0.0000001971,0.0000001942,0.0000001949,0.0000001915,0.0000002012,0.0000001949,0.0000001933,0.0000001950,0.0000001870,0.0000002092,0.0000001972,0.0000002008,0.0000001207,0.0000001319,0.0000001393,0.0000001402,0.0000001440,0.0000001387,0.0000001345,0.0000001448,0.0000001434,0.0000001482,0.0000001478,0.0000001475,0.0000001580,0.0000001491,0.0000001557,0.0000001538,0.0000001510,0.0000001589,0.0000001633,0.0000001635,0.0000001566,0.0000001621,0.0000001564,0.0000001650,0.0000001614,0.0000001559,0.0000001545,0.0000001660,0.0000001660,0.0000001586,0.0000001702,0.0000001668,0.0000001705,0.0000001625,0.0000001636,0.0000001701,0.0000001685,0.0000001729,0.0000001662,0.0000001664,0.0000001810,0.0000001741,0.0000001724,0.0000001717,0.0000001722,0.0000001758,0.0000001674,0.0000001696,0.0000001794,0.0000001720,0.0000001706,0.0000001841,0.0000001720,0.0000001773,0.0000001790,0.0000001747,0.0000001832,0.0000001806,0.0000001816,0.0000001817,0.0000001799,0.0000001808,0.0000001742,0.0000001886,0.0000001771,0.0000001897,0.0000001771,0.0000001833,0.0000001805,0.0000001751,0.0000001802,0.0000001790,0.0000001841,0.0000001852,0.0000001911,0.0000001809,0.0000001870,0.0000001818,0.0000001886,0.0000001807,0.0000001894,0.0000001846,0.0000001867,0.0000001864,0.0000001868,0.0000001876,0.0000001846,0.0000001798,0.0000001869,0.0000001835,0.0000001860,0.0000001824,0.0000001898,0.0000001850,0.0000001837,0.0000001851,0.0000001796,0.0000002007,0.0000001876,0.0000001910,0.0000001159,0.0000001265,0.0000001322,0.0000001331,0.0000001363,0.0000001319,0.0000001270,0.0000001378,0.0000001344,0.0000001399,0.0000001388,0.0000001400,0.0000001492,0.0000001412,0.0000001459,0.0000001471,0.0000001436,0.0000001502,0.0000001538,0.0000001548,0.0000001493,0.0000001543,0.0000001476,0.0000001555,0.0000001528,0.0000001494,0.0000001462,0.0000001558,0.0000001563,0.0000001506,0.0000001606,0.0000001586,0.0000001619,0.0000001554,0.0000001544,0.0000001623,0.0000001602,0.0000001639,0.0000001580,0.0000001569,0.0000001737,0.0000001657,0.0000001631,0.0000001629,0.0000001628,0.0000001675,0.0000001596,0.0000001614,0.0000001714,0.0000001639,0.0000001623,0.0000001743,0.0000001648,0.0000001678,0.0000001695,0.0000001672,0.0000001718,0.0000001718,0.0000001732,0.0000001744,0.0000001711,0.0000001725,0.0000001670,0.0000001799,0.0000001679,0.0000001804,0.0000001686,0.0000001749,0.0000001709,0.0000001653,0.0000001717,0.0000001683,0.0000001759,0.0000001773,0.0000001827,0.0000001724,0.0000001775,0.0000001718,0.0000001796,0.0000001699,0.0000001798,0.0000001762,0.0000001773,0.0000001766,0.0000001754,0.0000001778,0.0000001752,0.0000001709,0.0000001793,0.0000001743,0.0000001761,0.0000001736,0.0000001810,0.0000001749,0.0000001739,0.0000001772,0.0000001719,0.0000001906,0.0000001790,0.0000001819,0.0000001088,0.0000001205,0.0000001259,0.0000001274,0.0000001300,0.0000001256,0.0000001197,0.0000001316,0.0000001275,0.0000001324,0.0000001318,0.0000001338,0.0000001404,0.0000001340,0.0000001388,0.0000001399,0.0000001351,0.0000001439,0.0000001465,0.0000001484,0.0000001432,0.0000001462,0.0000001402,0.0000001488,0.0000001450,0.0000001415,0.0000001392,0.0000001464,0.0000001491,0.0000001427,0.0000001528,0.0000001517,0.0000001544,0.0000001478,0.0000001481,0.0000001553,0.0000001513,0.0000001567,0.0000001496,0.0000001503,0.0000001641,0.0000001580,0.0000001562,0.0000001545,0.0000001552,0.0000001592,0.0000001508,0.0000001526,0.0000001646,0.0000001552,0.0000001533,0.0000001643,0.0000001567,0.0000001594,0.0000001614,0.0000001587,0.0000001618,0.0000001627,0.0000001652,0.0000001651,0.0000001619,0.0000001628,0.0000001586,0.0000001714,0.0000001601,0.0000001715,0.0000001607,0.0000001641,0.0000001623,0.0000001556,0.0000001624,0.0000001597,0.0000001664,0.0000001680,0.0000001729,0.0000001634,0.0000001679,0.0000001643,0.0000001702,0.0000001613,0.0000001745,0.0000001687,0.0000001684,0.0000001687,0.0000001662,0.0000001693,0.0000001652,0.0000001623,0.0000001716,0.0000001670,0.0000001659,0.0000001656,0.0000001722,0.0000001651,0.0000001633,0.0000001679,0.0000001646,0.0000001822,0.0000001706,0.0000001726,0.0000001037,0.0000001144,0.0000001208,0.0000001210,0.0000001236,0.0000001201,0.0000001150,0.0000001258,0.0000001199,0.0000001263,0.0000001254,0.0000001275,0.0000001341,0.0000001275,0.0000001316,0.0000001320,0.0000001288,0.0000001363,0.0000001384,0.0000001411,0.0000001363,0.0000001389,0.0000001335,0.0000001427,0.0000001370,0.0000001343,0.0000001308,0.0000001401,0.0000001409,0.0000001362,0.0000001440,0.0000001434,0.0000001478,0.0000001410,0.0000001408,0.0000001456,0.0000001438,0.0000001493,0.0000001415,0.0000001428,0.0000001554,0.0000001493,0.0000001470,0.0000001469,0.0000001484,0.0000001495,0.0000001439,0.0000001432,0.0000001562,0.0000001482,0.0000001462,0.0000001566,0.0000001488,0.0000001503,0.0000001541,0.0000001516,0.0000001539,0.0000001541,0.0000001562,0.0000001579,0.0000001512,0.0000001533,0.0000001527,0.0000001628,0.0000001513,0.0000001633,0.0000001534,0.0000001557,0.0000001541,0.0000001476,0.0000001545,0.0000001512,0.0000001590,0.0000001590,0.0000001637,0.0000001545,0.0000001595,0.0000001557,0.0000001617,0.0000001525,0.0000001661,0.0000001603,0.0000001597,0.0000001613,0.0000001567,0.0000001611,0.0000001575,0.0000001547,0.0000001628,0.0000001574,0.0000001592,0.0000001562,0.0000001640,0.0000001575,0.0000001562,0.0000001603,0.0000001570,0.0000001727,0.0000001621,0.0000001658,0.0000000980,0.0000001087,0.0000001146,0.0000001150,0.0000001182,0.0000001132,0.0000001098,0.0000001194,0.0000001144,0.0000001197,0.0000001210,0.0000001210,0.0000001288,0.0000001206,0.0000001242,0.0000001259,0.0000001216,0.0000001298,0.0000001309,0.0000001325,0.0000001297,0.0000001306,0.0000001261,0.0000001352,0.0000001315,0.0000001268,0.0000001236,0.0000001343,0.0000001349,0.0000001301,0.0000001370,0.0000001360,0.0000001418,0.0000001337,0.0000001336,0.0000001408,0.0000001374,0.0000001412,0.0000001350,0.0000001363,0.0000001481,0.0000001426,0.0000001400,0.0000001408,0.0000001417,0.0000001422,0.0000001349,0.0000001363,0.0000001495,0.0000001398,0.0000001383,0.0000001485,0.0000001411,0.0000001434,0.0000001475,0.0000001443,0.0000001464,0.0000001465,0.0000001497,0.0000001500,0.0000001436,0.0000001461,0.0000001459,0.0000001557,0.0000001421,0.0000001541,0.0000001457,0.0000001480,0.0000001454,0.0000001405,0.0000001460,0.0000001435,0.0000001513,0.0000001518,0.0000001555,0.0000001461,0.0000001519,0.0000001483,0.0000001534,0.0000001445,0.0000001598,0.0000001523,0.0000001524,0.0000001556,0.0000001484,0.0000001524,0.0000001496,0.0000001466,0.0000001555,0.0000001488,0.0000001516,0.0000001477,0.0000001559,0.0000001502,0.0000001491,0.0000001524,0.0000001497,0.0000001642,0.0000001547,0.0000001575,0.0000000936,0.0000001026,0.0000001092,0.0000001097,0.0000001122,0.0000001071,0.0000001051,0.0000001123,0.0000001095,0.0000001139,0.0000001138,0.0000001149,0.0000001221,0.0000001146,0.0000001185,0.0000001186,0.0000001164,0.0000001236,0.0000001239,0.0000001264,0.0000001236,0.0000001251,0.0000001195,0.0000001282,0.0000001252,0.0000001208,0.0000001172,0.0000001276,0.0000001271,0.0000001247,0.0000001307,0.0000001291,0.0000001352,0.0000001266,0.0000001279,0.0000001338,0.0000001286,0.0000001342,0.0000001292,0.0000001297,0.0000001408,0.0000001360,0.0000001338,0.0000001338,0.0000001341,0.0000001345,0.0000001289,0.0000001280,0.0000001427,0.0000001334,0.0000001327,0.0000001418,0.0000001326,0.0000001367,0.0000001392,0.0000001370,0.0000001396,0.0000001399,0.0000001424,0.0000001419,0.0000001361,0.0000001386,0.0000001392,0.0000001480,0.0000001355,0.0000001463,0.0000001388,0.0000001402,0.0000001400,0.0000001343,0.0000001389,0.0000001360,0.0000001422,0.0000001445,0.0000001480,0.0000001392,0.0000001457,0.0000001417,0.0000001460,0.0000001363,0.0000001523,0.0000001425,0.0000001449,0.0000001473,0.0000001412,0.0000001459,0.0000001424,0.0000001413,0.0000001472,0.0000001400,0.0000001436,0.0000001400,0.0000001474,0.0000001413,0.0000001424,0.0000001457,0.0000001428,0.0000001563,0.0000001473,0.0000001489,0.0000000898,0.0000000977,0.0000001038,0.0000001048,0.0000001065,0.0000001013,0.0000000991,0.0000001068,0.0000001045,0.0000001082,0.0000001078,0.0000001091,0.0000001149,0.0000001081,0.0000001131,0.0000001110,0.0000001111,0.0000001171,0.0000001176,0.0000001208,0.0000001162,0.0000001188,0.0000001122,0.0000001225,0.0000001186,0.0000001150,0.0000001120,0.0000001209,0.0000001215,0.0000001196,0.0000001257,0.0000001222,0.0000001268,0.0000001201,0.0000001215,0.0000001277,0.0000001227,0.0000001273,0.0000001226,0.0000001236,0.0000001346,0.0000001309,0.0000001264,0.0000001269,0.0000001287,0.0000001283,0.0000001222,0.0000001206,0.0000001350,0.0000001273,0.0000001258,0.0000001353,0.0000001264,0.0000001287,0.0000001337,0.0000001294,0.0000001333,0.0000001340,0.0000001360,0.0000001358,0.0000001297,0.0000001324,0.0000001323,0.0000001412,0.0000001285,0.0000001392,0.0000001324,0.0000001329,0.0000001328,0.0000001283,0.0000001325,0.0000001302,0.0000001352,0.0000001374,0.0000001408,0.0000001323,0.0000001366,0.0000001354,0.0000001370,0.0000001287,0.0000001445,0.0000001348,0.0000001375,0.0000001405,0.0000001345,0.0000001378,0.0000001337,0.0000001340,0.0000001404,0.0000001330,0.0000001368,0.0000001319,0.0000001416,0.0000001350,0.0000001364,0.0000001390,0.0000001347,0.0000001483,0.0000001394,0.0000001420,0.0000000856,0.0000000925,0.0000000979,0.0000001001,0.0000001005,0.0000000967,0.0000000954,0.0000001024,0.0000000989,0.0000001029,0.0000001037,0.0000001031,0.0000001086,0.0000001028,0.0000001079,0.0000001045,0.0000001050,0.0000001114,0.0000001111,0.0000001137,0.0000001094,0.0000001140,0.0000001068,0.0000001168,0.0000001134,0.0000001095,0.0000001072,0.0000001155,0.0000001156,0.0000001137,0.0000001188,0.0000001152,0.0000001207,0.0000001136,0.0000001156,0.0000001218,0.0000001167,0.0000001193,0.0000001167,0.0000001174,0.0000001273,0.0000001241,0.0000001198,0.0000001200,0.0000001222,0.0000001217,0.0000001146,0.0000001138,0.0000001276,0.0000001225,0.0000001181,0.0000001296,0.0000001185,0.0000001208,0.0000001282,0.0000001239,0.0000001283,0.0000001271,0.0000001299,0.0000001289,0.0000001230,0.0000001259,0.0000001253,0.0000001351,0.0000001219,0.0000001328,0.0000001257,0.0000001246,0.0000001270,0.0000001209,0.0000001241,0.0000001231,0.0000001270,0.0000001310,0.0000001342,0.0000001269,0.0000001300,0.0000001283,0.0000001301,0.0000001228,0.0000001356,0.0000001283,0.0000001307,0.0000001330,0.0000001288,0.0000001304,0.0000001287,0.0000001280,0.0000001344,0.0000001264,0.0000001288,0.0000001257,0.0000001347,0.0000001272,0.0000001295,0.0000001337,0.0000001282,0.0000001410,0.0000001324,0.0000001353,0.0000000816,0.0000000881,0.0000000930,0.0000000951,0.0000000961,0.0000000915,0.0000000898,0.0000000979,0.0000000943,0.0000000983,0.0000000986,0.0000000976,0.0000001036,0.0000000975,0.0000001024,0.0000000989,0.0000000999,0.0000001054,0.0000001059,0.0000001074,0.0000001044,0.0000001074,0.0000001006,0.0000001103,0.0000001063,0.0000001038,0.0000001019,0.0000001096,0.0000001101,0.0000001100,0.0000001120,0.0000001104,0.0000001155,0.0000001085,0.0000001100,0.0000001176,0.0000001106,0.0000001133,0.0000001111,0.0000001112,0.0000001209,0.0000001165,0.0000001148,0.0000001139,0.0000001174,0.0000001159,0.0000001070,0.0000001097,0.0000001215,0.0000001166,0.0000001115,0.0000001224,0.0000001133,0.0000001151,0.0000001214,0.0000001177,0.0000001223,0.0000001213,0.0000001237,0.0000001226,0.0000001166,0.0000001200,0.0000001206,0.0000001294,0.0000001158,0.0000001266,0.0000001186,0.0000001189,0.0000001206,0.0000001156,0.0000001180,0.0000001169,0.0000001200,0.0000001247,0.0000001280,0.0000001213,0.0000001224,0.0000001230,0.0000001240,0.0000001176,0.0000001310,0.0000001239,0.0000001246,0.0000001274,0.0000001218,0.0000001231,0.0000001222,0.0000001216,0.0000001274,0.0000001219,0.0000001227,0.0000001194,0.0000001276,0.0000001200,0.0000001232,0.0000001268,0.0000001210,0.0000001340,0.0000001248,0.0000001285,0.0000000771,0.0000000840,0.0000000877,0.0000000893,0.0000000925,0.0000000858,0.0000000850,0.0000000929,0.0000000897,0.0000000950,0.0000000940,0.0000000916,0.0000000985,0.0000000918,0.0000000977,0.0000000940,0.0000000941,0.0000001001,0.0000001004,0.0000001012,0.0000001004,0.0000001032,0.0000000954,0.0000001056,0.0000001011,0.0000000990,0.0000000973,0.0000001041,0.0000001047,0.0000001044,0.0000001054,0.0000001057,0.0000001102,0.0000001029,0.0000001027,0.0000001120,0.0000001048,0.0000001083,0.0000001057,0.0000001065,0.0000001144,0.0000001110,0.0000001085,0.0000001101,0.0000001117,0.0000001100,0.0000001021,0.0000001045,0.0000001145,0.0000001098,0.0000001060,0.0000001178,0.0000001081,0.0000001098,0.0000001164,0.0000001123,0.0000001173,0.0000001157,0.0000001173,0.0000001166,0.0000001104,0.0000001143,0.0000001144,0.0000001222,0.0000001095,0.0000001196,0.0000001120,0.0000001134,0.0000001132,0.0000001098,0.0000001123,0.0000001107,0.0000001149,0.0000001184,0.0000001219,0.0000001134,0.0000001152,0.0000001161,0.0000001175,0.0000001107,0.0000001240,0.0000001179,0.0000001183,0.0000001211,0.0000001161,0.0000001169,0.0000001152,0.0000001151,0.0000001209,0.0000001153,0.0000001158,0.0000001138,0.0000001207,0.0000001150,0.0000001171,0.0000001190,0.0000001157,0.0000001265,0.0000001196,0.0000001221,0.0000000733,0.0000000799,0.0000000826,0.0000000855,0.0000000881,0.0000000806,0.0000000811,0.0000000879,0.0000000845,0.0000000907,0.0000000895,0.0000000882,0.0000000943,0.0000000873,0.0000000935,0.0000000899,0.0000000898,0.0000000955,0.0000000953,0.0000000953,0.0000000947,0.0000000985,0.0000000902,0.0000001003,0.0000000959,0.0000000934,0.0000000925,0.0000000985,0.0000000996,0.0000000989,0.0000000988,0.0000001008,0.0000001043,0.0000000976,0.0000000986,0.0000001059,0.0000000998,0.0000001027,0.0000001003,0.0000001015,0.0000001079,0.0000001055,0.0000001032,0.0000001049,0.0000001068,0.0000001058,0.0000000972,0.0000000996,0.0000001101,0.0000001039,0.0000001006,0.0000001119,0.0000001027,0.0000001050,0.0000001108,0.0000001074,0.0000001111,0.0000001101,0.0000001109,0.0000001105,0.0000001061,0.0000001080,0.0000001097,0.0000001170,0.0000001058,0.0000001139,0.0000001068,0.0000001082,0.0000001079,0.0000001050,0.0000001060,0.0000001060,0.0000001103,0.0000001126,0.0000001158,0.0000001086,0.0000001087,0.0000001110,0.0000001116,0.0000001044,0.0000001165,0.0000001111,0.0000001124,0.0000001163,0.0000001114,0.0000001111,0.0000001099,0.0000001081,0.0000001149,0.0000001090,0.0000001102,0.0000001074,0.0000001150,0.0000001078,0.0000001127,0.0000001143,0.0000001106,0.0000001202,0.0000001132,0.0000001152,0.0000000690,0.0000000753,0.0000000772,0.0000000815,0.0000000822,0.0000000759,0.0000000761,0.0000000837,0.0000000793,0.0000000858,0.0000000856,0.0000000831,0.0000000904,0.0000000838,0.0000000891,0.0000000867,0.0000000863,0.0000000918,0.0000000911,0.0000000910,0.0000000904,0.0000000934,0.0000000857,0.0000000953,0.0000000912,0.0000000883,0.0000000882,0.0000000941,0.0000000948,0.0000000938,0.0000000943,0.0000000956,0.0000000982,0.0000000927,0.0000000945,0.0000001010,0.0000000939,0.0000000988,0.0000000941,0.0000000972,0.0000001039,0.0000001002,0.0000000982,0.0000000993,0.0000001019,0.0000001004,0.0000000911,0.0000000949,0.0000001044,0.0000000993,0.0000000963,0.0000001055,0.0000000966,0.0000000996,0.0000001057,0.0000001014,0.0000001064,0.0000001043,0.0000001056,0.0000001057,0.0000001014,0.0000001014,0.0000001043,0.0000001100,0.0000001013,0.0000001087,0.0000001025,0.0000001032,0.0000001023,0.0000001007,0.0000001008,0.0000001006,0.0000001048,0.0000001069,0.0000001095,0.0000001032,0.0000001015,0.0000001063,0.0000001069,0.0000000993,0.0000001094,0.0000001055,0.0000001076,0.0000001101,0.0000001063,0.0000001056,0.0000001039,0.0000001029,0.0000001096,0.0000001031,0.0000001059,0.0000001029,0.0000001086,0.0000001030,0.0000001071,0.0000001082,0.0000001050,0.0000001137,0.0000001077,0.0000001093,0.0000000659,0.0000000709,0.0000000733,0.0000000783,0.0000000781,0.0000000710,0.0000000728,0.0000000797,0.0000000747,0.0000000818,0.0000000826,0.0000000787,0.0000000849,0.0000000792,0.0000000847,0.0000000815,0.0000000816,0.0000000861,0.0000000865,0.0000000872,0.0000000867,0.0000000897,0.0000000815,0.0000000913,0.0000000861,0.0000000834,0.0000000824,0.0000000893,0.0000000908,0.0000000898,0.0000000896,0.0000000916,0.0000000937,0.0000000865,0.0000000899,0.0000000964,0.0000000890,0.0000000937,0.0000000878,0.0000000921,0.0000000996,0.0000000951,0.0000000922,0.0000000928,0.0000000969,0.0000000959,0.0000000857,0.0000000902,0.0000000998,0.0000000946,0.0000000907,0.0000001005,0.0000000928,0.0000000958,0.0000001001,0.0000000951,0.0000001010,0.0000000994,0.0000001004,0.0000001006,0.0000000960,0.0000000968,0.0000001003,0.0000001052,0.0000000963,0.0000001038,0.0000000970,0.0000000987,0.0000000976,0.0000000944,0.0000000967,0.0000000945,0.0000000994,0.0000001009,0.0000001046,0.0000000987,0.0000000970,0.0000001008,0.0000001016,0.0000000944,0.0000001039,0.0000001009,0.0000001027,0.0000001052,0.0000001022,0.0000001004,0.0000000975,0.0000000976,0.0000001037,0.0000000980,0.0000001013,0.0000000992,0.0000001040,0.0000000967,0.0000001018,0.0000001027,0.0000000992,0.0000001084,0.0000001021,0.0000001043,0.0000000630,0.0000000673,0.0000000684,0.0000000741,0.0000000754,0.0000000684,0.0000000683,0.0000000752,0.0000000712,0.0000000780,0.0000000783,0.0000000734,0.0000000805,0.0000000752,0.0000000799,0.0000000777,0.0000000774,0.0000000826,0.0000000822,0.0000000837,0.0000000817,0.0000000851,0.0000000777,0.0000000874,0.0000000834,0.0000000800,0.0000000780,0.0000000842,0.0000000876,0.0000000847,0.0000000851,0.0000000872,0.0000000879,0.0000000815,0.0000000853,0.0000000922,0.0000000850,0.0000000889,0.0000000838,0.0000000878,0.0000000941,0.0000000897,0.0000000882,0.0000000890,0.0000000926,0.0000000917,0.0000000823,0.0000000854,0.0000000950,0.0000000896,0.0000000860,0.0000000964,0.0000000884,0.0000000899,0.0000000946,0.0000000898,0.0000000959,0.0000000950,0.0000000945,0.0000000963,0.0000000913,0.0000000917,0.0000000955,0.0000001002,0.0000000914,0.0000000992,0.0000000923,0.0000000936,0.0000000937,0.0000000895,0.0000000931,0.0000000894,0.0000000953,0.0000000954,0.0000000992,0.0000000947,0.0000000931,0.0000000956,0.0000000967,0.0000000894,0.0000000979,0.0000000966,0.0000000978,0.0000001006,0.0000000967,0.0000000948,0.0000000904,0.0000000925,0.0000000983,0.0000000933,0.0000000952,0.0000000947,0.0000000985,0.0000000917,0.0000000974,0.0000000982,0.0000000939,0.0000001051,0.0000000979,0.0000000983,0.0000000594,0.0000000634,0.0000000654,0.0000000705,0.0000000719,0.0000000644,0.0000000651,0.0000000715,0.0000000677,0.0000000747,0.0000000752,0.0000000703,0.0000000779,0.0000000715,0.0000000756,0.0000000737,0.0000000738,0.0000000795,0.0000000784,0.0000000801,0.0000000776,0.0000000810,0.0000000751,0.0000000840,0.0000000801,0.0000000767,0.0000000734,0.0000000800,0.0000000826,0.0000000818,0.0000000817,0.0000000833,0.0000000839,0.0000000773,0.0000000820,0.0000000873,0.0000000811,0.0000000843,0.0000000803,0.0000000835,0.0000000895,0.0000000862,0.0000000842,0.0000000840,0.0000000876,0.0000000875,0.0000000782,0.0000000803,0.0000000914,0.0000000853,0.0000000830,0.0000000918,0.0000000834,0.0000000850,0.0000000902,0.0000000848,0.0000000924,0.0000000909,0.0000000900,0.0000000910,0.0000000855,0.0000000868,0.0000000915,0.0000000951,0.0000000874,0.0000000942,0.0000000871,0.0000000891,0.0000000892,0.0000000841,0.0000000885,0.0000000850,0.0000000902,0.0000000915,0.0000000941,0.0000000904,0.0000000881,0.0000000910,0.0000000912,0.0000000861,0.0000000938,0.0000000922,0.0000000923,0.0000000958,0.0000000923,0.0000000901,0.0000000858,0.0000000880,0.0000000926,0.0000000894,0.0000000903,0.0000000894,0.0000000928,0.0000000871,0.0000000923,0.0000000940,0.0000000888,0.0000000988,0.0000000931,0.0000000936,0.0000000564,0.0000000604,0.0000000623,0.0000000668,0.0000000682,0.0000000610,0.0000000622,0.0000000673,0.0000000641,0.0000000708,0.0000000708,0.0000000672,0.0000000744,0.0000000681,0.0000000714,0.0000000695,0.0000000704,0.0000000750,0.0000000728,0.0000000766,0.0000000741,0.0000000761,0.0000000706,0.0000000804,0.0000000764,0.0000000732,0.0000000694,0.0000000761,0.0000000787,0.0000000774,0.0000000754,0.0000000793,0.0000000795,0.0000000736,0.0000000772,0.0000000834,0.0000000763,0.0000000810,0.0000000756,0.0000000785,0.0000000842,0.0000000828,0.0000000805,0.0000000800,0.0000000826,0.0000000827,0.0000000749,0.0000000755,0.0000000880,0.0000000816,0.0000000774,0.0000000866,0.0000000787,0.0000000810,0.0000000851,0.0000000811,0.0000000886,0.0000000864,0.0000000860,0.0000000858,0.0000000803,0.0000000822,0.0000000864,0.0000000891,0.0000000838,0.0000000890,0.0000000824,0.0000000843,0.0000000856,0.0000000799,0.0000000849,0.0000000804,0.0000000855,0.0000000880,0.0000000897,0.0000000865,0.0000000832,0.0000000863,0.0000000868,0.0000000815,0.0000000877,0.0000000870,0.0000000872,0.0000000913,0.0000000870,0.0000000857,0.0000000815,0.0000000833,0.0000000878,0.0000000853,0.0000000851,0.0000000854,0.0000000885,0.0000000817,0.0000000879,0.0000000896,0.0000000846,0.0000000942,0.0000000881,0.0000000886,0.0000000547,0.0000000574,0.0000000597,0.0000000628,0.0000000644,0.0000000574,0.0000000592,0.0000000638,0.0000000614,0.0000000666,0.0000000671,0.0000000636,0.0000000709,0.0000000629,0.0000000676,0.0000000672,0.0000000661,0.0000000711,0.0000000693,0.0000000725,0.0000000705,0.0000000715,0.0000000667,0.0000000761,0.0000000728,0.0000000706,0.0000000666,0.0000000725,0.0000000740,0.0000000735,0.0000000721,0.0000000751,0.0000000757,0.0000000705,0.0000000725,0.0000000799,0.0000000726,0.0000000763,0.0000000719,0.0000000742,0.0000000805,0.0000000787,0.0000000774,0.0000000758,0.0000000793,0.0000000786,0.0000000719,0.0000000714,0.0000000841,0.0000000778,0.0000000734,0.0000000823,0.0000000753,0.0000000774,0.0000000802,0.0000000777,0.0000000843,0.0000000816,0.0000000820,0.0000000822,0.0000000766,0.0000000777,0.0000000824,0.0000000844,0.0000000799,0.0000000854,0.0000000787,0.0000000798,0.0000000812,0.0000000755,0.0000000815,0.0000000765,0.0000000798,0.0000000840,0.0000000855,0.0000000820,0.0000000788,0.0000000813,0.0000000829,0.0000000782,0.0000000833,0.0000000820,0.0000000826,0.0000000876,0.0000000825,0.0000000810,0.0000000781,0.0000000796,0.0000000839,0.0000000816,0.0000000819,0.0000000804,0.0000000833,0.0000000781,0.0000000836,0.0000000844,0.0000000805,0.0000000884,0.0000000836,0.0000000846,0.0000000528,0.0000000542,0.0000000568,0.0000000598,0.0000000608,0.0000000555,0.0000000567,0.0000000604,0.0000000575,0.0000000626,0.0000000635,0.0000000603,0.0000000681,0.0000000604,0.0000000644,0.0000000645,0.0000000629,0.0000000674,0.0000000657,0.0000000698,0.0000000669,0.0000000685,0.0000000630,0.0000000727,0.0000000703,0.0000000669,0.0000000640,0.0000000689,0.0000000710,0.0000000698,0.0000000681,0.0000000715,0.0000000714,0.0000000663,0.0000000694,0.0000000745,0.0000000691,0.0000000723,0.0000000688,0.0000000706,0.0000000751,0.0000000750,0.0000000746,0.0000000728,0.0000000755,0.0000000751,0.0000000687,0.0000000684,0.0000000806,0.0000000735,0.0000000687,0.0000000788,0.0000000718,0.0000000747,0.0000000765,0.0000000745,0.0000000803,0.0000000775,0.0000000789,0.0000000779,0.0000000730,0.0000000733,0.0000000788,0.0000000804,0.0000000760,0.0000000811,0.0000000755,0.0000000764,0.0000000782,0.0000000724,0.0000000781,0.0000000725,0.0000000760,0.0000000785,0.0000000819,0.0000000778,0.0000000761,0.0000000775,0.0000000782,0.0000000745,0.0000000794,0.0000000780,0.0000000780,0.0000000822,0.0000000795,0.0000000766,0.0000000739,0.0000000765,0.0000000797,0.0000000754,0.0000000780,0.0000000765,0.0000000792,0.0000000751,0.0000000792,0.0000000810,0.0000000777,0.0000000842,0.0000000789,0.0000000809,0.0000000501,0.0000000519,0.0000000532,0.0000000563,0.0000000583,0.0000000530,0.0000000538,0.0000000580,0.0000000551,0.0000000599,0.0000000596,0.0000000565,0.0000000645,0.0000000569,0.0000000618,0.0000000612,0.0000000605,0.0000000637,0.0000000625,0.0000000673,0.0000000637,0.0000000643,0.0000000598,0.0000000690,0.0000000662,0.0000000637,0.0000000609,0.0000000649,0.0000000672,0.0000000670,0.0000000643,0.0000000677,0.0000000679,0.0000000619,0.0000000656,0.0000000696,0.0000000662,0.0000000684,0.0000000652,0.0000000677,0.0000000718,0.0000000715,0.0000000717,0.0000000696,0.0000000721,0.0000000704,0.0000000643,0.0000000649,0.0000000764,0.0000000690,0.0000000636,0.0000000748,0.0000000678,0.0000000702,0.0000000729,0.0000000716,0.0000000761,0.0000000736,0.0000000744,0.0000000743,0.0000000686,0.0000000693,0.0000000760,0.0000000766,0.0000000724,0.0000000781,0.0000000709,0.0000000728,0.0000000744,0.0000000689,0.0000000749,0.0000000690,0.0000000714,0.0000000743,0.0000000774,0.0000000745,0.0000000732,0.0000000737,0.0000000745,0.0000000709,0.0000000757,0.0000000746,0.0000000731,0.0000000775,0.0000000765,0.0000000723,0.0000000692,0.0000000729,0.0000000759,0.0000000716,0.0000000748,0.0000000727,0.0000000747,0.0000000709,0.0000000761,0.0000000780,0.0000000733,0.0000000805,0.0000000758,0.0000000771,0.0000000473,0.0000000496,0.0000000509,0.0000000537,0.0000000550,0.0000000507,0.0000000517,0.0000000557,0.0000000525,0.0000000559,0.0000000561,0.0000000539,0.0000000607,0.0000000540,0.0000000587,0.0000000582,0.0000000577,0.0000000610,0.0000000596,0.0000000643,0.0000000607,0.0000000624,0.0000000569,0.0000000665,0.0000000622,0.0000000602,0.0000000586,0.0000000619,0.0000000646,0.0000000644,0.0000000617,0.0000000647,0.0000000633,0.0000000599,0.0000000627,0.0000000666,0.0000000637,0.0000000647,0.0000000617,0.0000000643,0.0000000683,0.0000000683,0.0000000684,0.0000000667,0.0000000683,0.0000000680,0.0000000605,0.0000000610,0.0000000725,0.0000000654,0.0000000597,0.0000000711,0.0000000645,0.0000000668,0.0000000702,0.0000000684,0.0000000720,0.0000000707,0.0000000716,0.0000000704,0.0000000654,0.0000000654,0.0000000730,0.0000000726,0.0000000688,0.0000000743,0.0000000676,0.0000000699,0.0000000716,0.0000000658,0.0000000712,0.0000000653,0.0000000680,0.0000000713,0.0000000733,0.0000000709,0.0000000685,0.0000000690,0.0000000709,0.0000000676,0.0000000707,0.0000000710,0.0000000695,0.0000000737,0.0000000717,0.0000000693,0.0000000648,0.0000000695,0.0000000729,0.0000000682,0.0000000720,0.0000000688,0.0000000697,0.0000000683,0.0000000719,0.0000000739,0.0000000699,0.0000000777,0.0000000717,0.0000000734,0.0000000451,0.0000000466,0.0000000486,0.0000000513,0.0000000528,0.0000000477,0.0000000494,0.0000000534,0.0000000498,0.0000000534,0.0000000530,0.0000000520,0.0000000575,0.0000000511,0.0000000562,0.0000000553,0.0000000550,0.0000000578,0.0000000566,0.0000000613,0.0000000575,0.0000000589,0.0000000543,0.0000000636,0.0000000593,0.0000000566,0.0000000551,0.0000000583,0.0000000609,0.0000000611,0.0000000581,0.0000000616,0.0000000609,0.0000000570,0.0000000593,0.0000000632,0.0000000606,0.0000000617,0.0000000590,0.0000000605,0.0000000657,0.0000000651,0.0000000653,0.0000000633,0.0000000645,0.0000000644,0.0000000580,0.0000000574,0.0000000697,0.0000000621,0.0000000575,0.0000000672,0.0000000619,0.0000000640,0.0000000667,0.0000000646,0.0000000682,0.0000000664,0.0000000670,0.0000000671,0.0000000622,0.0000000628,0.0000000690,0.0000000689,0.0000000651,0.0000000708,0.0000000648,0.0000000662,0.0000000680,0.0000000628,0.0000000675,0.0000000629,0.0000000651,0.0000000675,0.0000000693,0.0000000673,0.0000000652,0.0000000664,0.0000000670,0.0000000644,0.0000000668,0.0000000667,0.0000000665,0.0000000705,0.0000000682,0.0000000668,0.0000000608,0.0000000654,0.0000000688,0.0000000649,0.0000000677,0.0000000650,0.0000000667,0.0000000647,0.0000000691,0.0000000701,0.0000000663,0.0000000730,0.0000000686,0.0000000699,0.0000000429,0.0000000445,0.0000000460,0.0000000490,0.0000000501,0.0000000460,0.0000000466,0.0000000504,0.0000000476,0.0000000500,0.0000000509,0.0000000488,0.0000000544,0.0000000493,0.0000000532,0.0000000527,0.0000000525,0.0000000542,0.0000000532,0.0000000574,0.0000000553,0.0000000554,0.0000000520,0.0000000599,0.0000000565,0.0000000531,0.0000000529,0.0000000557,0.0000000581,0.0000000589,0.0000000552,0.0000000578,0.0000000576,0.0000000539,0.0000000568,0.0000000592,0.0000000577,0.0000000580,0.0000000555,0.0000000583,0.0000000633,0.0000000614,0.0000000622,0.0000000593,0.0000000614,0.0000000608,0.0000000553,0.0000000542,0.0000000665,0.0000000588,0.0000000545,0.0000000627,0.0000000583,0.0000000598,0.0000000625,0.0000000625,0.0000000648,0.0000000627,0.0000000633,0.0000000640,0.0000000584,0.0000000593,0.0000000657,0.0000000658,0.0000000613,0.0000000670,0.0000000621,0.0000000633,0.0000000642,0.0000000591,0.0000000652,0.0000000598,0.0000000610,0.0000000641,0.0000000660,0.0000000642,0.0000000623,0.0000000626,0.0000000633,0.0000000616,0.0000000633,0.0000000633,0.0000000637,0.0000000659,0.0000000643,0.0000000638,0.0000000577,0.0000000617,0.0000000657,0.0000000607,0.0000000640,0.0000000613,0.0000000634,0.0000000617,0.0000000646,0.0000000664,0.0000000621,0.0000000683,0.0000000649,0.0000000647,0.0000000410,0.0000000415,0.0000000442,0.0000000463,0.0000000471,0.0000000439,0.0000000451,0.0000000486,0.0000000454,0.0000000476,0.0000000478,0.0000000463,0.0000000520,0.0000000464,0.0000000500,0.0000000502,0.0000000492,0.0000000517,0.0000000504,0.0000000553,0.0000000527,0.0000000527,0.0000000497,0.0000000569,0.0000000539,0.0000000508,0.0000000501,0.0000000520,0.0000000557,0.0000000569,0.0000000526,0.0000000548,0.0000000541,0.0000000503,0.0000000534,0.0000000570,0.0000000547,0.0000000547,0.0000000531,0.0000000559,0.0000000603,0.0000000576,0.0000000578,0.0000000551,0.0000000590,0.0000000578,0.0000000523,0.0000000518,0.0000000636,0.0000000562,0.0000000514,0.0000000588,0.0000000557,0.0000000552,0.0000000600,0.0000000591,0.0000000619,0.0000000589,0.0000000601,0.0000000609,0.0000000558,0.0000000569,0.0000000634,0.0000000615,0.0000000582,0.0000000642,0.0000000596,0.0000000612,0.0000000608,0.0000000560,0.0000000622,0.0000000571,0.0000000591,0.0000000599,0.0000000622,0.0000000610,0.0000000591,0.0000000603,0.0000000607,0.0000000586,0.0000000608,0.0000000598,0.0000000602,0.0000000632,0.0000000610,0.0000000601,0.0000000545,0.0000000580,0.0000000622,0.0000000580,0.0000000602,0.0000000589,0.0000000601,0.0000000587,0.0000000622,0.0000000638,0.0000000588,0.0000000659,0.0000000615,0.0000000611,0.0000000394,0.0000000387,0.0000000422,0.0000000437,0.0000000452,0.0000000420,0.0000000425,0.0000000460,0.0000000438,0.0000000447,0.0000000453,0.0000000442,0.0000000498,0.0000000444,0.0000000484,0.0000000484,0.0000000463,0.0000000489,0.0000000489,0.0000000524,0.0000000500,0.0000000501,0.0000000462,0.0000000540,0.0000000511,0.0000000486,0.0000000474,0.0000000478,0.0000000535,0.0000000544,0.0000000496,0.0000000524,0.0000000520,0.0000000470,0.0000000509,0.0000000552,0.0000000513,0.0000000512,0.0000000508,0.0000000529,0.0000000566,0.0000000547,0.0000000557,0.0000000526,0.0000000548,0.0000000547,0.0000000502,0.0000000492,0.0000000599,0.0000000537,0.0000000486,0.0000000558,0.0000000532,0.0000000532,0.0000000567,0.0000000564,0.0000000594,0.0000000562,0.0000000571,0.0000000579,0.0000000522,0.0000000541,0.0000000599,0.0000000589,0.0000000553,0.0000000606,0.0000000574,0.0000000577,0.0000000584,0.0000000533,0.0000000589,0.0000000550,0.0000000561,0.0000000570,0.0000000592,0.0000000572,0.0000000568,0.0000000570,0.0000000572,0.0000000555,0.0000000573,0.0000000560,0.0000000578,0.0000000598,0.0000000579,0.0000000575,0.0000000521,0.0000000537,0.0000000600,0.0000000555,0.0000000572,0.0000000563,0.0000000567,0.0000000562,0.0000000591,0.0000000596,0.0000000555,0.0000000628,0.0000000581,0.0000000584,0.0000000375,0.0000000370,0.0000000398,0.0000000409,0.0000000431,0.0000000398,0.0000000396,0.0000000436,0.0000000409,0.0000000425,0.0000000430,0.0000000421,0.0000000473,0.0000000426,0.0000000446,0.0000000467,0.0000000443,0.0000000463,0.0000000468,0.0000000498,0.0000000478,0.0000000479,0.0000000434,0.0000000515,0.0000000488,0.0000000460,0.0000000452,0.0000000450,0.0000000505,0.0000000518,0.0000000471,0.0000000494,0.0000000497,0.0000000439,0.0000000482,0.0000000530,0.0000000485,0.0000000491,0.0000000483,0.0000000501,0.0000000536,0.0000000521,0.0000000518,0.0000000504,0.0000000522,0.0000000523,0.0000000474,0.0000000465,0.0000000567,0.0000000516,0.0000000456,0.0000000535,0.0000000505,0.0000000505,0.0000000535,0.0000000532,0.0000000556,0.0000000534,0.0000000546,0.0000000546,0.0000000494,0.0000000516,0.0000000561,0.0000000566,0.0000000520,0.0000000576,0.0000000546,0.0000000547,0.0000000559,0.0000000510,0.0000000562,0.0000000523,0.0000000526,0.0000000548,0.0000000555,0.0000000547,0.0000000531,0.0000000545,0.0000000547,0.0000000523,0.0000000539,0.0000000532,0.0000000555,0.0000000568,0.0000000552,0.0000000540,0.0000000499,0.0000000517,0.0000000577,0.0000000530,0.0000000538,0.0000000532,0.0000000545,0.0000000534,0.0000000557,0.0000000568,0.0000000528,0.0000000595,0.0000000551,0.0000000565,0.0000000350,0.0000000349,0.0000000375,0.0000000392,0.0000000413,0.0000000375,0.0000000376,0.0000000420,0.0000000390,0.0000000399,0.0000000413,0.0000000402,0.0000000448,0.0000000402,0.0000000417,0.0000000439,0.0000000421,0.0000000436,0.0000000441,0.0000000479,0.0000000456,0.0000000449,0.0000000419,0.0000000484,0.0000000462,0.0000000445,0.0000000433,0.0000000436,0.0000000474,0.0000000495,0.0000000455,0.0000000472,0.0000000464,0.0000000410,0.0000000461,0.0000000508,0.0000000467,0.0000000474,0.0000000463,0.0000000480,0.0000000514,0.0000000485,0.0000000495,0.0000000467,0.0000000492,0.0000000501,0.0000000447,0.0000000439,0.0000000543,0.0000000493,0.0000000431,0.0000000500,0.0000000482,0.0000000482,0.0000000510,0.0000000509,0.0000000521,0.0000000507,0.0000000524,0.0000000522,0.0000000469,0.0000000490,0.0000000529,0.0000000541,0.0000000494,0.0000000545,0.0000000514,0.0000000515,0.0000000534,0.0000000488,0.0000000529,0.0000000492,0.0000000506,0.0000000516,0.0000000529,0.0000000525,0.0000000512,0.0000000517,0.0000000519,0.0000000500,0.0000000516,0.0000000501,0.0000000533,0.0000000536,0.0000000526,0.0000000507,0.0000000474,0.0000000491,0.0000000544,0.0000000507,0.0000000516,0.0000000496,0.0000000517,0.0000000506,0.0000000526,0.0000000535,0.0000000497,0.0000000560,0.0000000520,0.0000000536,0.0000000329,0.0000000332,0.0000000349,0.0000000376,0.0000000390,0.0000000354,0.0000000356,0.0000000389,0.0000000365,0.0000000380,0.0000000393,0.0000000383,0.0000000423,0.0000000386,0.0000000396,0.0000000420,0.0000000407,0.0000000416,0.0000000420,0.0000000459,0.0000000436,0.0000000427,0.0000000397,0.0000000459,0.0000000439,0.0000000422,0.0000000408,0.0000000415,0.0000000449,0.0000000464,0.0000000431,0.0000000454,0.0000000443,0.0000000387,0.0000000433,0.0000000481,0.0000000442,0.0000000453,0.0000000440,0.0000000454,0.0000000485,0.0000000457,0.0000000460,0.0000000447,0.0000000458,0.0000000480,0.0000000415,0.0000000416,0.0000000517,0.0000000469,0.0000000414,0.0000000482,0.0000000460,0.0000000455,0.0000000484,0.0000000487,0.0000000501,0.0000000484,0.0000000502,0.0000000492,0.0000000443,0.0000000476,0.0000000503,0.0000000519,0.0000000469,0.0000000524,0.0000000488,0.0000000488,0.0000000505,0.0000000469,0.0000000509,0.0000000478,0.0000000482,0.0000000491,0.0000000499,0.0000000502,0.0000000490,0.0000000494,0.0000000494,0.0000000477,0.0000000489,0.0000000474,0.0000000506,0.0000000500,0.0000000496,0.0000000486,0.0000000456,0.0000000470,0.0000000525,0.0000000478,0.0000000486,0.0000000468,0.0000000495,0.0000000479,0.0000000496,0.0000000517,0.0000000474,0.0000000523,0.0000000500,0.0000000512,0.0000000312,0.0000000320,0.0000000333,0.0000000363,0.0000000369,0.0000000338,0.0000000342,0.0000000370,0.0000000343,0.0000000365,0.0000000370,0.0000000359,0.0000000408,0.0000000365,0.0000000372,0.0000000399,0.0000000391,0.0000000390,0.0000000400,0.0000000433,0.0000000415,0.0000000411,0.0000000374,0.0000000425,0.0000000415,0.0000000399,0.0000000386,0.0000000387,0.0000000426,0.0000000434,0.0000000416,0.0000000420,0.0000000425,0.0000000371,0.0000000411,0.0000000452,0.0000000422,0.0000000435,0.0000000421,0.0000000432,0.0000000462,0.0000000435,0.0000000440,0.0000000427,0.0000000431,0.0000000454,0.0000000392,0.0000000395,0.0000000496,0.0000000440,0.0000000383,0.0000000459,0.0000000439,0.0000000434,0.0000000461,0.0000000466,0.0000000472,0.0000000458,0.0000000472,0.0000000465,0.0000000419,0.0000000455,0.0000000487,0.0000000490,0.0000000450,0.0000000494,0.0000000463,0.0000000459,0.0000000486,0.0000000447,0.0000000481,0.0000000455,0.0000000450,0.0000000469,0.0000000480,0.0000000467,0.0000000465,0.0000000471,0.0000000472,0.0000000458,0.0000000471,0.0000000453,0.0000000478,0.0000000472,0.0000000467,0.0000000465,0.0000000433,0.0000000450,0.0000000501,0.0000000449,0.0000000454,0.0000000444,0.0000000464,0.0000000459,0.0000000476,0.0000000484,0.0000000453,0.0000000496,0.0000000470,0.0000000481,0.0000000291,0.0000000309,0.0000000318,0.0000000338,0.0000000341,0.0000000321,0.0000000329,0.0000000347,0.0000000324,0.0000000351,0.0000000349,0.0000000345,0.0000000380,0.0000000352,0.0000000358,0.0000000383,0.0000000370,0.0000000377,0.0000000383,0.0000000418,0.0000000395,0.0000000394,0.0000000356,0.0000000405,0.0000000396,0.0000000383,0.0000000376,0.0000000366,0.0000000411,0.0000000408,0.0000000386,0.0000000399,0.0000000402,0.0000000361,0.0000000385,0.0000000426,0.0000000404,0.0000000422,0.0000000397,0.0000000406,0.0000000440,0.0000000414,0.0000000415,0.0000000405,0.0000000409,0.0000000428,0.0000000371,0.0000000370,0.0000000476,0.0000000416,0.0000000365,0.0000000436,0.0000000421,0.0000000416,0.0000000440,0.0000000448,0.0000000452,0.0000000432,0.0000000447,0.0000000446,0.0000000393,0.0000000432,0.0000000463,0.0000000474,0.0000000428,0.0000000457,0.0000000439,0.0000000438,0.0000000456,0.0000000423,0.0000000455,0.0000000433,0.0000000432,0.0000000438,0.0000000446,0.0000000435,0.0000000442,0.0000000452,0.0000000444,0.0000000431,0.0000000451,0.0000000432,0.0000000455,0.0000000452,0.0000000436,0.0000000439,0.0000000408,0.0000000430,0.0000000473,0.0000000433,0.0000000439,0.0000000420,0.0000000443,0.0000000438,0.0000000449,0.0000000456,0.0000000434,0.0000000458,0.0000000445,0.0000000457,0.0000000277,0.0000000292,0.0000000309,0.0000000322,0.0000000321,0.0000000308,0.0000000313,0.0000000331,0.0000000313,0.0000000335,0.0000000328,0.0000000331,0.0000000361,0.0000000331,0.0000000343,0.0000000363,0.0000000344,0.0000000357,0.0000000366,0.0000000401,0.0000000373,0.0000000380,0.0000000334,0.0000000388,0.0000000376,0.0000000365,0.0000000362,0.0000000348,0.0000000387,0.0000000393,0.0000000365,0.0000000383,0.0000000375,0.0000000345,0.0000000368,0.0000000412,0.0000000375,0.0000000403,0.0000000370,0.0000000391,0.0000000425,0.0000000397,0.0000000398,0.0000000384,0.0000000385,0.0000000403,0.0000000350,0.0000000350,0.0000000444,0.0000000394,0.0000000351,0.0000000411,0.0000000398,0.0000000389,0.0000000418,0.0000000417,0.0000000423,0.0000000402,0.0000000429,0.0000000422,0.0000000369,0.0000000413,0.0000000441,0.0000000452,0.0000000405,0.0000000434,0.0000000420,0.0000000425,0.0000000432,0.0000000398,0.0000000430,0.0000000411,0.0000000412,0.0000000417,0.0000000426,0.0000000419,0.0000000424,0.0000000432,0.0000000422,0.0000000412,0.0000000428,0.0000000418,0.0000000435,0.0000000424,0.0000000416,0.0000000427,0.0000000389,0.0000000411,0.0000000446,0.0000000402,0.0000000417,0.0000000403,0.0000000433,0.0000000410,0.0000000426,0.0000000435,0.0000000413,0.0000000429,0.0000000416,0.0000000437,0.0000000262,0.0000000272,0.0000000298,0.0000000308,0.0000000300,0.0000000286,0.0000000300,0.0000000310,0.0000000296,0.0000000316,0.0000000309,0.0000000313,0.0000000351,0.0000000314,0.0000000323,0.0000000344,0.0000000328,0.0000000339,0.0000000347,0.0000000382,0.0000000356,0.0000000362,0.0000000325,0.0000000373,0.0000000358,0.0000000341,0.0000000343,0.0000000329,0.0000000367,0.0000000376,0.0000000347,0.0000000360,0.0000000351,0.0000000331,0.0000000342,0.0000000396,0.0000000359,0.0000000387,0.0000000352,0.0000000366,0.0000000405,0.0000000378,0.0000000381,0.0000000358,0.0000000362,0.0000000384,0.0000000335,0.0000000329,0.0000000428,0.0000000372,0.0000000333,0.0000000393,0.0000000378,0.0000000370,0.0000000397,0.0000000397,0.0000000404,0.0000000384,0.0000000400,0.0000000402,0.0000000354,0.0000000384,0.0000000406,0.0000000435,0.0000000384,0.0000000415,0.0000000397,0.0000000406,0.0000000411,0.0000000382,0.0000000411,0.0000000396,0.0000000390,0.0000000393,0.0000000406,0.0000000402,0.0000000399,0.0000000412,0.0000000403,0.0000000386,0.0000000415,0.0000000399,0.0000000401,0.0000000403,0.0000000390,0.0000000403,0.0000000369,0.0000000387,0.0000000430,0.0000000384,0.0000000399,0.0000000385,0.0000000413,0.0000000391,0.0000000402,0.0000000417,0.0000000395,0.0000000404,0.0000000393,0.0000000409,0.0000000253,0.0000000254,0.0000000289,0.0000000292,0.0000000285,0.0000000266,0.0000000278,0.0000000296,0.0000000288,0.0000000302,0.0000000299,0.0000000298,0.0000000338,0.0000000300,0.0000000305,0.0000000329,0.0000000309,0.0000000321,0.0000000331,0.0000000355,0.0000000332,0.0000000347,0.0000000311,0.0000000357,0.0000000338,0.0000000322,0.0000000322,0.0000000311,0.0000000344,0.0000000354,0.0000000329,0.0000000335,0.0000000338,0.0000000314,0.0000000324,0.0000000378,0.0000000344,0.0000000360,0.0000000337,0.0000000352,0.0000000386,0.0000000361,0.0000000364,0.0000000345,0.0000000340,0.0000000370,0.0000000323,0.0000000307,0.0000000405,0.0000000354,0.0000000312,0.0000000384,0.0000000352,0.0000000345,0.0000000372,0.0000000377,0.0000000378,0.0000000366,0.0000000384,0.0000000377,0.0000000335,0.0000000365,0.0000000382,0.0000000412,0.0000000369,0.0000000391,0.0000000377,0.0000000389,0.0000000394,0.0000000369,0.0000000390,0.0000000387,0.0000000372,0.0000000381,0.0000000387,0.0000000380,0.0000000379,0.0000000397,0.0000000389,0.0000000367,0.0000000402,0.0000000382,0.0000000385,0.0000000385,0.0000000367,0.0000000380,0.0000000353,0.0000000372,0.0000000412,0.0000000362,0.0000000383,0.0000000365,0.0000000389,0.0000000377,0.0000000376,0.0000000394,0.0000000384,0.0000000386,0.0000000378,0.0000000381,0.0000000237,0.0000000241,0.0000000275,0.0000000279,0.0000000269,0.0000000253,0.0000000269,0.0000000281,0.0000000275,0.0000000286,0.0000000284,0.0000000279,0.0000000316,0.0000000288,0.0000000292,0.0000000309,0.0000000292,0.0000000304,0.0000000312,0.0000000336,0.0000000313,0.0000000331,0.0000000288,0.0000000340,0.0000000320,0.0000000299,0.0000000306,0.0000000288,0.0000000331,0.0000000333,0.0000000315,0.0000000320,0.0000000320,0.0000000300,0.0000000311,0.0000000354,0.0000000320,0.0000000343,0.0000000324,0.0000000335,0.0000000371,0.0000000341,0.0000000351,0.0000000325,0.0000000323,0.0000000353,0.0000000312,0.0000000296,0.0000000382,0.0000000339,0.0000000292,0.0000000368,0.0000000338,0.0000000326,0.0000000351,0.0000000362,0.0000000358,0.0000000352,0.0000000360,0.0000000353,0.0000000316,0.0000000338,0.0000000368,0.0000000397,0.0000000352,0.0000000376,0.0000000362,0.0000000366,0.0000000370,0.0000000351,0.0000000373,0.0000000368,0.0000000349,0.0000000360,0.0000000370,0.0000000353,0.0000000357,0.0000000380,0.0000000368,0.0000000341,0.0000000382,0.0000000363,0.0000000368,0.0000000369,0.0000000348,0.0000000359,0.0000000337,0.0000000356,0.0000000391,0.0000000345,0.0000000370,0.0000000351,0.0000000373,0.0000000353,0.0000000355,0.0000000374,0.0000000372,0.0000000367,0.0000000368,0.0000000367,0.0000000226,0.0000000226,0.0000000263,0.0000000266,0.0000000256,0.0000000239,0.0000000253,0.0000000270,0.0000000259,0.0000000271,0.0000000266,0.0000000268,0.0000000299,0.0000000273,0.0000000278,0.0000000294,0.0000000284,0.0000000292,0.0000000295,0.0000000316,0.0000000293,0.0000000317,0.0000000274,0.0000000324,0.0000000309,0.0000000273,0.0000000289,0.0000000269,0.0000000313,0.0000000317,0.0000000293,0.0000000308,0.0000000306,0.0000000291,0.0000000300,0.0000000331,0.0000000308,0.0000000322,0.0000000311,0.0000000321,0.0000000355,0.0000000327,0.0000000334,0.0000000308,0.0000000306,0.0000000338,0.0000000293,0.0000000287,0.0000000357,0.0000000319,0.0000000281,0.0000000350,0.0000000324,0.0000000307,0.0000000325,0.0000000346,0.0000000342,0.0000000337,0.0000000341,0.0000000345,0.0000000303,0.0000000326,0.0000000344,0.0000000377,0.0000000332,0.0000000351,0.0000000345,0.0000000353,0.0000000353,0.0000000329,0.0000000356,0.0000000353,0.0000000334,0.0000000344,0.0000000347,0.0000000334,0.0000000343,0.0000000364,0.0000000347,0.0000000327,0.0000000363,0.0000000343,0.0000000348,0.0000000352,0.0000000329,0.0000000340,0.0000000315,0.0000000332,0.0000000369,0.0000000323,0.0000000353,0.0000000335,0.0000000354,0.0000000334,0.0000000336,0.0000000351,0.0000000355,0.0000000353,0.0000000349,0.0000000342,0.0000000212,0.0000000217,0.0000000254,0.0000000254,0.0000000244,0.0000000224,0.0000000237,0.0000000255,0.0000000249,0.0000000260,0.0000000253,0.0000000250,0.0000000285,0.0000000264,0.0000000270,0.0000000279,0.0000000268,0.0000000276,0.0000000281,0.0000000296,0.0000000283,0.0000000301,0.0000000260,0.0000000310,0.0000000296,0.0000000260,0.0000000266,0.0000000256,0.0000000297,0.0000000302,0.0000000276,0.0000000289,0.0000000282,0.0000000275,0.0000000288,0.0000000317,0.0000000296,0.0000000309,0.0000000292,0.0000000306,0.0000000341,0.0000000310,0.0000000328,0.0000000286,0.0000000288,0.0000000321,0.0000000277,0.0000000269,0.0000000339,0.0000000306,0.0000000267,0.0000000333,0.0000000310,0.0000000299,0.0000000310,0.0000000330,0.0000000326,0.0000000320,0.0000000321,0.0000000330,0.0000000293,0.0000000313,0.0000000327,0.0000000358,0.0000000311,0.0000000330,0.0000000327,0.0000000330,0.0000000330,0.0000000310,0.0000000336,0.0000000337,0.0000000316,0.0000000327,0.0000000331,0.0000000319,0.0000000330,0.0000000350,0.0000000323,0.0000000313,0.0000000344,0.0000000324,0.0000000335,0.0000000334,0.0000000311,0.0000000310,0.0000000296,0.0000000310,0.0000000350,0.0000000309,0.0000000327,0.0000000320,0.0000000339,0.0000000322,0.0000000322,0.0000000333,0.0000000339,0.0000000328,0.0000000336,0.0000000323,0.0000000196,0.0000000201,0.0000000236,0.0000000241,0.0000000230,0.0000000213,0.0000000227,0.0000000242,0.0000000239,0.0000000244,0.0000000235,0.0000000236,0.0000000276,0.0000000254,0.0000000257,0.0000000260,0.0000000259,0.0000000261,0.0000000263,0.0000000284,0.0000000269,0.0000000286,0.0000000245,0.0000000296,0.0000000285,0.0000000247,0.0000000254,0.0000000246,0.0000000280,0.0000000289,0.0000000269,0.0000000269,0.0000000265,0.0000000257,0.0000000269,0.0000000298,0.0000000281,0.0000000287,0.0000000270,0.0000000292,0.0000000330,0.0000000290,0.0000000307,0.0000000276,0.0000000278,0.0000000302,0.0000000260,0.0000000256,0.0000000329,0.0000000295,0.0000000258,0.0000000321,0.0000000296,0.0000000284,0.0000000292,0.0000000309,0.0000000313,0.0000000306,0.0000000300,0.0000000318,0.0000000278,0.0000000298,0.0000000311,0.0000000348,0.0000000297,0.0000000316,0.0000000308,0.0000000306,0.0000000317,0.0000000297,0.0000000323,0.0000000321,0.0000000299,0.0000000311,0.0000000311,0.0000000303,0.0000000306,0.0000000336,0.0000000306,0.0000000294,0.0000000333,0.0000000303,0.0000000319,0.0000000318,0.0000000295,0.0000000294,0.0000000281,0.0000000295,0.0000000337,0.0000000295,0.0000000309,0.0000000307,0.0000000316,0.0000000311,0.0000000306,0.0000000312,0.0000000321,0.0000000312,0.0000000318,0.0000000294,0.0000000181,0.0000000193,0.0000000226,0.0000000227,0.0000000222,0.0000000207,0.0000000214,0.0000000224,0.0000000225,0.0000000234,0.0000000226,0.0000000218,0.0000000262,0.0000000244,0.0000000245,0.0000000250,0.0000000246,0.0000000254,0.0000000251,0.0000000273,0.0000000257,0.0000000271,0.0000000231,0.0000000273,0.0000000267,0.0000000234,0.0000000249,0.0000000236,0.0000000260,0.0000000279,0.0000000255,0.0000000259,0.0000000252,0.0000000247,0.0000000260,0.0000000282,0.0000000266,0.0000000275,0.0000000257,0.0000000276,0.0000000317,0.0000000279,0.0000000288,0.0000000257,0.0000000260,0.0000000288,0.0000000252,0.0000000240,0.0000000306,0.0000000281,0.0000000241,0.0000000306,0.0000000283,0.0000000272,0.0000000284,0.0000000291,0.0000000297,0.0000000285,0.0000000285,0.0000000305,0.0000000264,0.0000000284,0.0000000293,0.0000000328,0.0000000279,0.0000000296,0.0000000299,0.0000000295,0.0000000304,0.0000000275,0.0000000308,0.0000000307,0.0000000283,0.0000000296,0.0000000303,0.0000000292,0.0000000289,0.0000000319,0.0000000290,0.0000000276,0.0000000318,0.0000000285,0.0000000301,0.0000000299,0.0000000277,0.0000000282,0.0000000269,0.0000000280,0.0000000319,0.0000000284,0.0000000298,0.0000000292,0.0000000296,0.0000000301,0.0000000292,0.0000000302,0.0000000310,0.0000000302,0.0000000298,0.0000000283,0.0000000173,0.0000000182,0.0000000208,0.0000000213,0.0000000210,0.0000000200,0.0000000206,0.0000000216,0.0000000214,0.0000000220,0.0000000215,0.0000000211,0.0000000248,0.0000000233,0.0000000233,0.0000000236,0.0000000222,0.0000000241,0.0000000242,0.0000000250,0.0000000243,0.0000000256,0.0000000213,0.0000000254,0.0000000254,0.0000000225,0.0000000246,0.0000000225,0.0000000247,0.0000000270,0.0000000237,0.0000000246,0.0000000239,0.0000000237,0.0000000250,0.0000000267,0.0000000254,0.0000000265,0.0000000242,0.0000000261,0.0000000299,0.0000000261,0.0000000275,0.0000000241,0.0000000244,0.0000000273,0.0000000237,0.0000000230,0.0000000285,0.0000000265,0.0000000223,0.0000000293,0.0000000268,0.0000000266,0.0000000266,0.0000000278,0.0000000281,0.0000000272,0.0000000272,0.0000000293,0.0000000250,0.0000000275,0.0000000278,0.0000000311,0.0000000268,0.0000000276,0.0000000285,0.0000000279,0.0000000282,0.0000000261,0.0000000296,0.0000000283,0.0000000271,0.0000000279,0.0000000279,0.0000000275,0.0000000278,0.0000000300,0.0000000282,0.0000000267,0.0000000302,0.0000000266,0.0000000288,0.0000000281,0.0000000270,0.0000000258,0.0000000246,0.0000000270,0.0000000294,0.0000000266,0.0000000286,0.0000000284,0.0000000283,0.0000000286,0.0000000275,0.0000000284,0.0000000297,0.0000000291,0.0000000286,0.0000000265,0.0000000164,0.0000000172,0.0000000193,0.0000000199,0.0000000198,0.0000000189,0.0000000194,0.0000000204,0.0000000205,0.0000000211,0.0000000206,0.0000000202,0.0000000235,0.0000000219,0.0000000220,0.0000000227,0.0000000211,0.0000000232,0.0000000234,0.0000000241,0.0000000229,0.0000000244,0.0000000198,0.0000000239,0.0000000240,0.0000000216,0.0000000238,0.0000000214,0.0000000233,0.0000000252,0.0000000228,0.0000000238,0.0000000226,0.0000000224,0.0000000236,0.0000000255,0.0000000246,0.0000000254,0.0000000223,0.0000000252,0.0000000283,0.0000000256,0.0000000258,0.0000000229,0.0000000236,0.0000000260,0.0000000230,0.0000000219,0.0000000269,0.0000000250,0.0000000208,0.0000000273,0.0000000259,0.0000000255,0.0000000255,0.0000000269,0.0000000267,0.0000000256,0.0000000261,0.0000000280,0.0000000239,0.0000000262,0.0000000264,0.0000000289,0.0000000250,0.0000000265,0.0000000271,0.0000000269,0.0000000265,0.0000000252,0.0000000279,0.0000000265,0.0000000257,0.0000000263,0.0000000264,0.0000000263,0.0000000264,0.0000000285,0.0000000267,0.0000000246,0.0000000289,0.0000000248,0.0000000280,0.0000000272,0.0000000253,0.0000000245,0.0000000234,0.0000000260,0.0000000273,0.0000000252,0.0000000275,0.0000000275,0.0000000275,0.0000000274,0.0000000255,0.0000000275,0.0000000283,0.0000000284,0.0000000271,0.0000000252,0.0000000155,0.0000000157,0.0000000183,0.0000000187,0.0000000193,0.0000000180,0.0000000186,0.0000000193,0.0000000198,0.0000000203,0.0000000195,0.0000000196,0.0000000223,0.0000000208,0.0000000214,0.0000000213,0.0000000197,0.0000000222,0.0000000223,0.0000000224,0.0000000220,0.0000000228,0.0000000189,0.0000000225,0.0000000233,0.0000000205,0.0000000220,0.0000000210,0.0000000215,0.0000000242,0.0000000216,0.0000000224,0.0000000218,0.0000000213,0.0000000222,0.0000000241,0.0000000230,0.0000000240,0.0000000210,0.0000000236,0.0000000268,0.0000000235,0.0000000247,0.0000000215,0.0000000225,0.0000000245,0.0000000216,0.0000000209,0.0000000258,0.0000000229,0.0000000197,0.0000000250,0.0000000246,0.0000000241,0.0000000242,0.0000000253,0.0000000253,0.0000000243,0.0000000256,0.0000000267,0.0000000226,0.0000000248,0.0000000252,0.0000000271,0.0000000243,0.0000000254,0.0000000262,0.0000000257,0.0000000250,0.0000000241,0.0000000272,0.0000000252,0.0000000245,0.0000000252,0.0000000246,0.0000000245,0.0000000255,0.0000000272,0.0000000250,0.0000000235,0.0000000278,0.0000000231,0.0000000267,0.0000000258,0.0000000241,0.0000000237,0.0000000219,0.0000000243,0.0000000262,0.0000000238,0.0000000256,0.0000000260,0.0000000261,0.0000000261,0.0000000240,0.0000000260,0.0000000266,0.0000000273,0.0000000259,0.0000000242};

mat Y1(Han_data1);
mat Y2(Han_data2);

Y1.reshape(331,49);
Y2.reshape(331,100);


//cout<<"Y1 is: "<<endl;
//cout<<Y1.n_rows<<endl;
//cout<<Y1.n_cols<<endl;

//cout<<"Y2 is: "<<endl;


//cout<<Y2.n_rows<<endl;
//cout<<Y2.n_cols<<endl;


if(nstudy>50)
 {
   nstudy=50;
 }


int test1=int(floor(Z_score*10.0));
int test2=int(ceil(Z_score*10.0));

double adjust1=Y1(test1,nstudy);

double adjust2=Y1(test2,nstudy);

double statisticRandomEffects2_ = Z_score* Z_score;

double p = 0.5 * gsl_cdf_chisq_Q (statisticRandomEffects2_, 1.0) + 0.5 * gsl_cdf_chisq_Q (statisticRandomEffects2_, 2.0);

//cout<<"p is: "<<p<<endl;

double p_final =0.0;

if (test2 < 331) {
      int m = nstudy - 2;
      double d4 = 0.0;
      d4 = adjust1;
//      cout<<"I am here1"<<endl;
      double d5 = 0.5 * gsl_cdf_chisq_Q(test1 / 10.0, 1.0) + 0.5 * gsl_cdf_chisq_Q(test1 / 10.0, 2.0);
      double d6 = d4 / d5;
//      cout<<"I am here2"<<endl;
      double d7 = adjust2;
      double d8 = 0.5 * gsl_cdf_chisq_Q(test2 / 10.0, 1.0) + 0.5 * gsl_cdf_chisq_Q(test2 / 10.0, 2.0);
//      cout<<"I am here3"<<endl;
      double d9 = d7 / d8;
      double d10 = d6 + (d9 - d6) * (test2 - test1 / 10.0) / 0.1;
//      cout<<"I am here4"<<endl;
      p_final= d10 * p;
//      cout<<"p_final (1) here is: "<<p_final<<endl;
    } else
    {
    int k = nstudy - 2;
//    cout<<"I am here5"<<endl;
    double d1 = Y1(330,k);

    double d2 = 0.5 * gsl_cdf_chisq_Q(33.0, 1.0) + 0.5 * gsl_cdf_chisq_Q(33.0, 2.0);
    double d3 = d1 / d2;
    if(d3*p>1)
      {
        p_final= 1.0;
      } else
      {
        p_final= d3*p;
      }
//     cout<<"p_final (2) here is: "<<p_final<<endl;
    }
double D2=0.0;
int b1=1;
//cout<<"I am here6"<<endl;
for(int i=0; i<Cov1.n_rows; i++)
  {
  //  cout<<"I am here7"<<endl;
    for(int j=i; j<Cov1.n_rows; j++)
     {
       b1++;
       if(j==i)
        {
          D2+=1;
        } else
        {
          D2+=Cov1(i,j);
        }
     }
  }

//cout<<"I am here8"<<endl;

double D3 = D2 / b1;
    int j = int(floor(D3 * 100.0));
    if (j >= 100)
      j = 99; 
    if(j<0)
     {
       j=0;
     }
    double d4 = Y2(nstudy-2,j) / Y2(nstudy-2,0);
//    cout<<"d4 is: "<<d4<<endl;
    p_final=p_final*d4;
   if(p_final>1.0)
    {
      p_final=1.0;
    }
   return p_final;
   
}
