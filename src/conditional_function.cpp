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
      ss>>test;
      Y.push_back(test);
      ss>>test1;
      Sd.push_back(test1);
      ss>>test2;
      ss>>test2;
      Beta.push_back(test);
      
      Y.push_back(test2);
      test=test/test1;
      
      if(test<0)
       {
         test=-test;
       }

      Y.push_back(test);
      
      out_gemma.push_back(Y);
      Z_score.push_back(test);
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
    cout<<"Come into extractvariant function to extract genotype for tartet variants"<<endl; 
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



















double ConstructCovRandom(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > > &beta_,  vector<vector<vector<double> > > &sd_, int index, mat &Res_fixed, mat &Res_random, string meta_mode, string output, vector<vector<string> >  &variant_, map<string, bool> &pos, map<string, bool> &pos1, map <string, vector<double>> &Res, map<string, double > &MAF, vector<double> &VARIANCE, map<int,bool> &data_index, vector<double> &z_score, vector<double> &beta, vector<double> &sd, map<string, string> &ALLELE,string cor_file, bool write_cor)
 {

   mat cov_z_score = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_beta    = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_z_score_copy = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_sd      = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);
   mat cov_sd_copy      = mat(z_score_[0][index].size(),z_score_.size(), fill::zeros);

//   cout<<"cov_z_score_copy is: "<<cov_z_score_copy<<endl;

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
         if(abs(cov_z_score(i,j))>3.3 || abs(cov_z_score(i,j))<=0.00000001)
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
                    Cov(X,Y)=0.99*cov_sd(i,X)*cov_sd(i,Y);
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
















double MetaAnalysis(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > >&beta_, vector<vector<vector<double> > >&sd_,int index, mat &Res_fixed, mat &Res_random, string meta_mode, string output, vector<vector<string> > &variant_, map<string, bool> &pos, int colocalization, map<string, string> &GWAS_file, map<string, double > &MAF, vector<double > &VARIANCE, map<int,bool> &data_index,string gene, bool INDEX, bool finemap, map<string, string> & ALLELE, double p_cutoff,string cor_file, bool write_cor)
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


   double P=ConstructCovRandom(z_score_, beta_, sd_ ,index , Res_fixed, Res_random, meta_mode, output1, variant_, pos, pos1, Res, MAF, VARIANCE, data_index, test1, test2,test3, ALLELE, cor_file, write_cor);
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
















double MetaAnalysis(vector<vector<vector<double> > > &z_score_, vector<vector<vector<double> > >&beta_, vector<vector<vector<double> > >&sd_,int index, mat &Res_fixed, mat &Res_random, string meta_mode, string output, vector<vector<string> > &variant_, map<string, bool> &pos, int colocalization, map<string, string> &GWAS_file, map<string, double > &MAF, vector<double > &VARIANCE, map<int,bool> &data_index,string gene,vector<double> &test1_1, vector<double> &test2_1, vector<double> &test3_1, map<string, string> &ALLELE, string cor_file, bool write_cor)
 {
//   cout<<"Come into MetaAnalysis function to perform final meta-analysis"<<endl;
   string output1= output+"_statistical_signal";
   map <string, bool> pos1;
   map <string, vector<double>> Res;
   pos1=pos;
//   cout<<"Perform real meta-analysis"<<endl;
   
   double P=ConstructCovRandom(z_score_, beta_, sd_ ,index , Res_fixed, Res_random, meta_mode, output1, variant_, pos, pos1, Res, MAF, VARIANCE, data_index,test1_1,test2_1, test3_1, ALLELE, cor_file, write_cor);
   
   
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







int MM_conditional_analysis(vector<string> &file_, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,vector<string> &covariate_, vector <string> &grm_file_,string outputFile,vector<string> &geno_file_, string meta_mode, vector<string> &allele_file_, map<string, bool> &pos, int colocalization, map<string, string> &GWAS_file, map<string, double> &MAF,map<int,bool> &data_index, int nested, bool finemap, double p_cutoff, bool cis_summary, string cor_file, vector<string> &causal)
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
            MetaAnalysis(z_score_nested,  beta_nested, sd_nested, x, Res_fixed, Res_random, "fixed", Output1, variant_nested, pos, 0,GWAS_file, MAF, VARIANCE, data_index1, gene, test1_1, test2_1, test3_1, ALLELE, test_file,false);
           
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
            P=MetaAnalysis(z_score_nested2,  beta_nested2, sd_nested2, x, Res_fixed_2, Res_random_2, "random", Output2, variant_nested2, pos, 0,GWAS_file, MAF, VARIANCE2, data_index2, gene, INDEX, finemap, ALLELE,p_cutoff,cor_file,false);
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
          P1=MetaAnalysis(z_score_,  beta_, sd_, x, Res_fixed, Res_random, meta_mode, Output, variant_, pos, colocalization,GWAS_file, MAF, VARIANCE, data_index, gene, INDEX, finemap, ALLELE,p_cutoff,cor_file, write_cor);
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
                string PEAK=variant_[indx][int(Res_fixed_2(0,0))];
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
