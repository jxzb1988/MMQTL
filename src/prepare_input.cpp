
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <map>
#include <bitset>
using namespace std;
using namespace arma;



int Read_phe(string input, vector<vector<string> > &Ind,vector<vector<double> > &Phe, string gene,map<int,bool> &data_index);


int Readbim(string file)
 {
   int ibuf=0;
   string A;
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

int Readfam(string file)
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


int Read_geno(string geno, vector<string> &geno_,vector<vector<string> > &Ind,vector<vector<double> > &Phe ,int chr, int start, int end, map<string,bool> &pos, string gene, map<int, bool> &data_index, int colocalization, bool finemap, int cis_length, int cov_length,vector<string> &cov, int LD_index);


      int READInputFile(string input, vector<string> &output)
         {
           ifstream fin(input.c_str(), std::ifstream::in);
           vector<double> words;
           string line;
           while(fin && getline(fin, line))
             {
               stringstream ss(line);
               string Str;
               ss >> Str;
               output.push_back(Str);
             }
           return(0);
         }


      int READInputFile(string input, vector<string> &output, map<int, bool> &data_index)
         {
           cout<<"Come into function to extract genotype"<<endl;
           ifstream fin(input.c_str(), std::ifstream::in);
           vector<double> words;
           string line;
           int b=-1;
/*	   for(map<int,bool>::iterator it=data_index.begin();it!=data_index.end();it++)
//            {
//              if(it->second)
               {
                 cout<<"dataset "<<it->first<<" is available"<<endl;
               }
            }
*/
           while(fin && getline(fin, line))
             {  
               b++;
	       map<int,bool>::iterator Ite;
               Ite=data_index.find(b);
             
               if(Ite!=data_index.end())
                {
                  if(!Ite->second)
                   {
                     output.push_back(string("NA"));
                     continue;
                   }
                }
               stringstream ss(line);
               string Str;
               ss >> Str;
//	       cout<<"File is: "<<Str<<endl;
               output.push_back(Str);
             }
           return(0);
         } 


       int READInputFile(string input, vector<string> &output, vector<string> &allele_file)
         {
//           int index = 0;
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




void Readbim(string bimfile, map<int, vector<int> > &order,vector<string> &variant)
     {
       int ibuf=-1;
       string A;
//       string cbuf="0";
//       double dbuf=0.0;
       string str_buf;
       ifstream Bim(bimfile.c_str());
       if(!Bim) throw("Error: can not open the file ["+bimfile+"] to read.");
       while(Bim && getline(Bim, A))
        {

          vector<int> pos;
          ibuf++;
          stringstream ss(A);
          int waste;
          string var;
          ss>>waste;
          pos.push_back(waste);
          ss>>var;ss>>waste;
          ss>>waste;
          variant.push_back(var);
          pos.push_back(waste);
//          cout<<"One variant stored"<<endl;
          order.insert(pair<int,vector<int>>(ibuf,pos));
        }
   Bim.close();
     }
int Readbim(string bimfile, string output_file)
     {
       string A;
       string str_buf;
       ifstream Bim(bimfile.c_str());
       fstream  FHOU(output_file.c_str(),ios::out);
       if(!FHOU)
        {
       //   cout<<"Error for input or output"<<endl;
          return (1);
        }
       if(!Bim) throw("Error: can not open the file ["+bimfile+"] to read.");
       while(Bim && getline(Bim, A))
        {
          stringstream ss(A);
          string waste;
          string var;
          ss>>waste;
          ss>>var;ss>>waste;
          ss>>waste;
          ss>>waste;
          FHOU<<var<<"\t"<<waste<<endl;
        }
       Bim.clear();
       Bim.close();
       return 0;
     }
void Read_anno(string gene, string anno, int &chr, int &start, int &end)
         {
        //   cout<<"Come into read_anno function"<<endl;
           ifstream FHIN1(anno.c_str(),ifstream::in);
           string a;
//           int b=0;
           bool index=false;
           while(FHIN1 &&getline(FHIN1,a))
            {

                  stringstream ss(a);
                  string word;
                  int Chr;
                  int End;
                  int Start;
                  ss>>Chr;
                  ss>>Start;
                  ss>>End;
                  ss>>word;
          //        cout<<"gene is: "<<gene<<", and word is: "<<word<<endl;
                  if(gene == word)
                   {
                     chr=Chr;
                     start=Start;
                     end = End;
                     index=true;
                     break;
                   }
              }
           if(!index)
            {
             // cout<<"There is no annotation infor for "<<gene<<", please have a check, or change to a new gene ID "<<endl;
              exit(EXIT_FAILURE);
            }
           FHIN1.clear();
           FHIN1.close();
           
         }




void Match_ind(map<string, int> &ind_geno, vector<string> &ind_phe, vector<int> &index, map<string,int> &index2, vector<bool> &INDinc )
             {
               cout<<"Come to function Match_ind"<<endl;
               for(int i=0;i<ind_phe.size();i++)
                {
                  cout<<"Check for phenotype individual "<<ind_phe[i]<<endl;
                  map<string,int>::iterator Ite;
                  Ite=ind_geno.find(ind_phe[i]);
                  if(Ite!=ind_geno.end())
                   {
                     cout<<"Yes, genotype for the explored individuals is available"<<endl;
                     index.push_back(Ite->second);
                     index2.insert(pair<string, int> (ind_phe[i],i));
                     INDinc[Ite->second]=true;
                   }
                }
             }



bool Match_cis(vector<int> &Order, int chr, int start, int end, int dis)
               {
              //   cout<<"Order0 is: "<<Order[0]<<", Order1 is: "<<Order[1]<<", and Order2 is: "<<Order[2]<<endl;
//                 cout<<"Chr is: "<<chr<<", start is: "<<start<<", and end is: "<<end<<endl;
                 if(Order[0]==chr && (Order[1]>=start-dis && Order[1]<=end+dis))
                  {
                    return true;
                  } else
                  {
                    return false;
                  }
               }








void Readfam(string file,map<string, int> &ind, vector<string> &individuals)
     {
       ifstream Fam(file.c_str());

       string A;
       int ibuf=-1;
       if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
//       cout<<"Reading PLINK FAM file from ["+file+"]."<<endl;
       while(Fam && getline(Fam, A))
        {
          ibuf++;
          stringstream ss(A);
          string TT;
          ss>>TT;
          ss>>TT;
          individuals.push_back(TT);
          ind.insert(pair<string, int>(TT, ibuf));
       }
      Fam.clear();
      Fam.close();
    }




//int main(string input,vector<string> &input_, string geno_file, vector<string> &geno_file_, string gene, string anno)
//int prepare_input(int argc, char *argv[])


//int Prepare_input(string input, string geno_file, string gene, string anno, vector<string> &geno_file_, vector<string> &allele_file_, vector<string> &y, map<string, bool> & pos, map<int,bool> &data_index, int colocalization, bool finemap)


int Prepare_input(string input, string geno_file, string gene, string anno, vector<string> &geno_file_, vector<string> &allele_file_, vector<string> &y, map<string, bool> &pos, map<int, bool> &data_index, int colocalization, bool finemap,int cis_length, int cov_length, vector<string> &cov, int LD_index)
         {
//           cout<<"Come into Prepare_input function"<<endl;
  //         cout<<"input is: "<<input<<endl;
    //       cout<<"geno_file is: "<<geno_file<<endl;
      //     cout<<"annotation file is: "<<anno<<endl;
        //   cout<<"gene is: "<<gene<<endl;
//           if(argc<=0)
//            {
//              return 1;
//            }
           int chr=-9;
           int start=-9;
           int end=-9;
//           string input="data/phe_file";
//           string geno_file="data/geno_file";
  //         vector<string> input_;
  //         vector<string> geno_file_;
//           string gene="ENSG00000070010";
//           string anno="data/gene_annotation_chr22";
           vector<vector<string>> Ind;
           vector<vector<double>> Phe;
           if(gene!="")
            {
               Read_anno(gene, anno, chr, start, end);
            } else
            {
//              cout<<"Please provide the gene ID"<<endl;
              exit (EXIT_FAILURE);
            }
           if(chr==-9 || start==-9 || end==-9)
            {
  //            cout<<"NO position infor dound, EXIT"<<endl;
              exit(EXIT_FAILURE);
            }
   //        start=(start-3000000)>0?(start-3000000):0;
   //        end  =end  +3000000;
//           cout<<"chr is: "<<chr<<", start is: "<<start<<", and end is: "<<end<<endl; 
  //         cout<<"Come to extract phenotype"<<endl;
//           map<int, bool> data_index;         
           int state=Read_phe(input, Ind, Phe, gene, data_index);
           if(state!=0)
            {
              return -1;
            }
    //       cout<<"Extracting phenotype is over"<<endl;    
      //     cout<<"Get input file name"<<endl;
//           cout<<"extract phenotype database number is: "<<Ind.size()<<endl;
//           cout<<"Availability of phenotype infor "<<endl;
/*           
           for(map<int,bool>::iterator it=data_index.begin();it!=data_index.end();it++)
            {
              if(it->second)
               {
                 cout<<"dataset "<<it->first<<" is available"<<endl;
               } else
               {
                 cout<<"dataset "<<it->first<<" is not available"<<endl;
               }
            }
*/
           Read_geno(geno_file, geno_file_, Ind, Phe, chr, start, end, pos, gene, data_index,  colocalization,finemap,cis_length, cov_length, cov, LD_index);
           for(int i=0;i<data_index.size();i++)
//           for(int i=0;i<Ind.size();i++) 
            {
              if(data_index[i])
               {
                 y.push_back(gene+string("/merged_genotype_phenotype_tissue_")+to_string(i));
                 allele_file_.push_back(gene+string("/affected_allele_tissue_")+to_string(i));
               } else
               {
                 y.push_back(string("NA"));
                 allele_file_.push_back(string("NA"));
               }
          //     cout<<"Ouptut geno+phenotye and allele file name"<<endl;

      //        cout<<"Output individuals for data "<<i<<endl;
      //        for(unsigned int j=0;j<Ind[i].size();j++)
      //         {
      //           cout<<" "<<Ind[i][j];
      //         }
      //        cout<<endl;
            }
//           cout<<"Start to extract genotype"<<endl;
   //        return 0;
//           Read_geno(geno_file, geno_file_, Ind, Phe, chr, start, end);
           return 0;
         }


void outputcov(string cov,vector<bool> &INDinc)
 {
   ifstream FHIN1(cov.c_str(),ifstream::in);

    string a;
    double val;
    string cov_adjust=cov+string("_adjust");

    fstream  FHOU(cov_adjust.c_str(),ios::out);
              if(!FHOU)
               {
                 cout<<"Error for input or output"<<endl;
                 return;
    //             return (1);
               }


    int b=-1;
    cout<<"There are "<<INDinc.size()<<" individuals to extract"<<endl;
    cout<<"The input covariate file is: "<<cov<<endl;
    cout<<"The output covariate file is: "<<cov_adjust<<endl;
    for(int i=1;i<INDinc.size();i++)
     {
       cout<<"Useful individual "<<INDinc[i]<<endl;
     }
    while(FHIN1 &&getline(FHIN1,a))
     {
       b++;
       if(INDinc[b])
        {
          stringstream ss(a);
          ss>>val;
          FHOU<<val;
          while(ss && ss>>val)
           {
             //Changed by Biao Zeng, 10/26/2020, 12:54 PM.
             val+=0.001*randn();
             FHOU<<"\t"<<val;
           }
          FHOU<<endl;
        }
      }

    FHIN1.clear();
    FHIN1.close();
    FHOU.clear();
    FHOU.close();
 }


int Read_geno(string geno, vector<string> &geno_,vector<vector<string> > &Ind,vector<vector<double> > &Phe ,int chr, int start, int end, map<string,bool > &pos, string gene, map<int, bool> &data_index, int colocalization, bool finemap, int  cis_length, int cov_length, vector<string> &cov, int LD_index)
       {
            READInputFile(geno, geno_, data_index);
            string a;
            for(int I=0,indx=-1;I<geno_.size();I++)
             {
                if(geno_[I] =="NA")
                 {
                   continue;
                 }
                indx++;
//                cout<<"Prepare data for dataset"<<I<<endl;
                int i=0, j=0, k=0;
                string input=geno_[I];
//                cout<<"Genotype file name is: "<<input<<endl;
                string strfam=".fam";
                string strbim=".bim";
                string famfile=input+strfam;
                string bimfile=input+strbim;
                int nsnp=Readbim(bimfile);
                int nind=Readfam(famfile);
                string output_allele=gene+string("/affected_allele_tissue_")+to_string(I);
                cout<<"bimfile is: "<<bimfile<<endl;
                Readbim(bimfile,output_allele);
//                cout<<"Preparation for allele file is over"<<endl;
                vector<int> index;
                map<int,vector<int> > order;
                vector<string> variant;
                vector<string> variant1;
//                cout<<"bim file is "<<bimfile<<endl;
                Readbim(bimfile,order, variant);
                int zeng=0;
                int biao=0;
                int cis_variant=0;
                vector<bool> SNPinc;
                vector<int> pos1;
                for(int Zeng=0;Zeng<nsnp;Zeng++)
                 {
//                   cout<<"Infor is: "<<order[Zeng][0]<<order[Zeng][1]<<order[Zeng][2]<<endl;
//                   cout<<"chr is: "<<chr<<endl;
//                   cout<<"start is: "<<start<<endl;
//                   cout<<"end is: "<<end<<endl;
                   if(Match_cis(order[Zeng],chr, start, end, cov_length)) 
                    {
                      SNPinc.push_back(true);
                      zeng++;
                    } else
                    {
                      SNPinc.push_back(false);
                    }
                   if(Match_cis(order[Zeng],chr, start, end, cis_length))
                    {
                      pos.insert(pair<string, bool>(variant[Zeng],true));
                      pos1.push_back(zeng-1);
                      variant1.push_back(variant[Zeng]);
                      cis_variant++;
                    } else
                    {
                      pos.insert(pair<string, bool>(variant[Zeng],false));
                    }
                 }

                if(cis_variant==0)
                 {
                   cout<<"There is no genetic variant available for analysis, please have a check on the input genotype file "<<I<<". EXIT"<<endl;
                   exit(EXIT_FAILURE); 
                 }

//                cout<<"Here1_1"<<endl;
                map<string, int> ind_geno;
                vector<string> individuals;
                Readfam(famfile,ind_geno, individuals);
                vector<int> ind_index;
                map<string, int> ind_index2;
                
                vector<bool> INDinc;
//                cout<<"Here1_2"<<endl;
                for(int BIAO=0;BIAO<nind;BIAO++)
                 {
                   INDinc.push_back(false);
                 }
                cout<<"Here3_1"<<endl;
                Match_ind(ind_geno, Ind[indx], ind_index,ind_index2, INDinc); 
//                Match_ind(ind_geno, Ind[I], ind_index,ind_index2, INDinc);
                cout<<"Number of matched individuals is: "<<INDinc.size();
                cout<<"Here4_1"<<endl;
                mat _snp_2 =mat(zeng,ind_index2.size(), fill::zeros);
                mat _snp_1 =mat(zeng,ind_index2.size(), fill::zeros);
                char Chr[1];
                input=input+".bed";
                cout<<"Genotype file is "<<input<<endl;
                fstream BIT(input.c_str(), ios::in|ios::binary);
                if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
                for(i=0; i<3; i++) BIT.read(Chr,1); // skip the first three bytes
                int snp_indx=0, indi_indx=0;
                int snp_indx_1=0;
                for(j=0, snp_indx=0; j<nsnp; j++)
                 {
                   if(SNPinc[j]==false)
                    {
                      for(int i=0; i<nind;i+=4)
                       {
                         char ch[1];
                         BIT.read(ch,1);
                       }
                      continue;
                    }
                   for(int i=0, indi_indx=0; i<nind;)
                    {
                      char ch[1];
                      bitset<8> b;
                      BIT.read(ch,1);
                      if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
                      b=ch[0];
                      k=0;
                      while(k < 7 && i < nind)
                       {
                         if(INDinc[i]==false)
                          {
                            k+=2;
                            i++;
                            continue;
                          }
                         _snp_2(snp_indx,indi_indx)=(!b[k++]);
                         _snp_1(snp_indx,indi_indx)=(!b[k++]);
                         indi_indx++;
                         i++;  
                       }
                    }
                  snp_indx++;
                }
               BIT.clear();
               BIT.close();
               cout<<"Get dosage"<<endl;
               mat dosage=mat(_snp_2.n_rows,_snp_2.n_cols, fill::zeros);

              if(cov.size() !=0 && cov[I] !="NA")
               {
                 cout<<"There is covariate file available, and extract covariate information"<<endl;
                 outputcov(cov[I],INDinc);
               }
               dosage = dosage + 0.0001*mat(_snp_2.n_rows,_snp_2.n_cols, fill::randn);
   
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

     if(I==LD_index &&  ( colocalization>0  || finemap==true))
//     if(I==-1)  
      {
//        cout<<"Come to calculate LD r2"<<endl;
        string output=gene+string("/ldFile");
        fstream  FHOU(output.c_str(),ios::out);
              if(!FHOU)
               {
                 cout<<"Error for input or output"<<endl;
                 return (1);
               }
/*
       for(int indx1=0;indx1<pos1.size();indx1++)
         {
           cout<<"Variant is "<<variant1[i]<<endl;
           cout<<"Genotype is: "<<endl;
           cout<<dosage.row(pos1[indx1])<<endl;

         }
*/
        FHOU<<"variant";
       for(int i=0;i<pos1.size();i++)
        {
             FHOU<<"\t"<<variant1[i];
        }
       FHOU<<endl;
       mat dosage1=mat(cis_variant,_snp_2.n_cols, fill::zeros);
       for(int i=0;i<cis_variant;i++)
        {
          for(int j=0;j<_snp_2.n_cols;j++)
           {
             dosage1(i,j)=dosage(pos1[i],j);;
           }
        }
//        cout<<"Output dosage1"<<endl;
//        cout<<dosage1<<endl;
//        cout<<"Output is over"<<endl;
        if(cis_variant<=19000)
         {
           mat R = mat(cis_variant,cis_variant, fill::zeros);
           R = cor(dosage1.t());
           R.replace(datum::nan, 0);
           for(int indx1=0;indx1<pos1.size();)
            {
              FHOU<<variant1[indx1];

              for(int indx2=0;indx2<pos1.size();)
               {
                 FHOU<<"\t"<<R(indx1,indx2); //changed by Biao Zeng 8:45 am, 04/05/2020
                 indx2++;
               }
              FHOU<<endl;
              indx1++;
            }
         } else
         {
        cout<<"Oh, no, too many variants in the region, calculate the LD one by one"<<endl;
        for(int indx1=0;indx1<pos1.size();)
         {
           FHOU<<variant1[indx1];
           for(int indx2=0;indx2<pos1.size();)
            {
//               cout<<"Genotype for variant "<<indx1<<" is: "<<endl;
//               cout<<dosage.row(pos1[indx1])<<endl;
//               cout<<"Genotype for variant "<<indx2<<" is: "<<endl;
//               cout<<dosage.row(pos1[indx2])<<endl;

                 mat r2=cor(dosage.row(pos1[indx1]),dosage.row(pos1[indx2]));
                 FHOU<<"\t"<<r2(0,0);
              indx2++;
            }
           FHOU<<endl;
           indx1++;
          }
        }

      
     }



//     cout<<"Output dosage and phenotype"<<endl;
              string output_tissue=gene+string("/merged_genotype_phenotype_tissue_")+to_string(I);
              fstream  FHOU(output_tissue.c_str(),ios::out);
              if(!FHOU)
               {
                 cout<<"Error for input or output"<<endl;
                 return (1);
               }
//              cout<<"Start"<<endl;
              FHOU<<"ind";
              FHOU<<"     "<<"ILphe";
              for(int i=0;i<SNPinc.size();i++)
               {
                 if(SNPinc[i]==true)
                  {
                    FHOU<<"	"<<variant[i];
                  }
               }
              FHOU<<endl;
              for(int i=0, indx1=0;i<INDinc.size();)
               {
                 if(INDinc[i]==false)
                  { 
                    i++;
                    continue;
                  }
  //               cout<<individuals[ind_index[indx1]];
//                 string Ind=individuals[ind_index[indx1]];
                 string Ind=individuals[i];
//                 FHOU<<individuals[ind_index2[Ind]];
                 FHOU<<Ind;
                 FHOU<<"	"<<Phe[indx][ind_index2[Ind]];
//                 cout<<"	"<<Phe[I][ind_index2[Ind]];
                 for(int j=0, indx2=0;j<SNPinc.size();)
                  {
  //                  cout<<"i is: "<<i<<endl;
    //                cout<<"j is: "<<j<<endl;
                    if(SNPinc[j]==false)
                     {
                       j++;
                       continue;
                     }
                    FHOU<<"	"<<dosage(indx2,indx1);
      //              cout<<"	"<<dosage(indx2,indx1);
                    indx2++;
                    j++;
                  }
                 FHOU<<endl;
                 indx1++;
                 i++;
               }
              FHOU.clear();
              FHOU.close();
  //            cout<<"Process for file "<<I<<" is over"<<endl;
            }
         return 0;
       }

int Read_phe(string input, vector<vector<string> > &Ind,vector<vector<double> > &Phe, string gene, map<int, bool> &data_index)
  {
    vector<string> input_;
    READInputFile(input, input_);
    string a;
    for(int i=0;i<input_.size();i++)
     {
       int b=0;
       string data=input_[i];
       ifstream FHIN1(data.c_str(),ifstream::in);
       if(!FHIN1)
        {
          exit(EXIT_FAILURE);
        }
       vector<string> ind;
               vector<double> phe;
            int found=0;
            while(FHIN1 &&getline(FHIN1,a))
             {
               b++;
               if(b==1)
                {
                  stringstream ss(a);
                  string word;
                  ss>>word;
                  while(ss && ss>>word)
                   {
                     ind.push_back(word);
                   }
                 } else
                 {
                   double c=0;
                   stringstream ss(a);
                   string waste;

                   ss>>waste;
                   if(waste == gene)
                    {
                      vector<double> T;
                      while(ss && ss>>c)
                       {
                         phe.push_back(c);
                       }
                      Phe.push_back(phe);
                      Ind.push_back(ind);
                      found=1;
                      break;
                    }

                  }
              }
                     
            FHIN1.clear();
            FHIN1.close();
            if(found==1)
             {
               data_index.insert(pair<int, bool> (i, true));
             } else
             {
               data_index.insert(pair<int, bool> (i, false));
             }
     }
   if(Phe.size()==0)
    {
      cout<<"No phenotype data available for explored gene "<<gene<<endl;
   //   exit(EXIT_FAILURE);
      return -1;
    }
   for(int i=0;i<Phe.size();i++)
     {
       double Mean=0;
       int num=0;
       for(int j=0;j<Phe[0].size();j++)
        {
          if(Phe[i][j]!=-9)
           {
             Mean+=Phe[i][j];
             num++;
           }
        }
       Mean/=num;
       for(int j=0;j<Phe[0].size();j++)
        {
          if(Phe[i][j]==-9)
           {
             Phe[i][j]=Mean;
           }
        }
     }
    return 0;
  }
