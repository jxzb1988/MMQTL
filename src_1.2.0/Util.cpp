
//The part of the function in the code from adapted from CAVIAR developped by Hormozdiari, F et al.



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <map>

using namespace std;

/*
void readbim(string bimfile, mpa<int, vector<int>> &order,vector<string> &variant)
     {
       int ibuf=-1;
       string A;
       string cbuf="0";
       double dbuf=0.0;
       string str_buf;
       ifstream Bim(bimfile.c_str());
       if(!Bim) throw("Error: can not open the file ["+file+"] to read.");
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
          order.insert(pair<string,vector<int>>(ibuf,pos));
        }
   Bim.close();
     }


void readbim(string bimfile, string output_file)
     {
       string A;
       string str_buf;
       ifstream Bim(bimfile.c_str());
       fstream  FHOU(output_file.c_str(),ios::out);
       if(!FHOU)
        {
          cout<<"Error for input or output"<<endl;
          return (1);
        }
       if(!Bim) throw("Error: can not open the file ["+file+"] to read.");
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
       Bim.close();
     }





void read_anno(string gene, string anno, int chr, int start, int end)
         {

           ifstream FHIN1(anno.c_str(),ifstream::in);
           string a;
           int b=0;
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
              cout<<"There is no annotation infor for "<<gene<<", please have a check, or change to a new gene ID "<<endl;
              exit(EXIT_FAILURE);
            }
           FHIN1.clear();
           FHIN1.close();
         }




void match_ind(ind_geno, vector<string> &ind_phe, vector<int> &index)
             {
               for(int i=0;i<ind_phe.size();i++)
                {
                  map<string,int>::iterator Ite;
                  Ite=ind_geno.find(ind_phe[i]);
                  if(Ite!=order.end())
                   {
                     index.push_back(Ite->second) ;
                   }
                }
             }



bool match_cis(vector<int> &Order, int chr, int start, int end)
               {
                 if(Order[0]==chr && (Order[1]>=start && Order[1]<=end))
                  {
                    return true;
                  } else
                  {
                    return false;
                  }
               }








void readfam(string file,map<string, int> &ind, vector<string> &individuals)
     {
       ifstream Fam(file.c_str());

       string A;
       int ibuf=-1;
       if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
       cout<<"Reading PLINK FAM file from ["+file+"]."<<endl;
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




void prepare_input(string input,vector<string> &input_, string geno_file, vector<string> &geno_file_, string gene, string anno)
         {
           int chr;
           int start;
           int end;
           vector<vector<string>> Ind;
           vector<vector<double>> Phe;
           if(gene!="")
            {
               read_anno(gene, anno, chr, start, end);
            } else
            {
              cout<<"Please provide the gene ID"<<endl;
              exit (EXIT_FAILURE);
            }
           read_phe(input, Ind, Phe, gene);

           read_geno(geno_file, geno_file_, Ind, Phe, chr, start, end);

         }

void read_geno(string geno, vector<string> &geno_,vector<vector<string>> &Ind,vector<vector<double>> &Phe ,int chr, int start, int end)
       {
            ReadInputFile(geno, geno_);
            string a;
            int b=0;
            for(int I=0;I<geno_.size();I++)
             {
                int i=0, j=0, k=0;
                string input=geno_[i];
                string strfam=".fam";
                string strbim=".bim";
                string famfile=input+strfam;
                string bimfile=input+strbim;
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
               BIT.clear();
               BIT.close();
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
              map<int,vector<int>> order;
              vector<string> variant;
              readbim(bimfile,order, variant);
              string output_allele="affected_allele_tissue_"+to_string(I);
              readbim(bimfile,output_allele);
              vector<int> index;
              for(int x=0;x<dosage.n_rows;x++)
               {
                 if(match_cis(order[x],chr, start, end))
                  {
                    index.push_back(x);
                  }
               }
              map<string, int> ind_geno;
              vector<string> individuals;
              readfam(famfile,ind_geno, individuals);
              vector<int> ind_index;

              match_ind(ind_geno, Ind[I], ind_index);
              if(ind_index.size()<=0)
               {
                 cout<<"There is no sample overlapped between genotype and phenotype, please have a check on the input files. Exit"<<endl;
                 exit(EXIT_FAILURE);
               }
              mat target = dosage.submat(index,ind_index);
              string output_tissue="merged_genotype_phenotype_tissue_"+to_string(I);
              fstream  FHOU(output_tissue.c_str(),ios::out);
              if(!FHOU)
               {
                 cout<<"Error for input or output"<<endl;
                 return (1);
               }
              FHOU<<"ind";
              FHOU<<"     "<<"ILphe";
              for(int i=0;i<index.size();i++)
               {
                 FHOU<<"  "<<variant[index[i]];
               }
              FHOU<<endl;
              for(int i=0;i<ind_index.size();i++)
               {
                 FHOU<<individuals[ind_index[i]];
                 FHOU<<Phe[I][ind_index[i]];
                 for(int j=0;j<index.size();j++)
                  {
                    FHOU<<"       "<<target(ind_index[i],index[j]);
                  }
                 FHOU<<endl;
               }
              FHOU.clear();
              FHOU.close();
            }
       }






void read_phe(string input, vector<vector<string>> &Ind,vector<vector<double>> &Phe, string gene)
  {
            vector<string> input_;
            ReadInputFile(input, input_);
            string a;
            int b=0;
            for(int i=0;i<input_.size();i++)
             {
               string data=input_[i];
               ifstream FHIN1(data.c_str(),ifstream::in);
               vector<string> ind;
               vector<double> phe;
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
                      break;
                    }

                  }
              }
            FHIN1.clear();
            FHIN1.close();
         }
  }





*/







long int fact(int n) {
        if(n==0)
                return 1;
        return n* fact(n-1);
}

void copyConfigure(double *dest, double *src, int size) {
	for(int i = 0; i < size; i++) 
		dest[i] = src[i];
}

double min(double a, double b) {
	if(a>b)
		return b;
	else
		return a;
}

long int nCr(int n, int r) {
        long int result = 1;
        for(int i = n; i > n-r; i--)
                result *= i;
        return result/fact(r);
}

void printVector(char * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%c, ", data[i]);
}

void printVector(int * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%d, ", (int)data[i]);
}

void printVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%lf, ", data[i]);
}

void diffVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
		result[i] = data1[i] - data2[i];
}

void sumVector(double * data1, double * data2, int size, double * result) {
        for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] + data2[i];
}

double multVector(double * data1, double * data2, int size) {
	double res = 0;
	for(int i = 0; i < size; i++ ) 
                res += data1[i] * data2[i];
	return res;
}

void dotVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] * data2[i];
}

void multVectorMatrix(double *vector, double * matrix, int size, double * result) {
	double total_row = 0;
	for(int i = 0; i < size; i++) {
		total_row = 0;
		for(int j = 0; j < size; j++) {
			total_row += vector[j] * matrix[i + j * size];
		}
		result[i]= total_row;
	}
}




void import_first_row(string data, map<string, int> &variant)
  {
    ifstream FHIN1(data.c_str(),ifstream::in);
    string a;
    int b=0;
    if(FHIN1 &&getline(FHIN1,a))
     {

          stringstream ss(a);
          string word;
          ss>>word;
          ss>>word;
          while(ss && ss>>word)
           {
//             variant.push_back(word);
//             variant.insert();
             variant.insert(pair<string,int>(word,1));
           }
      }

    FHIN1.clear();
    FHIN1.close();
  }


void import_first_col(string data, map<string, int> &variant)
  {
    ifstream FHIN1(data.c_str(),ifstream::in);
    string a;
    int b=0;
    while(FHIN1 &&getline(FHIN1,a))
     {

          stringstream ss(a);
          string word;
          ss>>word;
             variant.insert(pair<string,int>(word,1));
      }

    FHIN1.clear();
    FHIN1.close();
  }



/*
void importData(string fileName, double * vector) {
	int index = 0;
	double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        fin >> data;
        while (fin.good()) {
                vector[index] = data;
                index++;
                fin >> data;
        }
        fin.close();
}

void importData(string fileName, int * vector) {
        int index = 0;
        double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( fin.good()  ){
                fin >> data;
                vector[index] = (int)data;
                index++;
        }
        fin.close();
}

*/

void importData(string fileName, double * vector) {
        int index = 0;
        double data = 0;
        string dataS = "";
        std::stringstream dataSS;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        fin >> dataS;
        stringstream(dataS) >> data;
        while ( fin.good()) {
                vector[index] = data;
                index++;
                fin >> dataS;
                if (dataS.find("nan") == string::npos)
                        stringstream(dataS) >> data;
                else
                        data = 0;
        }
        fin.close();
}

void importData(string fileName, int * vector) {
        int index = 0;
        double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( fin.good()  ){
                fin >> data;
                vector[index] = (int)data;
                index++;
        }
        fin.close();
}








/*
	The column index starts by 1 in this implemenation
*/



void importDataSecondColumn(string fileName, double * vector) {
	int index = 0;
	string line = "";
	string dataS = "";
	double data = 0.0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( getline(fin, line) ){
		istringstream iss(line);
		iss >> dataS;
		iss >> data;
	        vector[index] = (double)data;
                index++;
        }
//	cout << "reach=" << index << endl;
        fin.close();
}

/*
	The column index starts by 1 in this implemenation
*/

/*
void importDataNthColumn(string fileName, double * vector, int colNum, int ignore=0) {
        int index = 0;
        string line = "";
        string dataS = "";
        double data = 0.0;
	ifstream fin(fileName.c_str(), std::ifstream::in);
        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> dataS;
                for(int i = 0; i < colNum-1;i++)
			iss >> data;
                vector[index] = (double)data;
                index++;
        }
 //       cout << "reach=" << index << endl;
        fin.close();
}

*/



void importDataNthColumn(string fileName, double * vector, int colNum, int ignore=0) {
        int index = 0;
        string line = "";
        string dataS = "";
        double data = 0.0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        for(int i = 0; i < ignore; i++)
                getline(fin, line);

        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> dataS;
                for(int i = 0; i < colNum-1;i++)
                        iss >> data;
                vector[index] = (double)data;
                index++;
        }
        fin.close();
}



void importDataFirstColumn(string fileName, string * list, int ignore=0) {
        int index = 0;
        string data = "";
        string line = "";
        ifstream fin(fileName.c_str(), std::ifstream::in);
        for(int i = 0; i < ignore; i++)
                getline(fin, line);

        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> data;
                list[index] = data;
                index++;
        }
        fin.close();
}

void importDataFirstColumn(string fileName, map<string,bool> &MAP, int ignore=0) {
        int index = 0;
        string data = "";
        string line = "";
        ifstream fin(fileName.c_str(), std::ifstream::in);
        for(int i = 0; i < ignore; i++)
                getline(fin, line);

        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> data;
                MAP.insert(pair<string, bool>(data,true));
//                list[index] = data;
//                index++;
        }
        fin.close();
}

void importDataFirstColumn(string fileName, double * list) {
   //     cout<<"Come into to extract variants names"<<endl;
 	int index = 0;
        double data = 0.0;
        string line = "";
	ifstream fin(fileName.c_str(), std::ifstream::in);
   //     cout<<"I am safe in extracting variants"<<endl;
        while( getline(fin, line) ){
   //             cout<<"Come into"<<endl;
   //             cout<<"line is: "<<line<<endl;
//		istringstream iss(line);
                istringstream iss;
   //             string test_b="Biao Zeng is the best";
          //      istringstream iss_test;
   //             cout<<"Tested string is: "<<test_b<<endl;
        //        iss_test.str(test_b);
   //             cout<<"It was changed"<<endl;
                iss.str(line);
   //             cout<<"Still safe 1"<<endl;
//                cout<<"Come into"<<endl;
                iss >> data;
   //             cout<<"Still safe 2"<<endl;
                list[index] = (double) data ;
   //             cout<<"Still safe 3"<<endl;
		index++;
          }
        double mean=0;
        for(int i=0;i<index;i++)
         {
           mean+=list[i];
         }
        mean/=index;
        for(int i=0;i<index;i++)
         {
           list[i]-=mean;
         }
     //   double mean=0;
     //   for(int i=0;i<index;i++)
     //    {
     //      mean+=list[index]; 
     //    }
     //   mean=mean/index;
     //   for(int i=0;i<index;i++)
     //    {
     //      list[index]-=mean;
     //    }
//	cout << "FINISH" << endl;
        fin.close();
}

void fileSize(string fileName, int & size) {
    //    cout<<"Come into to extract variant number"<<endl;
	size = 0;
	double data = 0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( fin.good()  ){

		fin >> data;
      //          cout<<" "<<data;
		size++;
	}
     //   cout<<endl;
     //   cout<<"size is: "<<size<<endl;
	fin.close();
}

string convertInt(int number) {
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void resetVector(char *data, int size){
	for(int i = 0; i < size; i++)
		data[i] = '0';
}

void resetVector(int * data, int size) {
	for(int i = 0; i < size; i++)
		data[i] = 0;
}

void resetVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                data[i] = 0;
}

void exportVector2File(string fileName, char * data, int size) {
	ofstream outfile(fileName.c_str(), ios::out );
	for (int i = 0; i < size; i++)
		outfile << data[i] << " ";
	outfile << endl;
	outfile.close();
}

void exportVector2File(string fileName, double * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out );
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        outfile << endl;
        outfile.close();
}

void exportVector2File(string fileName, int * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out );
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        outfile << endl;
        outfile.close();
}

void export2File(string fileName, int data) {
        ofstream outfile(fileName.c_str(), ios::out );
        outfile << data << endl;
        outfile.close();
}

int snp2Gene(int * G, int snpId, int snpCount, int geneCount) {
	for(int i = 0; i < geneCount; i++) {
		if(G[snpId*geneCount + i] == 1)
			return i;
	}
	return -1;
}

void setIdentitymatrix(int * G, int snpCount, int geneCount) {
	for(int i = 0; i < snpCount; i++) {
		for(int j = 0; j < geneCount; j++) {
			G[i*geneCount + j] = 0;
		}
		G[i*geneCount + (i/(snpCount/geneCount))] = 1;
	}
	  for(int i = 0; i < snpCount; i++) {
                for(int j = 0; j < geneCount; j++) {
                        printf("%d ", G[i*geneCount+j]);
                }
		printf("\n");
        }
}

void makeSigmaPositiveSemiDefinite(double * sigma, int size) {
	int gsl_tmp = 0;
	double matDet  = 0;
        double addDiag = 0;
        bool positive = false;
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
 //       cout<<"Value of size is: "<<size<<endl;	
	//gsl_set_error_handler_off();	
	gsl_matrix * tmpResultMatrix = gsl_matrix_calloc (size, size);	
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
	gsl_permutation *p = gsl_permutation_alloc(size);
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
        do{
		for(int i = 0; i < size; i++) {
                	for (int j = 0; j < size; j++) {
                      		if(i==j)
					gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]+addDiag);
				else
					gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]);
			}
        	}
	
		gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
       		matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);	
//		cout << matDet << "\t" << addDiag << endl;
		if(matDet > 0 ) 
			positive = true;
		else {
//			cout << "add" << endl;
			addDiag+=0.1;		
		}
	} while(!positive);
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
	for(int i = 0; i < size*size; i++){
                if(i%(size+1) == 0)
                        sigma[i] = sigma[i] + addDiag;
        }
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
}
