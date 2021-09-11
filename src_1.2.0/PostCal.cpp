#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>
#include <omp.h>
#include "Util.h"
#include "PostCal.h"
using namespace arma;

//void printGSLPrint(mat &A, int row, int col) {
//	for(int i = 0; i < row; i++) {
//		for(int j = 0; j < col; j++)
//			printf("%g ", A(i, j));
//		printf("\n");
//	}	
//}

string PostCal::convertConfig2String(int * config, int size) {
	string result = "0";
	for(int i = 0; i < size; i++)
		if(config[i]==1)
			result+= "_" + convertInt(i);
	return result;
}

//Calculate the likelihood of each causal state.
double PostCal::likelihood(int * configure_x) {
	int causalCount = 0;
	int index_C = 0;
        double matDet = 0;
	double res    = 0;

	for(int i = 0; i < snpCount; i++)
         { 
           causalCount += configure_x[i];
         }
	if(causalCount == 0){
                mat tmpResultMatrix11 = trans(stat)*stat;
		res = tmpResultMatrix11(0,0);	
		baseValue = res;
                matDet =1;
		res = res - baseValue;
		return( exp(-res/2)/sqrt(abs(matDet)) );
	}

        vector<int> noz;
	for(int i = 0; i < snpCount; i++) {
                if (configure_x[i] == 0)	continue;
                else {
                        noz.push_back(i);
        
                }
        }
        mat X_noz=mat(number,noz.size(),fill::zeros);
        for(int i=0;i<noz.size();i++)
         {
           for(int j=0;j<number;j++)
            {
              X_noz(j,i)=X(j,noz[i]);
            }
         }

        double inv_ve=10.0;
        vector<double> test;
        vector<double> test_tmplogDet;
        double start=0.2;
        double end=0.8;
        double win=0.1;
        for(double x=start;x<=end;x+=win)
          {
            mat XtX=trans(X_noz)*X_noz;
  //             cout<<"XtX is: "<<XtX(0,0)<<endl;
               mat inv_XtX=pinv(XtX);
    //           cout<<"inv_XtX is: "<<inv_XtX(0,0)<<endl;
               mat xy=trans(X_noz)*stat;
      //         cout<<"xy is: "<<xy(0,0)<<endl;
               mat res_test=inv_XtX*trans(X_noz)*stat;
               mat resi=stat-X_noz*res_test;
           //    double  phe_variance=stddev(resi.col(0));
               colvec residuals = stat - X_noz * res_test;
               double s2=0;
               
               for(int i_x=0;i_x<number;i_x++)
                {
                  s2+=residuals(i_x,0)*residuals(i_x,0);
                }
               s2/=number;
        //       cout<<"s2 is: "<<s2<<endl;
             mat eVb=xy/s2;
        //    cout<<"eVb is: "<<eVb<<endl;
            mat eVg_inv=trans(X_noz)*X_noz/s2;
          //  cout<<"eVg_inv is: "<<eVg_inv(0,0)<<endl;
            mat covariance2=mat(noz.size(),noz.size(),fill::eye)+eVg_inv*(x*x)*mat(noz.size(),noz.size(),fill::eye);
            matDet =log(abs(det(covariance2)));
            mat covariance =(x*x)*mat(noz.size(),noz.size(),fill::eye)*pinv(covariance2);
          //  cout<<"covariance is: "<<covariance<<endl;
            mat tmp_AA=trans(eVb)*(covariance)*eVb;
            res =tmp_AA(0,0);
            test.push_back(res);
            test_tmplogDet.push_back(matDet);
	    if(matDet==0) 
             {
		exit(0);
       	     }
           }
         res=0;
         for(int i_test=0;i_test<test.size();i_test++)
          {
            res+=exp(test[i_test]/2-test_tmplogDet[i_test]/2);
          }
      double  test_t=res/int(1+(end-start)/win);
   //   cout<<"res is: "<<test_t<<endl;	
      return(res/int(1+(end-start)/win));
}


//Scan for each causal state, and store them in a vecotr, used in parallel manner
int PostCal::nextBinary(int * data, int size) {
	int i = 0;
	int total_one = 0;	
	int index = size-1;
        int one_countinus_in_end = 0;

        while(index >= 0 && data[index] == 1) {
                index = index - 1;
                one_countinus_in_end = one_countinus_in_end + 1;
	}
	if(index >= 0) {
        	while(index >= 0 && data[index] == 0) {
               	 index = index - 1;	
		}
	}
        if(index == -1) {
                while(i <  one_countinus_in_end+1 && i < size) {
                        data[i] = 1;
                        i=i+1;
		}
                i = 0;
                while(i < size-one_countinus_in_end-1) {
                        data[i+one_countinus_in_end+1] = 0;
                        i=i+1;
		}
	}
        else if(one_countinus_in_end == 0) {
                data[index] = 0;
                data[index+1] = 1;
	} else {
                data[index] = 0;
                while(i < one_countinus_in_end + 1) {
                        data[i+index+1] = 1;
			if(i+index+1 >= size)
				printf("ERROR3 %d\n", i+index+1);
                        i=i+1;
		}
                i = 0;
                while(i < size - index - one_countinus_in_end - 2) {
                        data[i+index+one_countinus_in_end+2] = 0;
			if(i+index+one_countinus_in_end+2 >= size) {
				printf("ERROR4 %d\n", i+index+one_countinus_in_end+2);
			}
                        i=i+1;
		}
	}
	i = 0;
	total_one = 0;
	for(i = 0; i < size; i++)
		if(data[i] == 1)
			total_one = total_one + 1;
	
	return(total_one);		
}
int PostCal::decomp(vector<int> &str, int *data, int num)
  {
           int init_start = 0;
           int init_end = 0;
           int test=str.size();
           if(str[0]!=-1)
            {
              for (int idx = 0; idx < str.size(); idx++)
               {
     //            cout<<" "<<str[idx]<<endl;
                for (int i_x = init_start; i_x<init_start + str[idx]; i_x++)
                 {
                   data[i_x] = 0;
                 }
                int last = init_start + str[idx];
                data[last] = 1;
                init_start = init_start + str[idx] + 1;
              }
             for (int x = init_start + 1; x<num; x++)
              {
                data[x] = 0;
              }
            } else
            {
              for (int y =0; y<num; y++)
               {
                 data[y] = 0;
               }
            }
          int num_x=0;
          for (int x = 0; x<num; x++)
           {
             if(data[x]==1)
              {
                num_x++;
              }
           }
          return num_x;
  }
string PostCal::convert_symbol(int  *data, int num,vector<int> &output)
  {
    string str;
    string sym = "";
    int a = 0;
    int index=0; 
    for (int n = 0; n<num; n++)
     {
       if (data[n] == 1)
        {
          ostringstream temp;  //temp as in temporary
          temp<<a;
          string x=temp.str();      //str is temp as string
          output.push_back(a);
          str += x;
          index=1;
          a = 0;
        }
       else
        {
          a++;
        }
     }
    if(index==0)
     {
       output.push_back(-1);
     }
    return str;
  }
double PostCal::computeTotalLikelihood(double * stat,map<int,double>& Weight, int nthread) {	
	double sumLikelihood = 0;
	double tmp_likelihood = 0;
	long int total_iteration = 0 ;
	int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data	

	for(long int i = 0; i <= maxCausalSNP; i++)
		total_iteration = total_iteration + nCr(snpCount, i);
//	cout << snpCount << endl;
//	cout << "Max Causal=" << maxCausalSNP << endl;
//	cout << "Total="      << total_iteration << endl;
	for(long int i = 0; i < snpCount; i++) 
		configure[i] = 0;
        vector<vector<int> > test;
        vector<int> res;
        convert_symbol(configure, snpCount,res);
        test.push_back(res);
        for(long int i = 1; i < total_iteration; i++) 
         {
          nextBinary(configure, snpCount);
          vector<int> res;
          convert_symbol(configure, snpCount,res);
          test.push_back(res);
         }
        for(long int i=0;i<snpCount;i++)
         {
           configure[i]=0;
         }
        double  prior_x=1;
        for(long int zeng = 0; zeng < snpCount; zeng++)
         {
           double test_x=0.01;
           if(configure[zeng]==1)
            {
              test_x=test_x;
            } else
            {
              test_x=1-test_x;
            }
           prior_x=prior_x*test_x;
         }
        int test_x=0;
        tmp_likelihood = likelihood(configure) * prior_x;
//        cout<<"tmp_likelihood is: "<<tmp_likelihood<<endl;
        sumLikelihood += tmp_likelihood;
        for(int j = 0; j < snpCount; j++)
         {
           postValues[j] = postValues[j] + tmp_likelihood * configure[j];
         }
        int confi[nthread][snpCount];
        for(long int i=0;i<nthread;i++)
         {
           for(long int j_x=0;j_x<snpCount;j_x++)
            {
               confi[i][j_x]=configure[j_x];
            }
         }
        int tid;
        omp_set_num_threads(nthread);
        prior_x=1;
        int num;
        int nloops=0;
        #pragma omp parallel   private (tid,prior_x,num,tmp_likelihood,nloops)
          {
             #pragma omp for
	     for(long int i = 1; i < total_iteration; i++) {
                double prior=1;
                tid=omp_get_thread_num();
                num=snpCount;
               #pragma omp critical
               {
               int init_start = 0;
               int init_end = 0;
               vector<int> str=test[i];
               int test=str.size();
               if(str[0]!=-1)
                {
                  for (int idx = 0; idx < str.size(); idx++)
                   {
                     for (int i_x = init_start; i_x<init_start + str[idx]; i_x++)
                       {
                   confi[tid][i_x] = 0;
                 }
                int last = init_start + str[idx];
                confi[tid][last] = 1;
                init_start = init_start + str[idx] + 1;
              }
             for (int x = init_start; x<num; x++)
              {
                confi[tid][x] = 0;
              }
            } else
            {
              for (int y =0; y<num; y++)
               {
                 confi[tid][y] = 0;
               }
            }
          for (int x = 0; x<num; x++)
           {
             if(confi[tid][x]==1)
              {
              }
           }
           }
                nloops++;
                prior_x=1;
                #pragma omp critical
                 {
                    for(long int zeng = 0; zeng < snpCount; zeng++)
         {
           double test_x=1;
           if(confi[tid][zeng]==1)
            {
              test_x=0.01;
            } else
            {
              test_x=1-0.01;
            }
           prior_x=prior_x*test_x;
         }
                 }
                tmp_likelihood = likelihood(confi[tid]) * prior_x;
//                cout<<"tmp_likelihood is: "<<tmp_likelihood<<", prior_x is: "<<prior_x<<endl;
                #pragma omp critical
                 {
                   sumLikelihood += tmp_likelihood;
		   for(int j = 0; j < snpCount; j++) 
                    {
                        postValues[j] = postValues[j] + tmp_likelihood * confi[tid][j];
                    }
                 }
           }
           tid = omp_get_thread_num();
         }
        free(configure);
        return(sumLikelihood);
}
double PostCal::findOptimalSetGreedy(double * stat, char * configure, int *rank,  double inputRho, map<int, double>& Weight, int nthread) {
	int index = 0;
        double rho = 0;
        double total_post = 0;
        totalLikeLihood = computeTotalLikelihood(stat, Weight, nthread);
	for(int i = 0; i < snpCount; i++)
		total_post += postValues[i];
	printf("Total Likelihood= %e SNP=%d \n", total_post, snpCount);
        std::vector<data> items;
        std::set<int>::iterator it;
        for(int i = 0; i < snpCount; i++) {
             items.push_back(data(postValues[i]/total_post, i, 0));
        }
        printf("\n");
        std::sort(items.begin(), items.end(), by_number());
        for(int i = 0; i < snpCount; i++)
                rank[i] = items[i].index1;
        for(int i = 0; i < snpCount; i++)
                configure[i] = '0';
        do{
                rho += postValues[rank[index]]/total_post;
                configure[rank[index]] = '1';
                printf("%d %e\n", rank[index], rho);
                index++;
        } while( rho < inputRho);

        printf("\n");
	return(0);
}
