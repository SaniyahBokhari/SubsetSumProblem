// Saniyah Bokhari
// THIS CODE IMPLEMENTED FOR CRAY XMT
// Compile w/o Printing: gcc -O4 -o SubsetSumCRAY SubsetSumCRAY.c -lm
// Compile w Printing: gcc -DPRINTTIME -O4 -o SubsetSumCRAY SubsetSumCRAY.c -lm
// run:
// ./SubsetSumCRAY 100 1000 500 1000 127 1000
// ./SubsetSumCRAY numObjLo numObjHi sizeLo sizeHi seed runs

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#ifdef MTA
#include <sys/mta_task.h>
#include <machine/runtime.h>
#include <luc/luc_exported.h>
#include <snapshot/client.h>
#define UINT unsigned
#else
#define UINT unsigned long
#include "festubs.h"
#include <sys/time.h>
//#include <omp.h>
#endif

void bitStringIntInArray(UINT d[],int numWords,int index, UINT x[64]);
void setBitInArray(UINT d[], int numWords, int index);
void clrBitInArray(UINT d[], int numWords, int index);
int tstBitInArray(UINT d[], int numWords, int index);

int main(int argc, char *argv[]){

   //variables for timing
   #ifdef MTA
   int t0, t1,
       tMem1, tMem2, tMem3, tMem4, tTab1, tTab2,
       tFree1, tFree2, tBac1, tBac2,
       elapseMem, elapseTab, 
       elapseFree, elapseBac, delta;
   #else
   struct timeval tvMem1, tvMem2, tvMem3, tvMem4, tvTab1, tvTab2,
                  tvFree1, tvFree2, tvBac1, tvBac2;
   long tMem, tTab, tFree, tBac;
   #endif

   double secMem, secTab, secFree, secBac, secPer;

   int numObj, seed;
   int numObjLo, numObjHi, sizeLo, sizeHi;

   int *vals;
   int *ans;
   UINT **d; //array to store binary 1's and 0's
   
   int count = 0;
   int runs;

   assert(argc==7);
   numObjLo = atoi(argv[1]); 
   numObjHi = atoi(argv[2]);
   sizeLo = atoi(argv[3]);
   sizeHi = atoi(argv[4]);
   seed = atoi(argv[5]);
   runs = atoi(argv[6]);
   
   srand(seed); //seed random number generator
   do{ 
      //printf("count= %d runs= %d\n",count,runs);
      UINT sumSize = 0;
      int k;
      numObj = numObjLo+(int)((numObjHi-numObjLo)*((double)rand()/(double)(RAND_MAX)));

      #ifdef PRINTIT
          printf("(n %d %d %d)(s %d)\n",
          numObjLo, numObjHi, numObj, seed);
      #endif

      #ifdef MTA
      #pragma mta fence
      t0 =mta_get_clock(0);
      #pragma mta fence
      t1 =mta_get_clock(0);
      delta=t1-t0;
      #pragma mta fence
      tMem1 = mta_get_clock(0);
      #else
      gettimeofday(&tvMem1, NULL);
      #endif

      vals = (int *)malloc((1+numObj)*sizeof(int));  
      assert(vals!=NULL);
      ans = (int *)malloc((1+numObj)*sizeof(int)); 
      assert(ans!=NULL);

      #ifdef MTA
      #pragma mta fence
      tMem2 = mta_get_clock(0);
      #else
      gettimeofday(&tvMem2, NULL);
      #endif
      
      vals[0] = 0;
      
      #ifdef PRINTIT
         printf("the problem:\n");
      #endif
          
      for(k=1;k<=numObj;k++){
         vals[k]= sizeLo+(int)(sizeHi-sizeLo)*(((double)rand()/(double)(RAND_MAX)));
         #ifdef PRINTIT
         printf("%d  ",vals[k]);
         #endif
         sumSize = sumSize + vals[k];
      }    
      #ifdef PRINTIT
         printf("\n");
      #endif

     UINT halfSize = sumSize/2; 
     int i,j;
   
     int words; //number of words needed
     words = (int)(ceil(((double)(1+sumSize))/64.0));
   
//If sumSize is Odd then no subset is possible

   if((sumSize % 2)!= 0){
      #ifdef PRINTTABLE
      printf(" -- No subset was Found as Sumsize is odd -- \n");
      #endif
      free(vals);
      free(ans);
      continue;
   }

      #ifdef MTA
      #pragma mta fence
      tMem3 = mta_get_clock(0);
      #else
      gettimeofday(&tvMem3, NULL);
      #endif

      //Allocate memory for Binary table
      //i.e. 1's and 0's
      d = (UINT **)malloc((1+numObj)*sizeof(UINT *));  
      assert(d!=NULL);
      for(i=0;i<=numObj;i++){
         d[i]=  (UINT *)malloc(words*sizeof(UINT));
         assert(d[i]!=NULL);    
      }

      #ifdef MTA
      #pragma mta fence
      tMem4 = mta_get_clock(0);
      #else
      gettimeofday(&tvMem4, NULL);
      #endif
 
      #ifdef MTA
      #pragma mta fence
      tTab1 = mta_get_clock(0);
      #else
      gettimeofday(&tvTab1, NULL);
      #endif

      //Create Table of 1's and 0's

      //Initialization
      #pragma mta assert parallel
      for(i=0; i<=numObj; i++){
         writexf(&d[i][0],1);
      }
      #pragma mta assert parallel
      for(i=0; i<=numObj; i++){
         #pragma mta assert parallel
         for(j=1;j<words;j++){
            writexf(&d[i][j],0);
         }
      }
      //Fill the table
      int w;
      for(i=1; i<=numObj; i++){
         int pos = vals[i];
            #pragma mta assert parallel
            for(j=0; j<words; j++){
               d[i][j] = d[i-1][j];
            } 
            #pragma mta assert parallel
            for(w=0;w<words;w++){
               int thisW=w;
               int POS=pos;
               #pragma mta serial
               for(j=0;j<=63;j++){
                  int thisPos = thisW*64+j;
                  int thatPos=POS+thisPos;
                  if(((tstBitInArray(d[i-1],words,thisPos))== 1) && ((thatPos) <= halfSize)){
                        setBitInArray(d[i],words,thatPos);
                  }
               }
            }   
       }

      #ifdef MTA
      #pragma mta fence
      tTab2 = mta_get_clock(0);
      #else
      gettimeofday(&tvTab2, NULL);
      #endif
      //Last element of table must be equal to halfSize
      //Otherwise no subset is possible

      if((tstBitInArray(d[numObj],words,halfSize))!= 1){
         #ifdef PRINTTABLE
         printf("++ No subset was Found as last element != halfSize ++\n");
         #endif
         for(i=0;i<=numObj;i++){ 
            free(d[i]);
         }
         free(d);
         free(vals);
         free(ans);
         continue;
      }

      #ifdef PRINTTABLE
      //PRINTING
      UINT sumSize2 = 0; 
      printf(" S  W    ");  
      for(j=0;j<=halfSize;j++){
         printf(" %2d",j);
      }
      printf("\n");
      for(i=0;i<=numObj;i++){
         sumSize2 = sumSize2 + vals[i];
         printf("%2d %2d  %2d ",sumSize2,vals[i],i);
         for(j=0;j<=halfSize;j++){
            if((tstBitInArray(d[i],words,j))== 1){
               printf(" 1 ");
            }
            if((tstBitInArray(d[i],words,j))== 0){
               printf(" 0 ");
            }
           }
         printf("\n");
      }
      #endif

      #ifdef MTA
      #pragma mta fence
      tBac1 = mta_get_clock(0);
      #else
      gettimeofday(&tvBac1, NULL);
      #endif
      int testSum = 0;
      int index = 0;
      i=numObj;
      j = halfSize;
      int j2;

      if((tstBitInArray(d[i],words,j))== 1){
         for(i=numObj;i>0;i--){ 
            if(((tstBitInArray(d[i-1],words,j))==1)&&((tstBitInArray(d[i],words,j))==1)){
               //go up
            }
            else{
               for(k=j-1;k>=0;k--){
                  j2 = j;
                  int diff;
                  if((tstBitInArray(d[i],words,k))==1){
                     diff = j2-k;
                     if(diff == vals[i]){
                        j = k;
                        ans[index] = vals[i]; 
                        //break;
                        testSum +=ans[index];
                        #ifdef PRINTSUBSET  
                        printf(" %2d %2d %2d %2d\n",ans[index],vals[i],halfSize,testSum);
                        #endif
                        index++;
                        break;
                     }
                  }
               }
            }
         } 
      } 
      assert(testSum==halfSize);   
      #ifdef PRINTSUBSET 
      if(testSum == halfSize){
         printf("Sum Correct\n");
      }
      #endif
      #ifdef MTA
      #pragma mta fence
      tBac2 = mta_get_clock(0);
      #else
      gettimeofday(&tvBac2, NULL);
      #endif

      #ifdef MTA
      #pragma mta fence
      tFree1 = mta_get_clock(0);
      #else   
      gettimeofday(&tvFree1, NULL);
      #endif
      //Free memory
      free(vals);
      free(ans);
      for(i=0;i<=numObj;i++){
        free(d[i]);
      }
      free(d);
      #ifdef MTA
      #pragma mta fence
      tFree2 = mta_get_clock(0);
      #else   
      gettimeofday(&tvFree2, NULL);
      #endif

      count++;
      #ifdef PRINTTIME

      //Calculating the timings
    
      #ifdef MTA
      elapseMem=tMem2-tMem1-delta;
      elapseMem=elapseMem+tMem4-tMem3-delta;
      secMem= (double)(elapseMem)/(double)mta_clock_freq();
      elapseTab=tTab2-tTab1-delta;
      double cPer=(double)elapseTab/(double)(numObj*halfSize);
      secTab= (double)(elapseTab)/(double)mta_clock_freq();
      elapseBac=tBac2-tBac1-delta;
      secBac= (double)(elapseBac)/(double)mta_clock_freq();
      elapseFree=tFree2-tFree1-delta;
      secFree= (double)(elapseFree)/(double)mta_clock_freq();
      
      #else
      //Timing for mallocing
      tMem = (tvMem2.tv_sec - tvMem1.tv_sec)*1000000 + tvMem2.tv_usec - tvMem1.tv_usec;
      tMem = tMem + (tvMem4.tv_sec - tvMem3.tv_sec)*1000000 + tvMem4.tv_usec - tvMem3.tv_usec;
      secMem = (double)tMem/1000000.0;

      //Timing for creating the table
      tTab = (tvTab2.tv_sec - tvTab1.tv_sec)*1000000 + tvTab2.tv_usec - tvTab1.tv_usec;
      secTab = (double)tTab/1000000.0;

      //Timing for backtracking
      tBac = (tvBac2.tv_sec - tvBac1.tv_sec)*1000000 + tvBac2.tv_usec - tvBac1.tv_usec;
      secBac = (double)tBac/1000000.0;

     //Timing for freeing
      tFree = (tvFree2.tv_sec - tvFree1.tv_sec)*1000000 + tvFree2.tv_usec - tvFree1.tv_usec;
      secFree = (double)tFree/1000000.0;
      #endif
      secPer=(secMem+secTab+secBac+secFree)/(double)(numObj*halfSize);
      long area = (long)numObj*(long)halfSize;
      printf("%4d Obj %4d a %10ld ",
             count,numObj,area);
      printf("sMem %10g sTab %10g sBac %10g sFree %10g sPer %10g",
             secMem,secTab,secBac,secFree,secPer);
      #ifdef MTA
      printf(" cPer %10g\n",cPer);
      #else
      printf("\n");
      #endif 
      #endif

     }while (count < runs);
   return 0;
}

void bitStringIntInArray(UINT d[],int numWords,int index, UINT x[64])
{
// converts a 64 bit integer into a string of 1s and 0s
    UINT q;
    UINT r;
    int bit;
    bit = index/64;  
    //if(bit>numWords){
    //  printf("1 ERROR: P IS OUT OF RANGE");
    //  exit(0);
    //}
    for(q=0,r=1; r>0; r<<=1,q++){
       x[q] = (((d[bit] & r)==r) ? 1 : 0);
    }
}
#pragma mta inline
void setBitInArray(UINT d[], int numWords, int index){
   //  in array d[numWords] set bit index to 1
   int word;
   word = index/64;  
   int bit;
   bit = index % 64;      
   UINT t=(1ULL<<bit);
   //UINT w=d[word];
   //w= w | t;        
   //d[word]=w;
   d[word]= d[word] | t;
}
#pragma mta inline
void clrBitInArray(UINT d[], int numWords, int index){
   //   in array d[numWords] set bit index to 0
   int word;
   word = index/64;
   int bit;
   bit = index % 64;
   //UINT w=d[word];
   //w= w & ~(1ULL<<bit);        
   //d[word]=w;
   d[word]=d[word] & ~(1ULL<<bit);
}

#pragma mta inline
int tstBitInArray(UINT d[], int numWords, int index){
   //   in array d[numWords] return the value of bit index
   int bit;
   bit = index % 64;
   int word;
   word = index/64;
   UINT t=(1ULL<<bit);
   if((d[word] & t)>0)
      return(1);
   else
      return(0);    
}
