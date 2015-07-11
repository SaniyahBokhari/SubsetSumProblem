//Saniyah Bokhari
//word oriented code started JAN 1 2011
//Compile w/o Printing: gcc -O4 -o SubsetSumOMP SubsetSumOMP.c -lm
//Compile w Printing: gcc -DPRINTTIME -O4 -o SubsetSumOMP SubsetSumOMP.c -lm
//run:
//./SubsetSumOMP 100 1000 500 1000 127 1000
//./SubsetSumOMP numObjLo numObjHi sizeLo sizeHi seed runs

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
#define UINT uint64_t
#include <sys/time.h>
#include <omp.h>
#define BPW 4
#endif
#define MAXOBJ  1000
#define MAXWORD 10000
   UINT d[MAXOBJ][MAXWORD]; //array to store binary 1's and 0's
   int vals[MAXOBJ+1];
   int answ[MAXOBJ+1];

char* bitString(UINT p);
void bitStringIntInArray(UINT d[],int numWords,int index, UINT x[64]);
void setBitInArray(UINT d[], int numWords, int index);
void clrBitInArray(UINT d[], int numWords, int index);
int tstBitInArray(UINT d[], int numWords, int index);

unsigned getsize(void)
{
/* function to return current memory use size                       JRS 8/10  */
/* return is type 'unsigned'; in kilowords, 4B/word in Linux                  */
/* return 0 means error opening process file structure with print to 'stderr' */
    unsigned size;
    char buf[64];

    snprintf(buf, 64, "/proc/%u/statm", (unsigned)getpid());
    FILE* pf = fopen(buf, "r");
    if (pf) {
        fscanf(pf, "%u", &size);
        fclose(pf);
    } else {
        fprintf(stderr, "getsize: error opening %s\n", buf);
        return(0);
    }
    return(size);
}


int main(int argc, char *argv[]){

   //variables for timing
   printf("start main %u\n",getsize()*BPW);
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
      printf("start do loop %u\n",getsize()*BPW);
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

      printf("after 1d arrays %u\n",getsize()*BPW);

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
     words = (int)(ceil(((double)(halfSize))/64.0));
   
     //If sumSize is Odd then no subset is possible

     if((sumSize % 2)!= 0){
        #ifdef PRINTTABLE
        printf(" -- No subset was Found as Sumsize is odd -- \n");
        #endif
        printf("abort, odd; size now: %u\n",getsize()*BPW);
        continue;
     }
     double warea=(double)words*(double)numObj;
     double thebytes=8.0*warea;
     printf("words %d numObj %d w-area %10lg bytes %lg\n",words,numObj,warea,thebytes); 
     if(thebytes > 1.29E10){
       printf(" -- Aborted, > 12G  -- \n");
       printf("abort, size; size now: %u\n",getsize()*BPW);
      
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

      printf("after 2d arrays %u\n",getsize()*BPW);
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
      #pragma omp parallel for 
      for(i=0; i<=numObj; i++){
         d[i][0]=1ULL;
      }
      #pragma omp parallel for 
      for(i=0; i<=numObj; i++){
         //#pragma omp parallel for
         for(j=1;j<words;j++){
            d[i][j]=0ULL;
         }
      }
      //Fill the table
      int w;
      for(i=1; i<=numObj; i++){
         int size = vals[i];
         int bit, wNew, start;
         UINT upper, temp11, temp12, temp21,temp22;
         #pragma omp parallel for  
         for(j=0; j<words; j++){
            d[i][j] = d[i-1][j];
         } 
         for(start=0;start<2;start++){
            #pragma omp parallel for private(bit, wNew,upper, temp11,temp12,temp21,temp22) 
            for(w=start;w<words;w=w+2){
               bit = (size % 64);
               wNew = w + (size/64);
               upper=d[i-1][w];
               if(((w+1)*64+size)<=(halfSize+63)){
                  temp11=d[i][wNew];
                  temp12=temp11|((upper)<<bit);
                  d[i][wNew]=temp12;
                  if(bit!=0){
                     temp21 = d[i][wNew+1];
                     temp22 = temp21|((upper)>>(64-bit));
                     d[i][wNew+1]=temp22;
                  }
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

         printf("abort, no solution; size now: %u\n",getsize()*BPW);
         continue;

      }

      

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
                        answ[index] = vals[i]; 
                        //break;
                        testSum +=answ[index];
                        #ifdef PRINTSUBSET  
                        printf(" %2d %2d %2d %2d\n",answ[index],vals[i],halfSize,testSum);
                        #endif
                        index++;
                        break;
                     }
                  }
               }
            }
         } 
      } 
      //assert(testSum==halfSize);   
      #ifdef PRINTSUBSET 
      if(testSum == halfSize){
         printf("Sum Correct\n");
      }
      else
         printf("Sum incorrect\n");
      #endif
      #ifdef MTA
      #pragma mta fence
      tBac2 = mta_get_clock(0);
      #else
      gettimeofday(&tvBac2, NULL);
      #endif

      #ifdef PRINTTABLE
      //PRINTING
      UINT sumSize2 = 0; 
/*
      printf("           ");  
      for(j=0;j<=halfSize;j++){
         printf(" %2d",j/100);
      }
      printf("\n");
*/
      printf("  S   W    ");  
      for(j=0;j<=halfSize;j++){
         printf(" %2d",j%100);
      }
      printf("\n");
      for(i=0;i<=numObj;i++){
         sumSize2 = sumSize2 + vals[i];
         printf("%3d %3d  %2d ",sumSize2,vals[i],i);
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
      tFree1 = mta_get_clock(0);
      #else   
      gettimeofday(&tvFree1, NULL);
      #endif
      //Free memory
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
      double cPer=(double)elapseTab/(double)(numObj*words);
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
      long area = (long)numObj*(long)words;
      long bitarea = (long)numObj*(long)halfSize;
      secPer=(secTab)/(double)(area);
      printf("%3d Ob %4d w %7d a %10ld ba %10ld ",
             count,numObj,words,area,bitarea);
      printf("sMem %10g sTab %8g sBac %7g sFree %10g sPer %10g",
             secMem,secTab,secBac,secFree,secPer);
      #ifdef MTA
      printf(" cPer %5g\n",cPer);
      #else
      printf("\n");
      #endif 
      #endif

      printf("end do loop %u\n",getsize()*BPW);
     }while (count < runs);
   printf("end prog  %u\n",getsize()*BPW);
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
//#pragma mta inline
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
//#pragma mta inline
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

//#pragma mta inline
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
char* bitString(UINT p)
{
// converts a 64 bit integer into a string of 1s and 0s

    static char str[65] = {0};
    UINT q;
    UINT r;

    for(q=0ULL,r=1ULL<<63; r>0ULL; r>>=1ULL,q++){
       str[q] = (((p & r)==r) ? '1' : '0');
    }
    return str;
}

