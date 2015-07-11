// Saniyah Bokhari
// April 15, 2011
// SubsetSumNvidia.cu; starts from test10.cu
//            fresh blocking code
// nvcc  -DPRINTIT -DPRINTTIME -DPRINTSUBSET -DPRINTTABLE -o SubsetSumNvidia SubsetSumNvidia.cu


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <sys/time.h>
#define UINT unsigned
#ifdef OSC
#define mainMemory 2.14E10
#define deviceMemory 3.0E9
#else
#define mainMemory 2.0E9
#define deviceMemory 3.0E8
#endif

void checkCUDAError(const char *msg) {
   cudaError_t err = cudaGetLastError();
   if( cudaSuccess != err) {
      fprintf(stderr, "Cuda error: %s: %s.\n", msg,
                       cudaGetErrorString( err));
      exit(EXIT_FAILURE);
   }
}
#define MIN(x,y)  ((x)<(y)?(x):(y))

__device__ __constant__ int valsDev[16384];//put in constant memory

int tstBitInArray(UINT d[], int numWords, int index);

__global__ void kernel0(UINT* devPtr){
   int tid = blockIdx.x * blockDim.x + threadIdx.x;
   UINT* row = (UINT*)((char*)devPtr);
   if(tid==0)
      row[tid]=1U;
   else
      row[tid]=0U;
}

__global__ void kernel1(UINT* devPtr, UINT* devPtr2, int pitch, int i){
   UINT* row = (UINT*)((char*)devPtr + i * pitch);
   UINT* row2;
   if(i==0)
      row2 = devPtr2;
   else
      row2 = (UINT*)((char*)devPtr + (i-1) * pitch);
   int tid = blockDim.x * blockIdx.x + threadIdx.x;
   row[tid] = row2[tid];
}

__global__ void kernel2(UINT* devPtr, UINT* devPtr2, int pitch, int width, UINT
halfSize, int i, int si, int start){
   int bit, wNew;
   UINT upper, temp11, temp12, temp21,temp22;
   int size = valsDev[si];
   UINT* row = (UINT*)((char*)devPtr + i * pitch);
   UINT* row2;
   if(i==0)
      row2 = devPtr2;
   else
      row2 = (UINT*)((char*)devPtr + (i-1) * pitch);
   int w = blockDim.x * blockIdx.x + threadIdx.x;
   if(((w & 1)==start)&&(w<width+1)){
      bit = (size % 32);
      wNew = w + (size/32);
      upper=row2[w];
      if(((w+1)*32+size)<=(halfSize+31)){
         temp11=row[wNew];
         temp12=temp11|((upper)<<bit);
         row[wNew]=temp12;
         if(bit!=0){
            temp21 = row[wNew+1];
            temp22 = temp21|((upper)>>(32-bit));
            row[wNew+1]=temp22;
         }
      }   
   }
}

int main(int argc, char *argv[])
{
   int *vals;
   UINT *d; //linearized array to store binary 1's and 0's
   UINT *dDev;
   UINT *eDev;
   int *ans;
   
   int numObj, seed;
   int numObjLo, numObjHi, sizeLo, sizeHi;
   int count = 0;
   int runs;

   assert(argc==8);
   numObjLo = atoi(argv[1]); 
   numObjHi = atoi(argv[2]);
   sizeLo = atoi(argv[3]);
   sizeHi = atoi(argv[4]);
   seed = atoi(argv[5]);
   int block_size = atoi(argv[6]);
   runs = atoi(argv[7]);

   srand(seed); //seed random number generator

   struct timeval tvTab1, tvTab2,
                  tvFree1, tvFree2, tvBac1, tvBac2;
   long tTab, tFree, tBac;
   double secTab, secFree, secBac, secPer;

do{
  
   UINT sumSize = 0;
   int k;
   numObj = numObjLo+(int)((numObjHi-numObjLo)*((double)rand()/(double)(RAND_MAX)));
   
   #ifdef PRINTIT
      printf("(n %d %d %d)(s %d)\n",
      numObjLo, numObjHi, numObj, seed);
   #endif

   vals = (int *)malloc((1+numObj)*sizeof(int));  
   assert(vals!=NULL);
   ans = (int *)malloc((1+numObj)*sizeof(int)); 
   assert(ans!=NULL);
   
   int rows=1+numObj;

   vals[0] = 0;
      
   #ifdef PRINTIT
      printf("the problem:\n");
   #endif
          
   for(k=1;k<=numObj;k++){
      vals[k]= sizeLo+(int)((sizeHi-sizeLo)*(((double)rand()/(double)(RAND_MAX))));
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
   words = (int)(ceil(((double)(halfSize))/32.0));
   
   
   //If sumSize is Odd then no subset is possible

   if((sumSize % 2)!= 0){
      #ifdef PRINTTABLE
      printf(" -- No subset was Found as Sumsize is odd -- \n");
      #endif
      free(vals);
      free(ans);
      continue;
   }
   double warea=(double)words*(double)numObj;
   double thebytes=(double)(sizeof(UINT))*warea;
   double gig = thebytes/(double)(1.073E9);
   printf("words %d numObj %d w-area %10lg bytes %lg gig %lg\n",words,numObj,warea,thebytes,gig); 
   if(thebytes > mainMemory){
      fprintf(stderr," -- Aborted, required main Memory %lg  > %lg -- \n",thebytes,mainMemory);fflush(stderr);
      free(vals);
      free(ans);
      exit(1);
   }

   int num_elements = (words+1);
   int grid_size = (int)ceil((double)(num_elements) / (double)block_size);
   if(grid_size>65535){
      fprintf(stderr,"problem too large for grid\n");fflush(stderr);
      exit(1);
   }

   d = (UINT *)malloc((1+numObj)*(words+1)*sizeof(UINT));  
   assert(d!=NULL);
   
   size_t pitch;
   
   int rowsPerChunk = MIN((int)ceil(deviceMemory/(sizeof(UINT)*(double)(words+1))),(1+numObj));
   int chunks = (int)ceil((double)rows/(double)rowsPerChunk);
   fprintf(stderr,"rows = %d ",rows);
   fprintf(stderr," chunks = %d  rowsPerChunk = %d\n",chunks, rowsPerChunk);
   gettimeofday(&tvTab1, NULL);
   cudaMallocPitch((void**)&dDev,&pitch,(words+1)*sizeof(UINT),rowsPerChunk);  
   assert(dDev!=NULL);
   cudaMalloc((void**)&eDev,pitch);  
   assert(eDev!=NULL);

   cudaMemcpyToSymbol("valsDev", vals, rows*sizeof(int));
   int ch, startRow, endRow;
   #ifdef PRINTERROR
   int num_bytes = num_elements * sizeof(UINT);
   printf("num_elements (words+1) = %d (%lu bytes)",num_elements,num_elements*sizeof(UINT));
   printf("pitch = %lu ",pitch);
   printf("numObj = %d\n",numObj);
   printf("d is = %d (%lg MB) ",(1+numObj)*(words+1),(double)((1+numObj)*(words+1)*sizeof(UINT))/1000000.0);
   printf("dDev is = %lu (%lg MB)\n",(pitch)*(1+numObj),(double)((pitch)*(1+numObj))/1000000.0);
   printf("num_bytes = %d, ",num_bytes);
   printf("grid_size %d, ",grid_size);
   printf("block_size %d\n",block_size);
   for(ch=0;ch<chunks;ch++){
      startRow=ch*rowsPerChunk;
      endRow=MIN(((ch+1)*rowsPerChunk),rows);
      printf("chunk %d, row %d to row %d\n", ch, startRow, endRow);
   }
   #endif
 
   int thisChunkSize;
   gettimeofday(&tvTab1, NULL);
   kernel0<<<grid_size,block_size>>>(dDev); //init first row
   for(ch=0;ch<chunks;ch++){                  //for each chunk
      startRow=ch*rowsPerChunk;               //starting row == 0th chunk row
      endRow=MIN(((ch+1)*rowsPerChunk),rows); //ending row == last chunk row + 1
      thisChunkSize = endRow - startRow;
      for(i=0; i<thisChunkSize; i++){     //iterate through chunk
         if(((ch==0)&&(i==0))){ //don't do this for 0th row of 0th chunk
            //printf("ch = %d i = %d ; do nothing\n", ch, i);
         }
         else{
            int start;
            int si=ch*rowsPerChunk+i;//the real index
            kernel1<<<grid_size,block_size>>>(dDev, eDev, pitch, i);//copy previous row
            for(start=0;start<2;start++){
               kernel2<<<grid_size,block_size>>>(dDev, eDev, pitch, words+1, halfSize, i, si, start);
            } 
         }
      }
      int lastRowInChunk = rowsPerChunk;
      //this is incorrect for last row in last chunk, but that is never used anyway
      cudaMemcpy(eDev, dDev+(lastRowInChunk-1)*(pitch/sizeof(UINT)), pitch, cudaMemcpyDeviceToDevice);
      int rowsToTransfer = thisChunkSize;
      cudaMemcpy2D(d+ch*rowsPerChunk*(words+1), (words+1)*sizeof(UINT),
                   dDev, pitch,
                   (words+1)*sizeof(UINT),rowsToTransfer,cudaMemcpyDeviceToHost);
                  //copy this chunk out to host in the appropriate location
   }
   gettimeofday(&tvTab2, NULL);


    if((tstBitInArray(&d[numObj*(words+1)],words,halfSize))!= 1){
       #ifdef PRINTTABLE
       printf("++ No subset was Found as last element != halfSize ++\n");
       #endif
    }    
    
    gettimeofday(&tvBac1, NULL);
    int testSum = 0;
    int index = 0;
    i = numObj;
    j = halfSize;
    int j2;
    if((tstBitInArray(&d[i*(words+1)],words,j))== 1){
       for(i=numObj;i>0;i--){
          if(((tstBitInArray(&d[(i-1)*(words+1)],words,j))== 1)&&((tstBitInArray(&d[i*(words+1)],words,j))== 1)){
             //go up
          }
          else{
            j2 = j;
            int diff;
            diff = j2-vals[i];
            j = diff;
            ans[index] = vals[i];
            //break;
            testSum +=ans[index];
            #ifdef PRINTSUBSET  
            printf(" %2d %2d %u %2d\n",ans[index],vals[i],halfSize,testSum);
            #endif
            index++;
         }
      }
    }
    #ifdef PRINTSUBSET 
    if(testSum == halfSize)
      printf("Sum Correct\n");
    else
       printf("Sum incorrect\n");
    #endif
    gettimeofday(&tvBac2, NULL);


    #ifdef PRINTTABLE
    UINT sumSize2 = 0; 
    printf("  S   W    ");  
    for(j=0;j<=halfSize;j++){
       printf(" %2d",j%100);
    }
    printf("\n");
    for(i=0;i<=numObj;i++){
       sumSize2 = sumSize2 + vals[i];
       printf("%3u %3d  %2d ",sumSize2,vals[i],i);
       for(j=0;j<=halfSize;j++){
          if((tstBitInArray(&d[i*(words+1)],words,j))== 1){
             printf(" 1 ");
          }
          if((tstBitInArray(&d[i*(words+1)],words,j))== 0){
             printf(" 0 ");
          }
       }
       printf("\n");
    }
    #endif
    
    gettimeofday(&tvFree1, NULL);
    free(vals);
    free(ans);
    free(d);
    cudaFree(dDev);
    cudaFree(eDev);
    gettimeofday(&tvFree2, NULL);

      #ifdef PRINTTIME

      //Timing for creating the table
      tTab = (tvTab2.tv_sec - tvTab1.tv_sec)*1000000 + tvTab2.tv_usec - tvTab1.tv_usec;
      secTab = (double)tTab/1000000.0;

      //Timing for backtracking
      tBac = (tvBac2.tv_sec - tvBac1.tv_sec)*1000000 + tvBac2.tv_usec - tvBac1.tv_usec;
      secBac = (double)tBac/1000000.0;

     //Timing for freeing
      tFree = (tvFree2.tv_sec - tvFree1.tv_sec)*1000000 + tvFree2.tv_usec - tvFree1.tv_usec;
      secFree = (double)tFree/1000000.0;

     unsigned long bitarea = (unsigned long)numObj*(unsigned long)halfSize;
      secPer=(secTab)/(double)(warea);
      printf("%3d Ob %4d w %7d a %10lg gig %5lg ba %10lu ",
             count,numObj,words,warea,gig,bitarea);
      printf("sTab %8g sBac %7g sFree %10g sPer %10g N",
             secTab,secBac,secFree,secPer);
      printf("\n");fflush(stdout);
      #endif

    count++;
 
    }while (count < runs);

    return 0;

}

int tstBitInArray(UINT d[], int numWords, int index){
   //   in array d[numWords] return the value of bit index
   int bit;
   bit = index % 32;
   int word;
   word = index/32;
   UINT t=(1U<<bit);
   if((d[word] & t)>0)
      return(1);
   else
      return(0);    
}
