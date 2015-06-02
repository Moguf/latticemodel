#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define SIZE 4000

int main(void){
    double a[SIZE][SIZE];
    double b[SIZE][SIZE];
    double c[SIZE][SIZE];
    int i,j,k;
    time_t start,end;
    start=clock();

#pragma omp parallel for
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++){
            a[i][j]=0;
            b[i][j]=1;
            c[i][j]=0;
        }
    for(i=0;i<SIZE;i++)
        a[i][i]=1;
        
#pragma omp parallel for
    for(i=0;i<SIZE;i++){
        for(j=0;j<SIZE;j++){
            for(k=0;k<SIZE;k++){
                c[i][j]+=a[i][k]*b[k][j];
            }
        }
    }

    end=clock();

    printf("%lf\n",a[0][0]);    
    printf("time=%lf\n",(double)(end-start)/CLOCKS_PER_SEC);
    return 0;
}
