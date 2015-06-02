#include <stdio.h>
#include <stdlib.h>

typedef struct atom{
    int x[10];
    int y[10];
}atoms;

int main(void){
    atoms xy[10];
    int i,j;
    
    for(i=0;i<10;i++)
        for(j=0;j<10;j++){
            xy[i].x[j]=0;
            xy[i].y[j]=1;
        }

    for(i=0;i<10;i++)
        for(j=0;j<10;j++)
            printf("%d,%d\n",xy[i].x[j],xy[i].y[j]);
    
    
    return 0;
}
