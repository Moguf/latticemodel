#include <stdio.h>
#include <stdlib.h>


typedef struct unit_struct{
    int state;
    int xy[2];
} unit;

typedef struct const_units{
    double k;
    double T;
} constant;


void init(unit units[],int state[],int size){
    int i;
    for(i=0;i<size;i++){
        units[i].state=state[i];
        units[i].xy[0]=0;
        units[i].xy[1]=1;
    }
}

void main_loop(unit units[],int total_step,constant constants){
    int istep;
    total_step=10;
    double MAX=2147483647.0;
    for(istep=0;istep<total_step;istep++){
        printf("%lf\n",MAX/(rand()+1));
    }
}


int main(void){
    int seq[]={1,1,1,1,0,0,0,0,0};
    int total_step=1000;
    int seq_size=sizeof(seq)/sizeof(int);
    constant constants;
    unit *units;
    srand(0);
    constants.k=1;
    constants.T=1;
    
    units=(unit *)malloc(sizeof(unit)*seq_size);

    init(units,seq,seq_size);
    main_loop(units,total_step,constants);
    

    return 0;
}
