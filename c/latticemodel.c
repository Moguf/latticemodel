#include <stdio.h>
#include <stdlib.h>

//----------------------------------------//
#define TRUE 1
#define FALSE 0


typedef struct unit_struct{
    int state;
    int local;
    int vector;
    int xy[2];
} unit;

typedef struct const_values{
    double k;
    double T;
} constant;

//----------------------------------------//
int init(unit peptide[],int state[],int size);
int main_loop(unit peptide[],int total_step,constant constants,int size);

int test(unit peptide[],int state[],int size,int total_step,constant constants);

//----------------------------------------//
int main(void){
    int seq[]={1,1,1,1,0,0,0,0,0};
    int total_step=10;
    int seq_size=sizeof(seq)/sizeof(int);
    constant constants;
    unit *peptide;
    srand(0);
    constants.k=1;
    constants.T=1;
    
    peptide=(unit *)malloc(sizeof(unit)*seq_size);
    
    test(peptide,seq,seq_size,total_step,constants);
    /*
    init(peptide,seq,seq_size);
    main_loop(peptide,total_step,constants,seq_size);
    */

    return 0;
}


//----------------------------------------//
int init(unit peptide[],int state[],int size){
    int i;
    for(i=0;i<size;i++){
        peptide[i].state=state[i];
        peptide[i].xy[0]=0;
        peptide[i].xy[1]=i;
    }
    // local i=2 to i=size
    for(i=2;i<size;i++){
        peptide[i].local=0;
    }
    // vector i=1 to i=size
    for(i=0;i<size;i++){
        peptide[i].vector=0;
    }
    return TRUE;
}

int main_loop(unit peptide[],int total_step,constant constants,int size){
    int istep;
    double MAX=2147483647.0;
    for(istep=0;istep<total_step;istep++){
        printf("%lf\n",MAX/(rand()+1));
    }
    return TRUE;
}

int local2vector(){
    return TRUE;
}

int vector2xy(){
    return TRUE;
}

int move(){
    return TRUE;
}

int conerflip(){
    return TRUE;
}

int pointflip(){
    return TRUE;
}

int rigidrotation(){
    return TRUE;    
}

int test(unit peptide[],int seq[],int seq_size,int total_step,constant constants){
    int i;
    init(peptide,seq,seq_size);
    for(i=0;i<seq_size;i++){
        if(peptide[i].state!=seq[i]){
            printf("in init state.\n");
            printf("%dth error!!\n",i);
            return FALSE;
        }
        if(peptide[i].local!=0){
            printf("in init local.\n");
            printf("%dth error!!\n",i);
            return FALSE;
        }
        if(peptide[i].vector!=0){
            printf("in init vector.\n");
            printf("%dth error!!\n",i);
            return FALSE;
        }
        if(peptide[i].xy[0]!=0 &&peptide[i].xy[1]!=i &&){
            printf("in init xy.\n");
            printf("%dth error!!\n",i);
            return FALSE;
        }
    }
        
    main_loop(peptide,total_step,constants,seq_size);

    return TRUE;
}

