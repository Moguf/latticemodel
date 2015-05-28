#include <stdio.h>
#include <stdlib.h>

//----------------------------------------//
#define TRUE 1
#define FALSE 0


typedef struct unit_struct{
    int state;
    int local;
    // Local is position from i-2th i-1th beads.
    // ex) Position of ith bead is left from line between i-2th i-1th beads of view.
    // ex) In that case, Local is -1.
    // ex) In 'right' case, Local is 1.
    // ex) In 'straght' case, Local is 0.
    int vector;
    // Vector is vector of ith bead.
    // The vector base is (i-1)th bead and toward ith bead.
    // The vector has one of {0,3,6,9}.
    int xy[2];
} unit;

typedef struct const_values{
    double k;
    double T;
} constant;

//-------------------- define func --------------------//
int init(unit peptide[],int state[],int size);
int main_loop(unit peptide[],int total_step,constant constants,int size);
int move(unit peptide[],int seq_size);

int local2vector(unit peptide[],int seq_size);
int vector2xy(unit peptide[],int seq_size);
int setlocal(unit peptide[],int point,int value,int seq_size);

//---------- move sets ----------//
int flip(unit peptide[],int point,int type,int seq_size);
int conerflip(unit peptide[],int point,int seq_size);

int test(unit peptide[],int state[],int size,int total_step,constant constants);
int show(unit peptide[],int seq_size);
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
    for(i=0;i<size;i++){
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

int local2vector(unit peptide[],int seq_size){
    int i;
    for(i=2;i<seq_size;i++){
        if(peptide[i].local==0){
            peptide[i].vector=peptide[i-1].vector;
        }else if(peptide[i].local==1){
            peptide[i].vector=(peptide[i-1].vector+3)%12;
        }else{
            peptide[i].vector=(peptide[i-1].vector+9)%12;
        }
    }
    return TRUE;
}

int vector2xy(unit peptide[],int seq_size){
    int i;
    for(i=2;i<seq_size;i++){
        if(peptide[i].vector==0){
            peptide[i].xy[0]=peptide[i-1].xy[0];
            peptide[i].xy[1]=peptide[i-1].xy[1]+1;
        }else if(peptide[i].vector==3){
            peptide[i].xy[0]=peptide[i-1].xy[0]+1;
            peptide[i].xy[1]=peptide[i-1].xy[1];
        }else if(peptide[i].vector==6){
            peptide[i].xy[0]=peptide[i-1].xy[0];
            peptide[i].xy[1]=peptide[i-1].xy[1]-1;
        }else{
            peptide[i].xy[0]=peptide[i-1].xy[0]-1;
            peptide[i].xy[1]=peptide[i-1].xy[1];
        }
    }
    return TRUE;
    
}

int setlocal(unit peptide[],int point,int value,int seq_size){
    peptide[point].local+=value;
    if(peptide[point].local<=-2)
        peptide[point].local=-1;
    else if(peptide[point].local>=2)
        peptide[point].local=1;
    return TRUE;
}

int move(unit peptide[],int seq_size){
    return TRUE;
}


//-------------------- move sets --------------------//

int conerflip(unit peptide[],int point,int seq_size){

    return TRUE;
}

int flip(unit peptide[],int point,int type,int seq_size){
    // type=1(right),type=0(left)
    if(type==1){ 
        setlocal(peptide,point,1,seq_size);
        local2vector(peptide,seq_size);
        vector2xy(peptide,seq_size);
    }else{
        setlocal(peptide,point,-1,seq_size);
        local2vector(peptide,seq_size);
        vector2xy(peptide,seq_size);
    }
    return TRUE;
}



//-------------------- test --------------------//
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
        if(peptide[i].xy[0]!=0 &&peptide[i].xy[1]!=i){
            printf("in init xy.\n");
            printf("%dth error!!\n",i);
            return FALSE;
        }
    }
        
    main_loop(peptide,total_step,constants,seq_size);
    
    move(peptide,seq_size);

    flip(peptide,2,0,seq_size);
    printf("flip\n");
    show(peptide,seq_size);
    printf("\
---\n\
  0th,local=  0,vector=  0,xy=[  0,  0]\n\
  1th,local=  0,vector=  0,xy=[  0,  1]\n\
  2th,local= -1,vector=  9,xy=[ -1,  1]\n\
  3th,local=  0,vector=  9,xy=[ -2,  1]\n\
  4th,local=  0,vector=  9,xy=[ -3,  1]\n\
  5th,local=  0,vector=  9,xy=[ -4,  1]\n\
  6th,local=  0,vector=  9,xy=[ -5,  1]\n\
  7th,local=  0,vector=  9,xy=[ -6,  1]\n\
  8th,local=  0,vector=  9,xy=[ -7,  1]\n\
----------------------------------------\n");

    init(peptide,seq,seq_size);
    flip(peptide,3,1,seq_size);
    flip(peptide,5,1,seq_size);
    show(peptide,seq_size);
    printf("\
---\n\
  0th,local=  0,vector=  0,xy=[  0,  0]\n\
  1th,local=  0,vector=  0,xy=[  0,  1]\n\
  2th,local=  0,vector=  0,xy=[  0,  2]\n\
  3th,local=  1,vector=  3,xy=[  1,  2]\n\
  4th,local=  0,vector=  3,xy=[  2,  2]\n\
  5th,local=  1,vector=  6,xy=[  2,  1]\n\
  6th,local=  0,vector=  6,xy=[  2,  0]\n\
  7th,local=  0,vector=  6,xy=[  2, -1]\n\
  8th,local=  0,vector=  6,xy=[  2, -2]\n\
----------------------------------------\n");


    printf("conerflip\n");
    printf("before---\n");
    init(peptide,seq,seq_size);
    flip(peptide,3,0,seq_size);
    show(peptide,seq_size);
    printf("after----\n");
    conerflip(peptide,3,seq_size);    
    show(peptide,seq_size);
    printf("\
answer---\n\
  0th,local=  0,vector=  0,xy=[  0,  0]\n\
  1th,local=  0,vector=  0,xy=[  0,  1]\n\
  2th,local= -1,vector=  9,xy=[ -1,  1]\n\
  3th,local=  1,vector=  0,xy=[ -1,  2]\n\
  4th,local= -1,vector=  9,xy=[ -2,  2]\n\
  5th,local=  0,vector=  9,xy=[ -3,  2]\n\
  6th,local=  0,vector=  9,xy=[ -4,  2]\n\
  7th,local=  0,vector=  9,xy=[ -5,  2]\n\
  8th,local=  0,vector=  9,xy=[ -6,  2]\n\
----------------------------------------\n");


    return TRUE;
}



int show(unit peptide[],int seq_size){
    int i;
    for(i=0;i<seq_size;i++){
        printf("%3dth,local=%3d,vector=%3d,xy=[%3d,%3d]\n",i,peptide[i].local,peptide[i].vector,peptide[i].xy[0],peptide[i].xy[1]);
    }
    return TRUE;
}

