#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "MT.h"
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
int init(unit peptide[],int state[],int seq_size);
int main_loop(unit peptide[],int total_step,constant constants,int seq_size);
int move(unit peptide[],int seq_size);
int compare_states(unit peptide1[],unit peptide2[],int seq_size);


int local2vector(unit peptide[],int seq_size);
int vector2xy(unit peptide[],int seq_size);
int setlocal(unit peptide[],int point,int value,int seq_size);

int copy(unit tmp_peptide[],unit peptide[],int seq_size);
int mc(unit peptide1[],unit peptide2[],constant constants,int seq_size);
int calc_energy(unit peptide[],int seq_size);
//---------- move sets ----------//
int flip(unit peptide[],int point,int type,int seq_size);
int cornerflip(unit peptide[],int point,int seq_size);
int searchcorner(unit peptide[],int seq_size);

int flipflop(unit peptide[],int point,int seq_size);
int searchcorner(unit peptide[],int seq_size);
//------------------------------//


int show(unit peptide[],int seq_size);
//----------------------------------------//
int main(void){
    //int seq[]={1,0,0,1,0,0,1,1};
    int seq[]={1,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,1,1};
    int total_step=1000*1000*10;
    int seq_size=sizeof(seq)/sizeof(int);
    constant constants;
    unit *peptide;
    clock_t start,end;
    start=clock();

    
    srand(0);
    constants.k=1;
    constants.T=0.4;
    
    peptide=(unit *)malloc(sizeof(unit)*seq_size);
    
    //test(peptide,seq,seq_size,total_step,constants);

    init(peptide,seq,seq_size);
    main_loop(peptide,total_step,constants,seq_size);
    printf("Energy =%4d\n",calc_energy(peptide,seq_size));
    end=clock();
    printf("%8.3lf\n",(double)(end-start)/CLOCKS_PER_SEC);
    
    return 0;
}


//----------------------------------------//
int init(unit peptide[],int state[],int seq_size){
    int i;
    for(i=0;i<seq_size;i++){
        peptide[i].state=state[i];
        peptide[i].xy[0]=0;
        peptide[i].xy[1]=i;
    }

    for(i=0;i<seq_size;i++){
        peptide[i].local=0;
    }

    for(i=0;i<seq_size;i++){
        peptide[i].vector=0;
    }
    return TRUE;
}

int main_loop(unit peptide[],int total_step,constant constants,int seq_size){
    int istep;
    unit tmp_peptide[seq_size];
    init_genrand(10);
    
    for(istep=0;istep<total_step;istep++){
        copy(tmp_peptide,peptide,seq_size);
        move(tmp_peptide,seq_size);
        if(TRUE!=compare_states(tmp_peptide,peptide,seq_size))
            mc(peptide,tmp_peptide,constants,seq_size);
    }
    
    show(peptide,seq_size);
    return TRUE;
}

int compare_states(unit peptide1[],unit peptide2[],int seq_size){
    int i;

    for(i=2;i<seq_size;i++)
        if(peptide1[i].local!=peptide2[i].local)
            return FALSE;
    return TRUE;
}

int mc(unit peptide1[],unit peptide2[],constant constants,int seq_size){
    int E1,E2;
    double prob;
    E1=calc_energy(peptide1,seq_size);
    E2=calc_energy(peptide2,seq_size);
    
    if(E2<E1){
        copy(peptide1,peptide2,seq_size);
    }else{
        prob=genrand_real3();
        if(prob<exp(-(E2-E1)/(constants.k*constants.T))){
            copy(peptide1,peptide2,seq_size);
        }
    }

    return TRUE;
}

int calc_energy(unit peptide[],int seq_size){
    int energy=0;
    int i,j;
    int x,y;

    for(i=0;i<seq_size;i++)
        for(j=i+1;j<seq_size;j++)
            if(peptide[i].xy[0]==peptide[j].xy[0] && peptide[i].xy[1]==peptide[j].xy[1])
                return 1000;

    for(i=0;i<seq_size;i++){
        if(peptide[i].state==1){
            for(j=i+3;j<seq_size;j++){
                if(peptide[j].state==1){
                    x=peptide[i].xy[0]-peptide[j].xy[0];
                    y=peptide[i].xy[1]-peptide[j].xy[1];
                    if(1==(x*x+y*y)){
                        energy--;
                    }
                }
            }
        }
    }
    return energy;
}

int copy(unit tmp_peptide[],unit peptide[],int seq_size){
    int i=0;
    for(i=0;i<seq_size;i++){
        tmp_peptide[i].local=peptide[i].local;
        tmp_peptide[i].vector=peptide[i].vector;
        tmp_peptide[i].xy[0]=peptide[i].xy[0];
        tmp_peptide[i].xy[1]=peptide[i].xy[1];
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
    int switcher=genrand_int32()%2;
    int point;
    int leftright;
    
    if(switcher==0){
        // flip
        // point = 0,1,2,...,seq_size-1
        // leftright = 0(left),1(right);
        point=genrand_int32()%(seq_size-3)+2;
        leftright=genrand_int32()%2;
        if(leftright==peptide[point].local)
            return TRUE;
        flip(peptide,point,leftright,seq_size);
        return TRUE;
        
    }else if(switcher==1){
        // cornerflip
        // searche corner.
        point=searchcorner(peptide,seq_size);
        if(point==0)
            return TRUE;
        cornerflip(peptide,point,seq_size);

    }
    return TRUE;
}


//-------------------- move sets --------------------//
int searchcorner(unit peptide[],int seq_size){
    int point=0;
    int i;
    int cornernumber=0;
    int select;
    int start=3;
    
    for(i=start;i<seq_size-1;i++)
        if(peptide[i].local!=0){
            cornernumber++;
            point=i;
        }

    if(cornernumber==0)
        return point;
    else if(cornernumber==1)
        return point;
    else{
        select=genrand_int32()%cornernumber;
        for(i=start;i<seq_size-1;i++){
            if(peptide[i].local!=0){
                if(select==0)
                    return i;
                select--;
            }
        }
        return point;
    }
    return point;
}

int cornerflip(unit peptide[],int point,int seq_size){
    if(peptide[point].local==1){
        setlocal(peptide,point-1,1,seq_size);
        setlocal(peptide,point,-2,seq_size);
        setlocal(peptide,point+1,1,seq_size);
        local2vector(peptide,seq_size);
        vector2xy(peptide,seq_size);
    }else{
        setlocal(peptide,point-1,-1,seq_size);
        setlocal(peptide,point,2,seq_size);
        setlocal(peptide,point+1,-1,seq_size);
        local2vector(peptide,seq_size);
        vector2xy(peptide,seq_size);
    }
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


int show(unit peptide[],int seq_size){
    int i;
    for(i=0;i<seq_size;i++){
        printf("%3dth,local=%3d,vector=%3d,xy=[%3d,%3d]\n",i,peptide[i].local,peptide[i].vector,peptide[i].xy[0],peptide[i].xy[1]);
    }
    return TRUE;
}

