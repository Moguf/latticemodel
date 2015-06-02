#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "MT.h"
//----------------------------------------//
#define TRUE 1
#define FALSE 0


typedef struct unit_struct{
    int *state;
    int *local;
    int *vector;
    int *xy[2];
    double E;
} unit;

typedef struct const_values{
    double k;
    double T;
    double ikT;
} constant;

//-------------------- define func --------------------//
int init(unit *peptide,int state[],int seq_size);
int main_loop(unit *peptide,int total_step,constant constants,int seq_size);
int move(unit *peptide,int seq_size);
int compare_states(unit *peptide1,unit *peptide2,int seq_size);


int local2vector(unit *peptide,int seq_size);
int vector2xy(unit *peptide,int seq_size);
int setlocal(unit *peptide,int point,int value,int seq_size);

int copy(unit *tmp_peptide,unit *peptide,int seq_size);
int mc(unit *peptide1,unit *peptide2,constant constants,int seq_size);
int calc_energy(unit *peptide,int seq_size);

void myfree(unit *peptide);
//---------- move sets ----------//
int flip(unit *peptide,int point,int type,int seq_size);
int cornerflip(unit *peptide,int point,int seq_size);
int searchcorner(unit *peptide,int seq_size);

int show(unit peptide[],int seq_size);
//----------------------------------------//
int main(int argv,char *argc[]){
    //int seq[]={1,0,0,1,0,0,1,1};
    int seq[]={1,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,1,1};
    //int seq[]={1,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,1};
    /* int seq[]={0,0,1,1,1,0,1,1,1,1, */
    /*            1,1,1,1,0,0,0,1,1,1, */
    /*            1,1,1,1,1,1,1,0,1,0, */
    /*            0,0,1,1,1,1,1,1,1,1, */
    /*            1,1,1,1,0,0,0,0,1,1, */
    /*            1,1,1,1,0,1,1,0,1,0}; */
    
    int total_step=atoi(argc[1]);
    int seq_size=sizeof(seq)/sizeof(int);
    constant constants;
    unit *peptide;
    clock_t start,end;
    
    start=clock();
    
    srand(2);

    constants.k=1;
    constants.T=2.0;
    constants.ikT=1/(constants.k*constants.T);
    peptide=(unit *)malloc(sizeof(unit));
    
    
    //test(peptide,seq,seq_size,total_step,constants);

    init(peptide,seq,seq_size);
    main_loop(peptide,total_step,constants,seq_size);
    printf(">> %4d\n",calc_energy(peptide,seq_size));
    end=clock();
    printf("%8.3lf\n",(double)(end-start)/CLOCKS_PER_SEC);

    myfree(peptide);
    return 0;
}

void myfree(unit *peptide){
    free(peptide->state);
    free(peptide->local);
    free(peptide->vector);
    free(peptide->xy[0]);
    free(peptide->xy[1]);
    free(peptide);
}
//----------------------------------------//
int init(unit *peptide,int state[],int seq_size){
    int i;
    peptide->state=(int *)malloc(sizeof(int)*seq_size);
    peptide->local=(int *)malloc(sizeof(int)*seq_size);
    peptide->vector=(int *)malloc(sizeof(int)*seq_size);
    peptide->xy[0]=(int *)malloc(sizeof(int)*seq_size);
    peptide->xy[1]=(int *)malloc(sizeof(int)*seq_size);
    peptide->E=0;
    
    for(i=0;i<seq_size;i++){
        peptide->state[i]=state[i];
        peptide->xy[0][i]=0;
        peptide->xy[1][i]=i;
        peptide->local[i]=0;
        peptide->vector[i]=0;
    }

    return TRUE;
}

int main_loop(unit *peptide,int total_step,constant constants,int seq_size){
    int istep;
    unit *tmp_peptide;
    init_genrand(10);
    
    tmp_peptide=(unit *)malloc(sizeof(unit));
    init(tmp_peptide,peptide->state,seq_size);

    for(istep=0;istep<total_step;istep++){
        copy(tmp_peptide,peptide,seq_size);
        move(tmp_peptide,seq_size);
        if(TRUE!=compare_states(tmp_peptide,peptide,seq_size))
            mc(peptide,tmp_peptide,constants,seq_size);
        if(istep%(total_step/100)==0){
            constants.T=exp(-3*constants.T*istep/total_step);
            constants.ikT=1/(constants.T*constants.k);
        }
    }

    myfree(tmp_peptide);
    show(peptide,seq_size);
    return TRUE;
}

int compare_states(unit *peptide1,unit *peptide2,int seq_size){
    int i;

    for(i=2;i<seq_size;i++)
        if(peptide1->local[i]!=peptide2->local[i])
            return FALSE;
    return TRUE;
}

int mc(unit *peptide1,unit *peptide2,constant constants,int seq_size){
    int E1,E2;
    double prob;
    E1=calc_energy(peptide1,seq_size);
    E2=calc_energy(peptide2,seq_size);
    
    if(E2<E1){
        copy(peptide1,peptide2,seq_size);
    }else{
        prob=genrand_real3();
        if(prob<exp(-(E2-E1)*constants.ikT)){
            copy(peptide1,peptide2,seq_size);
        }
    }

    return TRUE;
}

int calc_energy(unit *peptide,int seq_size){
    int energy=0;
    int i,j;
    int x,y;

    for(i=0;i<seq_size;i++)
        for(j=i+1;j<seq_size;j++)
            if(peptide->xy[0][i]==peptide->xy[0][j] && peptide->xy[1][i]==peptide->xy[1][j])
                return 10000;

    for(i=0;i<seq_size;i++){
        if(peptide->state[i]==1){
            for(j=i+3;j<seq_size;j++){
                if(peptide->state[j]==1){
                    x=peptide->xy[0][i]-peptide->xy[0][j];
                    y=peptide->xy[1][i]-peptide->xy[1][j];
                    if(1==(x*x+y*y)){
                        energy--;
                    }
                }
            }
        }
    }
    return energy;
}

int copy(unit *tmp_peptide,unit *peptide,int seq_size){
    int i=0;
    for(i=0;i<seq_size;i++)
        tmp_peptide->state[i]=peptide->state[i];
    for(i=0;i<seq_size;i++)    
        tmp_peptide->local[i]=peptide->local[i];
    for(i=0;i<seq_size;i++)
        tmp_peptide->vector[i]=peptide->vector[i];
    for(i=0;i<seq_size;i++)
        tmp_peptide->xy[0][i]=peptide->xy[0][i];
    for(i=0;i<seq_size;i++)
        tmp_peptide->xy[1][i]=peptide->xy[1][i];
    return TRUE;
}


int local2vector(unit *peptide,int seq_size){
    int i;
    for(i=2;i<seq_size;i++){
        if(peptide->local[i]==0){
            peptide->vector[i-1]=peptide->vector[i-1];
        }else if(peptide->local[i]==1){
            peptide->vector[i]=(peptide->vector[i-1]+3)%12;
        }else{
            peptide->vector[i]=(peptide->vector[i-1]+9)%12;
        }
    }
    return TRUE;
}

int vector2xy(unit *peptide,int seq_size){
    int i;
    for(i=2;i<seq_size;i++){
        if(peptide->vector[i]==0){
            peptide->xy[0][i]=peptide->xy[0][i-1];
            peptide->xy[1][i]=peptide->xy[1][i-1]+1;
        }else if(peptide->vector[i]==3){
            peptide->xy[0][i]=peptide->xy[0][i-1]+1;
            peptide->xy[1][i]=peptide->xy[1][i-1];
        }else if(peptide->vector[i]==6){
            peptide->xy[0][i]=peptide->xy[0][i-1];
            peptide->xy[1][i]=peptide->xy[1][i-1]-1;
        }else{
            peptide->xy[0][i]=peptide->xy[0][i-1]-1;
            peptide->xy[1][i]=peptide->xy[1][i-1];
        }
    }
    return TRUE;
    
}

int setlocal(unit *peptide,int point,int value,int seq_size){
    peptide->local[point]+=value;
    if(peptide->local[point]<=-2)
        peptide->local[point]=-1;
    else if(peptide->local[point]>=2)
        peptide->local[point]=1;
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
        if(leftright==peptide->local[point])
            leftright*=-1;
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
int searchcorner(unit *peptide,int seq_size){
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

int cornerflip(unit *peptide,int point,int seq_size){
    if(peptide->local[point]==1){
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



int flip(unit *peptide,int point,int type,int seq_size){
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


int show(unit *peptide,int seq_size){
    int i;
    for(i=0;i<seq_size;i++){
        printf(">>> %3d %3d %3d\n",peptide->state[i],peptide->xy[0][i],peptide->xy[1][i]);
    }
    return TRUE;
}

