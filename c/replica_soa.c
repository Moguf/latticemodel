#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include "MT.h"
//----------------------------------------//
#define TRUE 1
#define FALSE 0

#define REPLICA_SIZE 10

typedef struct unit_struct{
    int *state;
    int *local;
    int *vector;
    int *xy[2];
    double ikT;
    double E;
} unit;

typedef struct const_values{
    int replica_size;
    int replica_step;
} constant;

//-------------------- define func --------------------//
int init(unit *peptides,int state[],int seq_size,int replica_size);

int main_loop(unit *peptides,int total_step,constant *constants,int seq_size);
int move(unit *peptide,int seq_size);
int compare_states(unit *peptide1,unit *peptide2,int seq_size);

int local2vector(unit *peptide,int seq_size);
int vector2xy(unit *peptide,int seq_size);
int setlocal(unit *peptide,int point,int value,int seq_size);

int copy(unit *tmp_peptide,unit *peptide,int seq_size);
int mc(unit *peptide1,unit *peptide2,constant *constants,int seq_size);
int calc_energy(unit *peptide,int seq_size);

void peptidesfree(unit *peptide,constant *constant);
void peptidefree(unit *peptide);

void test(unit *peptide,int total_step,constant *constants,int seq_size);
//---------- move sets ----------//
int flip(unit *peptide,int point,int type,int seq_size);
int cornerflip(unit *peptide,int point,int seq_size);
int searchcorner(unit *peptide,int seq_size);

int replica_exchange(unit *peptides,constant *constants,int seq_size);

int show(unit *peptide,int seq_size);
//----------------------------------------//
int main(int argv,char *argc[]){
    //int seq[]={1,0,0,1,0,0,1,1};
    //int seq[]={1,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,1,1};
    //int seq[]={1,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,1};
    int seq[]={0,0,1,1,1,0,1,1,1,1,
               1,1,1,1,0,0,0,1,1,1,
               1,1,1,1,1,1,1,0,1,0,
               0,0,1,1,1,1,1,1,1,1,
               1,1,1,1,0,0,0,0,1,1,
               1,1,1,1,0,1,1,0,1,0};
    
    int total_step=1000*1000*10;
    int seq_size=sizeof(seq)/sizeof(int);
    constant *constants;
    unit *peptides;

    //for time
    #ifdef _OPENMP
    clock_t start,end;
    start=clock();
    #else
    double start,end;
    start=omp_get_wtime();
    #endif

    //init_genrand(atoi(argc[2]));
    init_genrand(0);

    //initialize constants
    constants=(constant *)malloc(sizeof(constant));
    constants->replica_size=REPLICA_SIZE;
    constants->replica_step=100;

    //initialize peptides 
    peptides=(unit *)malloc(sizeof(unit)*constants->replica_size);
    init(peptides,seq,seq_size,constants->replica_size);

    main_loop(peptides,total_step,constants,seq_size);
    
    //for time
    #ifdef _OPENMP
    end=clock();
    printf(">%8.3lf\n",(double)(end-start)/CLOCKS_PER_SEC);
    #else
    end=omp_get_wtime();
    printf(">%8.3lf\n",end-start);
    #endif
    
    peptidesfree(peptides,constants);
    free(constants);

    return 0;
}

void peptidesfree(unit *peptides,constant *constants){
    int ireplica;    

    for(ireplica=0;ireplica<constants->replica_size;ireplica++){
        free(peptides[ireplica].state);
        free(peptides[ireplica].local);
        free(peptides[ireplica].vector);
        free(peptides[ireplica].xy[0]);
        free(peptides[ireplica].xy[1]);
    }
    free(peptides);
}

void peptidefree(unit *peptide){

    free(peptide->state);
    free(peptide->local);
    free(peptide->vector);
    free(peptide->xy[0]);
    free(peptide->xy[1]);
    free(peptide);
}


//----------------------------------------//
int init(unit *peptides,int state[],int seq_size,int replica_size){
    int i,ireplica;

    for(ireplica=0;ireplica<replica_size;ireplica++){
        peptides[ireplica].state=(int *)malloc(sizeof(int)*seq_size);
        peptides[ireplica].local=(int *)malloc(sizeof(int)*seq_size);
        peptides[ireplica].vector=(int *)malloc(sizeof(int)*seq_size);
        peptides[ireplica].xy[0]=(int *)malloc(sizeof(int)*seq_size);
        peptides[ireplica].xy[1]=(int *)malloc(sizeof(int)*seq_size);
        peptides[ireplica].ikT=1/(0.1*ireplica+0.1);
        peptides[ireplica].E=0;

        for(i=0;i<seq_size;i++){
            peptides[ireplica].state[i]=state[i];
            peptides[ireplica].xy[0][i]=0;
            peptides[ireplica].xy[1][i]=i;
            peptides[ireplica].local[i]=0;
            peptides[ireplica].vector[i]=0;
        }

    }
    return TRUE;
}

int main_loop(unit *peptides,int total_step,constant *constants,int seq_size){
    int istep;
    int ireplica;
    int i;
    int replica_size=constants->replica_size;
    unit *tmp_peptides,*stable_peptides;
    int minE=0;
    int minE_index=0;

    tmp_peptides=(unit *)malloc(sizeof(unit)*replica_size);
    stable_peptides=(unit *)malloc(sizeof(unit)*replica_size);
    init(tmp_peptides,peptides[0].state,seq_size,replica_size);
    init(stable_peptides,peptides[0].state,seq_size,replica_size);

    #pragma omp parallel private(istep,ireplica) num_threads(replica_size)
    {
        ireplica=omp_get_thread_num();
        for(istep=0;istep<total_step;istep++){
            copy(&tmp_peptides[ireplica],&peptides[ireplica],seq_size);
            move(&tmp_peptides[ireplica],seq_size);

            if(TRUE!=compare_states(&tmp_peptides[ireplica],&peptides[ireplica],seq_size))
                mc(&peptides[ireplica],&tmp_peptides[ireplica],constants,seq_size);

            if(istep%constants->replica_step==0){
                #pragma omp barrier
                #pragma omp master
                {
                    replica_exchange(peptides,constants,seq_size);
                    for(i=0;i<replica_size;i++){
                        if(stable_peptides[i].E>peptides[i].E)
                            copy(&stable_peptides[i],&peptides[i],seq_size);
                    }
                }
                #pragma omp barrier
            }
        }
    }

    minE=stable_peptides[0].E;
    
    for(i=0;i<replica_size;i++){
        printf("%dth E=%lf\n",i,stable_peptides[i].E);
        if(stable_peptides[i].E<minE){
            minE=stable_peptides[i].E;
            minE_index=i;
        }
    }

    show(&stable_peptides[minE_index],seq_size);
    peptidesfree(stable_peptides,constants);

    return TRUE;
}
int replica_exchange(unit *peptides,constant *constants,int seq_size){
    int i;
    double prob;
    double tmp_ikT=0;

    for(i=0;i<constants->replica_size-1;i++){
        tmp_ikT=0;
        prob=genrand_real3();
        if(prob<(peptides[i].E-peptides[i+1].E)*(peptides[i].ikT-peptides[i+1].ikT)){
            tmp_ikT=peptides[i].ikT;
            peptides[i].ikT=peptides[i+1].ikT;
            peptides[i+1].ikT=tmp_ikT;
        }
    }

    return TRUE;
}


int compare_states(unit *peptide1,unit *peptide2,int seq_size){
    int i;

    for(i=2;i<seq_size;i++)
        if(peptide1->local[i]!=peptide2->local[i])
            return FALSE;
    
    return TRUE;
}

int mc(unit *peptide1,unit *peptide2,constant *constants,int seq_size){
    int E1,E2;
    double prob;
    E1=calc_energy(peptide1,seq_size);
    E2=calc_energy(peptide2,seq_size);
    
    if(E2<E1){
        copy(peptide1,peptide2,seq_size);
    }else{
        prob=genrand_real3();
        if(prob<exp(-(E2-E1)*peptide1->ikT)){
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
    
    for(i=0;i<seq_size-3;i++){
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
    peptide->E=energy;

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

    tmp_peptide->ikT=peptide->ikT;
    tmp_peptide->E=peptide->E;

    return TRUE;
}


int local2vector(unit *peptide,int seq_size){
    int i;
    int j;

    for(i=2;i<seq_size;i++)
        peptide->vector[i]=0;
    
    for(i=2;i<seq_size;i++){
        if(peptide->local[i]==1){
            for(j=i;j<seq_size;j++)
                peptide->vector[j]=(peptide->vector[j]+3)%12;
        }else if(peptide->local[i]==-1){
            for(j=i;j<seq_size;j++)
                peptide->vector[j]=(peptide->vector[j]+9)%12;
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

int move(unit *peptide,int seq_size){
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
        // search corner.
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
    
    for(i=start;i<seq_size;i++)
        if(peptide->local[i]!=0){
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
            if(peptide->local[i]!=0){
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
void test(unit *peptide,int total_step,constant *constants,int seq_size){

    flip(peptide,3,1,seq_size);
    flip(peptide,4,0,seq_size);
    flip(peptide,4,1,seq_size);
    flip(peptide,4,1,seq_size);
    flip(peptide,7,1,seq_size);
    show(peptide,seq_size);
}

int show(unit *peptide,int seq_size){
    int i;

    for(i=0;i<seq_size;i++){
        //printf(">>> %3d %3d %3d\n",peptide->state[i],peptide->xy[0][i],peptide->xy[1][i]);
        printf(">>>%3d %3d %3d %3d %3d\n",peptide->state[i],peptide->local[i],peptide->vector[i],peptide->xy[0][i],peptide->xy[1][i]);
    }
    printf(">> %4d,T=%lf\n",calc_energy(peptide,seq_size),1/peptide->ikT);

    return TRUE;
}

