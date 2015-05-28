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

int local2vector(unit peptide[],int point,int seq_size);
int vector2xy(unit peptide[],int point,int seq_size);
int setlocal(unit peptide[],int point,int value,int seq_size);

//---------- move sets ----------//
int rigidrotation(unit peptide[],int point,int type,int seq_size);
int pointflip(unit peptide[],int point,int type,int seq_size);
int conerflip(unit peptide[],int point,int type,int seq_size);


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

int local2vector(unit peptide[],int point,int seq_size){
    int i;
    if(peptide[point].local==1)
        for(i=point;i<seq_size;i++)
            peptide[i].vector=(peptide[i].vector+3)%12;
    else
        for(i=point;i<seq_size;i++)
            peptide[i].vector+=(peptide[i].vector+9)%12;
    return TRUE;
}

int vector2xy(unit peptide[],int point,int seq_size){
    int i;
    int tmpxy[2];
    int basexy[2];
    basexy[0]=peptide[point-1].xy[0];
    basexy[1]=peptide[point-1].xy[1];
    printf("%d,%d\n",basexy[0],basexy[1]);
    if(peptide[point].vector==3){
        if(peptide[point-1].vector==0){
            for(i=point;i<seq_size;i++){
                //    _    
                //   |     [cos(-90) -sin(-90)][x] _ [x*cos(-90)-y*sin(-90)]
                //   |     [sin(-90)  cos(-90)][y] - [x*sin(-90)+y*cos(-90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]+tmpxy[1];
                peptide[i].xy[1] = basexy[1]-tmpxy[0];
            }
        }else if(peptide[point-1].vector==6){
            for(i=point;i<seq_size;i++){
                //         
                //   |     [cos(90) -sin(90)][x] _ [x*cos(90)-y*sin(90)]
                //   L     [sin(90)  cos(90)][y] - [x*sin(90)+y*cos(90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]-tmpxy[1];
                peptide[i].xy[1] = basexy[1]+tmpxy[0];
            }
        }
        return TRUE;
    }else if(peptide[point].vector == 6){
        if(peptide[point-1].vector==3){
            for(i=point;i<seq_size;i++){
                //        
                //  ____   [cos(-90) -sin(-90)][x] _ [x*cos(-90)-y*sin(-90)]
                //      |  [sin(-90)  cos(-90)][y] - [x*sin(-90)+y*cos(-90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]+tmpxy[1];
                peptide[i].xy[1] = basexy[1]-tmpxy[0];
            }
        }else if(peptide[point-1].vector==9){
            for(i=point;i<seq_size;i++){
                //         
                //  _____  [cos(90) -sin(90)][x] _ [x*cos(90)-y*sin(90)]
                // |       [sin(90)  cos(90)][y] - [x*sin(90)+y*cos(90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]-tmpxy[1];
                peptide[i].xy[1] = basexy[1]+tmpxy[0];
            }
        }
        return TRUE;
    }else if(peptide[point].vector == 9){
        if(peptide[point-1].vector==0){
            for(i=point;i<seq_size;i++){
                //  _      
                //   |     [cos(90) -sin(90)][x] _ [x*cos(90)-y*sin(90)]
                //   |     [sin(90)  cos(90)][y] - [x*sin(90)+y*cos(90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]-tmpxy[1];
                peptide[i].xy[1] = basexy[1]+tmpxy[0];
            }
        }else if(peptide[point-1].vector==6){
            for(i=point;i<seq_size;i++){
                //   |      
                //   |     [cos(-90) -sin(-90)][x] _ [x*cos(-90)-y*sin(-90)]
                //  -      [sin(-90)  cos(-90)][y] - [x*sin(-90)+y*cos(-90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]+tmpxy[1];
                peptide[i].xy[1] = basexy[1]-tmpxy[0];
            }
        }
        return TRUE;
    }else{
        if(peptide[point-1].vector==0){
            for(i=point;i<seq_size;i++){
                //        
                //        [cos(90) -sin(90)][x] _ [x*cos(90)-y*sin(90)]
                //  ____| [sin(90)  cos(90)][y] - [x*sin(90)+y*cos(90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]-tmpxy[1];
                peptide[i].xy[1] = basexy[1]+tmpxy[0];
            }
        }else if(peptide[point-1].vector==6){
            for(i=point;i<seq_size;i++){
                //         
                //        [cos(-90) -sin(-90)][x] _ [x*cos(-90)-y*sin(-90)]
                //  |____ [sin(-90)  cos(-90)][y] - [x*sin(-90)+y*cos(-90)]
                //
                tmpxy[0] = peptide[i].xy[0]-basexy[0];
                tmpxy[1] = peptide[i].xy[1]-basexy[1];
                peptide[i].xy[0] = basexy[0]+tmpxy[1];
                peptide[i].xy[1] = basexy[1]-tmpxy[0];
            }
        }
        return TRUE;
    }
    return TRUE;
}

int setlocal(unit peptide[],int point,int value,int seq_size){
    peptide[point].local=value;
    return TRUE;
}

int move(unit peptide[],int seq_size){
    return TRUE;
}


//-------------------- move sets --------------------//

int conerflip(unit peptide[],int point,int type,int seq_size){
    return TRUE;
}

int pointflip(unit peptide[],int point,int type,int seq_size){
    // type=1(right),type=0(left)
    if(type==1){ 
        setlocal(peptide,point,1,seq_size);
        local2vector(peptide,point,seq_size);
        vector2xy(peptide,point,seq_size);
    }else{
        setlocal(peptide,point,-1,seq_size);
        local2vector(peptide,point,seq_size);
        vector2xy(peptide,point,seq_size);
    }
    return TRUE;
}

int rigidrotation(unit peptide[],int point,int type,int seq_size){

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

    pointflip(peptide,2,0,seq_size);
    printf("pointflip\n");
    show(peptide,seq_size);
    
    init(peptide,seq,seq_size);
    pointflip(peptide,3,1,seq_size);
    pointflip(peptide,5,1,seq_size);
    printf("conerflip\n");
    show(peptide,seq_size);


    
    
    return TRUE;
}



int show(unit peptide[],int seq_size){
    int i;
    for(i=0;i<seq_size;i++){
        printf("%3dth,local=%3d,vector=%3d,xy=[%3d,%3d]\n",i,peptide[i].local,peptide[i].vector,peptide[i].xy[0],peptide[i].xy[1]);
    }
    return TRUE;
}

