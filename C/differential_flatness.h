#include <math.h>
#define MAX_POLY_ORDER 7 //max poly order of y=x0+x1*t+x2*t^2+xn*t^n is n+1

struct poly_coeff{
    float x[MAX_POLY_ORDER];
    float dx[MAX_POLY_ORDER];
    float ddx[MAX_POLY_ORDER];
    float dddx[MAX_POLY_ORDER]; 
    float ddddx[MAX_POLY_ORDER];
};

float t_m = 15;
float dt =0.01;

#define DISCSIZE 1500

struct discretized_poly{
    float x[DISCSIZE];
    float dx[DISCSIZE];
    float ddx[DISCSIZE];
    float dddx[DISCSIZE]; 
    float ddddx[DISCSIZE];
};

struct droneparameters{
    float m;
    float beta;
    float g; 
    float L;
    float I[3][3];
    float k_F;
    float k_M;
    float tau_g[3];
    float D[3][3];
    float A[3][3];
};

struct droneparameters data = {
    0.38905,
    0.5,
    9.81,
    0.2,
    {{0.0049,0,0},{0,0.0049,0},{0,0,0.0069}},
    1.91e-6,
    2.7e-7,
    {0,0,0},
    {{0,0,0},{0,0,0},{0,0,0}},
    {{0,0,0},{0,0,0},{0,0,0}}
};

void discretize_time(float t0, float t_vector[DISCSIZE],float dt){
    t_vector[0]=t0;
    int i=1;
    for(i=1;i<DISCSIZE;i++){
        t_vector[i]=t_vector[i-1]+dt;
    }
}
void discretize_poly(float t_vector[DISCSIZE], float y[DISCSIZE], float poly_coeff[MAX_POLY_ORDER]){
    
    int i ;
    int j ;
    for (i=0;i<DISCSIZE;i++){
        y[i]=0;
        for (j=0;j<MAX_POLY_ORDER;j++){
            y[i]=y[i]+poly_coeff[j]*pow(t_vector[i],j);            
        }
    }
}

void cross(float A[3],float B[3],float out[3]){
    out[0]=A[1]*B[2]-A[2]*B[1];
    out[1]=-1*(A[0]*B[2]-A[2]*B[0]);
    out[2]=A[0]*B[1]-A[1]*B[0];
}

float norm2(float A[3]){
    int i;
    float N=0;
    for(i=0;i<3;i++){
        N+=pow(A[i],2);
    }
    N=sqrtf(N);
    return N;
}

float dotproduct(float a[3], float b[3]){
    float c=0;
    int i;
    for(i=0;i<3;i++){
        c=c+a[i]*b[i];
    }
    return c; 
}

void vectordividescalar(float a[3],float b,float out[3]){
    int i;
    for(i=0;i<3;i++){
        out[i]=a[i]/b;
    }
}

void vectormultiplyscalar(float a[3],float b,float out[3]){
    int i;
    for(i=0;i<3;i++){
        out[i]=a[i]*b;
    }
}

void matrixvectormultiply(float A[3][3],float B[3],float out[3]){
    int i;
    int j;
    out[0]=0;
    out[1]=0;
    out[2]=0;

    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            out[i]=out[i]+A[i][j]*B[j];
        }
    }

}