
#include<stdio.h>
#include "differential_flatness.h"


struct poly_coeff poly_coeff_x = {
    {0,0,0,0.296296296445966,-0.0296296296495693,0.000790123457588112,-5.88608758630694e-15}, //x
    {0,0,0.888888889337897,-0.118518518598277,0.00395061728794056,-3.53165255178416e-14,0},  //dx
    {0,1.77777777867579,-0.355555555794831,0.0158024691517622,-1.76582627589208e-13,0,0},   //ddx
    {1.77777777867579,-0.711111111589663,0.0474074074552867,-7.06330510356833e-13,0,0,0},   //dddx
    {-0.711111111589663,0.0948148149105735,-2.11899153107050e-12,0,0,0,0}                   // ddddx
} ;

struct poly_coeff poly_coeff_y = {
    {0,0,0,0.148148148223096,1-0.0148148148247847,0.000395061728794061,-2.94313291698298e-15}, //y
    {0,0,0.444444444669289,-0.0592592592991389,0.00197530864397030,-1.76587975018979e-14,0},  //dy
    {0,0.888888889338579,-0.177777777897417,0.00790123457588121,-8.82939875094893e-14,0,0},   //ddy
    {0.888888889338579,-0.355555555794833,0.0237037037276436,-3.53175950037957e-13,0,0,0},   //dddy
    {-0.355555555794833,0.0474074074552873,-1.05952785011387e-12,0,0,0,0}                   // ddddy,
} ;

struct poly_coeff poly_coeff_z = {
    {0,0,0,0.0296296296296264,-0.00296296296296231,7.90123456789687e-05,9.68196804005007e-19}, //z
    {0,0,0.0888888888888791,-0.0118518518518492,0.000395061728394844,5.80918082403004e-18,0},  //dz
    {0,0.177777777777758,-0.0355555555555477,0.00158024691357937,2.90459041201502e-17,0,0},   //ddz
    {0.177777777777758,-0.0711111111110954,0.00474074074073812,1.16183616480601e-16,0,0,0},   //dddz
    {-0.0711111111110954,0.00948148148147625,3.48550849441802e-16,0,0,0,0}                   // ddddz
} ;

struct poly_coeff poly_coeff_psi = {
    {0,0,0,0,0,0,0}, //psi
   {0,0,0,0,0,0,0},  //dpsi
    {0,0,0,0,0,0,0},   //ddpsi
    {0,0,0,0,0,0,0},   //dddpsi
   {0,0,0,0,0,0,0}               // ddddpsi
} ;
float t_vector[DISCSIZE];
struct discretized_poly x_disc;
struct discretized_poly y_disc;
struct discretized_poly z_disc;
struct discretized_poly psi_disc;

float t[3];
float tt;
float B1;
float C1;
float D1;
float A2;
float C2;
float D2;
float B3;
float C3;
float D3; 
float E1;
float E2;
float E3;
// float omega_x;
// float omega_y;
// float omega_z; 

// float omega_x_dot;
// float omega_y_dot;
// float omega_z_dot; 

float omega[3];
float omega_dot[3];

float term1[3];
float term2[3];
float term3[3];
float term4[3];

float a[3];
float j[3];
float s[3];
float c;
float c_dot;
float tempvec[3];
float z_b[3];
float x_b[3];
float x_c[3];
float y_b[3];
float y_b_negative[3];
float y_c_negative[3];
float y_c[3];
float u[4][DISCSIZE];
float psi;
int main(){
    
    //------Discretize polynomials-----------: There must be better way to do this? 
    discretize_time(0,t_vector,dt);   

    discretize_poly(t_vector,x_disc.x,      poly_coeff_x.x); 
    discretize_poly(t_vector,x_disc.dx,     poly_coeff_x.dx); 
    discretize_poly(t_vector,x_disc.ddx,    poly_coeff_x.ddx); 
    discretize_poly(t_vector,x_disc.dddx,   poly_coeff_x.dddx);
    discretize_poly(t_vector,x_disc.ddddx,  poly_coeff_x.ddddx); 

    discretize_poly(t_vector,y_disc.x,      poly_coeff_y.x); 
    discretize_poly(t_vector,y_disc.dx,     poly_coeff_y.dx); 
    discretize_poly(t_vector,y_disc.ddx,    poly_coeff_y.ddx); 
    discretize_poly(t_vector,y_disc.dddx,   poly_coeff_y.dddx);
    discretize_poly(t_vector,y_disc.ddddx,  poly_coeff_y.ddddx); 

    discretize_poly(t_vector,z_disc.x,      poly_coeff_y.x); 
    discretize_poly(t_vector,z_disc.dx,     poly_coeff_y.dx); 
    discretize_poly(t_vector,z_disc.ddx,    poly_coeff_y.ddx); 
    discretize_poly(t_vector,z_disc.dddx,   poly_coeff_y.dddx);
    discretize_poly(t_vector,z_disc.ddddx,  poly_coeff_y.ddddx); 

    discretize_poly(t_vector,psi_disc.x,      poly_coeff_psi.x); 
    discretize_poly(t_vector,psi_disc.dx,     poly_coeff_psi.dx); 
    discretize_poly(t_vector,psi_disc.ddx,    poly_coeff_psi.ddx); 
    discretize_poly(t_vector,psi_disc.dddx,   poly_coeff_psi.dddx);
    discretize_poly(t_vector,psi_disc.ddddx,  poly_coeff_psi.ddddx); 

    int i;
    for(i=0;i<DISCSIZE;i++){
            t[0]=x_disc.ddx[i]; //why can't i do t={x_disc.ddx[i],y_disc.ddx[i],z_disc.ddx[i]}?
            t[1]=y_disc.ddx[i];
            t[2]=z_disc.ddx[i]+data.g;
            c=norm2(t); //magnitude of acceleration (=mass normalized thrust when no drag is included)

            //calculate z_b
            z_b[0]=t[0]/c;
            z_b[1]=t[1]/c;
            z_b[2]=t[2]/c; //body z-axis

            psi=psi_disc.x[i];
            x_c[0]=cosf(psi);
            x_c[1]=sinf(psi);
            x_c[2]=0;
            y_c[0]=-sinf(psi);
            y_c[1]=cosf(psi); 
            y_c[2]=0;
            cross(z_b,x_c,y_b); //calculate y_b
            vectordividescalar(y_b,norm2(y_b),y_b);

            vectormultiplyscalar(y_b,-1,y_b_negative);
            vectormultiplyscalar(y_c,-1,y_c_negative);
            cross(y_b,z_b,x_b);//calculate x_b;


            u[0][i]=data.m*c; //desired thrust
            j[0]=x_disc.dddx[i]; //jerk
            j[1]=y_disc.dddx[i];
            j[2]=z_disc.dddx[i];

            s[0]=x_disc.ddddx[i]; //snap
            s[1]=y_disc.ddddx[i];
            s[2]=z_disc.ddddx[i];

            //working toward body moments as described in http://arxiv.org/abs/1712.02402
            B1=c;
            C1=0;
            D1=dotproduct(x_b,j);
            A2=c;
            C2=C1;
            D2=dotproduct(y_b_negative,j);            
            B3=dotproduct(y_c_negative,z_b);
            cross(y_c,z_b,tempvec);
            C3=norm2(tempvec);
            D3=psi_disc.dx[i]*dotproduct(x_c,x_b);

            omega[0]=(-B1*C2*D3+B1*C3*D2-B3*C1*D2+B3*C2*D1)/(A2*(B1*C3-B3*C1));
            omega[1]=(-C1*D3+C3*D1)/(B1*C3-B3*C1);
            omega[2]=(B1*D3-B3*D1)/(B1*C3-B3*C1);
            c_dot=dotproduct(z_b,j);

            E1=dotproduct(x_b,s)-2*c_dot*omega[1]-c*omega[0]*omega[2];
            E2=dotproduct(y_b_negative,s)-2*c_dot*omega[0]+c*omega[1]*omega[2];
            E3=psi_disc.ddx[i]*dotproduct(x_c,x_b)+2*psi_disc.dx[i]*omega[2]*dotproduct(x_c,y_b)-2*psi_disc.dx[i]*omega[1]*dotproduct(x_c,z_b)-omega[0]*omega[1]*dotproduct(y_c,y_b)-omega[0]*omega[2]*dotproduct(y_c,z_b);

            omega_dot[0]=(-B1*C2*E3+B1*C3*E2-B3*C1*E2+B3*C2*E1)/(A2*(B1*C3-B3*C1));
            omega_dot[1]=(-C1*E3+C3*E1)/(B1*C3-B3*C1);
            omega_dot[2]=(B1*E3-B3*E1)/(B1*C3-B3*C1);

            matrixvectormultiply(data.I,omega_dot,term1);  //calculate term1: J*omega_dot (not dot product)
            matrixvectormultiply(data.I,omega,term2);      //calculate  term2: J*omega
            cross(omega,term2,term3);                      //calculate term3: omega x J*omega (cross)
            u[1][i]=term1[0]+term3[0]; // moment around body x-axis
            u[2][i]=term1[1]+term3[1];// moment around body y-axis
            u[3][i]=term1[2]+term3[2];// moment around body z-axis
    }

   int dummy=2; //for debug brakepoint
}

