
//Precession Period came out to be 340 years with my RK4.
//RK4 function is working but the acceleration part is wrong.




#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define NBODY 2
#define ME 5.96e24
#define MS 1.99e30
#define G  0.0000000000667
//#define G  6.67e-11
#define I1 8.01e37
#define I3 8.03e37
#define theta0 23.45
#define C0 2.0*3.14/86164.1 //not sure about this value

double pi = 3.141592653589793;
double mass[NBODY];

///res[] made here
//F function to update velocities and accerlations to compute K_2
//to be put in as a parameter to intergrate_rk2
void F(double *y,double t, double *res)
{

    double dx = y[5]-y[0];
    double dy = y[6]-y[1];
    double dz = y[7]-y[2];
    double d  = sqrt( dx*dx + dy*dy + dz*dz );
    double d3 = ( dx*sin(y[3]) - dy*cos(y[3]) )*sin(y[4]) + dz*cos(y[4]);
//    double d3 = ( dx - dy ) + dz;


  for(int i=0; i<5*NBODY; ++i) res[i] = y[5*NBODY+i];//store velocities in res[]-first 10 elements
  for(int i=0; i<NBODY; ++i){//store accelerations in res[]-last 10 elements
/*
      //ax
      res[5*NBODY+5*i+0] =  -(dx/d*d*d)*G*MS + (1/(ME*d*d*d*d*d))*3*G*MS*(I1-I3)*(  dx/2 - (5*d3*d3*dx)/(2*d*d) + d3*sin(y[3])*sin(y[4])  );

      //ay
      res[5*NBODY+5*i+1] = -(dy/d*d*d)*G*MS + (1/(ME*d*d*d*d*d)) *  3*G*MS*(I1-I3) * (  dy/2 - (5*d3*d3*dy)/(2*d*d)  + d3*cos(y[3])*sin(y[4])  );

      //az
      res[5*NBODY+5*i+2] = -(dz/d*d*d)*G*MS + (1/(ME*d*d*d*d*d)) *   3*G*MS*(I1-I3) * (  dz/2 - (5*d3*d3*dz)/(2*d*d)  + d3*cos(y[4])  );

      //aphi
      res[5*NBODY+5*i+3] = (1/sin(y[4]))*( -2*y[13]*y[14]*cos(y[4]) + (I3/I1)*C0*y[14] + (3*G*MS/pow(d,5))*(I1-I3)/(I1)*d3*( dx*cos(y[3]) + dy*sin(y[3])));

      //atheta
      res[5*NBODY+5*i+4] = y[13]*y[13]*sin(y[4])*cos(y[4]) - (I3/I1)*C0*y[13]*sin(y[4]) + (3*G*MS/d*d*d*d*d)*((I1-I3)/I1)*d3*(  (dx*sin(y[3]) - dy*cos(y[3]) )*cos(y[4]) - dz*sin(y[4]));
*/

      //ax
      res[5*NBODY+5*i+0] =  -(dx/d*d*d) ;

      //ay
      res[5*NBODY+5*i+1] = -(dy/d*d*d);

      //az
      res[5*NBODY+5*i+2] = -(dz/d*d*d) ;

      //aphi
      res[5*NBODY+5*i+3] = (1/sin(y[4]))*( -2*y[13]*y[14]*cos(y[4]) + (I3/I1)*C0*y[14] + (3*G*MS/pow(d,5))*(I1-I3)/(I1)*d3*( dx*cos(y[3]) + dy*sin(y[3])));

      //atheta
      res[5*NBODY+5*i+4] = y[13]*y[13]*sin(y[4])*cos(y[4]) ;

  }
}

void integrate_rk2(double *yold, int n, double t, double dt, double *ynew, void (*F)(double*, double,double*))
{
        //double deriv1[n]; dont use this notation!!it works but just avoid it
        double *deriv1 = (double*) malloc(n*sizeof(double)); //(double*) means whatever is made, it returns a double result
        double *yhalf = (double*) malloc(n*sizeof(double));

        F(yold,t,deriv1); //K_1
        for(int i = 0; i < n; ++i) yhalf[i] = yold[i]+deriv1[i]*dt/2;

        F(yhalf, t+dt/2, deriv1);//K_2
        for(int i = 0; i < n; ++i) yhalf[i] = yold[i]+deriv1[i]*dt/2;

        F(yhalf, t+dt/2, deriv1);//K_3
        for(int i = 0; i < n; ++i) yhalf[i] = yold[i]+deriv1[i]*dt/2;

        F(yhalf, t+dt/2, deriv1);//K_4
        for(int i = 0; i < n; ++i) ynew[i] = yold[i]+deriv1[i]*dt;

        free(deriv1);//release the memory you allocated. For your later conveiniences. When the code crashes, we want to see $
        free(yhalf);
}

int main(int argc, char** argv){
//earth initial position and velocity
  double r0= 149600000000 *(1+0.0167112303531389)*MS/(ME+MS);
//double r0= 0;
  double ke= (G*MS*ME/149600000000)*(1/(1+0.0167112303531389) - 1/2);
  double v0= sqrt((2*ke) / ME*(1+ME/MS));
//double v0= 0;

//store initial values for x,y,z,phi,theta, vx,vy,vz,vphi,vtheta for both sun and earth
  double y[10*NBODY] = {/*eartch*/r0, 0, 0, 0, 23.45, /*Sun*/ 0, 0, 0, 0, 0, /*earthspeed*/
                        v0, 0, 0, 2.0*3.14/86164.1, 0,
                        /*sunspeed*/ 0,0,0,0,0  };
  mass[0] = ME;
  mass[1] = MS;
  double t = 0;
  double dt =1;
  double tstop = 340;

/* Conservation of Angular momentum part. Do it LATER
//Getting the Velocities Store it to y[]
        double mom[2]={0,0};
        double mtot = 0;
        for(int i = 0; i<NBODY; i++)
        {
        mom[0] += mass[i]*y[2*NBODY+2*i+0];
        mom[1] += mass[i]*y[2*NBODY+2*i+1];
        mtot += mass[i];
        }
        for(int i = 0; i<NBODY; ++i)
        {
        y[2*NBODY+2*i+0] -= mom[0]/mtot;
        y[2*NBODY+2*i+1] -= mom[1]/mtot;
        }
*/

//PRINT INITIAL positions and velocities
  printf("%f ", t);
  for(int i = 0; i<NBODY*10; ++i) printf("%f ", y[i]);printf("\n");

//  printf("%.13f \n", G);
//  printf("%f \n", ME);
//  printf("%f \n", MS);
//  printf("%f \n", I3);
//  printf("%f \n", I1);

  int step = 0;
  for(; t<tstop-dt/2; t+=dt)
  {
    double ynew[10*NBODY];
    integrate_rk2(y, 10*NBODY, t, dt,ynew, F);
    if((++step)%10 == 0)
    {
      printf("%f ", t);
      for(int i = 0; i<NBODY*10; ++i)
      printf("%f ", y[i]);printf("\n");
    }
    for(int i=0; i<10*NBODY; ++i) y[i] = ynew[i];
  }
  return 0;
}



