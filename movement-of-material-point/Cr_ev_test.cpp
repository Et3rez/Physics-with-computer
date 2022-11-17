#include <iostream>
#include <math.h>
#include <stdio.h>
const int maxi=10;
using namespace std;
typedef double T_FL;
typedef T_FL C_T_A1[maxi+1]; // element [0] not used!

typedef struct 
{ 
    int neq; // number of equations
    bool l_rak, l_op; // true: rocket, drag on
    T_FL g,sm,fm,w,c,b,f_rak,am,t,deg;

    //  g [m^2]: gravitational acceleration
    //  sm, fm [kg]: starting mass, final mass
    //  w [m/s]: rocket gases speed
    //  c [kg/s]: fuel burning speed
    //  drag: -b*v^2, b [kg/m]
    //  f_rak=w*c [N]: rocket engine force
    //  am=sm-c*t [kg]: actual mass
    //  t [s]: time
} C_T_PAR;

typedef struct 
{
    int isel, nrec; // case selector, file recording parameter
    T_FL beta, h, x0, y0, v0x, v0y; // numerical & starting data
    bool l_elr, l_ver, l_up, l_dn; // flags
    FILE *nfl; // pointer to logical file
} C_T_START;

typedef struct 
{
    T_FL ro, x_anal, y_anal, vx_anal, vy_anal, erx, ery, evx, evy;
}C_T_ERR_ANAL; // analytical values & errors

// Cr-init

void C_par_init(C_T_PAR &PAR)
{ T_FL p;
PAR.neq=2;
PAR.deg=M_PI/180.0;
do 
    {cout << "gravitational acceleration [m/s^2] (>0; 0-9.81):";
    cin >> PAR.g;
    } while (PAR.g<0.0);

}


int main ()
{
    std::cout << "Audi dulta, dic pauca!\n";
    std::cout << maxi;
    return 0;
}


