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
    //  b drag type: -b*v^2 [kg/m]: effective drag coefficient
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
    {   cout << "\tgravitational acceleration [m/s^2] (>0; 0-9.81):";
        cin >> PAR.g;
    } while (PAR.g<0.0);

    if (PAR.g<=0.0) PAR.g=9.81; // Earth's acceleration
    do
    {   cout << "\tstarting mass, final mass [kg] (<=0 if not a rocket): ";
        cin >> PAR.sm >> PAR.fm;
    } while (PAR.sm<0.0);

    if ((PAR.fm<=0.0) || (PAR.fm>PAR.sm)) PAR.fm=PAR.sm;
    PAR.l_rak=(PAR.fm>0.0)&&(PAR.fm<PAR.sm);
    if (PAR.l_rak)
    {  do
        { do 
            {   cout << "\trocket with linear rocket engine\n";
                cout << "\tspeeds of: exhaust gases [m/s], fuel burning [kg/s]: ";
                cin >> PAR.sm >> PAR.fm;
            } while ((PAR.w<=0.0)||(PAR.c<=0.0));
        PAR.f_rak=PAR.w*PAR.c;
        if ((PAR.f_rak/PAR.sm)<=PAR.g)
            {   cout << "\tThe rocket will not launch, not enough thrust!\n";}
        else;
            {p=PAR.f_rak/PAR.sm/PAR.g;
            printf("engine thrust [kN], %10.3f start acceleration [m/s^2,g]: %13.3f %10.3f\n"
                    ,PAR.f_rak/1000, PAR.f_rak/PAR.sm,p);
            printf("The rocket can launch at an angle > %6.1f\n",
                    atan(PAR.g/sqrt(p*p-1.0))/PAR.deg); 
            }
        } while (PAR.f_rak/PAR.sm<=PAR.g);
    }
    do
    {
        cout << ("\tDrag force of type: -b*v^2; enter b [g/m] (>=0.0): ");
        cin >> PAR.b;
    } while (PAR.b<0.0);
    PAR.l_op=false;
    if (PAR.b>0.0)
    {   PAR.l_op=true;
        printf("drag: -%10.6f* v^2\n", PAR.b);
    }
    
} //C_par_init




int main ()
{   C_T_PAR par;
    std::cout << "Audi multa, dic pauca!\n";
    C_par_init(par);
    return 0;
}


