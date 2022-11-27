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

//Cr-accel
void C_acc1(C_T_PAR &PAR, C_T_A1 A1, C_T_A1 R1, C_T_A1 V1)
// PAR - parameters structure
// R1  - position vector [m]
// V1  - velocity vector [m/s]
// A1  - returns acceleration vector [m/s^2]
{   T_FL v,zx,zy;

    A1[1]=0.0; //ax=0
    A1[2]=-PAR.g;
    PAR.am=PAR.sm;

    //end of free throw acceleration formula. 
    //Next step will be performed only for rocket and/or for motion with drag.

    if ((PAR.l_rak)||(PAR.l_op))
    {   v=sqrt(V1[1]*V1[1]+V1[2]*V1[2]); // velocity
        zx=V1[1]/v; zy=V1[2]/2; // velocity versors
        if (PAR.l_rak)
        {   PAR.am=PAR.sm-PAR.c*PAR.t; //actual mass
            if (PAR.am>=PAR.fm) // add rocket acceleration
            {   A1[1]=A1[1]+PAR.f_rak/PAR.am*zx;
                A1[2]=A1[2]+PAR.f_rak/PAR.am*zy;
            } 
            else;
            {   cout << "End of engine operation\n";
                PAR.l_rak=false;
            }
        } 
        else PAR.am=PAR.fm; // actual mass = final mass

        // current mass am for rocket has been calculated.
        // next step will be performed if drag is on.

        if (PAR.l_op) //subtract deceleration due to drag
        {   v=v*v*PAR.b/PAR.am;
            A1[1]=A1[1]-v*zx;
            A1[2]=A1[2]-v*zy;
        }
    } 
} //C_acc1

//C-ev Euler and Verlet integrators

void C_e_step(C_T_PAR &PAR,C_T_A1 ACC1,C_T_A1 ROLD1, C_T_A1 RNEW1,
                            C_T_A1 VOLD1,C_T_A1 VNEW1, T_FL H, T_FL BETA)

// Euler step: Rnew=r(t+h), Rold=r(t), Vnew=v(t+h), Vold=v(t)
{   C_acc1(PAR,ACC1,ROLD1,VOLD1); // acceleration defined by ROLD1,VOLD1
        for (int k=1;k<=PAR.neq;k++)
        {   VNEW1[k]=VOLD1[k]+ACC1[k]*H; // new velocity
            RNEW1[k]=ROLD1[k]+(BETA*VOLD1[k]+(1.0-BETA)*VNEW1[k])*H+0.5*ACC1[k]*H*H;
        }   // new position by weighted sum of olf and new velocities
} // C_e_step

void C_v_step(C_T_PAR &PAR,C_T_A1 ACC1,C_T_A1 ROLD1,C_T_A1 RNEW1,
                C_T_A1 ROLDER1,C_T_A1 VOLD1,C_T_A1 VNEV1, T_FL H, T_FL BETA)

// Verlet step: Rnew=r(t+h),Rold=r(t),Rolder=r(t-h),Vold=v(t)
// Vnew=V(t+h); Verlet starts in 2nd step; first do single Euler step
{   C_acc1(PAR,ACC1,ROLD1,VOLD1); // acceleration defined by ROLD1, VOLD1
    for (int k=1;k<=PAR.neq;k++)
    {   RNEW1[k]=(1.0+BETA)*ROLD1[k]-BETA*ROLDER1[k]+ACC1[k]*H*H; 
            // new position, with a possibility of damping
        VOLD1[k]=0.5*(RNEW1[k]-ROLDER1[k])/H; //old velocity
    }
    VNEW1[k]=VOLD1[k]+ACC1[k]*H; // new velocity
} // C_v_step


int main ()
{ C_T_PAR par;
    std::cout << "Audi multa, dic pauca!\n";
    C_par_init(par);
    return 0;
}


