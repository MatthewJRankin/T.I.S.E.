#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Energy(double*, double*);
void Potential(double*,double*);

int main(void) {

double psi,psi_old,psi_squared,psi_Norm,omega,E,E_Analytical,V,x,A_psi,B_psi,C_psi,D_psi,A_ome,B_ome,C_ome,D_ome;
double k=20000,Sum,Simpson,A_Norm,xtemp;
const double pi =3.14159;
const double h = 6.62607e-34;
const double h_bar = h/(2*pi);
double m=9.10938e-31;
double dx=1e-13;

double n;
int flag;

FILE*file_pointer;

file_pointer=fopen("ElectronHarmonicWellNorm.txt", "w");

printf("Enter desired number of energy levels:\n");
scanf("%lf",&n);

for(double i=0;i<=n;i++)
    {
                double E_max=1e-14;
                double E_min=1e-24;
                Energy(&E_Analytical,&i);

           do{
                    int n_crossings=0;
                    omega=1;
                    psi=0;
                    E=0.5*(E_max+E_min);
                    flag=0;
                    Sum=0;
                    double L=sqrt(2*(E+1.6e-16)/k);

                    for(x=-L;x<=L;x+=dx)
                     {
                              int j=0;
                              psi_old=psi;
                              Potential(&V,&x);

                              A_psi = omega;
                              A_ome = (-2*m/(h_bar*h_bar))*(E-V)*psi;

                              xtemp=x+dx/2;
                              Potential(&V,&xtemp);

                              B_psi = omega + (dx/2)*A_ome;
                              B_ome = (-2*m/(h_bar*h_bar))*(E-V)*(psi+(dx/2)*A_psi);

                              C_psi = omega + (dx/2)*B_ome;
                              C_ome = (-2*m/(h_bar*h_bar))*(E-V)*(psi+(dx/2)*B_psi);

                              xtemp=x+dx;
                              Potential(&V,&xtemp);

                              D_psi = (omega + dx*C_ome);
                              D_ome = (-2*m/(h_bar*h_bar))*(E-V)*(psi+dx*C_psi);

                              psi = psi +(dx/6)*(A_psi+2*B_psi+2*C_psi+D_psi);
                              omega = omega+ (dx/6)*(A_ome+2*B_ome+2*C_ome+D_ome);

                              if ((psi_old*psi)<0) n_crossings++;

                              if (n_crossings>i)
                                {
                                 flag=1;
                                }
                                psi_squared=psi*psi;
                                if(j%2==0)
                                    {
                                        if (j<(2*L/dx)) Sum+=2*psi_squared;
                                        else Sum+=psi_squared;
                                    }
                                else Sum+=4*psi_squared;
                                j++;
                     }
            if (flag==1) {E_max=E;}
                    else {E_min=E;}
                    //printf("E = %g\t",E);
            }
        while (E_max-E_min>1e-26);

    Simpson=(dx/3)*Sum;
    A_Norm = sqrt(1/Simpson);
    printf("Simpson 1 = %g\t A = %g\n",Simpson,A_Norm);
    psi=0;
    omega=1;
    Sum=0;
    double L=sqrt(2*(E+1.6e-16)/k);

    for(x=-L;x<=L;x+=dx)
        {
                int j=0;
                Potential(&V,&x);
                A_psi = omega;
                A_ome = (-2*m/(h_bar*h_bar))*(E-V)*psi;

                xtemp=x+dx/2;
                Potential(&V,&xtemp);

                B_psi = omega + (dx/2)*A_ome;
                B_ome = (-2*m/(h_bar*h_bar))*(E-V)*(psi+(dx/2)*A_psi);

                C_psi = omega + (dx/2)*B_ome;
                C_ome = (-2*m/(h_bar*h_bar))*(E-V)*(psi+(dx/2)*B_psi);

                xtemp=x+dx;
                Potential(&V,&xtemp);

                D_psi = (omega + dx*C_ome);
                D_ome = (-2*m/(h_bar*h_bar))*(E-V)*(psi+dx*C_psi);

                psi = psi +(dx/6)*(A_psi+2*B_psi+2*C_psi+D_psi);
                omega = omega+ (dx/6)*(A_ome+2*B_ome+2*C_ome+D_ome);

                psi_Norm=A_Norm*psi;

                printf("x= %g\t n=%g\t Psi = %g\t E Analytical = %g\t E= %g\t Error=%g\n",x,i,psi_Norm,E_Analytical,E,(E-E_Analytical)*100/E_Analytical);

                fprintf(file_pointer,"%g\t %g\t %g\t %g\n",x,psi_Norm,V,E);

                psi_squared=psi_Norm*psi_Norm;
                if(j%2==0)
                    {
                        if (j<(2*L/dx)) Sum+=2*psi_squared;
                        else Sum+=psi_squared;
                    }
                else Sum+=4*psi_squared;
                j++;
     }

    Simpson=(dx/3)*Sum;
    printf("Simpson = %g\n",Simpson);
        }
    fclose(file_pointer);
    return EXIT_SUCCESS;

}

void Energy(double* a ,double* b){

double pi=3.1415926535897;
double m= 9.10938e-31;
double k=20000;

double h_bar=6.62607e-34/(2*pi);

*a=(2.0*(*b)+1.0)*(h_bar/2.0)*(sqrt(k/m));
return;
}
void Potential(double* a ,double* b)
{

double k=20000;
*a=0.5*k*(*b)*(*b);
return;
}

