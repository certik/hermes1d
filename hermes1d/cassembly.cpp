#include "stdio.h"
#include "math.h"

#include "cassembly.h"

void System::print_info()
{
    printf("# nodes: %d\n", this->nmesh);
}

void System::set_dof_A(int i, int j, double value)
{
    this->A->set_value(i, j, value);
}

void System::set_dof_B(int i, int j, double value)
{
    this->B->set_value(i, j, value);
}

double System::int_grad_u_grad_v(int i, int j)
{
    double a1 = 1.0;
    if (j == i - 1) {
        double hi = this->h(i);
        return -a1/hi;
    } else if (j == i) {
        double hi = this->h(i);
        double hi_plus_1 = this->h(i+1);
        return a1*(1./hi + 1./hi_plus_1);
    } else if (j == i + 1) {
        double hi_plus_1 = this->h(i+1);
        return -a1/hi_plus_1;
    } else
        return 0.;
}

double System::int_u_v(int i, int j)
{
    double a0 = 1.0;
    if (j == i - 1) {
        double hi = this->h(i);
        return a0*hi/6;
    } else if (j == i) {
        double hi = this->h(i);
        double hi_plus_1 = this->h(i+1);
        return a0*(hi/3.+hi_plus_1/3.);
    } else if (j == i + 1) {
        double hi_plus_1 = this->h(i+1);
        return a0*hi_plus_1/6;
    } else
        return 0.;
}

double System::int_grad_u_v_over_x(int i, int j)
{
    double a = this->mesh[j+1];
    if (j == i - 1) {
        double h = this->h(i);
        return (-h + (h-a)*(log(a - h) - log(a)))/pow(h, 2);
    } else if (j == i) {
        double h = this->h(i);
        double h2 = this->h(i+1);
        return (h + (h-a)*(log(a) - log(a - h)))/pow(h, 2) + (h2 + (a+h2)*(log(a) - log(a + h2)))/pow(h2, 2);
    } else if (j == i + 1) {
        double h2 = this->h(i+1);
        return (-h2 + a*log(a + h2) - a*log(a))/pow(h2, 2);
    } else
        return 0.;
}

double System::int_u_v_over_x(int i, int j)
{


    double a = this->mesh[j+1];
    if (j == i - 1) {
        double h = this->h(i);
        return (-2*pow(a,2)*log(a) + 2*a*h + 2*pow(a,2)*log(a - h) - 2*a*h*log(a - h) + 2*a*h*log(a) - pow(h,2))/(2*pow(h,2));
    } else if (j == i) {
        double h = this->h(i);
        double h2 = this->h(i+1);
        //printf("H: %d %d %f %f\n", i, j, a, h);
        return (-2*a*h - 2*pow(a,2)*log(a - h) - 2*pow(h,2)*log(a - h) + 2*pow(a,2)*log(a) + 2*pow(h,2)*log(a) - 4*a*h*log(a) + 4*a*h*log(a - h) + 3*pow(h,2))/(2*pow(h,2)) + (-2*a*h2 - 2*pow(a,2)*log(a) - 2*pow(h2,2)*log(a) + 2*pow(a,2)*log(a + h2) + 2*pow(h2,2)*log(a + h2) - 4*a*h2*log(a) + 4*a*h2*log(a + h2) - 3*pow(h2,2))/(2*pow(h2,2));
    } else if (j == i + 1) {
        double h2 = this->h(i+1);
        return (-2*pow(a,2)*log(a + h2) + 2*a*h2 + 2*pow(a,2)*log(a) - 2*a*h2*log(a + h2) + 2*a*h2*log(a) + pow(h2,2))/(2*pow(h2,2));
    } else
        return 0.;
}

double System::int_u_v_over_x2(int i, int j)
{


    double a = this->mesh[j+1];
    if (j == i - 1) {
        double h = this->h(i);
        return (a*pow(h,4)*log(a) - a*pow(h,4)*log(a - h) - 3*pow(a,2)*pow(h,3)*log(a) - 2*pow(a,3)*pow(h,2)*log(a - h) + 2*pow(a,3)*pow(h,2)*log(a) + 3*pow(a,2)*pow(h,3)*log(a - h) - 2*pow(a,2)*pow(h,3) + 2*a*pow(h,4))/(pow(a,2)*pow(h,4) - a*pow(h,5));
    } else if (j == i) {
        double h = this->h(i);
        double h2 = this->h(i+1);
        return (-4*pow(a,2)*pow(h2,3)*log(a + h2) - 2*a*pow(h2,4)*log(a + h2) - 2*pow(a,3)*pow(h2,2)*log(a + h2) + 2*a*pow(h2,4)*log(a) + 2*pow(a,3)*pow(h2,2)*log(a) + 4*pow(a,2)*pow(h2,3)*log(a) + 2*pow(a,2)*pow(h2,3) + 3*a*pow(h2,4) + pow(h2,5))/(pow(a,2)*pow(h2,4) + a*pow(h2,5)) + (-4*pow(a,2)*pow(h,3)*log(a - h) - 2*a*pow(h,4)*log(a) - 2*pow(a,3)*pow(h,2)*log(a) + 2*a*pow(h,4)*log(a - h) + 2*pow(a,3)*pow(h,2)*log(a - h) + 4*pow(a,2)*pow(h,3)*log(a) + 2*pow(a,2)*pow(h,3) - 3*a*pow(h,4) + pow(h,5))/(pow(a,2)*pow(h,4) - a*pow(h,5));
    } else if (j == i + 1) {
        double h2 = this->h(i+1);
        return (a*pow(h2,4)*log(a + h2) - a*pow(h2,4)*log(a) - 3*pow(a,2)*pow(h2,3)*log(a) - 2*pow(a,3)*pow(h2,2)*log(a) + 2*pow(a,3)*pow(h2,2)*log(a + h2) + 3*pow(a,2)*pow(h2,3)*log(a + h2) - 2*pow(a,2)*pow(h2,3) - 2*a*pow(h2,4))/(pow(a,2)*pow(h2,4) + a*pow(h2,5));
    } else
        return 0.;
}

double System::bilinear_form_A(int i, int j)
{
    double l = 1.0;
    return int_grad_u_grad_v(i, j)/2
        - int_grad_u_v_over_x(i, j)
        - int_u_v_over_x(i, j)
        + l*(l+1)*int_u_v_over_x2(i, j)/2
        ;
}

double System::bilinear_form_B(int i, int j)
{
    return int_u_v(i, j);
}

void System::assemble()
{
    printf("assembling...\n");
    int i;
    for (i = 0; i < this->nmesh-2; i++) {
        if (i>0) set_dof_A(i, i-1, bilinear_form_A(i, i-1));
        set_dof_A(i, i, bilinear_form_A(i, i));
        set_dof_A(i, i+1, bilinear_form_A(i, i+1));
        if (i>0) set_dof_B(i, i-1, bilinear_form_B(i, i-1));
        set_dof_B(i, i, bilinear_form_B(i, i));
        set_dof_B(i, i+1, bilinear_form_B(i, i+1));
    }
    i = this->nmesh-2;
    set_dof_A(i, i-1, bilinear_form_A(i-1, i-1-1));
    set_dof_A(i, i, bilinear_form_A(i-1, i-1));
    set_dof_B(i, i-1, bilinear_form_B(i-1, i-1-1));
    set_dof_B(i, i, bilinear_form_B(i-1, i-1));
}
