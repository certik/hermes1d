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
    double a = this->mesh[j];
    if (j == i - 1) {
        double h = this->h(i);
        return (-h + (h-a)*(log(a - h) - log(a)))/pow(h, 2);
    } else if (j == i) {
        double h = this->h(i);
        double h2 = this->h(i+1);
        return (h + a*log(a - h) + h*log(a) - a*log(a) - h*log(a - h))/pow(h, 2) + (h2 + a*log(a) + h2*log(a) - a*log(a + h2) - h2*log(a + h2))/pow(h2, 2);
    } else if (j == i + 1) {
        double h2 = this->h(i+1);
        return (-h2 + a*log(a + h2) - a*log(a))/pow(h2, 2);
    } else
        return 0.;
}

double System::bilinear_form_A(int i, int j)
{
    return int_grad_u_grad_v(i, j);
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
