#include "stdio.h"

#include "cassembly.h"

void System::print_info()
{
    printf("# nodes: %d\n", this->nmesh);
}

void System::set_dof_A(int i, int j, double value)
{
    printf("A: ");
    this->A->set_value(i, j, value);
}

void System::set_dof_B(int i, int j, double value)
{
    printf("B: ");
    this->B->set_value(i, j, value);
}

void System::assemble()
{
    printf("assembling...\n");
    double h0, h1, a0 = 1.0, a1 = 1.0;
    h0 = this->h(0);
    h1 = this->h(1);
    set_dof_A(0, 0, a1*(1./h0+1./h1));
    set_dof_A(0, 1, -a1/h1);
    set_dof_B(0, 0, a0*(h0/3. + h1/3.));
    set_dof_B(0, 1, a0*h1/6.);
    int i;
    for (i=1; i<this->nmesh-2; i++) {
        h0 = this->h(i);
        h1 = this->h(i+1);
        set_dof_A(i, i-1, -a1/h0);
        set_dof_A(i, i, a1*(1./h0 + 1./h1));
        set_dof_A(i, i+1, -a1/h1);
        set_dof_B(i, i-1, a0*h0/6);
        set_dof_B(i, i, a0*(h0/3.+h1/3.));
        set_dof_B(i, i+1, a0*h1/6);
    }
    i = this->nmesh-2;
    h0 = this->h(i-1);
    h1 = this->h(i);
    set_dof_A(i, i-1, -a1/h0);
    set_dof_A(i, i, a1*(1./h0 + 1./h1));
    set_dof_B(i, i-1, a0*h0/6);
    set_dof_B(i, i, a0*(h0/3.+h1/3.));
}
