#include "cassembly.h"
#include "stdio.h"

void A::print_info()
{
    printf("# nodes: %d\n", this->nmesh);
}

void set_dof(int i, int j, double value)
{
    printf("%d %d %f\n", i, j, value);
}

void A::assemble()
{
    printf("assembling...\n");
    double h0, h1, a0 = 1.0, a1 = 1.0;
    h0 = this->h(0);
    h1 = this->h(1);
    set_dof(0, 0, a1*(1./h0+1./h1)+a0*(h0/3. + h1/3.));
    set_dof(0, 1, -a1/h1 + a0*h1/6.);
    for (int i=1; i<this->nmesh-2; i++) {
        h0 = this->h(i);
        h1 = this->h(i+1);
        set_dof(i, i-1, -a1/h0 + a0*h0/6);
        set_dof(i, i+1, -a1/h1 + a0*h1/6);
    }
}
