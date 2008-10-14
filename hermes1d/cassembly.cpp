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
    int i;
    set_dof(1, 2, 3.4);
}
