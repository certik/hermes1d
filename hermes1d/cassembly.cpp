#include "cassembly.h"
#include "stdio.h"

void A::print_info()
{
    printf("nmesh: %d\n", this->nmesh);
    printf("%f %f\n", this->mesh[0], this->mesh[1]);
}
