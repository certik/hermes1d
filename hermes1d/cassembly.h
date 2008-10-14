#ifndef cassembly_h
#define cassembly_h

class A
{
public:
    void set_mesh(double *mesh, int nmesh)
    {
        this->mesh = mesh;
        this->nmesh = nmesh;
    }
    void print_info();
    void assemble();

    double h(int i)
    {
        return this->mesh[i+1] - this->mesh[i];
    }

    double *mesh;
    int nmesh;
};

#endif
