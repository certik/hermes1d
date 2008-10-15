#ifndef matrix_h
#define matrix_h

#include "common.h"

class SparseMatrix
{
public:
    SparseMatrix(int nvalues)
    {
        this->Ax = new double[nvalues];
        this->Ai = new int[nvalues];
        this->Aj = new int[nvalues];
        this->A_len = 0;
        this->A_max = nvalues;
    }

    ~SparseMatrix()
    {
        delete[] this->Ax;
        delete[] this->Ai;
        delete[] this->Aj;
    }

    void set_value(int i, int j, double val)
    {
        printf("setting value: %d %d %f\n", i, j, val);
        this->Ai[this->A_len] = i;
        this->Aj[this->A_len] = j;
        this->Ax[this->A_len] = val;
        this->A_len++;
        if (this->A_len >= this->A_max) error("A_len >= A_max");
    }
    int *Ai, *Aj;
    double *Ax;
    int A_len, A_max;
};

#endif
