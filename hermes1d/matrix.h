class SparseMatrix
{
public:
    SparseMatrix(int n)
    {

    }

    ~SparseMatrix()
    {
    }

    void set_value(int i, int j, double val)
    {
        printf("setting value: %d %d %f\n", i, j, val);
    }
};
