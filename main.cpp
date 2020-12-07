#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "simplex.h"

int main(int argc, char ** argv)
{
    matrix A({{  0, 0, 0.5, 4, 0, 1, 0, 0},
             {   0, 0, 0, 3, 1, 0, 1, 0},
             {   2, 1, 2, 0, 0, 0, 0, 1},});

    simplexAlgo algo(A,
                     {0, 1, 2, 17, 3, 0, 0, 0},
                     {29, 23, 8},
                     {4, 3, 5, 8, 4, 16, 13, 3},
                     {1, -1, 2, 3, 1, 0, 0, 0},
                     {1, -1, 2, 3, 1, 16, 13, 3},
                     {5, 6, 7}
    );
    algo.solveDual();
    /*
    matrix u({{2,3,1}});

    matrix delta = matrix({{0, 1, 2, 17, 3, 0, 0, 0}}) - u * A;
    delta.print();

    matrix b({{29, 23, 8}});
    matrix w({{0, 0, 0, 0, 0, 0, 0, 0}});
    matrix v({{0, 0, 0, 0, 0, 0, 0, 0}});
    for (int i = 0 ; i < w.getcolumnsNum(); ++i)
    {
        if (delta(0, i) >= 0)
        {
            w(0, i) = delta(0, i);
            v(0, i) = 0;
        }
        else
        {
            v(0, i) = -delta(0, i);
            w(0, i) = 0;
        }
    }
    matrix dl({{1, -1, 2, 3, 1, 0, 0, 0}});
    matrix dh({{4, 3, 5, 8, 4, 16, 13, 3}});

    matrix a = b * u.transpose() + dh * w.transpose() - dl * v.transpose();
    a.print();
    */
    /*
    matrix A({{  -1, 1, 1, 0, 0},
             {    1,-2, 0, 1, 0},
             {    4, 1, 0, 0, 1},});

    simplexAlgo algo(A,
                     {2, 1, 0, 0, 0},
                     {3, 6, 22},
                     {5, 4, 8, 14, 22},
                     {0, 0, 0, 0, 0},
                     {0, 0, 3, 6, 22},
                     {2, 3, 4}
    );
    algo.solveDual();
    */
    return 0;
}