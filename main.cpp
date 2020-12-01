#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "simplex.h"

int main(int argc, char ** argv)
{
    matrix A({{ 4, -1, 1, 0, 0},
             {  2,  1, 0, 1, 0},
             { 2, 1, 0, 0, -1},});

    simplexAlgo algo(A,
                     {5 , 3, 0, 0, 0},
                     {std::numeric_limits<double>::infinity(), 4, 4, 6, std::numeric_limits<double>::infinity()},
                     {0, 0, 0, 0, 2},
                     {0, 0, 0, 6, 2},
                     {2, 3, 4}
    );
    algo.solve();

    return 0;
}