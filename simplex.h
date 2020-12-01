#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "matrix.h"

namespace
{
    class simplexAlgo
    {
    public:
        simplexAlgo(matrix A, std::vector<double> C, std::vector<double> dh, std::vector<double> dl, std::vector<double> x, std::vector<int> J);

        void solve();
        std::vector<double> getX() { return x; }
        std::vector<int> getJ() { return J; }

        void setMaxStep(int i) { maxstep = i; }

    private:
        double e = 10e-7;
        int step = 0;
        int maxstep = 99;

    protected:
        matrix A;
        std::vector<double> C;
        std::vector<double> dh;
        std::vector<double> dl;
        std::vector<double> x;
        std::vector<int> J;
    };

    simplexAlgo::simplexAlgo(matrix A, std::vector<double> C, std::vector<double> dh, std::vector<double> dl, std::vector<double> x, std::vector<int> J) :
        A(A), C(C), dh(dh), dl(dl), x(x), J(J)
    {
    }

    void simplexAlgo::solve()
    {
        if (step == maxstep)
        {
            std::cout << "!timeout!\n";
            return;
        }
        std::cout << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n         STEP " << ++step << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n ";

        // GET A BASIS AND C BASIS
        std::vector<double> cb;
        matrix Ab = A;
        for (int i = A.getcolumnsNum() - 1 ; i >= 0 ; --i)
        {
            if (find(J.begin(), J.end(), i) == J.end())
                Ab = Ab.getMinor(-1, i);
        }
        for (auto const & j : J)
            cb.emplace_back(C[j]);

        //FIND U
        matrix u = matrix({cb}) * Ab.reverse();
        std::cout << "\n u: \n";
        u.print();

        //FIND DELTA
        std::vector<double> delta;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) != J.end())
                delta.emplace_back(0);
            else
                delta.emplace_back(C[i] - (u * A.column(i))(0, 0));
        }
        std::cout << "\n delta: \n";
        for (double const & delta_i : delta)
            std::cout << std::setprecision(2) << std::setw(5) << delta_i << " ";
        std::cout << "\n";

        //FIND l
        std::vector<double> L(A.getcolumnsNum());
        int j0 = -1;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) == J.end())
            {
                if (!((delta[i] >= 0 && std::abs((double) (x[i] - dh[i])) < e) || (delta[i] <= 0 && std::abs((double) (x[i] - dl[i])) < e)))
                {
                    if ((j0 == -1 || std::abs(delta[i]) > std::abs(delta[j0])) && delta[i] != 0)
                    {
                        L[i] = 0;
                        j0 = i;
                    }
                }
            }
        }
        //EVERY CRITERION HAS PASSED, FINISH
        if (j0 == -1)
        {
            std::cout << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n     Found answer!" << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n ";

            std::cout << "\n X: \n";
            for (double const & xi : x)
                std::cout << std::setprecision(2) << std::setw(5) << xi << " ";
            std::cout << "\n";

            std::cout << "\n Jb: \n";
            for (double const & Ji : J)
                std::cout << std::setprecision(2) << std::setw(5) << Ji << " ";
            std::cout << "\n";

            return;
        }
        std::cout << "\n j0: \n" << std::setprecision(2) << std::setw(5) << j0 << "\n";

        matrix l0 = Ab.reverse() * A.column(j0) * (delta[j0] < 0 ? 1 : -1);
        for (int i = 0; i < J.size() ; ++i)
            L[J[i]] = l0(i, 0);
        L[j0] = delta[j0] > 0 ? 1 : -1;

        std::cout << "\n L: \n";
        for (double const & l : L)
            std::cout << std::setprecision(2) << std::setw(5) << l << " ";
        std::cout << "\n";

        //FIND TETA
        std::vector<double> teta;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (i == j0)
                teta.emplace_back(dh[i]);
            else if (L[i] == 0)
                teta.emplace_back(std::numeric_limits<double>::infinity());
            else if (L[i] > 0)
                teta.emplace_back((dh[i] - x[i]) / L[i]);
            else if (L[i] < 0)
                teta.emplace_back((dl[i] - x[i]) / L[i]);
        }

        std::cout << "\n teta: \n";
        for (double const & t : teta)
            std::cout << std::setprecision(2) << std::setw(5) << t << " ";
        std::cout << "\n";

        //FIND TETA 0
        int js = -1;
        double teta0 = std::numeric_limits<double>::infinity();
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (teta0 > teta[i])
            {
                js = i;
                teta0 = teta[i];
            }
        }

        std::cout << "\n js: \n" << std::setprecision(2) << std::setw(5) << js << "\n";
        std::cout << "\n teta0: \n" << std::setprecision(2) << std::setw(5) << teta0 << "\n";

        //GET NEW J
        std::replace(J.begin(), J.end(), js, j0);
        std::sort(J.begin(), J.end());

        //GET NEW X
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
            x[i] += L[i] * teta0;

        std::cout << "\n X: \n";
        for (double const & xi : x)
            std::cout << std::setprecision(2) << std::setw(5) << xi << " ";
        std::cout << "\n";

        std::cout << "\n Jb: \n";
        for (double const & Ji : J)
            std::cout << std::setprecision(2) << std::setw(5) << Ji << " ";
        std::cout << "\n";

        //NEXT STEP
        solve();
    }
}

