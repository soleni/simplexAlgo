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
        simplexAlgo(matrix A, std::vector<double> C, std::vector<double> b, std::vector<double> dh, std::vector<double> dl, std::vector<double> x, std::vector<int> J);

        void solve();
        void solveDual();
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
        std::vector<double> b;
        std::vector<double> dh;
        std::vector<double> dl;
        std::vector<double> x;
        std::vector<int> J;
    };

    simplexAlgo::simplexAlgo(matrix A, std::vector<double> C, std::vector<double> b, std::vector<double> dh, std::vector<double> dl, std::vector<double> x, std::vector<int> J) :
        A(A), C(C), b(b), dh(dh), dl(dl), x(x), J(J)
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

    void simplexAlgo::solveDual()
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

        //FIND PLANS
        //H
        std::vector<double> pph;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) == J.end())
            {
                if (delta[i] >= 0)
                    pph.emplace_back(dh[i]);
                else
                    pph.emplace_back(dl[i]);
            }
        }
        //B
        matrix Ah = A;
        for (int i = A.getcolumnsNum() - 1 ; i >= 0 ; --i)
        {
            if (find(J.begin(), J.end(), i) != J.end())
                Ah = Ah.getMinor(-1, i);
        }
        matrix ppb = Ab.reverse() * (matrix({b}).transpose() - Ah * matrix({pph}).transpose());

        //FIND JS
        int ii = 0;
        int js = -1;
        int sign = 1;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) != J.end())
            {
                if (ppb(ii, 0) < dl[i])
                {
                    sign = 1;
                    js = i;
                }
                if (ppb(ii, 0) > dh[i])
                {
                    sign = -1;
                    js = i;
                }
                ++ii;
            }
        }
        //FINISH
        if (js == -1)
        {
            std::cout << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n     Found answer!" << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n\n";
            int i1 = 0;
            int i2 = 0;
            for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
            {
                if (find(J.begin(), J.end(), i) == J.end())
                {
                    std::cout << std::setprecision(2) << std::setw(5) << pph[i1] << " ";
                    ++i1;
                }
                else
                {
                    std::cout << std::setprecision(2) << std::setw(5) << ppb(i2, 0) << " ";
                    ++i2;
                }
            }
            return;
        }

        //PRINT PLAN
        std::cout << "\n plan: \n";
        int i1 = 0;
        int i2 = 0;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) == J.end())
            {
                std::cout << std::setprecision(2) << std::setw(5) << pph[i1] << " ";
                ++i1;
            }
            else
            {
                std::cout << std::setprecision(2) << std::setw(5) << ppb(i2, 0) << " ";
                ++i2;
            }
        }
        std::cout << "\n";

        //FIND DIRECTIONS
        std::vector<double> p_t;
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) != J.end())
            {
                if (i == js)
                    p_t.emplace_back(sign);
                else
                    p_t.emplace_back(0);
            }
        }
        matrix p_m = matrix({p_t}) * Ab.reverse();
        std::cout << "\n Direction: \n";
        p_m.print();
        std::cout << "\n";

        matrix pdh = p_m * Ah * -1;

        //FIND STEPS
        std::cout << "\n steps: \n";
        ii = 0;
        int j1 = -1;
        double step = std::numeric_limits<double>::infinity();
        for (int i = 0 ; i < A.getcolumnsNum() ; ++i)
        {
            if (find(J.begin(), J.end(), i) == J.end())
            {
                if (delta[i] * pdh(0, ii) < 0 && -delta[i] / pdh(0, ii) < step)
                {
                    j1 = i;
                    step = -delta[i] / pdh(0, ii);
                }
                std::cout << std::setprecision(2) << std::setw(5) << step << " ";
                ++ii;
            }
        }
        std::cout << "\n";
        if (step == std::numeric_limits<double>::infinity())
        {
            std::cout << "\nCANNOT SOLVE\n";
            return;
        }

        //GET NEW J
        std::replace(J.begin(), J.end(), js, j1);
        std::sort(J.begin(), J.end());

        std::cout << "\n Jb: \n";
        for (double const & Ji : J)
            std::cout << std::setprecision(2) << std::setw(5) << Ji << " ";
        std::cout << "\n";

        solveDual();
    }
}

