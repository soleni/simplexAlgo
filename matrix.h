#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>

namespace
{
    class matrix
    {
    public:
        matrix(int rows, int columns);
        matrix(std::vector<std::vector<double>> newMatrix);

        matrix makeI();
        matrix transpose();
        matrix reverse();

        matrix column(int i);
        matrix row(int i);

        double det();
        matrix getMinor(int i, int j);

        std::pair<int, int> getSize();
        int getRowsNum();
        int getcolumnsNum();

        double & operator()(int r, int c);
        matrix operator*(matrix m);
        matrix operator*(double val);
        matrix operator/(double val);

        void print();

    protected:
        std::vector<std::vector<double>> m_matrix;
    };

    matrix::matrix(int rows, int columns)
    {
        m_matrix = std::vector<std::vector<double>>(rows, std::vector<double>(columns, 0.));
    }

    matrix::matrix(std::vector<std::vector<double>> newMatrix)
    {
        if (newMatrix.empty())
            return;

        int rowsSize = newMatrix[0].size();
        for (auto const row : newMatrix)
        {
            if (row.size() != rowsSize)
                return;
        }

        m_matrix = newMatrix;
    }

    matrix matrix::makeI()
    {
        if (getcolumnsNum() != getRowsNum())
            return *this;

        for (int i = 0 ; i < getcolumnsNum() ; ++i)
            (*this)(i, i) = 1;

        return *this;
    }

    matrix matrix::transpose()
    {
        matrix newMatrix(getcolumnsNum(), getRowsNum());

        for (int i = 0 ; i < newMatrix.getRowsNum() ; ++i)
        {
            for (int j = 0 ; j < newMatrix.getcolumnsNum() ; ++j)
            {
                newMatrix(i, j) += (*this)(j, i);
            }
        }

        return newMatrix;
    }

    matrix matrix::reverse()
    {
        if (getcolumnsNum() != getRowsNum() || getRowsNum() == 1)
            return *this;

        matrix newMatrix(getRowsNum(), getcolumnsNum());
        for (int i = 0 ; i < getRowsNum() ; ++i)
        {
            for (int j = 0 ; j < getcolumnsNum() ; ++j)
            {
                newMatrix(i, j) = getMinor(i, j).det() * pow(-1, i + j);
            }
        }

        return newMatrix.transpose() / (*this).det();
    }

    matrix matrix::column(int i)
    {
        matrix newMatrix(getRowsNum(), 1);
        for (int j = 0 ; j < getRowsNum() ; ++j)
            newMatrix(j, 0) = (*this)(j, i);

        return newMatrix;
    }
    matrix matrix::row(int i)
    {
        matrix newMatrix(1, getcolumnsNum());
        for (int j = 0 ; j < getcolumnsNum() ; ++j)
            newMatrix(0, j) = (*this)(i, j);

        return newMatrix;
    }
    double matrix::det()
    {
        if (getcolumnsNum() != getRowsNum())
            return 0.;

        if (getRowsNum() == 1)
            return (*this)(0, 0);

        double value = 0;
        for (int i = 0 ; i < getcolumnsNum() ; ++i)
        {
            value += (*this)(0, i) * getMinor(0, i).det() * pow(-1, i);
        }

        return value;
    }

    matrix matrix::getMinor(int i, int j)
    {
        matrix newMatrix(getRowsNum() - (i >= 0 ? 1 : 0), getcolumnsNum() - (j >= 0 ? 1 : 0));
        int m_i = 0;
        int m_j = 0;
        for (int h = 0 ; h < getRowsNum() ; ++h)
        {
            if (h == i)
                continue;

            m_j = 0;
            for (int t = 0 ; t < getcolumnsNum() ; ++t)
            {
                if (t == j)
                    continue;

                newMatrix(m_i, m_j) = (*this)(h, t);
                ++m_j;
            }
            ++m_i;
        }

        return newMatrix;
    };

    std::pair<int, int> matrix::getSize()
    {
        return {getRowsNum(), getcolumnsNum()};
    }

    int matrix::getRowsNum()
    {
        if (m_matrix.empty())
            return 0;
        return m_matrix.size();
    }
    int matrix::getcolumnsNum()
    {
        if (m_matrix.empty())
            return 0;
        return m_matrix[0].size();
    }

    double & matrix::operator()(int r, int c)
    {
        if (0 <= r && r < getRowsNum() && 0 <= c && c < getcolumnsNum())
            return m_matrix[r][c];
        std::cout << "SIZE ERROR!!!\n";

        double dummy;
        return dummy;
    }

    matrix matrix::operator*(matrix m)
    {
        if (getcolumnsNum() != m.getRowsNum())
            return matrix(0, 0);

        matrix newMatrix(getRowsNum(), m.getcolumnsNum());
        for (int i = 0 ; i < newMatrix.getRowsNum() ; ++i)
        {
            for (int j = 0 ; j < newMatrix.getcolumnsNum() ; ++j)
            {
                for (int h = 0 ; h < getcolumnsNum() ; ++h)
                    newMatrix(i, j) += (*this)(i, h) * m(h, j);
            }
        }

        return newMatrix;
    }

    matrix matrix::operator*(double val)
    {
        matrix newMatrix = *this;

        for (int i = 0 ; i < newMatrix.getRowsNum() ; ++i)
        {
            for (int j = 0 ; j < newMatrix.getcolumnsNum() ; ++j)
            {
                newMatrix(i, j) *= val;
            }
        }

        return newMatrix;
    }

    matrix matrix::operator/(double val)
    {
        matrix newMatrix = *this;

        for (int i = 0 ; i < newMatrix.getRowsNum() ; ++i)
        {
            for (int j = 0 ; j < newMatrix.getcolumnsNum() ; ++j)
            {
                newMatrix(i, j) /= val;
            }
        }

        return newMatrix;
    }

    void matrix::print()
    {
        if (m_matrix.empty())
            return;

        for (auto & row : m_matrix)
        {
            for (auto & i : row)
                std::cout << std::setprecision(2) << std::setw(5) << i;
            std::cout << "\n";
        }
    }
}