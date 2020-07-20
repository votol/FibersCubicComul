#pragma once
#include <vector>
#include <set>
#include "APolA.h"
#include "AbOA.h"

namespace AbAl
{
    class ComulManager
    {
        unsigned int m_order;
        std::vector<Polynomial<std::vector<unsigned int>, int> > m_storage;
        Polynomial<std::vector<unsigned int>, int> m_build_comul_patern(const unsigned int&);
    public:
        ComulManager(const unsigned int& order):m_order(order)
        {

        }

        Polynomial<POperator, int> split(const POperator&);
        Polynomial<std::vector<unsigned int>, int>& get_commul(unsigned int in)
        {
            return m_storage[in - m_order - 1];
        }
    };
}
