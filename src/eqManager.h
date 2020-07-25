#pragma once
#include <map>
#include <vector>
#include <list>
#include <cstdint>
#include <memory>
#include <array>
#include "AbOA.h"
#include "AbComul.h"
#include "PolynomialUtils.h"

class IIndexManage
{
public:
    IIndexManage(){}
    virtual ~IIndexManage(){}
    virtual int GetIndex(const AbAl::POperator& ) = 0;
    virtual const std::map<AbAl::POperator, std::uint32_t>& GetAllData() = 0;
};


class EqManager
{
    class _IndexManager: public IIndexManage
    {
        std::map<AbAl::POperator, std::uint32_t> m_index;
        std::vector<AbAl::POperator> m_tasks;
        uint32_t m_current_ind;

        AbAl::POperator m_conj(const AbAl::POperator& in);
    public:
        _IndexManager();
        virtual int GetIndex(const AbAl::POperator& ) override;
        virtual const std::map<AbAl::POperator, std::uint32_t>& GetAllData() override;
        std::vector<AbAl::POperator>& getTasks();
    };

    AbAl::ComulManager m_com_manger;
    std::list<AbAl::GOperator> m_hamiltonians;

    std::unique_ptr<_IndexManager> m_ind_manag;

    AbAl::Polynomial<std::array<int,1>, std::complex<double>> m_converter(const AbAl::Polynomial<AbAl::POperator, int>&);
public:
    EqManager(const unsigned int& comul_order);
    IIndexManage* getIndexManeger(void);
    void addHamiltonian(AbAl::GOperator&&);
    clde::PolynomialC buildPoly(const AbAl::Polynomial<AbAl::POperator, int>&, unsigned int out_ind = 0 );
    std::vector<clde::PolynomialC> buildEquations();
};
