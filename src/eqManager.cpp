#include <array>
#include "eqManager.h"

using namespace clde;
using namespace AbAl;

POperator EqManager::_IndexManager::m_conj(const POperator& in)
{
    POperator result;
    for(auto&& particle: in)
    {
        result[-particle.first] = particle.second;
    }
    return result;
}

EqManager::_IndexManager::_IndexManager():m_current_ind(1)
{

}

int EqManager::_IndexManager::GetIndex(const POperator& in)
{
    if(in.size() == 0)
        return 0;
    auto find_it = m_index.find(in);
    if(find_it != m_index.end())
        return static_cast<int>(find_it->second);

    auto conj_input = m_conj(in);
    find_it = m_index.find(conj_input);
    if(find_it != m_index.end())
        return -static_cast<int>(find_it->second);

    m_index[in] = m_current_ind;
    m_tasks.push_back(in);
    m_current_ind ++;
    return static_cast<int>(m_current_ind) - 1;
}

const std::map<AbAl::POperator, std::uint32_t>& EqManager::_IndexManager::GetAllData()
{
    return m_index;
}

std::vector<AbAl::POperator>& EqManager::_IndexManager::getTasks()
{
    return m_tasks;
}

AbAl::Polynomial<std::array<int,1>, std::complex<double>> EqManager::m_converter(
        const AbAl::Polynomial<AbAl::POperator, int>& in)
{
    AbAl::Polynomial<std::array<int,1>, std::complex<double>> result;
    AbAl::Polynomial<std::array<int,1>, std::complex<double>> res_eq1;
    std::array<int,1> tmp_index;
    for(auto&& mon: in.data())
    {
        res_eq1 = AbAl::Polynomial<std::array<int,1>, std::complex<double>>(1);
        for(auto&& factor: mon.first)
        {
            tmp_index[0] = m_ind_manag->GetIndex(factor.first);
            if(tmp_index[0] == 0)
                continue;
            for(size_t ind = 0; ind < factor.second; ++ind)
                res_eq1 *= tmp_index;
        }
        result += res_eq1 * std::complex<double>(static_cast<double>(mon.second));
    }
    return result;
}


EqManager::EqManager(const unsigned int& comul_order):m_com_manger(comul_order)
{
    m_ind_manag = std::make_unique<_IndexManager>();
}

IIndexManage* EqManager::getIndexManeger(void)
{
    return m_ind_manag.get();
}

void EqManager::addHamiltonian(AbAl::GOperator&& in)
{
    m_hamiltonians.push_back(in);
}

namespace AbAl {
std::ostream& operator<<(std::ostream& os, const std::array<int,1>& in)
{
    os << in[0];
    return os;
}
}

clde::PolynomialC EqManager::buildPoly(const AbAl::Polynomial<POperator, int>& in, unsigned int out_ind)
{
    clde::PolynomialC result;
    auto poly = m_converter(in);
    for(auto&& mon: poly.data())
    {
        result.push_back(MonomialC());
        result.back().outInd = out_ind;
        result.back().coe = mon.second;
        for(auto&& factor: mon.first)
        {
            if(factor.first[0] > 0)
            {
                for(size_t ind = 0; ind < factor.second; ++ind)
                    result.back().inInds.push_back(static_cast<unsigned int>(factor.first[0] - 1));
            }
            else if(factor.first[0] < 0)
            {
                for(size_t ind = 0; ind < factor.second; ++ind)
                    result.back().inIndsC.push_back(static_cast<unsigned int>(-factor.first[0] - 1));
            }
        }
    }
    return result;
}

std::vector<PolynomialC> EqManager::buildEquations()
{
    std::vector<PolynomialC> result(m_hamiltonians.size());
    GOperator tmp_op;
    AbAl::Polynomial<std::array<int,1>, std::complex<double>> res_eq;
    unsigned int eq_ind = 0;
    unsigned int ham_ind;
    while(eq_ind != m_ind_manag->getTasks().size())
    {
        ham_ind = 0;
        for(auto&& ham : m_hamiltonians)
        {

            res_eq = AbAl::Polynomial<std::array<int,1>, std::complex<double>>();
            tmp_op = std::complex<double>(0.0, 1.0) * commute(ham, m_ind_manag->getTasks()[eq_ind]);
            for(auto&& op_summ: tmp_op)
            {
                res_eq += m_converter(m_com_manger.split(op_summ.first)) * op_summ.second;
            }
            for(auto&& mon: res_eq.data())
            {
                result[ham_ind].push_back(MonomialC());
                result[ham_ind].back().outInd = eq_ind;
                result[ham_ind].back().coe = mon.second;
                for(auto&& factor: mon.first)
                {
                    if(factor.first[0] > 0)
                    {
                        for(size_t ind = 0; ind < factor.second; ++ind)
                            result[ham_ind].back().inInds.push_back(static_cast<unsigned int>(factor.first[0] - 1));
                    }
                    else if(factor.first[0] < 0)
                    {
                        for(size_t ind = 0; ind < factor.second; ++ind)
                            result[ham_ind].back().inIndsC.push_back(static_cast<unsigned int>(-factor.first[0] - 1));
                    }
                }
            }

            ++ham_ind;
        }
        ++eq_ind;
    }


    return result;
}
