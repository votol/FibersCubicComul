#include <filesystem>
#include <fstream>
#include <cstdint>
#include "ModelBuilder.h"
#include "ParameterDefines.h"
#include "eqManager.h"

namespace fs = std::filesystem;
using namespace AbAl;

void ModelBuilder::m_build(const std::vector<double>& parameters,
                           const std::optional<std::ofstream*>& os = std::optional<std::ofstream*>())
{
    std::uint32_t n_fibs = static_cast<std::uint32_t>(Nfibs);
    //std::uint32_t n_fibs_calc = static_cast<std::uint32_t>(Nfibs_calc);
    EqManager manager(3);
    GOperator linear_op;
    GOperator nonlinear_op;
    for(std::uint32_t ind = 1; ind <= n_fibs; ++ind)
    {
        if(ind != n_fibs)
            linear_op += ak(ind) *a(ind + 1) + ak(ind + 1) *a(ind);
        nonlinear_op += ak(ind) * ak(ind) * a(ind) * a(ind);
    }
    manager.addHamiltonian(std::move(linear_op));
    manager.addHamiltonian(std::move(nonlinear_op));

    clde::PolynomialC intensity;
    AbAl::ComulManager com_manager(2);
    for(std::uint32_t ind = 1; ind <= n_fibs; ++ind)
    {
        clde::unitePoly(intensity,
                        manager.buildPoly(com_manager.split((ak(ind)*a(ind)).begin()->first), ind - 1),
                        std::complex<double>(1.0));
    }
    auto intensity_tmp = convertMonomials(intensity);
    for(auto&& mon: intensity_tmp)
    {
        if(mon.outInd %2 ==0)
        {
            m_intensity.push_back(mon);
            m_intensity.back().outInd /= 2;
        }
    }

    clde::PolynomialC third_order;
    clde::PolynomialC third_order_sp;
    AbAl::ComulManager com_manager1(3);
    for(std::uint32_t ind = 1; ind <= n_fibs; ++ind)
    {
        for(std::uint32_t ind1 = 1; ind1 <= n_fibs; ++ind1)
        {
            clde::unitePoly(third_order,
                            manager.buildPoly(com_manager1.split((ak(ind)*a(ind) * a(ind1)).begin()->first), (ind - 1) * n_fibs + ind1 - 1),
                            std::complex<double>(1.0));
            clde::unitePoly(third_order_sp,
                            manager.buildPoly(com_manager.split((ak(ind)*a(ind) * a(ind1)).begin()->first), (ind - 1) * n_fibs + ind1 - 1),
                            std::complex<double>(1.0));
        }
    }

    for(std::uint32_t ind = 1; ind <= n_fibs; ++ind)
    {
        for(std::uint32_t ind1 = 1; ind1 <= n_fibs; ++ind1)
        {
            clde::unitePoly(third_order,
                            manager.buildPoly(com_manager1.split((a(ind)*a(ind) * a(ind1)).begin()->first), n_fibs * n_fibs + (ind - 1) * n_fibs + ind1 - 1),
                            std::complex<double>(1.0));
            clde::unitePoly(third_order_sp,
                            manager.buildPoly(com_manager.split((a(ind)*a(ind) * a(ind1)).begin()->first), n_fibs * n_fibs + (ind - 1) * n_fibs + ind1 - 1),
                            std::complex<double>(1.0));
        }
    }

    m_third_order = convertMonomials(third_order);
    m_third_order_sp = convertMonomials(third_order_sp);

    auto ham_parts = manager.buildEquations();

    clde::Polynomial linear = convertMonomials(ham_parts[0]);
    clde::Polynomial nonlinear = convertMonomials(ham_parts[1]);

    if(os)
    {
        *(*os) << linear;
        *(*os) << nonlinear;
        *(*os) << m_intensity;
        *(*os) << m_third_order;
        *(*os) << m_third_order_sp;
    }

    clde::unitePoly(m_main, std::move(linear), -1.0);
    clde::unitePoly(m_main, std::move(nonlinear), -L/2.0);

    clde::PolynomialC init;
    for(auto&& el: manager.getIndexManeger()->GetAllData())
    {
        init.push_back(clde::MonomialC());
        init.back().coe = 1.0;
        init.back().outInd = el.second - 1;
        for(auto&& factor: el.first)
        {
            if(factor.first > 0)
            {
                for(size_t ind = 0; ind < factor.second; ++ind)
                    init.back().inInds.push_back(static_cast<unsigned int>(factor.first - 1));
            }
            else if(factor.first < 0)
            {
                for(size_t ind = 0; ind < factor.second; ++ind)
                    init.back().inIndsC.push_back(static_cast<unsigned int>(-factor.first - 1));
            }
        }
    }
    m_init = convertMonomials(init);

    if(os)
    {
        *(*os) << m_init;
    }
}

ModelBuilder::ModelBuilder(const std::vector<double>& parameters, std::optional<std::string> cache_file)
{
    std::uint32_t n_fibs = static_cast<std::uint32_t>(Nfibs);
    std::uint32_t n_fibs_calc = static_cast<std::uint32_t>(Nfibs_calc);
    if(cache_file)
    {
        if(fs::exists(fs::path(*cache_file)))
        {
            auto i_file = std::ifstream(*cache_file, std::ifstream::in | std::ifstream::binary);
            std::uint32_t fibs_read;
            std::uint32_t fibs_calc_read;
            i_file.read(reinterpret_cast <char*>(&fibs_read), sizeof (std::uint32_t));
            i_file.read(reinterpret_cast <char*>(&fibs_calc_read), sizeof (std::uint32_t));
            if(n_fibs != fibs_read || n_fibs_calc != fibs_calc_read)
            {
                i_file.close();
                auto o_file = std::ofstream(*cache_file, std::ofstream::out | std::ofstream::binary);
                o_file.write(reinterpret_cast <char*>(&n_fibs), sizeof (std::uint32_t));
                o_file.write(reinterpret_cast <char*>(&n_fibs_calc), sizeof (std::uint32_t));
                m_build(parameters, &o_file);
            }
            else
            {
                clde::Polynomial linear;
                clde::Polynomial nonlinear;
                i_file >> linear;
                i_file >> nonlinear;
                clde::unitePoly(m_main, std::move(linear), -1.0);
                clde::unitePoly(m_main, std::move(nonlinear), -L/2.0);
                i_file >> m_intensity;
                i_file >> m_third_order;
                i_file >> m_third_order_sp;
                i_file >> m_init;
            }
        }
        else
        {
            auto o_file = std::ofstream(*cache_file, std::ofstream::out | std::ofstream::binary);
            o_file.write(reinterpret_cast <char*>(&n_fibs), sizeof (std::uint32_t));
            o_file.write(reinterpret_cast <char*>(&n_fibs_calc), sizeof (std::uint32_t));
            m_build(parameters, &o_file);
        }
    }
    else
    {
        m_build(parameters);
    }
}

clde::Polynomial& ModelBuilder::intensity(void)
{
   return m_intensity;
}

clde::Polynomial& ModelBuilder::main_equations(void)
{
    return m_main;
}

clde::Polynomial& ModelBuilder::init_calculation(void)
{
    return m_init;
}

clde::Polynomial& ModelBuilder::third(void)
{
    return m_third_order;
}

clde::Polynomial& ModelBuilder::third_sp(void)
{
    return m_third_order_sp;
}
