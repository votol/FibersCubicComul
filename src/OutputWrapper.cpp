#include <iostream>
#include "OutputWrapper.h"

using namespace clde;

class PolynomialOutputWrapper: public PolynomialOutput
{
public:
    template<typename _InputIterator,
         typename = std::_RequireInputIter<_InputIterator>>
    PolynomialOutputWrapper(_InputIterator begin,  _InputIterator end,
                     const OperatorDimension& dim, const unsigned int& calc_count,
                     const std::shared_ptr<ICLmanager>& context): PolynomialOutput(begin, end, 1, dim, calc_count, context)
    {

    }

    const std::vector<double>& get_data(void)
    {
        return m_result;
    }
};

class OutputWrapper: public IOutput
{
    std::string m_name;
    std::shared_ptr<PolynomialOutputWrapper> m_getter;
    std::vector<size_t> m_dims;

public:
    OutputWrapper(const std::string& name, const size_t& dim,
                  const std::shared_ptr<PolynomialOutputWrapper>& getter):
        m_name(name), m_getter(getter)
    {
        m_dims.push_back(dim);
    }

    virtual const std::string& GetName() override
    {
        return m_name;
    }

    virtual const std::vector<double>& GetData() override
    {
        return m_getter->get_data();
    }

    virtual const std::vector<size_t>& GetDimensions() override
    {
        return m_dims;
    }
};

PolyOutputWrapper::PolyOutputWrapper(const clde::Polynomial& eqs,
                  const clde::OperatorDimension& oper_dim, const unsigned int& calc_count,
                  std::string name, const std::shared_ptr<clde::ICLmanager>& context)
{
    OperatorDimension outDim = PolynomialOperator::calculateDimension(eqs.begin(), eqs.end(), false);
    outDim.in_dim = oper_dim.out_dim;
    auto pointer = std::make_shared<PolynomialOutputWrapper>(eqs.begin(), eqs.end(), outDim, calc_count, context);
    m_calc_pointer = pointer;
    m_out_pointer = std::make_shared<OutputWrapper>(name, outDim.out_dim, pointer);
}

std::shared_ptr<clde::IDEOutput> PolyOutputWrapper::get_calc_pointer()
{
    return m_calc_pointer;
}

std::shared_ptr<IOutput> PolyOutputWrapper::get_out_pointer()
{
    return m_out_pointer;
}
