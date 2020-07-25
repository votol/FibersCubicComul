#pragma once
#include <memory>
#include "PolynomialOutput.h"
#include "OutputInterface.h"

class PolyOutputWrapper
{
    std::shared_ptr<clde::IDEOutput> m_calc_pointer;
    std::shared_ptr<IOutput> m_out_pointer;

public:
    PolyOutputWrapper(const clde::Polynomial& eqs,
                      const clde::OperatorDimension& oper_dim, const unsigned int& calc_count,
                      std::string name, const std::shared_ptr<clde::ICLmanager>& context);

    std::shared_ptr<clde::IDEOutput> get_calc_pointer();
    std::shared_ptr<IOutput> get_out_pointer();
};

