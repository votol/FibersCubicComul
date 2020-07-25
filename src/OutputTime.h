#pragma once
#include <memory>
#include "IDEOutput.h"
#include "OutputInterface.h"


class OutputTime
{
    std::shared_ptr<clde::IDEOutput> m_calc_pointer;
    std::shared_ptr<IOutput> m_out_pointer;

public:
    OutputTime(const unsigned int& calc_count,
               std::string name);

    std::shared_ptr<clde::IDEOutput> get_calc_pointer();
    std::shared_ptr<IOutput> get_out_pointer();
};

