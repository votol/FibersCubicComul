#include <iostream>
#include "OutputTime.h"

using namespace clde;

class TimeOutputCalculator: public IDEOutput
{
    std::vector<double> m_result;
    unsigned int m_current_output = 0;
public:
    TimeOutputCalculator(const unsigned int& calc_count):m_result(calc_count, 0.0)
    {

    }

    virtual void apply(const CLDataStorage<double>& , const std::vector<double>& params) override
    {
        if(m_current_output != m_result.size())
        {
            m_result[m_current_output] = params[0];
            m_current_output++;
        }
    }

    const std::vector<double>& get_data(void)
    {
        return m_result;
    }
};

class TimeOutputWrapper: public IOutput
{
    std::string m_name;
    std::shared_ptr<TimeOutputCalculator> m_getter;
    std::vector<size_t> m_dims;

public:
    TimeOutputWrapper(const std::string& name,
                  const std::shared_ptr<TimeOutputCalculator>& getter):
        m_name(name), m_getter(getter)
    {
        m_dims.push_back(1);
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

OutputTime::OutputTime(const unsigned int& calc_count,
           std::string name)
{
    auto pointer = std::make_shared<TimeOutputCalculator>(calc_count);
    m_calc_pointer = pointer;
    m_out_pointer = std::make_shared<TimeOutputWrapper>(name, pointer);
}

std::shared_ptr<clde::IDEOutput> OutputTime::get_calc_pointer()
{
    return m_calc_pointer;
}

std::shared_ptr<IOutput> OutputTime::get_out_pointer()
{
    return m_out_pointer;
}

