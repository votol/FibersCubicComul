#pragma once
#include <optional>
#include <string>
#include <vector>
#include "eqManager.h"

class ModelBuilder
{
    clde::Polynomial m_intensity;
    clde::Polynomial m_main;
    clde::Polynomial m_init;
    clde::Polynomial m_third_order;
    clde::Polynomial m_third_order_sp;
    void m_build(const std::vector<double>&, const std::optional<std::ofstream*>& os);
    void m_save(std::ifstream& os);
public:
    ModelBuilder(const std::vector<double>& params, std::optional<std::string> cache_file = std::optional<std::string>());
    clde::Polynomial& intensity(void);
    clde::Polynomial& main_equations(void);
    clde::Polynomial& init_calculation(void);
    clde::Polynomial& third(void);
    clde::Polynomial& third_sp(void);
};
