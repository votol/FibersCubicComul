#pragma once
#include <map>
#include <cstdint>
#include "AbOA.h"

class IIndexManage
{
public:
    virtual ~IIndexManage(){}
    virtual int GetIndex(const AbAl::POperator& ) = 0;
    virtual const std::map<AbAl::POperator, std::uint32_t> GetAllData(const AbAl::POperator& ) = 0;
};


class EqManager
{
    public:
};
