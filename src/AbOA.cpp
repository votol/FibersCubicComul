#include "AbOA.h"
#include <optional>
#include <queue>
#include <iostream>

using namespace AbAl;

struct _AncStorage
{
    POperator left;
    POperator right;
    int coe;
};

void clear_null(GOperator& op)
{
    for(auto it = op.begin(); it != op.end();)
    {
        if(it->second == 0.0)
        {
            it = op.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

std::map<POperator, int> mult_processor(std::queue<_AncStorage>& tasks)
{
    std::map<POperator, int> tmp_res;
    while(!tasks.empty())
    {
        auto current = tasks.front();
        tasks.pop();
        for(auto elem = current.right.begin();elem != current.right.end();)
        {
            if(elem->first < 0)
            {
                if(current.left.find(-elem->first) != current.left.end())
                {
                    for(uint8_t power = 0; power < elem->second; ++power)
                    {
                        _AncStorage new_summand;
                        new_summand.left = current.left;
                        new_summand.right = current.right;
                        new_summand.coe = current.coe;
                        new_summand.left[-elem->first] -= 1;
                        if(new_summand.left[-elem->first] == 0)
                            new_summand.left.erase(-elem->first);
                        new_summand.right[elem->first] -= power + 1;
                        if(new_summand.right[elem->first] == 0)
                            new_summand.right.erase(elem->first);
                        new_summand.coe *= current.left[-elem->first];
                        tasks.push(new_summand);
                        current.left[elem->first] += 1;
                    }
                }
                else
                {
                    auto find_it = current.left.find(elem->first);
                    if(find_it != current.left.end())
                    {
                        find_it -> second += elem->second;
                    }
                    else
                    {
                        current.left[elem->first] = elem->second;
                    }
                }
            }
            else if(elem->first > 0)
            {
                auto find_it = current.left.find(elem->first);
                if(find_it != current.left.end())
                {
                    find_it -> second += elem->second;
                }
                else
                {
                    current.left[elem->first] = elem->second;
                }
            }
            elem = current.right.erase(elem);
        }
        auto ins_it = tmp_res.find(current.left);
        if(ins_it != tmp_res.end())
            ins_it->second += current.coe;
        else
            tmp_res[current.left] = current.coe;
    }
    return tmp_res;
}


GOperator AbAl::operator*(const POperator& left, const POperator& right)
{
    std::queue<_AncStorage> tasks;
    _AncStorage init;
    init.left=left;
    init.right=right;
    init.coe = 1;
    tasks.push(init);

    std::map<POperator, int> tmp_res = mult_processor(tasks);

    GOperator result;
    for(auto&& elem: tmp_res)
    {
        result[elem.first] = std::complex<double>(static_cast<double>(elem.second));
    }
    return result;
}

GOperator AbAl::operator*(const std::complex<double>& coe, const POperator& op)
{
    GOperator result;
    result[op] = coe;
    return result;
}

GOperator AbAl::operator*(const POperator& op, const std::complex<double>& coe)
{
    return coe * op;
}

GOperator AbAl::operator*(const double& coe, const POperator& op)
{
    GOperator result;
    result[op] = std::complex<double>(coe);
    return result;
}

GOperator AbAl::operator*(const POperator& op, const double& coe)
{
    return coe * op;
}

GOperator AbAl::operator*(const int& coe, const POperator& op)
{
    GOperator result;
    result[op] = std::complex<double>(static_cast<double>(coe));
    return result;
}

GOperator AbAl::operator*(const POperator& op, const int& coe)
{
    return coe * op;
}

GOperator operator/(const POperator& op, const std::complex<double>& coe)
{
    GOperator result;
    result[op] = 1.0 / coe;
    return result;
}

GOperator operator/(const POperator& op, const double& coe)
{
    GOperator result;
    result[op] = 1.0 / std::complex<double>(coe);
    return result;
}

GOperator operator/(const POperator& op, const int& coe)
{
    GOperator result;
    result[op] = 1.0 / std::complex<double>(static_cast<double>(coe));
    return result;
}

GOperator AbAl::operator*(const GOperator& left, const GOperator& right)
{
    GOperator result;
    for(auto&& left_it: left)
    {
        for(auto&& right_it: right)
        {
            auto mul_res = left_it.first * right_it.first;
            for(auto&& mul_it: mul_res)
            {
                auto res_coe = left_it.second * right_it.second * mul_it.second;
                auto ins_it = result.find(mul_it.first);
                if(ins_it != result.end())
                {
                    ins_it->second += res_coe;
                }
                else
                {
                    result[mul_it.first] = res_coe;
                }
            }
        }
    }

    clear_null(result);
    return result;
}

GOperator AbAl::operator*(const std::complex<double>& coe, const GOperator& op)
{
    GOperator result(op);
    for(auto&& el: result)
    {
        el.second *= coe;
    }
    return result;
}

GOperator AbAl::operator*(const GOperator& op, const std::complex<double>& coe)
{
    return coe * op;
}

GOperator AbAl::operator*(const double& coe, const GOperator& op)
{
    GOperator result(op);
    for(auto&& el: result)
    {
        el.second *= std::complex<double>(coe);
    }
    return result;
}

GOperator AbAl::operator*(const GOperator& op, const double& coe)
{
    return coe * op;
}

GOperator AbAl::operator*(const int& coe, const GOperator& op)
{
    GOperator result(op);
    for(auto&& el: result)
    {
        el.second *= std::complex<double>(static_cast<double>(coe));
    }
    return result;
}

GOperator AbAl::operator*(const GOperator& op, const int& coe)
{
    return coe * op;
}

GOperator AbAl::operator/(const GOperator& op, const std::complex<double>& coe)
{
    GOperator result(op);
    for(auto&& el: result)
    {
        el.second /= coe;
    }
    return result;
}

GOperator AbAl::operator/(const GOperator& op, const double& coe)
{
    GOperator result(op);
    for(auto&& el: result)
    {
        el.second /= std::complex<double>(coe);
    }
    return result;
}

GOperator AbAl::operator/(const GOperator& op, const int& coe)
{
    GOperator result(op);
    for(auto&& el: result)
    {
        el.second /= std::complex<double>(static_cast<double>(coe));
    }
    return result;
}

GOperator AbAl::operator*(const GOperator& left, const POperator& right)
{
    GOperator result;
    for(auto&& elem: left)
    {
        GOperator tmp_res = elem.first * right;
        for(auto&& res_elem: tmp_res)
        {
            auto coe = res_elem.second * elem.second;
            auto find_it = result.find(res_elem.first);
            if(find_it != result.end())
            {
                find_it->second += coe;
            }
            else
            {
                result[res_elem.first] = coe;
            }
        }
    }
    clear_null(result);
    return result;
}

GOperator AbAl::operator*(const POperator& left, const GOperator& right)
{
    GOperator result;
    for(auto&& elem: right)
    {
        GOperator tmp_res = left * elem.first;
        for(auto&& res_elem: tmp_res)
        {
            auto coe = res_elem.second * elem.second;
            auto find_it = result.find(res_elem.first);
            if(find_it != result.end())
            {
                find_it->second += coe;
            }
            else
            {
                result[res_elem.first] = coe;
            }
        }
    }
    clear_null(result);
    return result;
}

GOperator& operator*=(GOperator& dst, const POperator& src)
{
    auto tmp = dst * src;
    dst = tmp;
    return dst;
}

GOperator& operator*=(GOperator& dst, const GOperator& src)
{
    auto tmp = dst * src;
    dst = tmp;
    return dst;
}

GOperator& operator*=(GOperator& dst, const std::complex<double>& src)
{
    for(auto&& el: dst)
    {
        el.second *= src;
    }
    clear_null(dst);
    return dst;
}

GOperator& operator*=(GOperator& dst, const double& src)
{
    for(auto&& el: dst)
    {
        el.second *= std::complex<double>(src);
    }
    clear_null(dst);
    return dst;
}

GOperator& operator*=(GOperator& dst, const int& src)
{
    for(auto&& el: dst)
    {
        el.second *= std::complex<double>(static_cast<double>(src));
    }
    clear_null(dst);
    return dst;
}

GOperator AbAl::operator-(const POperator& in)
{
    GOperator result;
    result[in] = std::complex<double>(-1.0, 0.0);
    return result;
}

GOperator AbAl::operator-(const GOperator& in)
{
    GOperator result(in);
    for(auto&& el: result)
    {
        el.second *= -1.0;
    }
    return result;
}

GOperator AbAl::operator+(const POperator& left, const POperator& right)
{
    GOperator result;
    result[left] = std::complex<double>(1.0);
    auto find_it = result.find(right);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(1.0);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else
        result[right] = std::complex<double>(1.0);

    return result;
}

GOperator AbAl::operator+(const std::complex<double>& num, const POperator& op)
{
    GOperator result;
    POperator fake_op;
    result[op] = std::complex<double>(1.0);
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += num;
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = num;
    return result;
}

GOperator AbAl::operator+(const POperator& op, const std::complex<double>& num)
{
    return num + op;
}

GOperator AbAl::operator+(const double& num, const POperator& op)
{
    GOperator result;
    POperator fake_op;
    result[op] = std::complex<double>(1.0);
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(num);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(num);
    return result;
}

GOperator AbAl::operator+(const POperator& op, const double& num)
{
    return num + op;
}

GOperator AbAl::operator+(const int& num, const POperator& op)
{
    GOperator result;
    POperator fake_op;
    result[op] = std::complex<double>(1.0);
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(static_cast<double>(num));
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(static_cast<double>(num));
    return result;
}

GOperator AbAl::operator+(const POperator& op, const int& num)
{
    return num + op;
}

GOperator AbAl::operator-(const POperator& left, const POperator& right)
{
    GOperator result;
    result[left] = std::complex<double>(1.0);
    auto find_it = result.find(right);
    if(find_it != result.end())
    {
        find_it->second -= std::complex<double>(1.0);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else
        result[right] = std::complex<double>(-1.0);
    return result;
}

GOperator AbAl::operator-(const std::complex<double>& num, const POperator& op)
{
    GOperator result;
    POperator fake_op;
    result[op] = std::complex<double>(-1.0);
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += num;
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = num;
    return result;
}

GOperator AbAl::operator-(const POperator& op, const std::complex<double>& num)
{
    return -num + op;
}

GOperator AbAl::operator-(const double& num, const POperator& op)
{
    GOperator result;
    POperator fake_op;
    result[op] = std::complex<double>(-1.0);
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(num);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(num);
    return result;
}

GOperator AbAl::operator-(const POperator& op, const double& num)
{
    return -num + op;
}

GOperator AbAl::operator-(const int& num, const POperator& op)
{
    GOperator result;
    POperator fake_op;
    result[op] = std::complex<double>(-1.0);
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(static_cast<double>(num));
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(static_cast<double>(num));
    return result;
}

GOperator AbAl::operator-(const POperator& op, const int& num)
{
    return -num + op;
}

GOperator AbAl::operator+(const GOperator& left, const GOperator& right)
{
    GOperator result(left);
    for(auto&& el: right)
    {
        auto find_it = result.find(el.first);
        if(find_it != result.end())
            find_it->second += el.second;
        else
            result.insert(el);
    }
    clear_null(result);
    return result;
}

GOperator AbAl::operator+(const std::complex<double>& num, const GOperator& op)
{
    GOperator result(op);
    POperator fake_op;
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += num;
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = num;
    return result;
}

GOperator AbAl::operator+(const GOperator& op, const std::complex<double>& num)
{
    return num + op;
}

GOperator AbAl::operator+(const double& num, const GOperator& op)
{
    GOperator result(op);
    POperator fake_op;
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(num);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(num);
    return result;
}

GOperator AbAl::operator+(const GOperator& op, const double& num)
{
    return num + op;
}

GOperator AbAl::operator+(const int& num, const GOperator& op)
{
    GOperator result(op);
    POperator fake_op;
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(static_cast<double>(num));
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(static_cast<double>(num));
    return result;
}

GOperator AbAl::operator+(const GOperator& op, const int& num)
{
    return num + op;
}

GOperator AbAl::operator-(const GOperator& left, const GOperator& right)
{
    GOperator result(left);
    for(auto&& el: right)
    {
        auto find_it = result.find(el.first);
        if(find_it != result.end())
            find_it->second -= el.second;
        else
            result[el.first] = -el.second;
    }
    clear_null(result);
    return result;
}

GOperator AbAl::operator-(const std::complex<double>& num, const GOperator& op)
{
    GOperator result(op);
    POperator fake_op;
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second -= num;
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = -num;
    return result;
}

GOperator AbAl::operator-(const GOperator& op, const std::complex<double>& num)
{
    return -num + op;
}

GOperator AbAl::operator-(const double& num, const GOperator& op)
{
    GOperator result(op);
    POperator fake_op;
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second -= std::complex<double>(num);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(-num);
    return result;
}

GOperator AbAl::operator-(const GOperator& op, const double& num)
{
    return -num + op;
}

GOperator AbAl::operator-(const int& num, const GOperator& op)
{
    GOperator result(op);
    POperator fake_op;
    auto find_it = result.find(fake_op);
    if(find_it != result.end())
    {
        find_it->second -= std::complex<double>(static_cast<double>(num));
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else if(num != 0.0)
        result[fake_op] = std::complex<double>(static_cast<double>(-num));
    return result;
}

GOperator AbAl::operator-(const GOperator& op, const int& num)
{
    return -num + op;
}

GOperator AbAl::operator+(const POperator& left, const GOperator& right)
{
    GOperator result(right);
    auto find_it = result.find(left);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(1.0);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else
        result[left] = std::complex<double>(1.0);
    return result;
}

GOperator AbAl::operator+(const GOperator& left, const POperator& right)
{
    GOperator result(left);
    auto find_it = result.find(right);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(1.0);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else
        result[right] = std::complex<double>(1.0);
    return result;
}

GOperator AbAl::operator-(const POperator& left, const GOperator& right)
{
    GOperator result(-right);
    auto find_it = result.find(left);
    if(find_it != result.end())
    {
        find_it->second += std::complex<double>(1.0);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else
        result[left] = std::complex<double>(1.0);
    return result;
}

GOperator AbAl::operator-(const GOperator& left, const POperator& right)
{
    GOperator result(left);
    auto find_it = result.find(right);
    if(find_it != result.end())
    {
        find_it->second -= std::complex<double>(1.0);
        if(find_it->second == 0.0)
            result.erase(find_it);
    }
    else
        result[right] = std::complex<double>(-1.0);
    return result;
}

GOperator& AbAl::operator+=(GOperator& dst, const POperator& src)
{
    auto find_it = dst.find(src);
    if(find_it != dst.end())
    {
        find_it->second += std::complex<double>(1.0);
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else
        dst[src] = std::complex<double>(1.0);
    return dst;
}

GOperator& AbAl::operator+=(GOperator& dst, const GOperator& src)
{
    for(auto&& el: src)
    {
        auto find_it = dst.find(el.first);
        if(find_it != dst.end())
            find_it->second += el.second;
        else
            dst.insert(el);
    }
    clear_null(dst);
    return dst;
}

GOperator& AbAl::operator+=(GOperator& dst, const std::complex<double>& src)
{
    POperator fake_op;
    auto find_it = dst.find(fake_op);
    if(find_it != dst.end())
    {
        find_it->second += src;
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else if(src != std::complex<double>(0.0))
        dst[fake_op] = src;
    return dst;
}

GOperator& AbAl::operator+=(GOperator& dst, const double& src)
{
    POperator fake_op;
    auto find_it = dst.find(fake_op);
    if(find_it != dst.end())
    {
        find_it->second += std::complex<double>(src);
    }
    else if(src != 0.0)
        dst[fake_op] = std::complex<double>(src);
    return dst;
}

GOperator& AbAl::operator+=(GOperator& dst, const int& src)
{
    POperator fake_op;
    auto find_it = dst.find(fake_op);
    if(find_it != dst.end())
    {
        find_it->second += std::complex<double>(static_cast<double>(src));
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else if(src != 0)
        dst[fake_op] = std::complex<double>(static_cast<double>(src));
    return dst;
}

GOperator& AbAl::operator-=(GOperator& dst, const POperator& src)
{
    auto find_it = dst.find(src);
    if(find_it != dst.end())
    {
        find_it->second -= std::complex<double>(1.0);
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else
        dst[src] = std::complex<double>(-1.0);
    return dst;
}

GOperator& AbAl::operator-=(GOperator& dst, const GOperator& src)
{
    for(auto&& el: src)
    {
        auto find_it = dst.find(el.first);
        if(find_it != dst.end())
            find_it->second -= el.second;
        else
            dst[el.first] = - el.second;
    }
    clear_null(dst);
    return dst;
}

GOperator& AbAl::operator-=(GOperator& dst, const std::complex<double>& src)
{
    POperator fake_op;
    auto find_it = dst.find(fake_op);
    if(find_it != dst.end())
    {
        find_it->second -= src;
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else if(src != 0.0)
        dst[fake_op] = -src;
    return dst;
}

GOperator& AbAl::operator-=(GOperator& dst, const double& src)
{
    POperator fake_op;
    auto find_it = dst.find(fake_op);
    if(find_it != dst.end())
    {
        find_it->second -= std::complex<double>(src);
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else if(src != 0.0)
        dst[fake_op] = std::complex<double>(-src);
    return dst;
}

GOperator& AbAl::operator-=(GOperator& dst, const int& src)
{
    POperator fake_op;
    auto find_it = dst.find(fake_op);
    if(find_it != dst.end())
    {
        find_it->second -= std::complex<double>(static_cast<double>(src));
        if(find_it->second == 0.0)
            dst.erase(find_it);
    }
    else if(src != 0.0)
        dst[fake_op] = std::complex<double>(static_cast<double>(-src));
    return dst;
}

std::ostream& AbAl::operator<<(std::ostream& os, const POperator& oper)
{
    for(auto&& elem: oper)
    {
        os << "{" << elem.first << "}";
        if(elem.second != 1)
        {
            os << "^" << int(elem.second);
        }
        os << " ";
    }
    return os;
}

std::ostream& AbAl::operator<<(std::ostream& os, const GOperator& oper)
{
    for(auto&& elem: oper)
    {
        os << elem.first;
        os << ": ";
        auto coe = elem.second;
        if(coe.real() != 0.0 && coe.imag() != 0.0)
        {
            os << coe.real();
            if(coe.imag() > 0)
                os<< "+";

            if(coe.imag() == 1.0)
            {
                os << "i";
            }
            else if(coe.imag() == -1.0)
            {
                os << "-i";
            }
            else
            {
                os << coe.imag() << "i";
            }
        }
        else if(coe.real() == 0.0 && coe.imag() != 0.0)
        {
            if(coe.imag() == 1.0)
            {
                os << "i";
            }
            else if(coe.imag() == -1.0)
            {
                os << "-i";
            }
            else
            {
                os << coe.imag() << "i";
            }
        }
        else if(coe.real() != 0.0 && coe.imag() == 0.0)
        {
            os << coe.real();
        }
        else
        {
            os << "0";
        }
        os <<std::endl;
    }
    return os;
}

GOperator AbAl::commute(const GOperator& left, const GOperator& right)
{
    GOperator result;
    for(auto&& left_it: left)
    {
        for(auto&& right_it: right)
        {
            auto mul_res = commute(left_it.first, right_it.first);
            for(auto&& mul_it: mul_res)
            {
                auto res_coe = left_it.second * right_it.second * mul_it.second;
                auto ins_it = result.find(mul_it.first);
                if(ins_it != result.end())
                {
                    ins_it->second += res_coe;
                }
                else
                {
                    result[mul_it.first] = res_coe;
                }
            }
        }
    }

    clear_null(result);
    return result;
}

GOperator AbAl::commute(const GOperator& left, const POperator& right)
{
    GOperator result;
    for(auto&& elem: left)
    {
        GOperator tmp_res = commute(elem.first, right);
        for(auto&& res_elem: tmp_res)
        {
            auto coe = res_elem.second * elem.second;
            auto find_it = result.find(res_elem.first);
            if(find_it != result.end())
            {
                find_it->second += coe;
            }
            else
            {
                result[res_elem.first] = coe;
            }
        }
    }
    clear_null(result);
    return result;
}

GOperator AbAl::commute(const POperator& left, const GOperator& right)
{
    GOperator result;
    for(auto&& elem: right)
    {
        GOperator tmp_res = commute(left, elem.first);
        for(auto&& res_elem: tmp_res)
        {
            auto coe = res_elem.second * elem.second;
            auto find_it = result.find(res_elem.first);
            if(find_it != result.end())
            {
                find_it->second += coe;
            }
            else
            {
                result[res_elem.first] = coe;
            }
        }
    }
    clear_null(result);
    return result;
}

GOperator AbAl::commute(const POperator& left, const POperator& right)
{
    std::queue<_AncStorage> tasks;
    _AncStorage init;
    init.left=left;
    init.right=right;
    init.coe = 1;
    tasks.push(init);
    init.left=right;
    init.right=left;
    init.coe = -1;
    tasks.push(init);

    std::map<POperator, int> tmp_res = mult_processor(tasks);

    GOperator result;
    for(auto&& elem: tmp_res)
    {
        result[elem.first] = std::complex<double>(static_cast<double>(elem.second));
    }

    return result;
}

POperator AbAl::a(const unsigned int& ind)
{
    POperator result;
    result[static_cast<int>(ind)] = 1;
    return result;
}

POperator AbAl::ak(const unsigned int& ind)
{
    POperator result;
    result[-static_cast<int>(ind)] = 1;
    return result;
}
