#pragma once
#include <map>
#include <iostream>
#include <cstdint>

namespace AbAl {
    template<typename Index, typename Field = double>
    class Polynomial
    {
    public:
        using Monomial = std::map<Index, std::uint16_t>;
    private:
        std::map<Monomial, Field> m_poly;

        Monomial mul_mon(const Monomial& left, const Monomial& right) const
        {
            Monomial result(left);
            for(auto&& el: right)
            {
                insert_power(result, el.first, el.second);
            }
            return result;
        }

        Monomial mul_mon(Monomial&& left, const Monomial& right) const
        {
            Monomial result(left);
            for(auto&& el: right)
            {
                insert_power(result, el.first, el.second);
            }
            return result;
        }

        Monomial mul_mon(const Monomial& left, Monomial&& right) const
        {
            Monomial result(left);
            for(auto&& el: right)
            {
                insert_power(result, el.first, el.second);
            }
            return result;
        }

        void insert_power(Monomial& mon, const Index& ind, const std::uint16_t& pw) const
        {
            if(ind.size() == 0)
                return;
            auto find_it = mon.find(ind);
            if(find_it != mon.end())
            {
                find_it->second += pw;
            }
            else
            {
                mon[ind] = pw;
            }
        }

        void clear_poly(std::map<Monomial, Field>& in) const
        {
            for(auto it = in.begin(); it!= in.end();)
            {
                if(it->second == static_cast<Field>(0.0))
                    it = in.erase(it);
                else
                    ++it;
            }
        }

    public:
        Polynomial(){}
        Polynomial(const Polynomial& in):m_poly(in.m_poly){}
        Polynomial(Polynomial&& in):m_poly(std::move(in.m_poly)){}
        Polynomial(const Index& in)
        {
            Monomial mon;
            if(in.size() != 0)
                mon[in] = 1;
            m_poly[mon] = static_cast<Field>(1.0);
        }

        Polynomial(Index&& in)
        {
            Monomial mon;
            if(in.size() != 0)
                mon[in] = 1;
            m_poly[mon] = static_cast<Field>(1.0);
        }

        Polynomial(const Field& in)
        {
            Monomial mon;
            m_poly[mon] = in;
        }
        Polynomial(Field&& in)
        {
            Monomial mon;
            m_poly[mon] = in;
        }

        Polynomial& operator=(const Polynomial& in)
        {
            m_poly = in.m_poly;
            return *this;
        }

        Polynomial& operator=(Polynomial&& in)
        {
            m_poly = std::move(in.m_poly);
            return *this;
        }

        Polynomial& operator=(const Index& in)
        {
            m_poly.clear();
            Monomial mon;
            mon[in] = 1;
            m_poly[mon] = static_cast<Field>(1.0);
            return *this;
        }

        Polynomial& operator=(Index&& in)
        {
            m_poly.clear();
            Monomial mon;
            mon[in] = 1;
            m_poly[mon] = static_cast<Field>(1.0);
            return *this;
        }

        Polynomial& operator=(const Field& in)
        {
            m_poly.clear();
            Monomial mon;
            m_poly[mon] = in;
            return *this;
        }

        Polynomial& operator=(Field&& in)
        {
            m_poly.clear();
            Monomial mon;
            m_poly[mon] = in;
            return *this;
        }

        Polynomial operator-()
        {
            Polynomial result(*this);
            for(auto&& el: result.m_poly)
            {
                el.second *= -static_cast<Field>(1.0);
            }
            return result;
        }


        Polynomial operator*(const Polynomial &in) const
        {
            Polynomial result;
            for(auto&& left: m_poly)
            {
                for(auto&& right: in.m_poly)
                {
                    auto mon_res = mul_mon(left.first, right.first);
                    auto find_it = result.m_poly.find(mon_res);
                    if(find_it != result.m_poly.end())
                    {
                        find_it->second += left.second * right.second;
                    }
                    else
                    {
                        result.m_poly[mon_res] = left.second * right.second;
                    }
                }
            }
            clear_poly(result.m_poly);
            return result;
        }

        Polynomial& operator*=(const Polynomial &in)
        {
            std::map<Monomial, Field> res_poly;
            for(auto&& left: m_poly)
            {
                for(auto&& right: in.m_poly)
                {
                    auto mon_res = mul_mon(left.first, right.first);
                    auto find_it = res_poly.find(mon_res);
                    if(find_it != res_poly.end())
                    {
                        find_it->second += left.second * right.second;

                    }
                    else
                    {
                       res_poly[mon_res] = left.second * right.second;
                    }
                }
            }
            m_poly = res_poly;
            clear_poly(m_poly);
            return *this;
        }

        Polynomial operator*(const Index &in) const
        {
            Polynomial result(*this);
            for(auto&& mon: result.m_poly)
            {
                insert_power(mon.first, in, 1);
            }
            return result;
        }

        Polynomial& operator*=(const Index &in)
        {
            std::map<Monomial, Field> result;
            for(auto&& mon: m_poly)
            {
                auto new_mon = mon.first;
                insert_power(new_mon, in, 1);
                result[new_mon] = mon.second;
            }
            m_poly = std::move(result);
            return *this;
        }

        Polynomial operator*(const Field &in) const
        {
            Polynomial result(*this);
            for(auto&& mon: result.m_poly)
            {
                mon.second *= in;
            }
            clear_poly(result.m_poly);
            return result;
        }

        Polynomial& operator*=(const Field &in)
        {
            for(auto&& mon: m_poly)
            {
                mon.second *= in;
            }
            clear_poly(m_poly);
            return *this;
        }

        Polynomial operator+(const Polynomial &in) const
        {
            Polynomial result(*this);
            for(auto&& summ : in.m_poly)
            {
                auto find_it = result.m_poly.find(summ.first);
                if(find_it != result.m_poly.end())
                {
                    find_it->second += summ.second;
                    if(find_it->second == static_cast<Field>(0.0))
                    {
                        result.m_poly.erase(find_it);
                    }
                }
                else
                {
                    result.m_poly[summ.first] = summ.second;
                }
            }
            return result;
        }

        Polynomial& operator+=(const Polynomial &in)
        {
            for(auto&& summ : in.m_poly)
            {
                auto find_it = m_poly.find(summ.first);
                if(find_it != m_poly.end())
                {
                    find_it->second += summ.second;
                    if(find_it->second == static_cast<Field>(0.0))
                    {
                        m_poly.erase(find_it);
                    }
                }
                else
                {
                    m_poly[summ.first] = summ.second;
                }
            }
            return *this;
        }

        Polynomial operator+(const Index &in) const
        {
            Polynomial result(*this);
            Monomial mon;
            mon[in] = 1;
            auto find_it = result.m_poly.find(mon);
            if(find_it != result.m_poly.end())
            {
                find_it->second += static_cast<Field>(1.0);
                if(find_it->second == static_cast<Field>(0.0))
                {
                    result.m_poly.erase(find_it);
                }
            }
            else
            {
                result.m_poly[mon] = static_cast<Field>(1.0);
            }
            return result;
        }

        Polynomial& operator+=(const Index &in)
        {
            Monomial mon;
            mon[in] = 1;
            auto find_it = m_poly.find(mon);
            if(find_it != m_poly.end())
            {
                find_it->second += static_cast<Field>(1.0);
                if(find_it->second == static_cast<Field>(0.0))
                {
                    m_poly.erase(find_it);
                }
            }
            else
            {
                m_poly[mon] = static_cast<Field>(1.0);
            }
            return *this;
        }

        Polynomial operator+(const Field &in) const
        {
            Polynomial result(*this);
            Monomial mon;
            auto find_it = result.m_poly.find(mon);
            if(find_it != result.m_poly.end())
            {
                find_it->second += in;
                if(find_it->second == static_cast<Field>(0.0))
                {
                    result.m_poly.erase(find_it);
                }
            }
            else
            {
                result.m_poly[mon] = in;
            }
            return result;
        }

        Polynomial& operator+=(const Field &in)
        {
            Monomial mon;
            auto find_it = m_poly.find(mon);
            if(find_it != m_poly.end())
            {
                find_it->second += in;
                if(find_it->second == static_cast<Field>(0.0))
                {
                    m_poly.erase(find_it);
                }
            }
            else
            {
                m_poly[mon] = in;
            }
            return *this;
        }

        std::map<Monomial, Field>& data()
        {
            return m_poly;
        }

        const std::map<Monomial, Field>& data()const
        {
            return m_poly;
        }
    };

    template<typename Index, typename Field = double>
    std::ostream& operator<<(std::ostream& os, const Polynomial<Index, Field>& dt)
    {
        for(auto&& mon : dt.data())
        {
            os << "+ " << mon.second << " ";
            for(auto&& factor: mon.first)
            {
                os << "<" << factor.first << ">";
                if(factor.second != 1)
                {
                    os << "^" << factor.second;
                }
                os << " ";
            }
        }
        os << std::endl;
        return os;
    }
}
