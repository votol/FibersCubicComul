#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <optional>
#include <memory>
#include "AbComul.h"

using namespace AbAl;

Polynomial<std::vector<unsigned int>, int> substitute(const Polynomial<std::vector<unsigned int>, int>& in, const std::vector<unsigned int>& values)
{
    Polynomial<std::vector<unsigned int>, int> result;
    Polynomial<std::vector<unsigned int>, int> res_sum;
    std::vector<unsigned int> tmp_index;
    for(auto&& in_sum : in.data())
    {
        res_sum = Polynomial<std::vector<unsigned int>, int>(1);
        for(auto&& factor: in_sum.first)
        {
            tmp_index.resize(factor.first.size());
            for(size_t ind = 0; ind < factor.first.size(); ++ind)
            {
                tmp_index[ind] = values[factor.first[ind] - 1];
            }
            for(size_t ind = 0; ind < factor.second; ++ind)
                res_sum *= tmp_index;
        }
        result += res_sum * in_sum.second;
    }
    return result;
}

void insert_part(POperator& in, const int& part)
{
    auto find_it = in.find(part);
    if(find_it != in.end())
    {
        find_it->second += 1;
    }
    else
    {
        in[part] = 1;
    }
}

Polynomial<POperator, int> substitute(const Polynomial<std::vector<unsigned int>, int>& in, const std::vector<int>& values)
{
    Polynomial<POperator, int> result;
    Polynomial<POperator, int> res_sum;
    POperator tmp_index;
    for(auto&& in_sum : in.data())
    {
        res_sum = Polynomial<POperator, int>(1);
        for(auto&& factor: in_sum.first)
        {
            tmp_index.clear();
            for(size_t ind = 0; ind < factor.first.size(); ++ind)
            {
                insert_part(tmp_index, values[factor.first[ind] - 1]);
            }
            for(size_t ind = 0; ind < factor.second; ++ind)
                res_sum *= tmp_index;
        }
        result += res_sum * in_sum.second;
    }
    return result;
}


class ComulBuilder
{
    class ISetManager
    {
    public:
        virtual ~ISetManager(){}
        virtual std::shared_ptr<std::vector<unsigned int>> GetState() = 0;
        virtual std::shared_ptr<std::vector<unsigned int>> GetRemainder() = 0;
        virtual bool next() = 0;
    };

    class CnkSetManager: public ISetManager
    {
        struct Task
        {
            unsigned int n;
            unsigned int k;
            unsigned int ind;
            std::vector<unsigned int> state;
            std::vector<unsigned int> remainder;
        };

        std::optional<ISetManager*> m_parent;
        std::shared_ptr<std::vector<unsigned int>> m_parent_list;
        std::shared_ptr<std::vector<unsigned int>> m_state;
        std::shared_ptr<std::vector<unsigned int>> m_remainder;
        unsigned int m_n;
        unsigned int m_k;

        std::queue<Task> tasks;

        bool m_build()
        {
            if(tasks.size() == 0)
                return false;
            while(true)
            {
                if(tasks.front().k == 0)
                {
                    for(unsigned int ins_ind = tasks.front().ind; ins_ind < m_n; ++ins_ind)
                    {
                        tasks.front().remainder[ins_ind - m_k] = (*m_parent_list)[ins_ind];
                    }
                    *m_state = std::move(tasks.front().state);
                    *m_remainder = std::move(tasks.front().remainder);
                    tasks.pop();
                    return true;
                }

                tasks.push(Task());
                tasks.back().n = tasks.front().n - 1;
                tasks.back().k = tasks.front().k - 1;
                tasks.back().ind = tasks.front().ind + 1;
                if(tasks.front().k >= tasks.front().n)
                {
                    tasks.back().state = std::move(tasks.front().state);
                    tasks.back().remainder = std::move(tasks.front().remainder);
                }
                else
                {
                    tasks.back().state = tasks.front().state;
                    tasks.back().remainder = tasks.front().remainder;
                }
                tasks.back().state[m_k - tasks.front().k] = (*m_parent_list)[tasks.front().ind];



                if(tasks.front().k < tasks.front().n)
                {
                    tasks.push(Task());
                    tasks.back().n = tasks.front().n - 1;
                    tasks.back().k = tasks.front().k;
                    tasks.back().ind = tasks.front().ind + 1;
                    tasks.back().state = std::move(tasks.front().state);
                    tasks.back().remainder = std::move(tasks.front().remainder);
                    tasks.back().remainder[m_n - m_k + tasks.front().k - tasks.front().n] = (*m_parent_list)[tasks.front().ind];
                }

                tasks.pop();
            }
        }

        void m_initialize()
        {
            tasks.push(Task());
            tasks.back().n = m_n;
            tasks.back().k = m_k;
            tasks.back().ind = 0;
            tasks.back().state.resize(m_k);
            tasks.back().remainder.resize(m_n - m_k);
            m_build();
        }

    public:
        CnkSetManager(unsigned int k, std::shared_ptr<std::vector<unsigned int>>& source,
                     std::optional<ISetManager*> parent = std::optional<ISetManager*>()):
            m_parent(parent), m_parent_list(source), m_n(static_cast<unsigned int>(source->size())), m_k(k)
        {
            m_state = std::make_shared<std::vector<unsigned int>>();
            m_remainder = std::make_shared<std::vector<unsigned int>>();
            m_initialize();
        }
        virtual std::shared_ptr<std::vector<unsigned int>> GetState() override
        {
            return m_state;
        }

        virtual std::shared_ptr<std::vector<unsigned int>> GetRemainder() override
        {
            return m_remainder;
        }

        virtual bool next() override
        {
            if(!m_build())
            {
                if(m_parent && (*m_parent)->next())
                {
                    m_initialize();
                    return true;
                }
                else
                    return false;
            }
            else
                return true;
        }
    };

    class FirstSetManager: public ISetManager
    {
        std::optional<ISetManager*> m_parent;
        std::shared_ptr<std::vector<unsigned int>> m_parent_list;
        std::shared_ptr<std::vector<unsigned int>> m_state;
        std::shared_ptr<std::vector<unsigned int>> m_remainder;

        void m_initialize()
        {
            (*m_state)[0] = (*m_parent_list)[0];

            m_remainder->clear();
            m_remainder->insert(m_remainder->end(), ++m_parent_list->begin(), m_parent_list->end());
        }

    public:
        FirstSetManager(std::shared_ptr<std::vector<unsigned int>> source, std::optional<ISetManager*> parent = std::optional<ISetManager*>()):
            m_parent(parent), m_parent_list(source)
        {
            m_state = std::make_shared<std::vector<unsigned int>>(1);
            m_remainder = std::make_shared<std::vector<unsigned int>>(source->size() - 1);
            m_initialize();
        }
        virtual std::shared_ptr<std::vector<unsigned int>> GetState() override
        {
            return m_state;
        }

        virtual std::shared_ptr<std::vector<unsigned int>> GetRemainder() override
        {
            return m_remainder;
        }

        virtual bool next() override
        {
            if(m_parent && (*m_parent)->next())
            {
                m_initialize();
                return true;
            }
            else
                return false;
        }
    };

    struct ProductTask
    {
        int n;
        int k;
        std::map<unsigned int, unsigned int> res;
    };

    std::list<std::map<unsigned int, unsigned int>> m_products;
    void gen_products(const unsigned int& order)
    {
        std::queue<ProductTask> tasks;
        tasks.push(ProductTask());
        tasks.back().k = static_cast<int>(order);
        tasks.back().n = static_cast<int>(order);

        while(tasks.size() != 0)
        {
            if ( tasks.front().n == 0 )
            {
                m_products.push_back(std::move(tasks.front().res));
            }
            else
            {
                if ( tasks.front().n - tasks.front().k >= 0)
                {
                    tasks.push(ProductTask());
                    tasks.back().n = tasks.front().n - tasks.front().k;
                    tasks.back().k = tasks.front().k;
                    if(tasks.front().k - 1 > 0)
                        tasks.back().res = tasks.front().res;
                    else
                        tasks.back().res = std::move(tasks.front().res);
                    auto inserting = static_cast<unsigned int>(tasks.front().k);
                    auto find_it = tasks.back().res.find(inserting);
                    if(find_it != tasks.back().res.end())
                        find_it->second += 1;
                    else
                        tasks.back().res[inserting] = 1;
                }

                if ( tasks.front().k - 1 > 0)
                {
                    tasks.push(ProductTask());
                    tasks.back().n = tasks.front().n;
                    tasks.back().k = tasks.front().k - 1;
                    tasks.back().res = std::move(tasks.front().res);
                }
            }
            tasks.pop();
        }
    }

    void fill_patern(const std::map<unsigned int, unsigned int>& schema,
                     std::list<std::list<std::vector<unsigned int>>>& storage
                     )
    {
        //make initial vector with all factors
        unsigned int part_count = 0;
        for(auto&& group: schema)
        {
            part_count += group.first * group.second;
        }
        std::shared_ptr<std::vector<unsigned int>> all_factors = std::make_shared<std::vector<unsigned int>>(part_count);
        for(unsigned int factor_ind = 1; factor_ind <= part_count; ++factor_ind)
            (*all_factors)[factor_ind - 1] = factor_ind;

        //building generation structure
        size_t output_number = 0;
        std::vector<std::unique_ptr<CnkSetManager> > groups_distributor(schema.size());
        unsigned int cnk_ind = 0;
        std::shared_ptr<std::vector<unsigned int>> & current_set = all_factors;
        std::optional<ISetManager*> gener_parent;
        for(auto&& group: schema)
        {
            output_number += group.second;
            groups_distributor[cnk_ind] = std::make_unique<CnkSetManager>(
                        group.second * group.first, current_set, gener_parent);
            current_set = groups_distributor[cnk_ind]->GetRemainder();
            gener_parent = groups_distributor[cnk_ind].get();
            cnk_ind++;
        }

        std::vector<std::unique_ptr<CnkSetManager> > groups_second_extr(output_number);
        std::vector<std::unique_ptr<FirstSetManager> > groups_first_extr(output_number);

        unsigned int process_ind = 0;
        cnk_ind = 0;
        for(auto&& group: schema)
        {
            current_set = groups_distributor[cnk_ind]->GetState();
            for(unsigned int in_gr_ind = 0; in_gr_ind < group.second; ++in_gr_ind)
            {
                groups_first_extr[process_ind] = std::make_unique<FirstSetManager>(
                            current_set, gener_parent);
                gener_parent = groups_first_extr[process_ind].get();
                current_set = groups_first_extr[process_ind]->GetRemainder();

                groups_second_extr[process_ind] = std::make_unique<CnkSetManager>(
                            group.first - 1, current_set, gener_parent);

                gener_parent = groups_second_extr[process_ind].get();
                current_set = groups_second_extr[process_ind]->GetRemainder();
                process_ind++;
            }
            cnk_ind++;
        }

        //generating
        do
        {
            storage.push_back(std::list<std::vector<unsigned int>>());
            for(size_t gener_ind = 0 ; gener_ind < groups_first_extr.size(); ++gener_ind)
            {
                storage.back().push_back(std::vector<unsigned int>((*groups_first_extr[gener_ind]->GetState())));
                storage.back().back().insert(storage.back().back().end(),
                                             groups_second_extr[gener_ind]->GetState()->begin(),
                                             groups_second_extr[gener_ind]->GetState()->end());
            }
        }while((*gener_parent)->next());
    }

public:
    static void print_element(const std::list<std::vector<unsigned int>>& in)
    {
        for(auto&& factor: in)
        {
            std::cout << "<";
            for(auto it = factor.begin(); it != factor.end(); ++it )
            {
                std::cout << *it;
                if(it != --factor.end())
                    std::cout << " ";
            }
            std::cout << "> ";
        }
        std::cout << std::endl;
    }

    ComulBuilder(){}
    std::list<std::list<std::vector<unsigned int>>> generate(const unsigned int& order)
    {
        std::list<std::list<std::vector<unsigned int>>> result;
        gen_products(order);
        for(auto&& pattern: m_products)
        {
            if(pattern.find(order) == pattern.end())
                fill_patern(pattern, result);
        }
        return result;
    }
};

void print(const std::set<std::set<unsigned int>>& in)
{
    for(auto&& row: in)
    {
        for(auto&& el: row)
        {
            std::cout << el << " ";
        }
        std::cout << std::endl;
    }
}

namespace AbAl
{
std::ostream& operator<<(std::ostream& os, const std::vector<unsigned int>& dt)
{
    for(auto it = dt.begin(); it != dt.end(); ++it)
    {
        os << *it;
        if(it != --dt.end())
            os << " ";
    }
    return os;
}
}

Polynomial<std::vector<unsigned int>, int> ComulManager::m_build_comul_patern(const unsigned int& order)
{
    std::vector<int> coes(order - 1);
    int cur_coe = 1;
    for(size_t ind = 0; ind < order -1; ++ind)
    {
        coes[ind] = cur_coe;
        cur_coe *= -static_cast<int>(ind + 2);
    }

    ComulBuilder builder;
    auto splittings = builder.generate(order);

    Polynomial<std::vector<unsigned int>, int> result;
    Polynomial<std::vector<unsigned int>, int> tmp_poly;

    for(auto&& splitting: splittings)
    {
        tmp_poly = Polynomial<std::vector<unsigned int>, int>(1);
        for(auto&& factor: splitting)
        {
            if(factor.size() > m_order)
            {
                auto tmp = substitute(m_storage[factor.size() - m_order - 1], factor);
                //std::cout << tmp;
                tmp_poly *= substitute(m_storage[factor.size() - m_order - 1], factor);
            }
            else
                tmp_poly *= factor;
        }
        result += tmp_poly * coes[splitting.size() - 2];
    }
    //std::cout << result;
    return result;
}



Polynomial<POperator, int> ComulManager::split(const POperator& in)
{
    Polynomial<POperator, int> result;
    unsigned int power = 0;
    for(auto&& particle: in)
    {
        power += particle.second;
    }
    std::vector<int> lin_elements(power);
    auto liner_it = lin_elements.begin();
    for(auto&& particle: in)
    {
        for(size_t ind = 0; ind < particle.second; ++ind)
        {
            *liner_it = particle.first;
            ++liner_it;
        }
    }

    if(power <= m_order)
    {
        result = in;
    }
    else
    {
        if(m_storage.size() < power - m_order)
        {
            size_t last_comul = m_storage.size();
            m_storage.resize(power - m_order);
            for(size_t cur_comul = last_comul; cur_comul < m_storage.size(); ++cur_comul)
            {
                m_storage[cur_comul] = m_build_comul_patern(static_cast<unsigned int>(cur_comul + m_order + 1));
            }
        }
        result = substitute(m_storage[power - m_order - 1], lin_elements);

    }
    return result;
}
