#include <iostream>
#include <map>
#include <set>
#include <list>
#include <complex>
#include <math.h>
#include "yaml-cpp/yaml.h"
#include "schema.h"
#include "NetCdfWriter.h"
#include "AbOA.h"
#include "AbComul.h"
#include "eqManager.h"
#include "ParametersInit.h"
#include "ParameterDefines.h"
#include "ModelBuilder.h"
#include "CLDataStorage.h"
#include "CLmanager.h"
#include "OutputWrapper.h"
#include "OutputTime.h"
#include "DERunge4.h"

using namespace AbAl;
using namespace clde;
using POperator_a = std::map<int, std::uint8_t, std::greater<int>>;

namespace AbAl
{
std::ostream& operator<<(std::ostream& os, const POperator_a& oper)
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
}

/*Polynomial<POperator_a, int> antinormilize(const POperator& in)
{
    Polynomial<POperator_a, int> result;
    GOperator tmp_in;
    tmp_in += in;
    while(tmp_in.size() != 0)
    {
        auto curr_prod = tmp_in.begin()->first;
        int coe = static_cast<int>(tmp_in.begin()->second.real());
        tmp_in.erase(tmp_in.begin());
        POperator positive, negative;
        for(auto&& el : curr_prod)
        {
            if(el.first > 0)
                positive[el.first] = el.second;
            else
                negative[el.first] = el.second;
        }
        tmp_in += commute(negative, positive) * coe;
        POperator_a summ;
        summ.insert(positive.begin(), positive.end());
        summ.insert(negative.begin(), negative.end());
        Polynomial<POperator_a, int> tmp(summ);
        tmp *= coe;
        result += tmp;
    }
    return result;
}

void insert_part(POperator_a& in, const int& part)
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

Polynomial<POperator_a, int> substitute_a(const Polynomial<std::vector<unsigned int>, int>& in, const std::vector<int>& values)
{
    Polynomial<POperator_a, int> result;
    Polynomial<POperator_a, int> res_sum;
    POperator_a tmp_index;
    for(auto&& in_sum : in.data())
    {
        res_sum = Polynomial<POperator_a, int>(1);
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


POperator conj(const POperator& in)
{
    POperator result;
    for(auto&& particle: in)
    {
        result[-particle.first] = particle.second;
    }
    return result;
}

GOperator res;
std::map<POperator, unsigned int> av_map;
std::list<GOperator> eqs;
unsigned int current_ind = 0;

void add_operator(const POperator& in)
{
    auto find_it = av_map.find(in);
    if(find_it != av_map.end())
        return;

    find_it = av_map.find(conj(in));
    if(find_it != av_map.end())
        return;

    av_map[in] = current_ind;
    ++ current_ind;
    eqs.push_back( std::complex<double>(0.0, 1.0) * commute(res, in));
}

void perform_test()
{
    unsigned int num = 64;
    for(unsigned int ind = 1; ind != num + 1; ++ind)
    {
        if(ind != num)
            res += ak(ind) * a(ind + 1) + ak(ind + 1) * a(ind);
        res += 0.1 * ak(ind) * ak(ind) * a(ind) * a(ind);
    }

    for(unsigned int ind = 1; ind != num + 1; ++ind)
    {
        add_operator(a(ind));

    }
    for(unsigned int ind = 1; ind != num + 1; ++ind)
    {
        for(unsigned int ind1 = 1; ind1 != num + 1; ++ind1)
        {
            for(unsigned int ind2 = 1; ind2 != num + 1; ++ind2)
            {
                add_operator((a(ind) * a(ind1) * a(ind2)).begin()->first);
                add_operator((ak(ind) * a(ind1) * a(ind2)).begin()->first);
                add_operator((ak(ind) * ak(ind1) * a(ind2)).begin()->first);
                add_operator((ak(ind) * ak(ind1) * ak(ind2)).begin()->first);
            }
        }
    }
    std::cout << eqs.size() << " ; " << av_map.size() << std::endl;
}

void perform_test1(std::vector<unsigned int>& params)
{
    using lIndex = std::set<unsigned int>;
    using Poly = Polynomial<lIndex, double>;
    auto u = Poly(lIndex{1,2,3});
    auto v = Poly(lIndex{4,5});
    auto res = (u + v) * (u +(-v));
    //res *= res;
    for(auto&& el:res.data())
    {
        for(auto&& ind: el.first)
        {
            std::cout << "{";
            for(auto&& part: ind.first)
            {
                std::cout << part << ",";
            }
            std::cout << "}^" << ind.second << " ";
        }
        std::cout << " : " <<el.second << std::endl;
    }
    unsigned int com_order = 3;
    ComulManager manager(com_order);
    POperator oper = (a(params[0]) * a(params[2]) * a(params[3]) * a(params[4])).begin()->first;
    auto poly1 = manager.split(oper);
    Polynomial<POperator_a, int> right;
    //std::cout << poly1;
    for(auto&& elem: poly1.data())
    {
        Polynomial<POperator_a, int> tmp(elem.second);
        for(auto&& factor: elem.first)
        {
            Polynomial<POperator_a, int> conv_tmp = antinormilize(factor.first);
            for(std::uint16_t ind =0 ;ind < factor.second; ++ind)
                tmp *= conv_tmp;
        }
        right += tmp;
    }

    auto left_tmp =  antinormilize(oper);
    //std::cout << left_tmp;
    Polynomial<POperator_a, int> left;
    for(auto&& summ: left_tmp.data())
    {
        Polynomial<POperator_a, int> tmp(summ.second);
        for(auto&& factor: summ.first)
        {
            Polynomial<POperator_a, int> new_factor;
            unsigned int power = 0;
            for(auto&& particle: factor.first)
            {
                power += particle.second;
            }

            if(power <= com_order)
                new_factor = factor.first;
            else
            {
                std::vector<int> lin_elements(power);
                auto liner_it = lin_elements.begin();
                for(auto&& particle: factor.first)
                {
                    for(size_t ind = 0; ind < particle.second; ++ind)
                    {
                        *liner_it = particle.first;
                        ++liner_it;
                    }
                }
                new_factor = substitute_a(manager.get_commul(power), lin_elements);
            }

            for(std::uint16_t ind =0 ;ind < factor.second; ++ind)
                tmp *= new_factor;
        }
        left += tmp;
    }
    //std::cout << right;
    left += right * (-1);
    std::cout << left;
    //auto poly = manager.split((ak(1) * a(1) * a(1) * a(1) * a(2)).begin()->first);
}*/

void print(const clde::Polynomial& in)
{
    for(auto it = in.begin();it != in.end();++it)
    {
        std::cout << *it;
        std::cout<<std::endl;
    }
    std::cout<<"----------------------"<<std::endl;
}

std::vector<double> calculate_init(const std::vector<double> in,
                                   const clde::Polynomial& eqs,
                                   const std::shared_ptr<ICLmanager>& context)
{
    CLDataStorage<double> in_vector(in, context);
    OperatorDimension operDim = PolynomialOperator::calculateDimension(
                eqs.begin(), eqs.end(), false);
    operDim.in_dim = in.size() - 1;
    CLDataStorage<double> out_vector(static_cast<unsigned int>(operDim.out_dim + 1), context);
    PolynomialOperator oper(eqs.begin(),
                            eqs.end(),
                            1,
                            operDim,
                            context);

    oper.apply(in_vector, out_vector, std::vector<double>());
    return out_vector.read();
}

std::vector<double> load_init(const std::vector<double>& parameters, const std::string& path)
{
    size_t nfibs = static_cast<size_t>(Nfibs);
    std::vector<double> result(1 + 2 * nfibs);
    std::ifstream ifs;
    double tmpd;

    ifs.open (path + "/init.bin", std::ifstream::in | std::ifstream::binary);

    for (unsigned int ind = 0; ind < nfibs; ind ++){
        ifs.read(reinterpret_cast<char *>(&tmpd), 8);
        result[1 + ind *2] = tmpd;
        //vec[1 + ind *2] = 0.0;
        result[2 + ind *2] = 0.0;
    }

    //vec[DI_r(0,0)] = 100.0;
    result[0] = 1.0;
    ifs.close();
    return result;
}


int main(int argc, char **argv)
{
    YAML::Node config = YAML::LoadFile(argv[1]);
    std::string output_dir = config["properties"]
								  [FibersCubicComulSchema::PROPERTY_output_path].as<std::string>();

    ParametersHolder parameters_holder_instance(config["parameters"]);
    std::shared_ptr<ICLmanager> manag = std::make_shared<CLmanager>(config["properties"]);

    ModelBuilder model(parameters_holder_instance.GetParameters(), argv[2]);
    //make main equations
    OperatorDimension operDim = PolynomialOperator::calculateDimension(
                model.main_equations().begin(), model.main_equations().end(), true);
    PolynomialOperator oper(model.main_equations().begin(),
                            model.main_equations().end(), 1, operDim, manag);

    //make outputs
    PolyOutputWrapper intens(model.intensity(),
                             operDim,
                             config["parameters"][FibersCubicComulSchema::PARAMETER_Nout].as<unsigned int>(),
                             FibersCubicComulSchema::OUTPUT_I,
                             manag);
    PolyOutputWrapper third(model.third(),
                             operDim,
                             config["parameters"][FibersCubicComulSchema::PARAMETER_Nout].as<unsigned int>(),
                             FibersCubicComulSchema::OUTPUT_TO,
                             manag);

    PolyOutputWrapper third_sp(model.third_sp(),
                             operDim,
                             config["parameters"][FibersCubicComulSchema::PARAMETER_Nout].as<unsigned int>(),
                             FibersCubicComulSchema::OUTPUT_TO_s,
                             manag);

    OutputTime time(config["parameters"][FibersCubicComulSchema::PARAMETER_Nout].as<unsigned int>(),
                    FibersCubicComulSchema::OUTPUT_time);

    std::list<IDEOutput*> c_outputs;
    c_outputs.push_back(intens.get_calc_pointer().get());
    c_outputs.push_back(third.get_calc_pointer().get());
    c_outputs.push_back(third_sp.get_calc_pointer().get());
    c_outputs.push_back(time.get_calc_pointer().get());

    // make init
    std::vector<double> alphas = load_init(parameters_holder_instance.GetParameters(), config["properties"]
            [FibersCubicComulSchema::PROPERTY_tmp_path].as<std::string>());
    std::vector<double> init_state = calculate_init(alphas,
                                                    model.init_calculation(),
                                                    manag);
    init_state[0] = 1.0;

    //std::cout << "Dimension: " << init_state.size() << std::endl;
    //std::cout << "Calc elements: " << model.main_equations().size() << std::endl;

    DERunge4 calc(manag, &oper);
    calc.SetTimeStep(config["parameters"][FibersCubicComulSchema::PARAMETER_dt].as<double>());
    calc.SetStepsNumber(config["parameters"][FibersCubicComulSchema::PARAMETER_Nsteps].as<unsigned int>());
    calc.SetOutputSteps(config["parameters"][FibersCubicComulSchema::PARAMETER_Nout].as<unsigned int>());
    calc.SetInitState(init_state);
    calc.SetOutputs(c_outputs);
    calc.calculate();

    std::vector<std::shared_ptr<IOutput> > outputs(4);
    outputs[0] = intens.get_out_pointer();
    outputs[1] = third.get_out_pointer();
    outputs[2] = third_sp.get_out_pointer();
    outputs[3] = time.get_out_pointer();
    NetCdfWriter netcdf_writer_instance(
                output_dir + "/output.nc",
                outputs,
                config["parameters"][FibersCubicComulSchema::PARAMETER_Nout].as<unsigned int>());
    return 0;
}
