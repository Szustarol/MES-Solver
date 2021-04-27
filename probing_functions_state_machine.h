#ifndef PROBING_FUNCTIONS_STATE_MACHINE
#define PROBING_FUNCTIONS_STATE_MACHINE

#include <vector>
#include <functional>
#include "gaussian_quadrature.h"

class probing_functions_state_machine{

    std::vector<std::function<double(double)>> functions;
    std::vector<std::function<double(double)>> derivatives;
    std::function<double(double)> k;

    int u_index, v_index;

    double h;

    std::vector<double> coefs;

    std::function<double(double)> bias;

    std::vector<std::pair<double, double>> defined_in_range;

    std::pair<double, double> find_smallest_integration_range(int u_idx, int v_idx, double from, double to){
        return std::make_pair(
            std::max({defined_in_range.at(u_idx).first, defined_in_range.at(v_idx).first, from}),
            std::min({defined_in_range.at(u_idx).second, defined_in_range.at(v_idx).second, to})
        );
    }


    double d_v(double x){
        return derivatives.at(v_index)(x);
    }

    double d_u(double x){
        return derivatives.at(u_index)(x);
    }

    double k_du_dv(double x){
        return k(x)*d_u(x)*d_v(x);
    }

public:

    probing_functions_state_machine(double from, double to, double n_splits){
        double dist = to-from;
        double h = dist/(n_splits-1);//e0

        this->h = h;

        functions.emplace_back([h](double x){
            if(x >= 0 and x <= h){
                return 1 - (x/h);
            }
            return (double)0;
        });

        derivatives.emplace_back([h](double x){
            if(x >= 0 and x <= h)
                return -1/h;
            return (double)0;
        });

        defined_in_range.emplace_back(std::make_pair(0, h));


        for(int i = 1; i < n_splits-1; i++){

            auto e = [h, i](double x){
                if(x >= h*(i-1) and x <= h*i){
                    return (x/h-i)+1;
                }
                else if(x >= h*i and x <= h*(i+1)){
                    return i-x/h+1;
                }
                else{
                    return (double)0;
                }
            };

            auto d_e = [h, i](double x){
                if(x == h*i)
                    return (double)0;
                if(x >= h*(i-1) and x <= h*i){
                    return 1/h;
                }
                else if(x >= h*i and x <= h*(i+1)){
                    return -1/h;
                }
                else{
                    return (double)0;
                }
            };

            derivatives.emplace_back(d_e);
            functions.emplace_back(e);

            defined_in_range.emplace_back(std::make_pair(h*(i-1), h*(i+1)));
        }
    }
    
    std::function<double(double)> get_function(int function){
        return functions.at(function);
    }

    void set_bias(std::function<double(double)> bias){
        this->bias = bias;
    }

    void set_coefs(std::vector<double> coefs){
        this->coefs = coefs;
    }

    double result_at(double x){
        double r = bias(x);
        for(size_t i = 0; i < coefs.size(); i++){
            r += functions.at(i)(x)*coefs.at(i);
        }
        return r;
    }

    void set_u(int index){
        u_index = index;
    }

    void set_v(int index){
        v_index = index;
    }

    void set_k(std::function<double(double)> k){
        this->k = k;
    }

    double probe(int function, double x){
        return functions.at(function)(x);
    }

    double u(double x){
        return functions.at(u_index)(x);
    }

    double v(double x){
        return functions.at(v_index)(x);
    }

    double B(double idx_u, double idx_v, double from, double to, int integration_points = 30){
        set_u(idx_u);
        set_v(idx_v);

        std::function<double(double)> f = [this](double x){
            return this->k_du_dv(x);
        };
        
        auto range = find_smallest_integration_range(idx_u, idx_v, from, to);
        
        return  d_u(1)*v(1) -u(0)*v(0) + gaussian_quadrature_integral<6>(f, range.first, range.second, integration_points);
    }


};

#endif