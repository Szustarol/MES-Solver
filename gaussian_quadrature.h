#ifndef GAUSSIAN_QUADRATURE
#define GAUSSIAN_QUADRATURE
#include <cmath>
#include <array>
#include <functional>

template<int number_of_quadrature_points>
constexpr std::array<double, number_of_quadrature_points> get_quadrature_factors(){
    std::array<double, number_of_quadrature_points> factors;
    if constexpr(number_of_quadrature_points == 2){
        factors = {1, 1};
    }
    else if constexpr(number_of_quadrature_points == 1){
        factors = {2};
    }
    else if constexpr(number_of_quadrature_points == 4){
        factors = {(18+sqrt(30.0))/36, (18-sqrt(30.0))/36, (18+sqrt(30.0))/36, (18-sqrt(30.0))/36};    
    }
    else if constexpr(number_of_quadrature_points == 5){
        factors = {
            128/225, (322+13*sqrt(70))/900, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900, (322-13*sqrt(70))/900
        };
    }
    else if constexpr(number_of_quadrature_points == 6){
        factors = {
            0.467913934572691,  0.467913934572691,
            0.360761573048139,  0.360761573048139,
            0.171324492379170,  0.171324492379170
        };
    }
    else{
        factors = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    }
    return factors;
}

template<int number_of_quadrature_points>
constexpr std::array<double, number_of_quadrature_points> get_quadrature_points(){
    std::array<double, number_of_quadrature_points> points;
    if constexpr(number_of_quadrature_points == 2){
        points = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
    }
    else if constexpr(number_of_quadrature_points == 1){
        points = {0};
    }
    else if constexpr(number_of_quadrature_points == 4){
        points = {sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)), sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)),
            -sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)), -sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0))};
    }
    else if constexpr(number_of_quadrature_points == 5){
        points = {
            0,
            1.0/3*sqrt(5-2*sqrt(10/7)),
            -1.0/3*sqrt(5-2*sqrt(10/7)),
            1.0/3*sqrt(5+2*sqrt(10/7)),
            -1.0/3*sqrt(5+2*sqrt(10/7))
        };
    }
    else if constexpr(number_of_quadrature_points == 6){
        points = {
            0.238619186083197, -0.238619186083197,
            0.661209386466265, -0.661209386466265,
            0.932469514203152, -0.932469514203152
        };
    }
    else{
        points = {-sqrt(3.0/5.0), 0, sqrt(3.0/5.0)};
    }
    return points;
}


template<int number_of_quadrature_points>
double gaussian_quadrature_integral_core(std::function<double(double)> f, double from, double to){
    double scaling_factor = (to-from)/2;   
    double average_factor = (from+to)/2;
    static const auto factors = get_quadrature_factors<number_of_quadrature_points>();
    static const auto points = get_quadrature_points<number_of_quadrature_points>();
    double sum = 0;

    for(int i = 0; i < number_of_quadrature_points; i++){
        sum += factors.at(i)*(f(scaling_factor*points.at(i)+average_factor));
    }

    sum *= scaling_factor;
    return sum;
}

template<int number_of_quadrature_points = 3>
double gaussian_quadrature_integral(std::function<double(double)> f, double from, double to, double n_points = 100){
    static_assert(number_of_quadrature_points >= 1 and number_of_quadrature_points <= 6);
    double interval_width = (to-from)/(n_points-1);
    double sum = 0;
    for(int i = 0; i < n_points-1; i++){
        double start = from + i*interval_width;
        sum += gaussian_quadrature_integral_core<number_of_quadrature_points>(f, start, start + interval_width);
    }
    return sum;
}


#endif