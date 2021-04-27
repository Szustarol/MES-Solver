/*

    Autor: Karol Szustakowski

    Implementacja używa następujących bibliotek:
    mathgl

    Program był kompilowany na systemie Fedora Linux, są potrzebne 
    następujące pakiety:
    - mathgl-devel (libgml-dev na ubuntu)

    Kompilacja

    g++ plik.cpp -std=c++17 -lmgl

*/

#include <iostream>
#include <cmath>
#include <array>
#include <functional>
#include <mgl2/mgl.h>
#include "gaussian_quadrature.h"
#include "linear_equations.h"
#include "probing_functions_state_machine.h"

double f(double x){
    return log(abs(x)+0.1);
}

void plot(std::function<double(double)> f, double from, double to, const std::string & output_file, double n_points,
    int alg_n, int alg_n_inrevrals, int alg_roots){
    mglGraph graph;
    graph.SetSize(1920, 1080);
    std::string title = std::to_string(alg_n) + std::string(" probing functions,\n ") + std::to_string(alg_n_inrevrals)
    + std::string(" integration points, \n") + std::to_string(alg_roots) + std::string(" Legendre polynomial roots");
    graph.LoadFont("adventor font");
    graph.SetFontSize(2.3);
    graph.Title(title.c_str());
    //graph.SetOrigin(0, 0);
    graph.Light(true);
    mglData graph_data_x(n_points);
    mglData graph_data_y(n_points);

    double interval_width = (to-from)/(n_points-1);

    double max_y, min_y;
    max_y = f(from);
    min_y = min_y;

    for(int i = 0; i < n_points; i++){
        double function_x = from + i * interval_width;
        double function_eval = f(function_x);
        graph_data_x.a[i] = function_x;
        graph_data_y.a[i] = function_eval;

        max_y = std::max(max_y, function_eval);
        min_y = std::min(min_y, function_eval);
    }

    double margin = 0.1;
    graph.SetRanges(from-margin, to+margin, min_y-margin, max_y+margin);
    graph.SetOrigin(0, 0);
    graph.Plot(graph_data_x, graph_data_y);
    graph.Axis();
    graph.WritePNG(output_file.c_str());

}

int main(){

    std::cout << "STARTING SOLVER" << std::endl;

    int n = 200;
    int number_of_integration_intervals = 150;
    std::string output_file = "test.png";

    linear_gauss_matrix_solver l_eq(n-1);

    std::cout << "CONSTRUCTING PROBING FUNCTIONS..." << std::endl;

    probing_functions_state_machine s(0, 2, n);
    s.set_k([](double x){
        if(x <= 1)
            return 1;
        return 2;
    });
    std::cout << "CALCULATING INTEGRALS..." << std::endl;

    //calculate B(e0, e0):
    double b_e0_e1 = s.B(0, 1, 0, 2, number_of_integration_intervals);
    double b_e1_e0 = s.B(1, 0, 0, 2, number_of_integration_intervals);

    l_eq.set(1, 0, b_e0_e1);
    l_eq.set(0, 1, b_e1_e0);

    for(int u = 0; u < n-2; u++){
        for(int v = 0; v < n-2; v++){
            if(abs(u-v) > 1)
                continue;
            //values of B aren't symmetric around x = 1
            double r = s.B(u, v, 0, 2, number_of_integration_intervals);
            l_eq.set(u, v, r);
        }
    }

    std::vector<double> R_side(n-1);

    for(int i = 0; i < n-1; i++){
        l_eq.set(i, i, s.B(i, i, 0, 2, number_of_integration_intervals));
        R_side.at(i) = 0;
    }

    R_side.at(0) = -20;

    l_eq.set_result_vector<std::vector<double>>(R_side);

    //l_eq._debug_print();

    std::cout << "SOLVING LINEAR EQUATIONS..." << std::endl;

    l_eq.solve();

    auto solution = l_eq.get_solution();

    s.set_bias([](double x){return 0;});

    s.set_coefs(solution);

    auto f = [&s](double x){
        return s.result_at(x);
    };
    std::cout << "PLOTTING.." << std::endl;

    plot(f, 0, 2, output_file, 500, n, number_of_integration_intervals, 6);

    std::cout << "DONE" << std::endl;
    std::cout << "OUTPUT SAVED TO: " << output_file << std::endl;
}