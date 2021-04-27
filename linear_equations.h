#ifndef LINEAR_EQUATIONS
#define LINEAR_EQUATIONS

#include <algorithm>
#include <vector>
#include <cmath>


class linear_solver_base{
protected:
    std::vector<std::vector<double>> coefficientsMatrix;
    std::vector<double> resultVector;

public:

    virtual void set(int row, int col, double value) = 0;
    virtual double get(int row, int col) = 0;

    template<typename T>
    void set_row(int row_index, const T & container){
        int idx = 0;
        for(auto & value : container){
            set(row_index, idx, value);
        }
    }

    template<typename T>
    void set_result_vector(const T & container){
        int idx = 0;
        for(auto & value : container){
            resultVector.at(idx++) = value;
        }
    }
};

class linear_sparse_laplace_solver : public linear_solver_base{

//maybe todo, maybe not - it will probably be very slow
private:
    typedef std::vector<std::vector<std::vector<std::vector<double>>>> dynamic_matrix;
    double determinant = -1;

    double det(dynamic_matrix & matrix, int from_x, int from_y, int width){
        //expand - if any of those is zero, no need to calculate det
        return 0;  
    };

public:

    void set(int row, int col, double value){
        coefficientsMatrix.at(row).at(col-row) = value;
    }

    double get(int row, int col){
        return coefficientsMatrix.at(row).at(col-row);
    }

    linear_sparse_laplace_solver(int number_of_equations){
        coefficientsMatrix.resize(number_of_equations);
        for(auto & row : coefficientsMatrix){
            row.resize(3, 0);
        }
        resultVector.resize(number_of_equations, 0);
    }

    void solve(){
        //minor [fromX, fromY, toX, toY]
        dynamic_matrix minor_values;
        minor_values.resize(coefficientsMatrix.size());
        for(auto & fromX : minor_values){
            fromX.resize(coefficientsMatrix.size());
            for(auto & fromY : fromX){
                fromY.resize(coefficientsMatrix.size());
                for(auto & toX : fromY){
                    toX.resize(coefficientsMatrix.size());
                }
            }
        }

        int matrix_size = coefficientsMatrix.size();

        for(int minor_size = 1; minor_size <= matrix_size; minor_size++){
            for(int start_x = 0; start_x < matrix_size-minor_size+1; start_x++){
                for(int start_y = 0; start_y < matrix_size- minor_size+1; start_y++){
                    minor_values.at(start_x).at(start_y).at(start_x + minor_size).at(start_y+minor_size)
                        = det(minor_values, start_x, start_y, minor_size);
                }
            }
        }
    }
};

class linear_gauss_matrix_solver : public linear_solver_base{


    void swap_rows(int row_1, int row_2){
        std::swap(
            coefficientsMatrix.at(row_1),
            coefficientsMatrix.at(row_2)
        );
        std::swap(
            resultVector.at(row_1),
            resultVector.at(row_2)
        );
    }

    void add_row_to_row(int adding_from, int adding_to, double factor = 1){
        auto & from = coefficientsMatrix.at(adding_from);
        auto & to = coefficientsMatrix.at(adding_to);
        for(size_t idx = 0; idx < coefficientsMatrix.size(); idx++){
            to.at(idx) += factor*from.at(idx);
        }
        resultVector.at(adding_to) += factor*resultVector.at(adding_from);
    }

public:
    void set(int row, int col, double value){
        coefficientsMatrix.at(row).at(col) = value;
    }

    double get(int row, int col){
        return coefficientsMatrix.at(row).at(col);
    }

    linear_gauss_matrix_solver(int number_of_equations){
        coefficientsMatrix.resize(number_of_equations);
        for(auto & row : coefficientsMatrix){
            row.resize(number_of_equations, 0);
        }
        resultVector.resize(number_of_equations, 0);
    }

    void _debug_print(){
        for(auto & row : coefficientsMatrix){
            for(auto & cell : row){
                std::cout << "\t" << cell;
            }
            std::cout << std::endl;
        }
    }

    void solve(){
        int n_rows = coefficientsMatrix.size();
        for(int eval_from = 0; eval_from < n_rows; eval_from++){
            int max_index = eval_from;
            double max_value = coefficientsMatrix.at(max_index).at(eval_from);

            for(int potential_idx = max_index+1; potential_idx < n_rows; potential_idx++){
                if(abs(coefficientsMatrix.at(potential_idx).at(eval_from)) > abs(max_value)){
                    max_value = coefficientsMatrix.at(potential_idx).at(eval_from);
                    max_index = potential_idx;
                }
            }

            if(max_value == 0)
                throw std::invalid_argument("Trying to solve a singular matrix");

            if(eval_from != max_index){
                swap_rows(eval_from, max_index);
            }

            for(int to_normalize = eval_from + 1; to_normalize < n_rows; to_normalize++){
                double ratio = coefficientsMatrix.at(to_normalize).at(eval_from)
                    / max_value;
                add_row_to_row(eval_from, to_normalize, -ratio);
                coefficientsMatrix.at(to_normalize).at(eval_from) = 0; // in case rounding errors make this non-zero
            }
        }
    }

    std::vector<double> get_solution(){
        int n_variables = coefficientsMatrix.size();

        std::vector<double> solution(n_variables);

        for(int i = n_variables-1; i >= 0; i--){
            double ratio = coefficientsMatrix.at(i).at(i);
            double sum = resultVector.at(i);
            for(int j = i+1; j < n_variables; j++){
                sum -= solution.at(j)*coefficientsMatrix.at(i).at(j);
            }
            solution.at(i) = sum / ratio;
        }

        return solution;
    }
};

#endif