#include "variable.hpp"

void LinVar::clear(){ 
    constant = 0.0;
    size = 0;
    coeffs.clear();
    vars.clear();
}

void LinVar::add_var(double coeff, const Variable &var) { 
    ++size;
    coeffs.push_back(coeff);
    vars.push_back(var);
}

void LinVar::view(){
    std::cout << std::endl;
    std::cout << constant << " + ";
    for(const auto &[coeff, var] : *this)
        std::cout << " + " << coeff << "*" << var.getName();
}

void LinVar::simplify() {
    std::unordered_map<std::string, Variable> var_map;
    std::unordered_map<std::string, double> coeff_map;
    for(const auto& [coeff, var] : *this){ 
        const std::string var_name = var.getName();
        var_map[var_name] = var;
        coeff_map[var_name] += coeff;
    }
    
    coeffs.clear();
    vars.clear();
    for (const auto& [var_name, var] : var_map) {
        double coeff = coeff_map[var_name];
        coeffs.push_back(coeff);
        vars.push_back(var);
    }
}

LinVar operator + (const LinVar &x, const LinVar &y){
    LinVar z = x; 
    z += y;
    return z;
}

LinVar operator - (const LinVar &x, const LinVar &y){
    LinVar z = x;
    z -= y;
    return z;
}

LinVar operator * (double c, const LinVar &x){ 
    LinVar z = x;
    z *= c;
    return z;
}

LinVar operator / (const LinVar &x, double c){
    return 1.0 / c * x;
}

void LinVar::operator += (const LinVar &x) { 
    constant += x.constant;
    for(const auto &[coeff, var] : x) add_var(coeff, var);
} 

void LinVar::operator -= (const LinVar &x) { 
    constant += -x.constant;
    for(const auto &[coeff, var] : x) add_var(-coeff, var);
} 

void LinVar::operator += (double c){
    constant += c;
}

void LinVar::operator -= (double c){
    constant -= c;
}

void LinVar::operator *= (double c){
    constant *= c;
    for (size_t i = 0; i < size; i++) coeffs[i] *= c; 
}

void LinVar::operator /= (double c){
    *this *= 1.0 / c;
}

std::ostream &operator << (std::ostream &stream, const LinVar &linvar){
    stream << linvar.constant;
    for(const auto &[coeff, var] : linvar)
        stream << " + " << coeff << "*" << var.getName();
    return stream;
}