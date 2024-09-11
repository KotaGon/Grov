#ifndef _VARIABLE_H_
#define _VARIABLE_H_

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "constant.hpp"
#include "iterator.hpp"

class Variable;
class LinVar;
class Grov;
class Constraint;

class Variable{
    private:
        double lb = 0.0, ub = 0.0;
        var_type type = CONTINUOUS;
        std::string name = "";      
        Variable *ptr_var = 0;  
        double solution = 0.0;
    public:
        Variable() = default;
        Variable(double lb, double ub, var_type type = CONTINUOUS, std::string name = "") :
            lb(lb), ub(ub), type(type), name(name){ 
            }
        const var_type getVarType() const { return type; }
        const double getLB() const { return lb; }
        const double getUB() const { return ub; }
        const std::string getName() const { return name; }
        const double getSolution() const { return ptr_var->solution; }

        friend class Grov;
};


LinVar operator + (const LinVar &x, const LinVar &y);
LinVar operator - (const LinVar &x, const LinVar &y);
LinVar operator * (double c, const LinVar &x);
LinVar operator / (const LinVar &x, double c);
std::ostream &operator << (std::ostream &stream, const LinVar &linvar);

class LinVar{ 
    private:
        // begin()とend()をテンプレートを使って実装
        using iterator = Iterator<std::vector<double>::iterator, std::vector<Variable>::iterator>;
        using const_iterator = Iterator<std::vector<double>::const_iterator, std::vector<Variable>::const_iterator>;

        double constant = 0.0;
        size_t size = 0;
        std::vector<double> coeffs;
        std::vector<Variable> vars;
    public:
        LinVar(double constant = 0.0) : constant(constant) { };
        LinVar(const Variable &var, double coeff = 1.0){
            add_var(coeff, var);
        }
        
        void clear();
        void add_var(double coeff, const Variable &var);
        void view();
        void simplify();
        
        iterator begin() { return iterator(coeffs.begin(), vars.begin()); }
        iterator end() { return iterator(coeffs.end(), vars.end()); }
        const_iterator begin() const { return const_iterator(coeffs.begin(), vars.begin()); }
        const_iterator end() const { return const_iterator(coeffs.end(), vars.end()); }

        void operator += (const LinVar &x);
        void operator -= (const LinVar &x);
        void operator *= (double c);
        void operator /= (double c);
        void operator += (double c);
        void operator -= (double c);
        friend LinVar operator + (const LinVar &x, const LinVar &y);
        friend LinVar operator - (const LinVar &x, const LinVar &y);
        friend LinVar operator * (double c, const LinVar &x);
        friend LinVar operator / (const LinVar &x, double c);
        friend std::ostream &operator << (std::ostream &stream, const LinVar &linvar);
        friend class Grov;
};

#endif 