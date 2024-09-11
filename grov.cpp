#include "grov.hpp"

void Grov::clear(){ 
    vars.clear();
    constraints.clear();
}

void Grov::setObjective(const LinVar &linvar, objective_sense optimize_sense){ 
    objective = linvar;
    objective.simplify();
    sense = optimize_sense;
}

Variable Grov::addVar(double lb, double ub, var_type type, std::string name){
    if(name == "") name = "var" + stringParse(vars.size());
    Variable var = Variable(lb, ub, type, name);
    vars.push_back(var);
    var.ptr_var = &vars.back();
    return var;
}

void Grov::addConstrant(const Constraint &constraint, std::string name){ 
    if(name == "") name = "cons" + stringParse(constraints.size());
    constraints.push_back(constraint);
    auto *ptr_const = &constraints.back();
    ptr_const->linvar.simplify();
    ptr_const->name = name;
    ptr_const->constant -= ptr_const->linvar.constant;
    ptr_const->linvar.constant = 0;
}

void Grov::view() { 
    std::cout << std::endl;
    std::cout << "Variables " << vars.size() << std::endl;
    std::cout << "Constraints " << constraints.size() << std::endl;
    for(auto &constr : constraints) std::cout << constr << std::endl;
}

status_code Grov::optimize(){



    return OPTIMAL;
}




