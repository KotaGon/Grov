#ifndef _GROV_H_
#define _GROV_H_

#include "constant.hpp"
#include "variable.hpp"
#include "constraint.hpp"
#include "utils.hpp"
#include "tensor.hpp"

class GrovSolver{
    private:
        
    public:
        GrovSolver() = default;
};

class Grov{

    private:
        std::vector<Variable> vars;
        std::vector<Constraint> constraints;
        LinVar objective;
        objective_sense sense;
    public:
        Grov(){ vars.reserve(100000000); constraints.reserve(100000000); }
        
        void clear();
        void setObjective(const LinVar &linvar, objective_sense optimize_sense);
        Variable addVar(double lb, double ub, var_type type, std::string name = "");
        void addConstrant(const Constraint &constraint, std::string name = "");
        void view();
        status_code optimize();
};

#endif 