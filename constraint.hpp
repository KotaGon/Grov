#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include "constant.hpp"
#include "variable.hpp"

class Grov;
class Constraint;

Constraint operator <= (const LinVar &linvar, double constant); 
Constraint operator == (const LinVar &linvar, double constant); 
std::ostream &operator << (std::ostream &stream, const Constraint &constraint);

class Constraint{ 
    private:
        double constant = 0;
        LinVar linvar;
        constraint_sense sense = constraint_sense::EQUAL;
        std::string name = "";
    public:  
        
        Constraint(double constant = 0.0) : constant(constant) { }
        Constraint(double constant, const LinVar &linvar, constraint_sense sense)
         : constant(constant), linvar(linvar), sense(sense) { }

        void clear();
        void view();
        
        friend Constraint operator <= (const LinVar &linvar, double constant); 
        friend Constraint operator == (const LinVar &linvar, double constant); 
        friend std::ostream &operator << (std::ostream &stream, const Constraint &constraint);

        friend class Grov;

};

#endif  