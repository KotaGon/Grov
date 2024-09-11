#include "constraint.hpp"

void Constraint::clear(){ constant = 0.0; linvar.clear(); }
void Constraint::view(){
    std::cout << std::endl;
    std::cout << name << " : " << linvar;
    if(constraint_sense::EQUAL == sense) std::cout << "==";
    else if(constraint_sense::LESS_EQUAL == sense) std::cout << "<=";
    else if(constraint_sense::GREATER_EQUAL == sense) std::cout << ">="; 
    std::cout << constant;
}

Constraint operator <= (const LinVar &linvar, double constant){ 
    return Constraint(constant, linvar, constraint_sense::LESS_EQUAL);
}

Constraint operator == (const LinVar &linvar, double constant){ 
    return Constraint(constant, linvar, constraint_sense::EQUAL);
}

std::ostream &operator << (std::ostream &stream, const Constraint &constraint){
    stream << constraint.name << " : " << constraint.linvar;
    if(constraint_sense::EQUAL == constraint.sense) stream << "==";
    else if(constraint_sense::LESS_EQUAL == constraint.sense) stream << "<=";
    else if(constraint_sense::GREATER_EQUAL == constraint.sense) stream << ">="; 
    stream << constraint.constant;
    return stream;
}

