#ifndef _CONSTANT_H_
#define _CONSTANT_H_

enum objective_sense{
    MINIMIZE,
    MAXIMIZE,
};
enum constraint_sense{
    EQUAL,
    GREATER_EQUAL,
    LESS_EQUAL,
};
enum var_type{ 
    CONTINUOUS,
    BINARY,
    INTERGER,
};

enum status_code{
    OPTIMAL,
};

#endif