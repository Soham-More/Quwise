#include "include/linalg.h"

void error_handler(const char* file, const char* function, size_t line_no)
{
    printf("in function %s, defined in file %s at line %zu:\n\t", function, file, line_no);
}


