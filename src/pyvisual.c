#include <include/pyvisual.h>

typedef struct PyViSectionData
{
    // section name
    const char* name;
    // all the Vec data, in this section
    DynStack/*Vec*/ data;
    PyViBase parameter;
} PyViSectionData;

PyVi pyviInitA(const char* filename)
{
    PyVi pyvi;
    pyvi.file = fopen(filename, "w");
    pyvi.sections = dynStackInit(sizeof(PyViSectionData));
    pyvi.parameters = dynStackInit(sizeof(PyViBase));

    return pyvi;
}

PyViSec pyviCreateSection(PyVi* pyvi, const char* section_name, PyViBase p)
{
    PyViSectionData section;
    section.name = section_name;
    section.data = dynStackInit(sizeof(Vec));
    section.parameter = p;

    dynStackPush(&pyvi->sections, &section);

    return (PyViSec){pyvi->sections.len - 1, pyvi};
}

PyViBase pyviCreateParameter(PyVi* pyvi, const char* param_name, Vec p)
{
    PyViBase param;
    param.axis = p;
    param.name = param_name;

    dynStackPush(&pyvi->parameters, &param);

    return param;
}

// push a vector fx varying with parameter x
void pyviSectionPush(PyViSec section, Vec fx)
{
    Vec tmp = vecCopyA(fx);

    PyViSectionData* sec = dynStackGet(section.pyvi->sections, section.id);

    // simply push to stack
    dynStackPush(&sec->data, &tmp);
}

void freePyVi(PyVi* pyvi)
{
    // close the file
    fclose(pyvi->file);
    pyvi->file = NULL;

    for(size_t i = 0; i < pyvi->sections.len; i++)
    {
        PyViSectionData* section = dynStackGet(pyvi->sections, i);
        for(size_t j = 0; j < section->data.len; j++)
        {
            freeVec(dynStackGet(section->data, j));
        }
        freeDynStack(&section->data);
    }

    freeDynStack(&pyvi->sections);
    freeDynStack(&pyvi->parameters);
}

void pyviPrintVec(PyVi pyvi, Vec x)
{
    for(size_t i = 0; i < x.len; i++)
    {
        fprintf(pyvi.file, "%.17g", *vecRef(x, i));
        if(i != x.len - 1) fprintf(pyvi.file, ",");
    }
}

// writes all the data to file
void pyviWrite(PyVi pyvi)
{
    // first add all parameters
    fprintf(pyvi.file, "[Parameters]\n");
    for(size_t i = 0; i < pyvi.parameters.len; i++)
    {
        PyViBase parameter = *(PyViBase*)dynStackGet(pyvi.parameters, i);
        
        fprintf(pyvi.file, "%s:", parameter.name);
        pyviPrintVec(pyvi, parameter.axis);
        fprintf(pyvi.file, "\n");
    }

    fprintf(pyvi.file, "[Sections]\n");

    // then add all sections
    for(size_t i = 0; i < pyvi.sections.len; i++)
    {
        PyViSectionData section = *(PyViSectionData*)dynStackGet(pyvi.sections, i);
        
        fprintf(pyvi.file, "(%s)->[%s]\n", section.name, section.parameter.name);
        
        for(size_t j = 0; j < section.data.len; j++)
        {
            fprintf(pyvi.file, "I[%zu]=", j); pyviPrintVec(pyvi, *(Vec*)dynStackGet(section.data, j));
            fprintf(pyvi.file, "\n");
        }
        fprintf(pyvi.file, "\n");
    }
}


