#include <include/stack.h>

#include <memory.h>
#include <stdlib.h>

#define __DSTACK_GET(dstack, x) dstack.data[dstack.element_size * x]
#define __DSTACK_PTR_GET(dstack, x) dstack->data[dstack->element_size * x]

// make a empty dyn stack
DynStack dynStackInit(size_t elem_len)
{
    return (DynStack){NULL, 0, 0, elem_len};
}
// Reserve N element space
//DynStack dynStackReserve(size_t new_len);

// push element_size data to stack, WARNING: assumes data is element_size large!
void dynStackPush(DynStack* dstack, void* data)
{
    // if we have enough reserve space then, simply copy the data onto stack
    if(dstack->reserved_size < dstack->len + 1)
    {
        // if not empty, try to realloc
        // Possible point of failiure
        dstack->data = realloc(dstack->data, 2 * dstack->element_size * (dstack->reserved_size + 1));

        if(!dstack->data) printf("[Dyn Stack] Fatal Error: realloc returned null!\n");

        dstack->reserved_size = 2 * (dstack->reserved_size + 1);
    }

    memcpy(&__DSTACK_PTR_GET(dstack, dstack->len), data, dstack->element_size);
    dstack->len++;
}

// access a element in dyn stack, safely
void* dynStackGet(DynStack dstack, size_t index)
{
    if (index < dstack.len) return &__DSTACK_GET(dstack, index);
    else
    {
        printf("[Dyn Stack] Warning: attempt to access index %zu in stack of size %zu!\n", index, dstack.len);
        return NULL;
    }
}

// free a dyn stack
void freeDynStack(DynStack* dstack)
{
    if(dstack->data) free(dstack->data);

    // empty out the dstack
    *dstack = dynStackInit(0);
}


