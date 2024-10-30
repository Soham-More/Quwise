#pragma once

#include <stdio.h>
#include <stdint.h>

typedef struct DynStack
{
    uint8_t* data;
    size_t len;
    size_t reserved_size;
    size_t element_size;
} DynStack;

// make a empty dyn stack
DynStack dynStackInit(size_t elem_len);
// Reserve N element space
DynStack dynStackReserve(size_t new_len);

// push element_size data to stack, WARNING: assumes data is element_size large!
void dynStackPush(DynStack* dstack, void* data);

// access a element in dyn stack, safely
void* dynStackGet(DynStack dstack, size_t index);

// free a dyn stack, WARNING: Free memory in the dstack before freeing the stack
void freeDynStack(DynStack* dstack);


