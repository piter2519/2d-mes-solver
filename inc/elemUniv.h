#pragma once
#include "gauss.h"

struct Surface
{
    double **N;
};

struct ElemUniv 
{
    double **dNdKsi;
    double **dNdEta;
    double **N;
    Surface s[4];
    ElemUniv(unsigned int npc);
    void funkcjaKsztaltu(unsigned int npc);
};