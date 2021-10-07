#ifndef PATCH_HPP
#define PATCH_HPP

#include "individual.h"
#include <vector>


// Patch struct
class Patch 
{
    public:
        std::vector <Individual> breeders;
        std::vector <double> payoffs;
        double total_payoff;

        Patch(int const Nbp);

        Patch(Patch const &other);

        void operator=(Patch const &other);
};


#endif
