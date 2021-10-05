#ifndef PATCH_HPP
#define PATCH_HPP

#include <vector>

struct Individual 
{
    bool is_hawk;
    double payoff;
    double prob_hawk;
};

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
