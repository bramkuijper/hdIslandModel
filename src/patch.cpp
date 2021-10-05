#include <vector>
#include <Rcpp.h>
#include "patch.h"

Patch::Patch(int const Nbp)
    :
        breeders{Nbp}
        ,payoffs{Nbp}
        ,total_payoff{0.0}
{
} //end patch constructor


Patch::Patch(Patch const &other) :
    breeders{other.breeders}
    ,payoffs{other.payoffs}
    ,total_payoff{0.0}
{
}

void Patch::operator=(Patch const &other)
{
    breeders = other.breeders;
    payoffs= other.payoffs;
    total_payoff = other.total_payoff;
} // end Patch::operator=(
