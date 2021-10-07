#include <vector>
#include <cassert>
#include "patch.h"

Patch::Patch(int const Nbp)
    :
        breeders(Nbp)
        ,payoffs(Nbp)
        ,total_payoff{0.0}
{
    assert(breeders.size() == payoffs.size());
} //end patch constructor


Patch::Patch(Patch const &other) :
    breeders{other.breeders}
    ,payoffs{other.payoffs}
    ,total_payoff{0.0}
{
    assert(breeders.size() == payoffs.size());
}

void Patch::operator=(Patch const &other)
{
    assert(other.breeders.size() == other.payoffs.size());
    breeders = other.breeders;
    payoffs= other.payoffs;
    total_payoff = other.total_payoff;
    
    assert(breeders.size() == payoffs.size());
} // end Patch::operator=(
