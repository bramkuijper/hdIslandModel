#include "individual.h"

Individual::Individual() :
    is_hawk{false}
    ,payoff{0.0}
    ,prob_hawk{0.0}
    ,disperse{{0.0,0.0}}
{
}

Individual::Individual(Individual const &other) :
    is_hawk{other.is_hawk}
    ,payoff{other.payoff}
    ,prob_hawk{other.prob_hawk}
    ,disperse{{other.disperse[0],other.disperse[1]}}
{
}

void Individual::operator=(Individual const &other)
{
    is_hawk = other.is_hawk;
    payoff = other.payoff;
    prob_hawk = other.prob_hawk;
    other.disperse[0] = other.disperse[0];
    other.disperse[1] = other.disperse[1];
}

