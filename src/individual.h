#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

class Individual
{
    public: 
        bool is_hawk;
        double payoff;
        double prob_hawk;
        double disperse[2];
        double id;

        Individual();

        Individual(Individual const &other);

        void operator=(Individual const &other);
};


#endif
