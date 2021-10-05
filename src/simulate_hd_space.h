#ifndef HD_SPACE_HPP
#define HD_SPACE_HPP

#include <Rcpp.h>
#include <random>
#include "patch.h"

class SimulateHDSpatial
{
    
    private:
        std::mt19937 rng_r; // random number generator
        std::uniform_real_distribution<double> uniform;
        std::vector <Patch> pop;

        // number of patches
        int Npatches;

        // breeders per patch
        int Nbp;

        double v;
        double c;

        double payoff_matrix[2][2];

        int generation;

        double init_pHawk;

        double d;

        // are strategies pure or not?
        bool is_pure;

        // per how many generations should we
        // output stats
        int output_nth_timestep;

        // mutation rate
        double mu;

        // max time the simulation runs
        int max_time;

        // stats function that returns by argument
        void stats(
            double &freq_Hawk
            ,double &sd_freq_Hawk
            ,double &mean_pHawk
            ,double &sd_pHawk);

        void interact_reproduce();

        double mutate(double const val);
        bool mutate(bool const val);

        void create_kid(Individual &parent, Individual &kid);

    public:
       SimulateHDSpatial(
            int const Npatches
            ,int const NbreedersPatch
            ,double const v
            ,double const c
            ,double const init_pHawk
            ,double const dispersal
            ,bool const is_pure
            ,int const output_nth_generation
            ,double const mu
            ,int const max_time);

        Rcpp::DataFrame run();

};

#endif
