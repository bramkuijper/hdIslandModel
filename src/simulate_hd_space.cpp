#include <Rcpp.h>
#include "patch.h"
#include "simulate_hd_space.h"

// assure number of breeders/patch is even
int _make_even(int const n)
{
    if (n < 1)
    {
        return(2);
    }

    if (n % 2 != 0)
    {
        return(n+1);
    }

    return(n);
}


SimulateHDSpatial::SimulateHDSpatial(
            int const Npatches
            ,int const NbreedersPatch
            ,double const v
            ,double const c
            ,double const init_pHawk
            ,double const dispersal
            ,bool const is_pure = true
            ,int const output_nth_timestep = 1
            ,double const mu = 0.01
            ,double const mu_d = 0.01
            ,int const max_time = 500) :
    rng_r((std::random_device())()) // init random number generator
    ,uniform(0.0, 1.0)
    ,pop(Npatches, Patch(_make_even(NbreedersPatch)))
    ,Npatches{Npatches}
    ,Nbp{_make_even(NbreedersPatch)}
    ,v{v}
    ,c{c}
    ,payoff_matrix{{v/2,0},{v,(v-c)/2}}
    ,generation{0}
    ,init_pHawk{init_pHawk}
    ,d{dispersal}
    ,is_pure{is_pure}
    ,output_nth_timestep{output_nth_timestep}
    ,mu{mu}
    ,mu_d{mu_d}
    ,max_time{max_time}
{
    // loop through all the patches
    for (size_t patch_idx = 0; patch_idx < Npatches; ++patch_idx)
    {
        assert(pop[patch_idx].breeders.size() == Nbp);
        assert(pop[patch_idx].payoffs.size() == Nbp);

        for (size_t ind_idx = 0; ind_idx < Nbp; ++ind_idx)
        {
            // initialize this population
            pop[patch_idx].breeders[ind_idx].is_hawk = 
                uniform(rng_r) < init_pHawk;
            
            pop[patch_idx].breeders[ind_idx].payoff = 0.0;

            pop[patch_idx].breeders[ind_idx].prob_hawk = init_pHawk;

            pop[patch_idx].breeders[ind_idx].disperse[0] = d;
            pop[patch_idx].breeders[ind_idx].disperse[1] = d;
        }
    }
} // SimulateHDSpatial::SimulateHDSpatial()



// interact within patches and reproduce
void SimulateHDSpatial::interact_reproduce()
{
    bool ind1_is_hawk,ind2_is_hawk;

    double baseline_fitness = c;

    double sum_payoffs = 0.0;
    double sum_d = 0.0;

    double local_sum_payoffs = 0.0;
    double local_sum_d = 0.0;

    double payoff1, payoff2;


    // make a distribution of payoffs from which we can sample
    // later on
    std::vector<double> payoffs_immigrant_parents(pop.size() * Nbp);
    
    // counter to keep track of payoff distribution()
    int payoff_idx = 0;

    for (size_t patch_idx = 0; 
            patch_idx < pop.size(); ++patch_idx)
    {
        // nonzero
        assert(pop[patch_idx].breeders.size() > 0);

        // even number?
        assert(pop[patch_idx].breeders.size() % 2 == 0);

        // shuffle stack of breeders so that 
        // we can have interactions
        // within random pairs
        if (pop[patch_idx].breeders.size() > 2)
        {

            // shuffle stack of breeders so that we can have interactions
            // within random pairs
            std::shuffle(pop[patch_idx].breeders.begin()
                    ,pop[patch_idx].breeders.end()
                    ,rng_r);
        }

        local_sum_payoffs = 0.0;
        local_sum_d = 0.0;

        assert(pop[patch_idx].breeders.size() == Nbp);
        assert(pop[patch_idx].payoffs.size() == Nbp);

        // loop through pairs
        for (size_t ind_idx = 0; 
                ind_idx < pop[patch_idx].breeders.size(); ind_idx+=2)
        {
            assert(ind_idx + 1 < pop[patch_idx].breeders.size());

            ind1_is_hawk = pop[patch_idx].breeders[ind_idx].is_hawk;
            ind2_is_hawk = pop[patch_idx].breeders[ind_idx+1].is_hawk;

            payoff1 = payoff2 = baseline_fitness;

            // if we have HH we need to calculate
            // payoff stochastically
            if (ind1_is_hawk && ind2_is_hawk)
            {
                if (uniform(rng_r) < 0.5)
                {
                    payoff1 += v;
                    payoff2 += -c;
                }
                else
                {
                    payoff1 += -c;
                    payoff2 += v;
                }
            }
            else
            {
                // payoff ind 1
                payoff1 += payoff_matrix[ind1_is_hawk][ind2_is_hawk];
                // payoff ind 2
                payoff2 += payoff_matrix[ind2_is_hawk][ind1_is_hawk];
            }

            sum_payoffs += d1 * payoff1 + d2 * payoff2;

            local_sum_payoffs += payoff1 + payoff2;

            assert(payoff_idx + 2 <= payoffs_immigrant_parents.size());

            payoffs_immigrant_parents[payoff_idx++] = payoff1;
            payoffs_immigrant_parents[payoff_idx++] = payoff2;
            

            pop[patch_idx].payoffs[ind_idx] = payoff1;
            pop[patch_idx].payoffs[ind_idx + 1] = payoff2;

            pop[patch_idx].breeders[ind_idx].payoff = payoff1;
            pop[patch_idx].breeders[ind_idx + 1].payoff = payoff2;

        } // size_t ind_idx = 0, ind_idx += 2

        pop[patch_idx].total_payoff = local_sum_payoffs;

    } // end for (size_t patch_idx = 0; 
            // patch_idx < pop.size(); ++patch_idx)


    // make a vector of individuals to use for recruiting purposes
    std::vector <Individual> recruits(Nbp);

    // make a probability distribution of global payoffs
    std::discrete_distribution<int> immigrant_payoffs(
        payoffs_immigrant_parents.begin()
        ,payoffs_immigrant_parents.end());

    // aux variable reflecting prob 
    // to sample local parent
    double prob_local;

    // total prob of immigration per patch
    double immigrants_per_patch = d * sum_payoffs / Npatches;
    double locals_per_patch;

    // now recruit
    for (int patch_idx = 0; patch_idx < pop.size(); ++patch_idx)
    {
        assert(pop[patch_idx].total_payoff > 0.0);

        assert(pop[patch_idx].total_payoff <= 
                Nbp * (baseline_fitness + v));

        // sample local parent
        locals_per_patch = pop[patch_idx].total_payoff * (1.0 - d);

        // competition function
        prob_local = locals_per_patch / 
            (locals_per_patch + immigrants_per_patch);

        // make distribution of local payoffs
        std::discrete_distribution<int> local_payoffs(
            pop[patch_idx].payoffs.begin()
            ,pop[patch_idx].payoffs.end());

        int sample, patch_sample_idx, parent_sample_idx;

        for (int ind_idx = 0; ind_idx < Nbp; ++ind_idx)
        {
            if (uniform(rng_r) < prob_local)
            {
                // sample local
                parent_sample_idx = local_payoffs(rng_r);
                patch_sample_idx = patch_idx;

            } // end prob_local
            else // sample a global individual
            {
                sample = immigrant_payoffs(rng_r);

                // get parent index
                parent_sample_idx = sample % Nbp;

                // get patch idx
                patch_sample_idx = int(floor((double) sample / Nbp)) % Npatches;
            }

            assert(parent_sample_idx >= 0);
            assert(parent_sample_idx < Nbp);
            assert(parent_sample_idx < pop[patch_sample_idx].breeders.size());
            assert(ind_idx < Nbp);
            assert(ind_idx >= 0);

            create_kid(
                    pop[patch_sample_idx].breeders[parent_sample_idx]
                    ,recruits[ind_idx]
                    );
        } //end for int ind_idx;

        assert(pop[patch_idx].breeders.size() == recruits.size());

        pop[patch_idx].breeders = recruits;
    } // end for (int patch_idx = 0; 
        //patch_idx < pop.size(); ++patch_idx)
} // end void SimulateHDSpatial::interact_reproduce()

void SimulateHDSpatial::create_kid(Individual &parent, Individual &kid)
{
    kid.payoff = 0.0;
    kid.prob_hawk = is_pure ? 
        0.0 
        : 
        mutate_prob(parent.prob_hawk, mu);


    kid.is_hawk = is_pure ? 
        mutate(parent.is_hawk)
        :
        uniform(rng_r) < kid.prob_hawk;

    kid.disperse[0] = mutate_prob(parent.disperse[0], mu_d);
    kid.disperse[1] = mutate_prob(parent.disperse[1], mu_d);

} // end create_kid()

double SimulateHDSpatial::mutate_prob(double const orig_val, double const mu)
{
    double val = orig_val;

    if (uniform(rng_r) < mu)
    {
        std::normal_distribution<double> mu_val(0.0, 0.1);

        val += mu_val(rng_r);

        if (val < 0.0)
        {
            val = 0.0;
        }
        else if (val > 1.0)
        {
            val = 1.0;
        }
    }

    return(val);
}

bool SimulateHDSpatial::mutate(bool const val)
{
    bool tval = val;

    if (uniform(rng_r) < mu)
    {
        tval = !tval;
    }

    return(tval);
}

Rcpp::DataFrame SimulateHDSpatial::run()
{
    // get a whole bunch of zero-filled vectors for the stats
    Rcpp::NumericVector freqHawk(floor((double)max_time / output_nth_timestep));
    Rcpp::NumericVector sd_freqHawk(floor((double)max_time / output_nth_timestep));
    Rcpp::NumericVector mean_pHawk(floor((double)max_time / output_nth_timestep));
    Rcpp::NumericVector sd_pHawk(floor((double)max_time / output_nth_timestep));
    Rcpp::NumericVector timestep_vec(floor((double)max_time / output_nth_timestep));

    double fHawk = 0.0;
    double sd_fHawk = 0.0; 
    double pHawk_x = 0.0;
    double sd_pHawk_x = 0.0;

    int stats_ctr = 0;

    for (generation = 0; generation < max_time; ++generation)
    {
        if (generation % output_nth_timestep == 0)
        {
            stats(fHawk, sd_fHawk, pHawk_x, sd_pHawk_x);

            assert(stats_ctr < floor((double)max_time / output_nth_timestep));
            freqHawk[stats_ctr] = fHawk;
            sd_freqHawk[stats_ctr] = sd_fHawk;
            mean_pHawk[stats_ctr] = pHawk_x;
            sd_pHawk[stats_ctr] = sd_pHawk_x;
            timestep_vec[stats_ctr] = generation;

            ++stats_ctr;
        }
        
        // check for user interruption
        Rcpp::checkUserInterrupt();

        interact_reproduce();

    }

    Rcpp::DataFrame simulation_data = Rcpp::DataFrame::create(
            Rcpp::Named("generation") = timestep_vec
            ,Rcpp::Named("freq_Hawk") = freqHawk
            ,Rcpp::Named("sd_freqHawk") = sd_freqHawk
            ,Rcpp::Named("mean_pHawk") = mean_pHawk 
            ,Rcpp::Named("sd_pHawk") = sd_pHawk);
            
    return(simulation_data);
}// end SimulateHDSpatial::run()

// write stats
void SimulateHDSpatial::stats(
        double &freq_Hawk
        ,double &sd_freq_Hawk
        ,double &mean_pHawk
        ,double &sd_pHawk)
{
    freq_Hawk = 0.0;
    sd_freq_Hawk = 0.0;
    mean_pHawk = 0.0;
    sd_pHawk = 0.0;

    double ss_freq_Hawk = 0.0;
    double ss_pHawk = 0.0;

    bool h;
    double p;

    for (size_t patch_idx = 0; patch_idx < pop.size(); ++patch_idx)
    {
        for (size_t ind_idx = 0;
                ind_idx < pop[patch_idx].breeders.size(); ++ind_idx)
        {
            h = pop[patch_idx].breeders[ind_idx].is_hawk;

            freq_Hawk += h;
            ss_freq_Hawk += h * h;

            p = pop[patch_idx].breeders[ind_idx].prob_hawk;

            assert(p >= 0);
            assert(p <= 1.0);

            mean_pHawk += p;
            ss_pHawk += p * p;
        }
    }

    int totalN =  pop.size() * Nbp;

    freq_Hawk /= totalN;
    mean_pHawk /= totalN;
    ss_pHawk /= totalN;
    ss_freq_Hawk /= totalN;

    sd_freq_Hawk = ss_freq_Hawk - freq_Hawk * freq_Hawk;
    sd_freq_Hawk = sd_freq_Hawk < 0.0 ? 0.0 : sqrt(sd_freq_Hawk);
    
    sd_pHawk = ss_pHawk - mean_pHawk * mean_pHawk;

    sd_pHawk = sd_pHawk < 0.0 ? 0.0 : sqrt(sd_pHawk);
} //end stats()
