#include <Rcpp.h>
#include "simulate_hd_space.h"

//' Run a single replicate individual-based simulation
//' of a hawk-dove
//' game in an island population.
//' @param Npatches Number of islands (integer)
//' @param Nbp Number of breeders per island
//'     (minimum 2 if you want interactions) (integer)
//' @param v Value of winning (floating-point value)
//' @param c Cost of losing a fight (floating-point value).
//' @param d Juvenile dispersal probability. If a juvenile disperses,
//'    it will move to a randomly chosen patch (floating-point value)
//' @param is_pure If \code{is_pure=TRUE}, Hawks and Doves are pure
//'    strategies, with mutation rates \code{mu} reflecting the probability
//'    that one pure strategy mutates into another. Alternatively, in case
//'    \code{is_pure=FALSE}, Hawks and Doves are the result of an evolving
//'    haploid locus that codes for probability \eqn{0 <= pHawk <= 1},
//'    which determines an individual's probability of playing Hawk.
//'    Here, the locus coding for \eqn{pHawk} gradually evolves.
//'    In case of a mutation event
//'    (ocurring at rate \code{mu}), the new value of the \eqn{pHwk}
//'    allele
//'    in individual \eqn{i} is given by \eqn{pHawk(i,t+1) = pHawk(i,t) +
//'    e}, where e is drawn from a normal distribution
//'    with mean 0 and standard deviation 0.01. \code{is_pure} is a boolean.
//' @param mu Mutation rate from Hawk to Dove, from Dove to Hawk
//'     (floating-point value between 0 and 1)
//' @param mu Mutation rate of dispersal of hawks and doves. If \code{mu_d=0}
//'     hawks and doves always disperse with probability \code{d}. If \code{mu_d>0}
//'     hawks and doves evolve their own dispersal rates
//'     (floating-point value between 0 and 1)
//' @param max_time Maximum number of generations after which
//'    the simulation ends -- in case of extinction simulation ends
//'    earlier (integer).
//' @param pHawk_init Initial sampling probability of Hawks.
//'    If \code{pHawk_init = 0.5}, the initial proportion of Hawks should be
//'    close to one half. It is not exactly one half, as the initial number
//'    of Hawks is not deterministic, but is sampled from a binomial
//'    distribution (floating-point value)
//' @param output_nth_generation Data is produced every
//'    \code{output_nth_generation} timesteps. Minimum number is 1 (output
//'    data every generation). This is done to save data in case you want to
//'    run the thing for
//'    say, values of \code{max_time = 50000}. Integer equal or larger than 1.
//' @return A \code{data.frame()} that contains all kinds of statistics
//'    collected during the running of the simulation.
//' @examples
//'    # test run of simulation (see parameters below), which is then
//'    # storted in the \code{data.frame} named \code{resulting.data}
//'    resulting.data <- runHDSpace(
//'                 Npatches=500, # 500 patches
//'                 Nbp=2, # 2 breeders per patch
//'                 v=1.0, # value of winning is 1
//'                 c=2.0, # cost of losing is 2
//'                 d=0.5, # 50% probability of dispersing as a juvenile
//'                 is_pure=TRUE, # pure strategy rather than mixed
//'                 mu=0.01, # mutation rate of a phenotype from H to D, D to H
//'                 mu_d=0, # dispersal does not evolve
//'                 max_time=10000, # maximum duration of the simulation
//'                 pHawk_init=0.5, # initial frequency of hawks
//'                 output_nth_generation=1 # output every generation
//'                 )
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame runHDSpace(
        int const Npatches=500
        ,int const Nbp=2
        ,double const v = 1.0
        ,double const c = 2.0
        ,double const d = 0.5
        ,bool const is_pure = true
        ,double const mu = 0.01
        ,double const mu_d = 0
        ,int max_time = 10000
        ,double pHawk_init = 0.5
        ,int output_nth_generation = 1
        ) {

    SimulateHDSpatial hds(
                Npatches
                ,Nbp
                ,v
                ,c
                ,pHawk_init
                ,d
                ,is_pure
                ,output_nth_generation
                ,mu
                ,mu_d
                ,max_time);

    Rcpp::DataFrame output = hds.run();

    return(output);
} // end load_simulation
