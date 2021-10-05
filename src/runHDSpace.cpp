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
//' @param d Juvenile dispersal probability. If juvenile disperses
//'    (floating-point value), it will move to a random patch
//' @param is_pure If \code{is_pure=TRUE}, Hawks and Doves are pure
//'    strategies, with mutation rates \code{mu} reflecting the probability
//'    that one pure strategy mutates into another. Alternatively, in case
//'    \code{is_pure=FALSE}, Hawks and Doves are the result of an evolving
//'    haploid locus that codes for probability \eqn{0 \leq pHawk \leq 1},
//'    which determines an individual's probability of playing Hawk.
//'    Here, \eqn{pHawk} gradually evolves. In case of a mutation event
//'    (ocurring at rate \code{mu}), the new value of the \eqn{pHwk}
//'    allele
//'    in individual \eqn{i} is given by \eqn{pHawk(i,t+1) = pHawk(i,t) +
//'    e}, where e is drawn from a normal distribution
//'    with mean 0 and standard deviation \code{sd_pHawkMixed} (see below)
//' @param mu Mutation rate from Hawk to Dove, from Dove to Hawk
//' @param max_time Maximum number of generations after which
//'    the simulation ends -- in case of extinction simulation ends
//'    earlier (integer).
//' @param pHawk_init Initial sampling probability of Hawks.
//'    If \code{pHawk_init = 0.5}, the initial proportion of Hawks should be
//'    close to one half. It is not exactly one half, as the initial number
//'    of Hawks is not deterministic, but is sampled from a binomial
//'    distribution
//' @param output_nth_generation Data is produced every
//'    \code{output_nth_generation} timesteps. Minimum number is 1 (output
//'    data every generation). This is done to save data in case you want to
//'    run the thing for
//'    say, values of \code{max_time = 50000}
//' @return A \code{data.frame()} that contains all kinds of statistics
//'    collected during the running of the simulation.
//' @examples
//' # quick test run of 10 generations and size of 10 individuals
//' # expecting ~100% hawks (but, yeah... sampling)
//' resulting.data <- runHDSpace(Npatches=500,
//'                 Nbp=2,
//'                 v=1.0,
//'                 c=2.0,
//'                 d=0.5,
//'                 is_pure=TRUE,
//'                 mu=0.01,
//'                 max_time=10000,
//'                 pHawk_init=0.5,
//'                 output_nth_generation=1)
//'
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
                ,max_time);

    Rcpp::DataFrame output = hds.run();

    return(output);
} // end load_simulation
