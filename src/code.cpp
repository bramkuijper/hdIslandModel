#include <Rcpp.h>

int Npatches = 500;
int Nbp = 10;

struct Individual {
  bool is_hawk;
  double phawk;
  double payoff;
};

struct Patch {
  Individual breeders[Nbp];
};

Patch Pop[Npatches];

init_population()
{
  for (int patch_idx = 0; patch_idx < Npatches; ++patch_idx)
  {
    for (int ind_idx = 0; ind_idx < Nbp; ++ind_idx)
    {
      Pop[patch_idx].breeders[ind_idx].is_hawk = init_hawk

    } // end for int ind_idx
  } // end for int patch_idx
}
