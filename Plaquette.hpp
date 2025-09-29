//
//

#ifndef Plaquette_hpp
#define Plaquette_hpp


#include <stdio.h>
#include <string>
#include <vector>
#include <set>

#include "Scaffold.hpp"
#include "Site.hpp"
#include "SystemType.hpp"
#include "head.hpp"


#include "SystemType.hpp"


class Plaquette{

    public:

    Plaquette(Site* center, int radius);
    float compute_1plaquette_energy(); 
    float compute_2plaquette_energy(Plaquette* otherPlaquette);

    bool intersectionEmpty(Plaquette* otherPlaquette);
    std::pair<bool,float> attemptSwapPlaquetteTo(Plaquette* destinationPlaquette, float Beta, std::vector <Site*> &full_to_empty_sites,  std::vector <Site*> &empty_to_full_sites);
    std::pair <bool, float> attemptFlipPlaquette(int rotation, float Beta, std::vector<int> corresponding_indices, std::vector<Site*> &full_to_empty_sites,  std::vector<Site*> &empty_to_full_sites, int N_Orientations);
    bool centerIsEmpty();
    void update_orientations();


    private:

    //Site* center; 

    std::vector<Site*> plaquette_sites ;
    std::set<Scaffold*> boundary_scaffolds;
    std::vector<int> orientations;
    int NPlaquetteSites ; 
    bool center_full; 

};
#endif 