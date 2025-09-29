//

#include "Plaquette.hpp"
#include "Site.hpp"
#include "Scaffold.hpp"
#include "head.hpp"
#include <ctime> // For random
#include <cstdlib> // For random
#include <iostream>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <set>

//#include "SystemType.hpp"
using namespace std ;

class Site;
class Scaffold;


Plaquette::Plaquette(Site* center, int radius):  plaquette_sites(vector<Site *> ()), 
boundary_scaffolds(set<Scaffold*> ()), orientations(vector<int>(0)), NPlaquetteSites(0), center_full(true)
{
    // Initialize the sites 
    set<Site*> set_plaquette_sites;
    (center)->add_orientation_recursive_flat(0, radius, set_plaquette_sites, plaquette_sites);

    NPlaquetteSites = plaquette_sites.size();
    // Get all the boundary scaffolds

    for (int i=0 ; i<NPlaquetteSites; i++){
        plaquette_sites[i]->scan_through_scaffolds(boundary_scaffolds);//,inner_scaffolds);
    }

    // Get the orientations
    for (int i=0 ; i<NPlaquetteSites; i++){
            orientations.push_back((plaquette_sites[i])->getOrientationRef());
    }

    if (center->isEmpty()){
        center_full = false;
    }

}



bool Plaquette::centerIsEmpty(){
    return (plaquette_sites[0]->isEmpty());
}


float Plaquette::compute_1plaquette_energy(){

    float energy = 0;

    for (auto it = boundary_scaffolds.begin(); it != boundary_scaffolds.end(); ++it) {
        energy+=(*it)->get_energy(); 
    }
    return energy;   
}

void Plaquette::update_orientations(){
    for (int i=0 ; i<NPlaquetteSites; i++){
            orientations[i]= (plaquette_sites[i])->getOrientationRef();
    }
}


float Plaquette::compute_2plaquette_energy(Plaquette* otherPlaquette){
    float energy = 0;

    set<Scaffold*> boundary_scaffolds_both ;
    for (auto it0 = (otherPlaquette->boundary_scaffolds).begin(); it0 != (otherPlaquette->boundary_scaffolds.end()); ++it0) {
        boundary_scaffolds_both.insert(*it0);
    }
    

    // Take the union of the two sets, with the intersection only one
    for (auto it = boundary_scaffolds.begin(); it != boundary_scaffolds.end(); ++it) {
        if ((otherPlaquette->boundary_scaffolds).find(*it) == (otherPlaquette->boundary_scaffolds).end()) {
            // Element not found in set1, add it
            boundary_scaffolds_both.insert(*it);}
    }

    //cout<<boundary_scaffolds.size()<<" boundary scaffolds"<<endl;
    for (auto it2 = boundary_scaffolds_both.begin(); it2 != boundary_scaffolds_both.end(); ++it2) {
        energy+=(*it2)->get_energy(); 
        /* if ((*it)->get_energy()!=0){
            cout<<(*it)->get_energy()<<" ";} */
    }
    return energy;   
}


bool Plaquette::intersectionEmpty(Plaquette* otherPlaquette){

    set<Site*> set1 ;
    set<Site*> set2 ;

    for (int i=0 ; i<NPlaquetteSites; i++){
        set1.insert((plaquette_sites[i]));
        set2.insert((otherPlaquette->plaquette_sites)[i]);
    }

    bool intersectionEmpty = true;
        for (Site* elem : set1) {
            if (set2.count(elem) > 0) {
                intersectionEmpty = false;
                break;
            } 
        }
    return intersectionEmpty;
}



pair<bool,float> Plaquette::attemptSwapPlaquetteTo(Plaquette* destinationPlaquette, float Beta, vector <Site*> &full_to_empty_sites,  vector <Site*> &empty_to_full_sites){
    
    //cout<< "swap plaquette "<<(this) <<" and "<<destinationPlaquette<< " of sizes "<<NPlaquetteSites<<" "<<destinationPlaquette->NPlaquetteSites<<endl;
    pair<bool, float> successAndDeltaE = {true, 0} ;//the boolean represents the success of swap, and the int is delta_energy
    
    this->update_orientations();
    destinationPlaquette->update_orientations();

    //if (true){//(m_orientationRef>0 && destinationSite.isEmpty()){
        if (this->intersectionEmpty(destinationPlaquette)){
            float old_energy = this->compute_2plaquette_energy(destinationPlaquette);
         
            for (int i=0 ; i<NPlaquetteSites; i++){
                (plaquette_sites[i])->setFullSite((destinationPlaquette->orientations)[i]);
                ((destinationPlaquette->plaquette_sites)[i])->setFullSite(orientations[i]);
            }

            float new_energy = this->compute_2plaquette_energy(destinationPlaquette);
           
            successAndDeltaE.second=(new_energy-old_energy);
            if (successAndDeltaE.second>0){
                
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

                if (exp (- static_cast <float> (successAndDeltaE.second*Beta)) < r){
                    // Do the reverse swap 
                    for (int i=0 ; i<NPlaquetteSites; i++){
                
                        (plaquette_sites[i])->setFullSite(orientations[i]);
                        ((destinationPlaquette->plaquette_sites)[i])->setFullSite((destinationPlaquette->orientations)[i]);
                    }   

                    successAndDeltaE.first=false;
                    successAndDeltaE.second=0;
                    
                    //cout<<"move refused "<<endl;
                }
                
            }
            if (successAndDeltaE.first){
                //cout<<"plaquette swap accepted "<<endl;

                // We find the sites that went from empty to full and vice versa
                for (int i=0 ; i<NPlaquetteSites; i++){
                    if (orientations[i]==0 && (destinationPlaquette->orientations)[i]!=0){
                        full_to_empty_sites.push_back((destinationPlaquette->plaquette_sites)[i]);
                        empty_to_full_sites.push_back(plaquette_sites[i]);
                    }
                    if ((destinationPlaquette->orientations)[i]==0 && orientations[i]!=0){
                        empty_to_full_sites.push_back((destinationPlaquette->plaquette_sites)[i]);
                        full_to_empty_sites.push_back(plaquette_sites[i]);
                    }
                
                    // We update the stored orientations
                    //orientations[i] = ((plaquette_sites[i])->getOrientationRef());
                }
                if (full_to_empty_sites.size()!=empty_to_full_sites.size()){
                cout<<"problem counting the full ("<<full_to_empty_sites.size()<<") and empty("<<empty_to_full_sites.size()<<")"<<endl;
                    /* for (auto it_origin = origin_plaquette_sites.begin(); it_origin != origin_plaquette_sites.end(); ++it_origin) {
                        cout<<"in origin "<< *it_origin<<" "<<endl;
                    }
                    for (auto it_dest = origin_plaquette_sites.begin(); it_dest != origin_plaquette_sites.end(); ++it_dest) {
                        cout<<"in dest "<< *it_dest<<" "<<endl;
    
                    } */
                }

            
            }
            //else{cout<<"plaquette swap refused"<<endl;}
        }
        else{
            // the two plaquettes overlap 
            successAndDeltaE.first=false;
        }
    
    //else{
    //    cout<< "The initial site was empty or the destination site was not. It can not be swapped" <<endl;
    //}
    return successAndDeltaE;
}






pair <bool, float> Plaquette::attemptFlipPlaquette(int rotation, float Beta, vector<int> corresponding_indices, vector<Site*> &full_to_empty_sites,  vector<Site*> &empty_to_full_sites, int N_Orientations){
        pair<bool, float> successAndDeltaE = {true, 0} ;

        this->update_orientations();
        
        float old_energy = this->compute_1plaquette_energy();
        cout<<"plaquette flip with rotation "<<rotation<<", Norientations "<< N_Orientations<<" on plaquette of size "<<NPlaquetteSites<<endl;

        for (int i=0 ; i<NPlaquetteSites; i++){
            
            int old_orientation_i = orientations[corresponding_indices[i]];
            //int current_orientation_i = orientations[i];
            //cout<<"site "<<i<<" had orientation "<<current_orientation_i<<" and now we change it"<<endl;
            /*cout<<"position "<<plaquette_sites_list[i]->get_real_space_position(10)[0]<<" "<<plaquette_sites_list[i]->get_real_space_position(10)[1]<<" ";
            cout<<"old position "<<plaquette_sites_list[corresponding_indices[i]]->get_real_space_position(10)[0]<<" "<<plaquette_sites_list[corresponding_indices[i]]->get_real_space_position(10)[1]<<endl;
            cout<<"index "<<i<<" is from "<<corresponding_indices[i]<<endl;
            cout<<endl; */
            

            if (old_orientation_i==0){
                (plaquette_sites[i])->setFullSite(old_orientation_i);
                //cout<<"0->?"<<endl;
            }
            else{
                //int new_orientation_i = (old_orientation_i-1)/6 * 6 + (old_orientation_i + rotation-1)%6 +1 ;// CAREFUL , this line is specific to the triangular lattice
                int new_orientation_i = (old_orientation_i-1)/N_Orientations * N_Orientations + (old_orientation_i + rotation-1)%N_Orientations +1 ;// CAREFUL , this line is specific to the triangular lattice
                //cout<<"old_orientation_i "<<old_orientation_i<<" -> new_orientation_i "<<new_orientation_i<<endl;
                //cout<<"new "<<new_orientation_i<<" old "<<old_orientation_i<<" rot "<<rotation<<endl;
                (plaquette_sites[i])->setFullSite(new_orientation_i); 
            }
            
        }

        //cout<<"compute new energy"<<endl;
        float new_energy = this->compute_1plaquette_energy();
        successAndDeltaE.second=(new_energy-old_energy);
        if (successAndDeltaE.second>0){
            
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

            if (exp (- static_cast <float> (successAndDeltaE.second*Beta)) < r){
                // we do the reverse swap 
                for (int i=0 ; i<NPlaquetteSites; i++){
                    (plaquette_sites[i])->setFullSite(orientations[i]);
                  
                }
                successAndDeltaE.first = false;
                //cout<<"plaquette flip refused "<<endl;
            }

        }
        if (successAndDeltaE.first){
            //cout<<"plaquette flip accepted "<<endl;
            //auto site_iterator3 = plaquette_sites.begin();
            for (int i=0 ; i<NPlaquetteSites; i++){
                if (orientations[i]==0 and plaquette_sites[i]->getOrientationRef()!=0){
                    empty_to_full_sites.push_back(plaquette_sites[i]);
                }
                if (orientations[i]!=0 and plaquette_sites[i]->getOrientationRef()==0){
                    full_to_empty_sites.push_back(plaquette_sites[i]);
                }
                // We update the stored orientations
                //orientations[i] = ((plaquette_sites[i])->getOrientationRef());
            }
            //cout<<"stored orientation updated"<<endl;

        if (full_to_empty_sites.size()!=empty_to_full_sites.size()){
                cout<<"problem counting the full ("<<full_to_empty_sites.size()<<") and empty("<<empty_to_full_sites.size()<<")"<<endl;
            } 

        }

    return successAndDeltaE;
}