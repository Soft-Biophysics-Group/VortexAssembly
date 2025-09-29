//
//  Site.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 22/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

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
#include <random>
#include "LogDuration.hpp"

//#include "SystemType.hpp"
using namespace std ;
class Scaffold ;

/*Site::Site(int orientationRef, vector<vector <Site*>> neighbourList) : m_orientationRef (orientationRef), m_neighbourList (vector<vector <Site *>>(0)), hashedPosition(0), belongingScaffolds(vector <Scaffold*> (0)), indexInBelongingScaffold(vector <int>()), belongingExtendedScaffolds(vector <Scaffold*> (0)), indexInBelongingExtendedScaffold(vector <int>())
{
    vector <Site *> :: iterator i;
    for (i=neighbourList.begin(); i != neighbourList.end() ; ++i){
        (*this).setNeighbourRelation(**i);
    }
}*/

Site::Site() : m_orientationRef(0), m_neighbourList(vector<vector <Site *>> (0)), hashedPosition(0), belongingScaffolds(vector<vector <Scaffold*>> (0)),  indexInBelongingScaffold(vector<vector <int>>())
{
    for (int i=0;i<maxDegreeOfNeighbourhood;i++){
        m_neighbourList.push_back(vector <Site*>(0));
        belongingScaffolds.push_back(vector<Scaffold*>(0));
        indexInBelongingScaffold.push_back(vector <int>());
    }
//    cout<< this << " has "<< m_neighbourList.size()<<" type of neighbors"<<endl;

}


/* vector <Site*> Site::getNeighbours (int degree){
    return m_neighbourList[degree];
}
 */

/*
Site*  Site::getAdress() const
{
    return this;
}
 */

void Site::setScaffoldAppartenance(Scaffold* s_adress, int indexInS, int degreeOfNeighbourhood){
    belongingScaffolds[degreeOfNeighbourhood].push_back(s_adress);
    indexInBelongingScaffold[degreeOfNeighbourhood].push_back(indexInS);
} 

  
                                  
                    
void Site::printAttributes(int const Nhash) const{
    cout << endl<< "Site description"<<endl;
    cout << "The adress of the site is : " <<this <<endl;
    cout << "The orienation is : " <<m_orientationRef <<endl;
    
    cout << "The hashed position is : "<<hashedPosition<<endl;
    vector <int> real_pos (this->get_real_space_position(Nhash));
    cout << "The real position is : "<<real_pos[0]<<" "<<real_pos[1]<<" "<<real_pos[2]<<endl;
    
    //Site::printVectorOfSite(m_neighbourList);
    cout << "Neighbours and their orientation : "<<endl;
    vector <Site *> ::const_iterator iter;
    for (iter = m_neighbourList[0].begin() ; iter!=m_neighbourList[0].end() ; iter++){
        cout << *iter << " -> "<< (**iter).m_orientationRef <<endl;
    }
    
    cout << "Second neighbours and their orientation : "<<endl;
    vector <Site *> ::const_iterator iter2;
    for (iter2 = m_neighbourList[1].begin() ; iter2!=m_neighbourList[0].end() ; iter2++){
        cout << *iter2 << " -> "<< (**iter2).m_orientationRef <<endl;
    }
    //printVectorOfAdress(m_neighbourList);
    
    cout << "Scaffolds and their structures : "<<endl;
    vector <Scaffold *> ::const_iterator it;
    for (it = belongingScaffolds[0].begin() ; it!=belongingScaffolds[0].end() ; it++){
        cout << *it << " -> "<< (*it)->getStructure() <<endl;
    }
    //printVectorOfAdress(belongingScaffolds);
    //Scaffold::printStructureRefs(belongingScaffolds);
    
    /*vector <Site *> :: const_iterator i ;
    for (i = m_neighbourList.begin(); i != m_neighbourList.end() ; ++i){
        cout << *i <<endl;}
    cout <<endl;*/
}

/*
void Site::printVectorOfSite(vector <Site *> listToPrint){
    //if listToPrint.site()!=0){
    vector <Site *> :: const_iterator i ;
    for (i = listToPrint.begin(); i != listToPrint.end() ; ++i){
        cout << *i <<endl;}
    cout <<endl;
        
}*/


bool Site::setNeighbourRelation(Site& neighbourSite, int degreeOfNeighbour){
    bool madeBond = false;
    bool madeReverseBond = false;
    if (!(this==&neighbourSite)){
    //vector <Site*> &listNeighbour (neighbourSite.m_neighbourList);
    if (! ( (*this).isNeighbourOf(neighbourSite,degreeOfNeighbour) ))
       { //if the current site is not already a neighbour of neighbourSite
            neighbourSite.m_neighbourList[degreeOfNeighbour].push_back(this); //add this to the list of neighbour adress for the neighbourSite
            madeBond = true;
    }
    if (! (neighbourSite.isNeighbourOf(*this, degreeOfNeighbour)) )
        { //if the neighbourSite is not already a neighbour of the current site
        m_neighbourList[degreeOfNeighbour].push_back(&neighbourSite); //add the neighbour adresss to the list of the current site neighbours adresses
            madeReverseBond=true;
    }
    }
    //else{cout << "found equality while setting neighbour relation"<<endl;}
    return(madeBond && madeReverseBond);
}


vector<vector<int>> Site::reorder_neighbors(int degreeOfNeighborhood, int NhashPosition, int Lx, int Ly, int Lz, vector<vector<int>> order_ref){
    
    vector <int> pos_ref = this->get_real_space_position(NhashPosition);
    vector <vector <int>> pos_neighbors;
    vector <Site *> ::const_iterator iter;
    for (iter = m_neighbourList[degreeOfNeighborhood].begin() ; iter!=m_neighbourList[degreeOfNeighborhood].end() ; iter++){
        pos_neighbors.push_back((*iter)->get_real_space_position(NhashPosition));
    }

    //cout<<"old order"<<endl;
    for (int i=0; i<pos_neighbors.size(); i++){
        pos_neighbors[i][0] = (pos_neighbors[i][0]-pos_ref[0]+Lx)%Lx;
        pos_neighbors[i][1] = (pos_neighbors[i][1]-pos_ref[1]+Ly)%Ly;
        pos_neighbors[i][2] = (pos_neighbors[i][2]-pos_ref[2]+Lz)%Lz;
        //cout<< pos_neighbors[i][0]<<" "<<pos_neighbors[i][1]<<" "<<pos_neighbors[i][2]<<endl;   
    }
    
    //cout<<endl;

    if (order_ref.size()!=0){
        vector <Site*> new_neighbor_list ;
        for (int j=0; j<order_ref.size(); j++){
            vector <int> delta_pos = order_ref[j];
            for (int k=0 ; k<pos_neighbors.size(); k++){
                if (pos_neighbors[k][0]==order_ref[j][0] && pos_neighbors[k][1]==order_ref[j][1] && pos_neighbors[k][2]==order_ref[j][2]){
                        new_neighbor_list.push_back(m_neighbourList[degreeOfNeighborhood][k]);
                }
            }
        }

        m_neighbourList[degreeOfNeighborhood] = new_neighbor_list;
    }
    

   /*  vector <vector <int>> pos_neighbors2;
    vector <Site *> ::const_iterator iter2;
    for (iter2 = m_neighbourList[degreeOfNeighborhood].begin() ; iter2!=m_neighbourList[degreeOfNeighborhood].end() ; iter2++){
        pos_neighbors2.push_back((*iter2)->get_real_space_position(NhashPosition));
    }
    cout<<"new order"<<endl;
    for (int l=0; l<pos_neighbors2.size(); l++){
        pos_neighbors2[l][0] = (pos_neighbors2[l][0]-pos_ref[0]+Lx)%Lx;
        pos_neighbors2[l][1] = (pos_neighbors2[l][1]-pos_ref[1]+Ly)%Ly;
        pos_neighbors2[l][2] = (pos_neighbors2[l][2]-pos_ref[2]+Lz)%Lz;
        cout<< pos_neighbors2[l][0]<<" "<<pos_neighbors2[l][1]<<" "<<pos_neighbors2[l][2]<<endl;   
    } */

    
    return(pos_neighbors);

}


void Site::set_position(int const i, int const  j,int  const  k, int const  N, int const N2){
    hashedPosition=hash_ints(i, j, k, N, N2);
    
}

vector<int> Site::get_real_space_position(int const N) const {
    return(unhash_ints(hashedPosition, N));
}

//void Site::fillWithRandomValue(int numberParticleType){
//    setFullSite( ( rand()%(numberParticleType) ) +1 );
//}
//
//void Site::fillWithNewRandomValue (int numberParticleType, int currentOrientation){
//    setFullSite(  1 + ( currentOrientation + rand()%(numberParticleType-1) )%numberParticleType );
//}

                                                                                                                       



void Site::setFullSite(int orientationRef){

    m_orientationRef = orientationRef ;

    for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){

        /* vector <Scaffold *> :: iterator it;
        vector <int>:: iterator indexIterator;
        indexIterator=indexInBelongingScaffold[degreeOfNeighbourhood].begin();
        for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
            (*it)->actualize_scaffold(m_orientationRef,*indexIterator);
            indexIterator++;
        }
 */
        const auto& scaffolds = belongingScaffolds[degreeOfNeighbourhood];
        const auto& indices = indexInBelongingScaffold[degreeOfNeighbourhood];
        for (size_t i = 0; i < scaffolds.size(); ++i) {
            scaffolds[i]->actualize_scaffold(m_orientationRef, indices[i]);
        }


}
}

void Site::changeScaffoldOrientation(){

     for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){

        vector <Scaffold *> :: iterator it;
//        vector <int>:: iterator indexIterator;
//        indexIterator=indexInBelongingScaffold[degreeOfNeighbourhood].begin();
        int indexIterator = 0;
        for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
            
//            (*it)->change_orientations(m_orientationRef,indexInBelongingScaffold[degreeOfNeighbourhood].at(*indexIterator));
            (*it)->change_orientations(m_orientationRef,indexInBelongingScaffold[degreeOfNeighbourhood][indexIterator]);
            indexIterator++;
        
        }
}
}

void Site::changeScaffoldStructure(){
//    if (m_orientationRef==0){
//          cout<<"New empty site ";
//      }
     for (int degreeOfNeighbourhood=0 ; degreeOfNeighbourhood<maxDegreeOfNeighbourhood ; degreeOfNeighbourhood++){

        vector <Scaffold *> :: iterator it;
        for (it = belongingScaffolds[degreeOfNeighbourhood].begin(); it!= belongingScaffolds[degreeOfNeighbourhood].end(); it++){
//            if (m_orientationRef==0){
//             cout << (*it)->getStructure()<<"->";
//            }
            
            (*it)->actualize_structure();

//            if (m_orientationRef==0){
//             cout << (*it)->getStructure()<<" ";
//            }
        }
}
}





int Site::getOrientationRef() const {
    return m_orientationRef;
}

int Site::getSiteStateTemporary() {
    return m_orientationRef;
}

float Site::getSiteEnergy(){
    float energy =0;
    vector <Scaffold *> :: iterator it;
    
//    cout << "Site energy ";
    for (it = belongingScaffolds[0].begin(); it!= belongingScaffolds[0].end(); it++){
        energy+=(*it)->get_energy();
//        cout<<energy<<", ";
    
    }
//    cout<<endl;
    //cout<<"Etot "<<energy<<endl;
    return energy;
}

// Now defined as inline in the hpp file 
/* void Site::setOrientation(int orientationValue){
    m_orientationRef=orientationValue;
} */

/* void Site::putParticleIn(Site &destinationSite, int new_state){
    // the destination site has the new_state
    // the original site has the state of the destination site (usually, empty)

    
    int orientationCopy = destinationSite.m_orientationRef;
    //Change both orientation at the same time
    //destinationSite.setOrientation(m_orientationRef);
    destinationSite.setOrientation(new_state);
    
    m_orientationRef = orientationCopy ;
    
    destinationSite.changeScaffoldOrientation();
    this->changeScaffoldOrientation();
    
    destinationSite.changeScaffoldStructure();
    this->changeScaffoldStructure();

} */
     

pair<bool,float> Site::attemptSwapTo(Site &destinationSite, float Beta, int new_state){
   
    
    pair<bool, float> successAndDeltaE = {true, 0} ;//the boolean represents the success of swap, and the int is delta_energy
    //if (m_orientationRef>0 && destinationSite.isEmpty()){
        
    int old_state = m_orientationRef ;
        
    float old_energy = getSiteEnergy() + destinationSite.getSiteEnergy();
    
    this->putParticleIn(destinationSite, new_state);
    float new_energy = getSiteEnergy() + destinationSite.getSiteEnergy();
    successAndDeltaE.second=new_energy-old_energy;
        if (successAndDeltaE.second>0){
            
            // Old way to draw r, slow 
            // float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

            // New way, faster 
            static thread_local std::mt19937 rng(std::random_device{}());
            static thread_local std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            float r = dist(rng);

            if (exp (- static_cast <float> (successAndDeltaE.second*Beta)) < r){
                destinationSite.putParticleIn(*this, old_state);
                successAndDeltaE.first=false;
                successAndDeltaE.second=0;
            }
        }
        
    return successAndDeltaE;
}

pair <bool, float> Site::attemptFlip(int new_state, float Beta){

    pair <bool, float> successAndDeltaE {true, 0} ;
    //if (!(this->isEmpty())){
        
    float old_orientationRef = m_orientationRef;
    float old_energy = this->getSiteEnergy();
    

//        cout << "Old energy flip: "<<old_energy;
    
//        this->fillWithNewRandomValue(numberParticleType, m_orientationRef);
    this->setFullSite(new_state);
    float new_energy = this->getSiteEnergy();
//        cout << " New energy flip: "<<new_energy<<endl;
    successAndDeltaE.second = new_energy - old_energy;
    if (successAndDeltaE.second>0){
        
        // Old way to draw r, slow 
        // float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        // New way, faster 
        static thread_local std::mt19937 rng(std::random_device{}());
        static thread_local std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        float r = dist(rng);


        if (exp (- static_cast <float> (successAndDeltaE.second*Beta)) < r){
            this->setFullSite(old_orientationRef);
            successAndDeltaE.second=0;
            successAndDeltaE.first =false;
        }
    }
 
        return successAndDeltaE;
}






void Site::get_plaquette_relative_positions(int level, int level_max, int NhashPosition, set<Site*> &plaquette_sites, vector<vector<int>> &relative_positions, vector <int> center_position, int Lx, int Ly, int Lz){
    // This function explores a plaquette and save the relative position compared to the center site for all the sites, with 
    // the same exploration algorithm that the other plaquette function, such that this list is the same for all the sites
    // Consequently, this function only needs to be used once at the beginning of the Monte-Carlo algorithm

    

    auto it = plaquette_sites.find(this);
    bool isNewPlaquetteSite =  (it == plaquette_sites.end());
    vector <int> sizes = {Lx, Ly, Lz};

    //cout << this << " new"<<isNewPlaquetteSite<<" level="<<level<<" plaquette size"<<plaquette_sites.size()<<endl;
    if (isNewPlaquetteSite){
         plaquette_sites.insert(this);
         vector <int> position = this->get_real_space_position(NhashPosition);
         for (int d=0 ; d<3 ; d++){
            position[d] = (position[d]-center_position[d]+sizes[d])%sizes[d];
         }
         relative_positions.push_back(position);
    }

    if (level!=level_max ){
        vector <Site *> ::const_iterator iter;
        for (iter = m_neighbourList[0].begin() ; iter!=m_neighbourList[0].end() ; iter++){
            (*iter)->get_plaquette_relative_positions(level+1, level_max, NhashPosition, plaquette_sites,relative_positions, center_position, Lx, Ly, Lz);
        } 
    }
}






void Site::add_orientation_recursive(int level, int level_max, vector <vector <int>> &old_orientations, set <Site*> &plaquette_sites){
    // This function builds the list of site that are within a radius level_max around the initial site, and save the orientations

    auto it = plaquette_sites.find(this);
    bool isNewPlaquetteSite =  (it == plaquette_sites.end());

    //cout << this << " new"<<isNewPlaquetteSite<<" level="<<level<<" plaquette size"<<plaquette_sites.size()<<endl;
    if (isNewPlaquetteSite){
        if (old_orientations.size()<level+1){
        old_orientations.push_back(vector<int>());
        }
        if (level_max>3){
            cout <<m_orientationRef<<" ";
        }
        old_orientations[level].push_back(m_orientationRef);
        //cout << "m" << m_orientationRef<<" ";
         plaquette_sites.insert(this);
    }
    if (level!=level_max ){
       
        vector <Site *> ::const_iterator iter;
        for (iter = m_neighbourList[0].begin() ; iter!=m_neighbourList[0].end() ; iter++){
            (*iter)->add_orientation_recursive(level+1, level_max, old_orientations, plaquette_sites);
        } 
    }
}

void Site::get_plaquette_relative_positions2( int level_max, int NhashPosition, vector<vector<int>> &relative_positions, vector <int> center_position, int Lx, int Ly, int Lz){
    std::vector <int> old_orientations;
    set<Site*> plaquette_sites ;
    vector<Site*> plaquette_sites_list;
    //this->add_orientation_recursive_flat(0, level_max, old_orientations, plaquette_sites, plaquette_sites_list);    
    this->add_orientation_recursive_flat(0, level_max, plaquette_sites, plaquette_sites_list);      
    vector <int> sizes = {Lx, Ly, Lz};

    // save the positions of the plaquette
    //auto site_iterator0 = plaquette_sites.begin();

    for (int i=0 ; i<plaquette_sites.size(); i++){

        vector <int> position = (plaquette_sites_list[i])->get_real_space_position(NhashPosition);
        //cout<<"site position "<<(*site_iterator0)->get_real_space_position(NhashPosition)[0]<<" "<<(*site_iterator0)->get_real_space_position(NhashPosition)[1]<<endl;
            
        for (int d=0 ; d<3 ; d++){
            position[d] = (position[d]-center_position[d]+sizes[d])%sizes[d];
        }
        relative_positions.push_back(position);
        //site_iterator0++;
        }
}



//void Site::add_orientation_recursive_flat(int level, int level_max, vector <int> &old_orientations, set <Site*> &plaquette_sites, vector <Site*> &plaquette_sites_list){
void Site::add_orientation_recursive_flat(int level, int level_max, set <Site*> &plaquette_sites, vector <Site*> &plaquette_sites_list){
 
     // This function builds the list of site that are within a radius level_max around the initial site, and save the orientations
    // The old orientations list is a one dimensional vector only 

    auto it = plaquette_sites.find(this);
    bool isNewPlaquetteSite =  (it == plaquette_sites.end());

    //cout << this << " new"<<isNewPlaquetteSite<<" level="<<level<<" plaquette size"<<plaquette_sites.size()<<endl;
    if (isNewPlaquetteSite){
        //old_orientations.push_back(m_orientationRef);
        plaquette_sites.insert(this);
        plaquette_sites_list.push_back(this);
        //cout<<"site in recursive position "<<(this)->get_real_space_position(10)[0]<<" "<<(this)->get_real_space_position(10)[1]<<endl;
            
    }
    if (level!=level_max ){
        vector <Site *> ::const_iterator iter;
        for (iter = m_neighbourList[0].begin() ; iter!=m_neighbourList[0].end() ; iter++){
            //(*iter)->add_orientation_recursive_flat(level+1, level_max, old_orientations, plaquette_sites, plaquette_sites_list);
            (*iter)->add_orientation_recursive_flat(level+1, level_max, plaquette_sites, plaquette_sites_list);
        } 
    }
}





void Site::set_new_orientations_plaquette(int level, int level_max, vector <int> &indices, vector<vector <int>> orientations, vector<Site*> &full_to_empty_sites, vector <Site*> &empty_to_full_sites, set <Site*> &explored_sites){

    //cout<< this << " old orientation "<<m_orientationRef <<" - new orientation "<<orientations[level][indices[level]] <<endl;

    auto it = explored_sites.find(this);
    bool isNewPlaquetteSite =  (it == explored_sites.end());
    

    if (isNewPlaquetteSite){

        if (m_orientationRef==0 && orientations[level][indices[level]]!=0){
        //cout << this << " empty to full"<<endl;
        empty_to_full_sites.push_back(this);
        }
        if (m_orientationRef!=0 && orientations[level][indices[level]]==0){
            //cout << this << " full to empty"<<endl;
            full_to_empty_sites.push_back(this);
        }
        if (level_max>3){
            cout<< orientations[level][indices[level]]<<" ";
        }
        this->setFullSite(orientations[level][indices[level]]);

        indices[level]+=1;
        explored_sites.insert(this);
    }

    if (level!=level_max ){        
        vector <Site *> ::const_iterator iter;
        for (iter = m_neighbourList[0].begin() ; iter!=m_neighbourList[0].end() ; iter++){
            (*iter)->set_new_orientations_plaquette(level+1, level_max, indices, orientations, full_to_empty_sites, empty_to_full_sites,explored_sites);
        } 
    }
}








void Site::scan_through_scaffolds(set< Scaffold *> &boundary_scaffolds){//, set <Scaffold *> &inner_scaffolds){
    vector <Scaffold *> :: iterator it;
    // We go through all the scaffolds, a scaffold belongs to two sites 
    // If a scaffold is found only once, it is a boundary scaffold and it is in the final set
    // If it is shared between two sites of the plaquette, it is added and then removed from the set 
    // and we won't count its energy in the end
    
    //cout<<"begin site "<<(this)<<" n boundaries "<<boundary_scaffolds.size()<<endl;
    for (it = belongingScaffolds[0].begin(); it!= belongingScaffolds[0].end(); it++){
            auto index = boundary_scaffolds.find(*it);
            
            if ((*it)->getEdge()>0){ // we only take into account the two particles scaffold
                if (index == boundary_scaffolds.end()) { // the scaffold is not in the set
                    boundary_scaffolds.insert(*it);
                    
                    
                } else {
                    boundary_scaffolds.erase(*it);
                    //inner_scaffolds.insert(*it);
                }
            }
    }

    //cout<<"end site "<<(this)<<" n boundaries "<<boundary_scaffolds.size()<<endl;
}


float Site::compute_2plaquette_energy(set<Site*> plaquette_sites1, set<Site*> plaquette_sites2){

    float energy = 0;

    // The boundaries of the first plaquette 
    set <Scaffold *> boundary_scaffolds1; 
    //set <Scaffold *> inner_scaffolds1; 
    for (Site* sitePtr1 : plaquette_sites1) {
        sitePtr1->scan_through_scaffolds(boundary_scaffolds1);//,inner_scaffolds1);
    }
    // The boundaries of the second plaquette 
    set <Scaffold *> boundary_scaffolds2; 
    //set <Scaffold *> inner_scaffolds2; 
    for (Site* sitePtr2 : plaquette_sites2) {
        sitePtr2->scan_through_scaffolds(boundary_scaffolds2);//,inner_scaffolds2);
    }
    
    // Take the union of the two sets, with the intersection only one
    for (auto it3 = boundary_scaffolds2.begin(); it3 != boundary_scaffolds2.end(); ++it3) {
        if (boundary_scaffolds1.find(*it3) == boundary_scaffolds1.end()) {
            // Element not found in set1, add it
            boundary_scaffolds1.insert(*it3);}
    }


    //cout<<boundary_scaffolds.size()<<" boundary scaffolds"<<endl;
    for (auto it = boundary_scaffolds1.begin(); it != boundary_scaffolds1.end(); ++it) {
        energy+=(*it)->get_energy(); 
        /* if ((*it)->get_energy()!=0){
            cout<<(*it)->get_energy()<<" ";} */
    }

    // This is to debug and calculate the energy of the bonds inside the plaquette 
    /* float energy_inner1 = 0;
    float energy_inner2 = 0;
    for (auto it_in1 = inner_scaffolds1.begin(); it_in1 != inner_scaffolds1.end(); ++it_in1) {
        energy_inner1+=(*it_in1)->get_energy(); 
         if ((*it)->get_energy()!=0){
            cout<<(*it)->get_energy()<<" ";} 
    }
    for (auto it_in2 = inner_scaffolds2.begin(); it_in2 != inner_scaffolds2.end(); ++it_in2) {
        energy_inner2+=(*it_in2)->get_energy(); 
         if ((*it)->get_energy()!=0){
            cout<<(*it)->get_energy()<<" ";} 
    }
    cout<< "energy_inner_1 "<<energy_inner1<<"(" <<inner_scaffolds1.size() <<"scaffolds) energy_inner_2 "<<energy_inner2<<"(" <<inner_scaffolds2.size() <<"scaffolds)"<<endl;
     */
    return energy;   

}




float Site::compute_1plaquette_energy(set<Site*> plaquette_sites){

    float energy = 0;

    // The boundaries of the first plaquette 
    set <Scaffold *> boundary_scaffolds; 
    //set <Scaffold *> inner_scaffolds; 
    for (Site* sitePtr : plaquette_sites) {
        sitePtr->scan_through_scaffolds(boundary_scaffolds);//,inner_scaffolds);
    }


    //cout<<boundary_scaffolds.size()<<" boundary scaffolds"<<endl;
    for (auto it = boundary_scaffolds.begin(); it != boundary_scaffolds.end(); ++it) {
        energy+=(*it)->get_energy(); 
    }

    /* // This is to debug and calculate the energy of the bonds inside the plaquette 
    float energy_inner = 0;
    for (auto it_in = inner_scaffolds.begin(); it_in != inner_scaffolds.end(); ++it_in) {
        cout<< (*it_in)->getStructure()<<" ";
        energy_inner+=(*it_in)->get_energy(); 
       
    }
    
    cout<< "energy_inner "<<energy_inner<<endl; */
    

  
    return energy;   

}











pair<bool,float> Site::attemptSwapPlaquetteTo(Site &destinationSite, float Beta, int plaquette_radius, vector <Site*> &full_to_empty_sites,  vector <Site*> &empty_to_full_sites){
    
    pair<bool, float> successAndDeltaE = {true, 0} ;//the boolean represents the success of swap, and the int is delta_energy
    if (m_orientationRef>0 && destinationSite.isEmpty()){
        
        // Save the old state and the old neighbors         
        //std::vector <std::vector <int>> origin_old_orientations;
        std::vector <int> origin_old_orientations;
        set <Site*> origin_plaquette_sites;
        vector<Site*> origin_plaquette_sites_list;

        //std::vector <std::vector <int>> destination_old_orientations;
        std::vector <int> destination_old_orientations;
        set <Site*> destination_plaquette_sites;
        vector<Site*> destination_plaquette_sites_list;
        
        //if (plaquette_radius>3){
        //cout<<"get orientations origin plaquette "<<endl;}
        //this->add_orientation_recursive(0, plaquette_radius, origin_old_orientations, origin_plaquette_sites);
        this->add_orientation_recursive_flat(0, plaquette_radius, origin_plaquette_sites, origin_plaquette_sites_list);
        
        //if (plaquette_radius>3){
        //cout<<"get orientations destination plaquette "<<endl;}
        //destinationSite.add_orientation_recursive(0, plaquette_radius, destination_old_orientations, destination_plaquette_sites);
        this->add_orientation_recursive_flat(0, plaquette_radius, destination_plaquette_sites, destination_plaquette_sites_list);
        
        // Check that the two plaquettes do not overlap 
        bool intersectionEmpty = true;
        for (Site* elem : origin_plaquette_sites) {
            if (destination_plaquette_sites.count(elem) > 0) {
                intersectionEmpty = false;
                successAndDeltaE.first=false;

                //cout<<"plaquette overlap"<<endl;   
                break;
            } 
        }
        //cout << "length plaquette "<<origin_plaquette_sites.size()<<endl;

        if (intersectionEmpty){
            //cout<<"plaquette don't overlap"<<endl;   
            // Calculate the old energy 
            //cout<<"old energy "<<endl;
            float old_energy = compute_2plaquette_energy(destination_plaquette_sites,origin_plaquette_sites);
            //cout <<"old energy "<<old_energy<<endl;

            for (int i=0 ; i<origin_plaquette_sites.size(); i++){
                origin_old_orientations.push_back((origin_plaquette_sites_list[i])->m_orientationRef);
                destination_old_orientations.push_back((destination_plaquette_sites_list[i])->m_orientationRef);

                (origin_plaquette_sites_list[i])->setFullSite(destination_old_orientations[i]);
                (destination_plaquette_sites_list[i])->setFullSite(origin_old_orientations[i]);
            }


            // Do the swap
            /* vector <int> indices_origin(plaquette_radius, 0); 
            vector <int> indices_destination(plaquette_radius, 0); 
            set <Site*> explored_sites01 ;
            set <Site*> explored_sites02 ;
            if (plaquette_radius>3){
            cout<<"new orientation origin "<<endl;}
            this->set_new_orientations_plaquette(0, plaquette_radius, indices_origin, destination_old_orientations, full_to_empty_sites, empty_to_full_sites,explored_sites01);
            if (plaquette_radius>3){
            cout<<"new orientation destination "<<endl;}
            destinationSite.set_new_orientations_plaquette(0, plaquette_radius, indices_destination, origin_old_orientations, full_to_empty_sites, empty_to_full_sites,explored_sites02);
     */
            // Calculate new energy 

            //cout<<"new energy "<<endl;
            float new_energy = compute_2plaquette_energy(destination_plaquette_sites,origin_plaquette_sites);
            //cout <<"new energy "<<new_energy<<endl;

            successAndDeltaE.second=(new_energy-old_energy);
            if (successAndDeltaE.second>0){
            //if (true){
                
                float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

                if (exp (- static_cast <float> (successAndDeltaE.second*Beta)) < r){
                //if (true){
                    // Do the reverse swap 
                    for (int i=0 ; i<origin_plaquette_sites.size(); i++){
                
                        (origin_plaquette_sites_list[i])->setFullSite(origin_old_orientations[i]);
                        (destination_plaquette_sites_list[i])->setFullSite(destination_old_orientations[i]);
                    }   

                    /* vector <int> indices_origin2(plaquette_radius, 0); 
                    vector <int> indices_destination2(plaquette_radius, 0);

                    set <Site*> explored_sites1 ;
                    set <Site*> explored_sites2 ;
                    //cout<<"reverse orientation origin"<<endl;
                    this->set_new_orientations_plaquette(0, plaquette_radius, indices_origin2, origin_old_orientations, full_to_empty_sites, empty_to_full_sites, explored_sites1);
                    //cout<<"reverse orientation destination"<<endl;
                    destinationSite.set_new_orientations_plaquette(0, plaquette_radius, indices_destination2, destination_old_orientations, full_to_empty_sites, empty_to_full_sites,explored_sites2);
                    */
                    successAndDeltaE.first=false;
                    successAndDeltaE.second=0;

                    // float last_energy = compute_2plaquette_energy(destination_plaquette_sites,origin_plaquette_sites);
                    //cout<<"before move "<<old_energy<<" after move "<<new_energy<<" after reverse "<<last_energy<<endl;
                }
                
            }
            if (successAndDeltaE.first){
            //auto site_iterator3 = plaquette_sites.begin();
            for (int i=0 ; i<origin_plaquette_sites.size(); i++){
                if (origin_old_orientations[i]==0 and destination_old_orientations[i]!=0){
                    full_to_empty_sites.push_back(destination_plaquette_sites_list[i]);
                    empty_to_full_sites.push_back(origin_plaquette_sites_list[i]);
                }
                if (destination_old_orientations[i]==0 and origin_old_orientations[i]!=0){
                    empty_to_full_sites.push_back(destination_plaquette_sites_list[i]);
                    full_to_empty_sites.push_back(origin_plaquette_sites_list[i]);
                }
                //site_iterator3++;
        }
        }

            
            if (full_to_empty_sites.size()!=empty_to_full_sites.size()){
                cout<<"problem counting the full ("<<full_to_empty_sites.size()<<") and empty("<<empty_to_full_sites.size()<<")"<<endl;
                    for (auto it_origin = origin_plaquette_sites.begin(); it_origin != origin_plaquette_sites.end(); ++it_origin) {
                        cout<<"in origin "<< *it_origin<<" "<<endl;
                    }
                    for (auto it_dest = origin_plaquette_sites.begin(); it_dest != origin_plaquette_sites.end(); ++it_dest) {
                        cout<<"in dest "<< *it_dest<<" "<<endl;
    
                    }

            }
            /* if (successAndDeltaE.first){
                cout<<"success!!! de="<<successAndDeltaE.second<<endl;
            } */
        }
        
    }
    else{
        cout<< "The initial site was empty or the destination site was not. It can not be swapped" <<endl;
    }
        
//    cout << " success DE "<< sucessAndDeltaE.first << " "<< sucessAndDeltaE.second<< endl;
    return successAndDeltaE;
}





pair <bool, float> Site::attemptFlipPlaquette(int rotation, float Beta, int plaquette_radius, vector<int> corresponding_indices, vector<Site*> &full_to_empty_sites,  vector<Site*> &empty_to_full_sites){
        pair<bool, float> successAndDeltaE = {true, 0} ;

        // get the plaquette 
        std::vector <int> old_orientations2;

        std::vector <int> old_orientations;
        std::vector <int> new_orientations;
        set<Site*> plaquette_sites;
        vector<Site*> plaquette_sites_list;
        //cout<<endl<<"before running flat"<<endl;
        //this->add_orientation_recursive_flat(0, plaquette_radius, old_orientations2, plaquette_sites,plaquette_sites_list);
        this->add_orientation_recursive_flat(0, plaquette_radius, plaquette_sites,plaquette_sites_list);

        //cout<<"compute old energy"<<endl;
        float old_energy = compute_1plaquette_energy(plaquette_sites);

        //cout<<plaquette_sites.size()<<" sites in plaquette of radius "<<plaquette_radius<<endl;

        // save the orientations of the plaquette
       
        //cout<<"save the orientations"<<endl;
        for (int j=0 ; j<plaquette_sites.size(); j++){

            old_orientations.push_back((plaquette_sites_list[j])->m_orientationRef);
            //cout<<"site position "<<(plaquette_sites_list[j])->get_real_space_position(10)[0]<<" "<<(plaquette_sites_list[j])->get_real_space_position(10)[1]<<endl;
            //cout<<"address "<<(*site_iterator0)<<" orientation "<<(*site_iterator0)->m_orientationRef<<endl;
            }


        // update all the sites 
        
        //cout<<"change the orientations"<<endl;
        for (int i=0 ; i<plaquette_sites.size(); i++){
            //cout<<"site position "<<(*site_iterator)->get_real_space_position(10)[0]<<" "<<(*site_iterator)->get_real_space_position(10)[1]<<endl;
            
            int old_orientation_i = old_orientations[corresponding_indices[i]];
            /* cout<<"position "<<plaquette_sites_list[i]->get_real_space_position(10)[0]<<" "<<plaquette_sites_list[i]->get_real_space_position(10)[1]<<" ";
            cout<<"old position "<<plaquette_sites_list[corresponding_indices[i]]->get_real_space_position(10)[0]<<" "<<plaquette_sites_list[corresponding_indices[i]]->get_real_space_position(10)[1]<<endl;
            cout<<"index "<<i<<" is from "<<corresponding_indices[i]<<endl;
            cout<<endl; */

            if (old_orientation_i==0){
                new_orientations.push_back(old_orientation_i);
                //cout << "old o new site:"<<old_orientations[i]<<" o old site:"<<old_orientation_i<<" new o new site:"<<(old_orientation_i ) <<endl;
            }
            else{
                new_orientations.push_back( (old_orientation_i + rotation-1)%6 +1 ); // CAREFUL , this line is specific to the triangular lattice
                //cout << "old o new site:"<<old_orientations[i]<<" o old site:"<<old_orientation_i<<" new o new site:"<<(old_orientation_i + rotation-1)%6 +1 <<endl;
            }
            
            (plaquette_sites_list[i])->setFullSite(new_orientations[i]);
            
        }

        //cout<<"compute new energy"<<endl;
        float new_energy = compute_1plaquette_energy(plaquette_sites);
        successAndDeltaE.second=(new_energy-old_energy);
        if (successAndDeltaE.second>0){
            
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

            if (exp (- static_cast <float> (successAndDeltaE.second*Beta)) < r){
                // we do the reverse swap 
                //auto site_iterator2 = plaquette_sites.begin();
                for (int i=0 ; i<plaquette_sites.size(); i++){
                    (plaquette_sites_list[i])->setFullSite(old_orientations[i]);
                    //site_iterator2++;
                }
                successAndDeltaE.first = false;
                
            }

        }
        if (successAndDeltaE.first){
            //auto site_iterator3 = plaquette_sites.begin();
            for (int i=0 ; i<plaquette_sites.size(); i++){
                if (new_orientations[i]==0 and old_orientations[i]!=0){
                    full_to_empty_sites.push_back(plaquette_sites_list[i]);
                }
                if (new_orientations[i]!=0 and old_orientations[i]==0){
                    empty_to_full_sites.push_back(plaquette_sites_list[i]);
                }
                //site_iterator3++;
        }
        }
        if (full_to_empty_sites.size()!=empty_to_full_sites.size()){
                cout<<"problem counting the full ("<<full_to_empty_sites.size()<<") and empty("<<empty_to_full_sites.size()<<")"<<endl;
                    /* for (auto it_origin = plaquette_sites.begin(); it_origin != plaquette_sites.end(); ++it_origin) {
                        cout<<"in origin "<< *it_origin<<" "<<endl;
                    }
                    for (auto it_dest = plaquette_sites.begin(); it_dest != origin_plaquette_sites.end(); ++it_dest) {
                        cout<<"in dest "<< *it_dest<<" "<<endl;
    
                    } */

            }

    return successAndDeltaE;
}










vector <Site*> Site::commonNeigbhours(Site &otherSite, int degreeOfNeighbourhood) const {
    vector <Site *> :: const_iterator it_this ;
    vector <Site *> :: const_iterator it_other ;
    vector <Site *> result;
    for (it_this = m_neighbourList[degreeOfNeighbourhood].begin(); it_this != m_neighbourList[degreeOfNeighbourhood].end() ; ++it_this){
        for (it_other = (otherSite.m_neighbourList[degreeOfNeighbourhood]).begin(); it_other !=(otherSite.m_neighbourList[degreeOfNeighbourhood]).end() ;
            it_other++){
            if ((*it_this)==(*it_other)){
                result.push_back(*it_this);
            }
        }
}
    return result;
}

vector <Scaffold*> Site::commonScaffold(Site &otherSite, int degreeOfNeighbourhood) const{
    vector <Scaffold *> :: const_iterator it_this ;
    vector <Scaffold *> :: const_iterator it_other ;
    vector <Scaffold *> result;
    for (it_this = belongingScaffolds[degreeOfNeighbourhood].begin(); it_this != belongingScaffolds[degreeOfNeighbourhood].end() ; ++it_this){
        for (it_other = (otherSite.belongingScaffolds)[degreeOfNeighbourhood].begin(); it_other !=(otherSite.belongingScaffolds)[degreeOfNeighbourhood].end() ;
            it_other++){
            if ((*it_this)==(*it_other)){
                result.push_back(*it_this);
            }
        }
}
    return result;
}

bool Site::isEmpty() const{
    return m_orientationRef==0;
}

bool Site::equals(Site siteToCompare){
    return (this == &siteToCompare);
}

bool operator==(Site &site1, Site &site2){
    return site1.equals(site2);
}

bool operator!=(Site &site1, Site &site2){
    return !site1.equals(site2);
}
/*
vector <int> Site::getNeighbourOrientationRef () const {
    long N_neighbour = m_neighbourList.size();
    vector <int> NeighbourOrientation (N_neighbour);
    for (int i=0 ; i<N_neighbour ; i++)
    {
        NeighbourOrientation[i]=(*(m_neighbourList[i])).m_orientationRef ;
    }
    return NeighbourOrientation ;
}*/


bool Site::isNeighbourOf(Site& site2, int degreeOfNeighbourhood) const {
    //vector <Site*> &listNeighbour (site2.m_neighbourList);
    bool answ ;
    if (site2.m_neighbourList[degreeOfNeighbourhood].size()==0){
        
        answ = false;
    }
    else{
        answ=(std::find(site2.m_neighbourList[degreeOfNeighbourhood].begin(), site2.m_neighbourList[degreeOfNeighbourhood].end(), this) != site2.m_neighbourList[degreeOfNeighbourhood].end()) ;
    }
    return answ;
}

bool Site::isSecondNeighbour(Site& site2, int degreeOfNeighbourhood) const{
    vector <Site *> :: const_iterator i ;
    bool answ (false); //default is false, whenever find a second neighbour relation, return true.
    for (i = (site2.m_neighbourList)[degreeOfNeighbourhood].begin(); i != (site2.m_neighbourList)[degreeOfNeighbourhood].end() ; ++i){
        if (this->isNeighbourOf((**i),degreeOfNeighbourhood) ){
            answ=true;
            break;
        }
    }
    return answ ;
}


