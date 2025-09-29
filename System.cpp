//
//  System.cpp
//  FrustratedSystems
//
//  Created by lara koehler on 24/09/2020.
//  Copyright © 2020 lara koehler. All rights reserved.
//

#include "System.hpp"
#include "Site.hpp"
#include "SystemType.hpp"
#include "Scaffold.hpp"
#include "head.hpp"
#include "small_classes.hpp"
#include "Plaquette.hpp"

//#include <sstream>

#include "json.hpp"
using json = nlohmann::json;

#include "LogDuration.hpp"

#include <ctime> // For random
#include <time.h>
#include <math.h>
#include <cstdlib> // For random
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <chrono>
#include <set>  
#include <random>

using namespace std;


System::System(json parameters_arg) :

scaffoldList(vector <Scaffold*>(0)),

extendedScaffoldList(vector <Scaffold*>(0)),
siteList(vector <Site*>(0)),
RealPositionMappedToSiteList(map<int,int>()),
fullSites(vector <Site*>(0)),
emptySites(vector <Site*>(0)),
Etot(0),
Emin(1000),
Lx(0), Ly(0), Lz(0),
NhashPosition(0),
Nsites(0), NparticleTypes(0), Nparticles(vector<int>(0)), NparticlesTot(0), N_Orientation(0),
latticeFilename(0),latticeTypeName(),N_Structures(0),LEL(vector <float>()),
Nsteps_per_T(0), Nmes_per_T(0),Ntemperatures(0),
parameters(parameters_arg),
directoryForResults(),
saveSnapshot(false), saveEnergies(false),saveTemperatures(false), saveMCsuccess(false),
oneGroundStateSaved(false),
AllEnergies(vector<double>(0)),AllMCsuccess(vector<float>(0)),
cVector(vector<int>(0)), dVector(vector<int>(0)),
//scaffoldsStats(vector <std::vector<int>>(0)),extendedScaffoldsStats(vector <std::vector<int>>(0)),
typeOfThisSystem(SystemType()),
cluster_probabilities(vector<float>(0)), index_correspondance_plaquette(vector<vector<vector<int>>>(((0)))),
plaquettes_list(vector<vector<Plaquette*>>()), compatible_plaquettes(vector<vector<pair<Plaquette*, Plaquette*>>>()), max_cluster_radius(0)
{
    
//    std::ifstream f(input_parameters);
//    json parameters = json::parse(f);
    
    unsigned long seed = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
    srand((int)seed);
    
    
    string const latticeType = parameters["lattice"];
    latticeTypeName = parameters["lattice"];
    Lx = parameters["Lx"];
    Ly = parameters["Ly"];
    Lz = parameters["Lz"];
    
        
    NparticleTypes = parameters["NparticleTypes"];
    Nparticles = parameters["Nparticles"].get<std::vector<int>>();
    for (int i=0 ; i<NparticleTypes ; i++){
        NparticlesTot+=Nparticles[i];
    }
    

    FlipRate = parameters["FlipRate"];
    SwapRate = parameters["SwapRate"];
    FlipSwapRate = parameters["FlipSwapRate"];
    MutationRate = parameters["MutationRate"];
    PlaquetteRate = 1-FlipRate -SwapRate - FlipSwapRate - MutationRate;

    // Precompute the rate of each possible transition
    T1 = SwapRate;
    T2 = T1 + FlipRate;
    T3 = T2 + FlipSwapRate;
    T4 = T3 + MutationRate;
    T5 = T4 + PlaquetteRate / 2.0f;
    
    Ntemperatures = parameters["Ntemperatures"];
    Nsteps_per_T = parameters["Nsteps_per_T"];
    Nmes_per_T = parameters["Nmes_per_T"];
    
   
    
    
    string const directoryInputs= parameters["directoryInputs"];
    directoryForResults = parameters["directoryResults"];
    string LELfileSuffix = parameters["LELfileSuffix"] ;
    saveTemperatures=parameters["saveTemperatures"] ;
    saveSnapshot= parameters["saveSnapshot"] ;
    saveEnergies = parameters["saveEnergies"] ;
    saveMCsuccess = parameters["saveMCsuccess"] ;
    
    
    
    //Initialize lattice, structure map, and energies in other class
    
    string const latticeFile = directoryInputs + latticeType + ".json" ;
    string const structuresFile = directoryInputs + latticeType + "_structures_"+to_string(NparticleTypes)+"particles.txt" ;
    string const LELfile = directoryForResults + LELfileSuffix ;
    
    typeOfThisSystem.initiate_values(structuresFile, latticeFile,  LELfile, NparticleTypes);
    typeOfThisSystem.initialize_mutation_and_flip();
        
    cVector = vector<int> (N_Structures, 0);
    dVector = vector<int> (N_Structures, 0);
        
    N_Orientation=typeOfThisSystem.getNOrientation();

    NhashPosition = max(Lx, Ly);
   // cout<<"NhashPosition"<<endl;
    NhashPosition = max(NhashPosition, Lz);
    N_Structures = typeOfThisSystem.get_N_Structures();

    //max_cluster_radius = int((min(Lx,Ly)-2)/3)+1;
    max_cluster_radius = min(int(sqrt(NparticlesTot)/3)+2,int((min(Lx,Ly)-2)/3)+1) ;
    max_cluster_radius = min(7, max_cluster_radius);
    float normalization_cluster_probabilities = 0;
    cout <<"Maximal radius of clusters of particle "<<max_cluster_radius<<endl;

    for (int i=1; i<max_cluster_radius; i++){
        normalization_cluster_probabilities += 1/(float(i)*float(i));
        cluster_probabilities.push_back(1/(float(i)*float(i)));
    }
    for (int i=0; i<max_cluster_radius-1; i++){
        cluster_probabilities[i]=cluster_probabilities[i]/normalization_cluster_probabilities ;
        if (i>0){
            cluster_probabilities[i]+=cluster_probabilities[i-1];
        }
        //cout<<"Probability of drawing cluster of radius "<<i+1<<" or lower is "<<cluster_probabilities[i]<<endl;
    }
    
    }




vector<Site*>  System::get_full_sites() const{
    return fullSites;
}

vector<Site*>  System::get_empty_sites() const{
    return emptySites;
}

Site* System::get_site(int index) const{
    return siteList[index] ;
    
}


string System::get_directoryResults(){
    return directoryForResults;
}

void System::build_lattice(int const Lx, int const Ly, int const Lz){
    
    LOG_DURATION(to_string(Lx)+"*"+to_string(Ly)+"*"+to_string(Lz)+" "+ latticeTypeName+ " lattice is initialized");

    int cursorInSiteList (0);
    for (int i=0 ; i<Lx ; i++){
        for(int j=0 ; j<Ly ; j++){
            for(int k=0 ; k<Lz ; k++){
//                cout<<"ijk"<<i<<j<<k<<endl;
                
                Site* current_site = new Site();
                (*current_site).set_position(i, j, k, NhashPosition, NhashPosition*NhashPosition);
                int pos ((*current_site).get_hashedPosition());
                siteList.push_back(current_site);
                //build 1-particle scaffold (orientations are 0 because sites are empty, edge is 0)
                Scaffold* currentScaffold = new Scaffold(0, typeOfThisSystem );
                current_site->setScaffoldAppartenance(currentScaffold, 0, 0);//current_site is member of current_scaffold, and is the 0th site of this scaffold(which will be a 1p scaffold anyway)
                scaffoldList.push_back(currentScaffold);
//                cout<<"scaffold ok"<<endl;
                Scaffold* currentExtendedScaffold = new Scaffold(0, typeOfThisSystem);
                current_site->setScaffoldAppartenance(currentExtendedScaffold, 0, 1);
                extendedScaffoldList.push_back(currentExtendedScaffold);                
                RealPositionMappedToSiteList[pos] =cursorInSiteList;
                
                //cout <<"ij "<<i<<" "<<j<<endl;
                //cout <<"pos "<<pos<<endl;
                //cout<< "cursorInSiteList "<<cursorInSiteList<<endl;
                cursorInSiteList++;

            }
        }
    }
    emptySites=siteList;
    Nsites=(int) siteList.size();
    //cout<<"A total of "<<Nsites<<" sites were build at the end of build_lattice function"<<endl;
    
    //create the file with the ref between adress and position of sites
    if (saveSnapshot){
        ofstream myStream ((directoryForResults+"_sites_positions.txt").c_str());
        vector <Site*> :: iterator it;
        if (!myStream){
            cout << endl << "CAREFUL, the folder "<<directoryForResults<<" does not exist !!" <<endl<<endl;
        }
        for (it=siteList.begin(); it!= siteList.end(); it++){
            myStream << (*it);
            vector<int> position = (*it)->get_real_space_position(NhashPosition);
            myStream << " "<< position[0] <<" "<<position[1]<<" "<<position[2]<<endl;
        }
        myStream<<"#END";
        myStream.flush();
        myStream.close();
        myStream.clear();
    }
    
}

void System::build_neighbour_relations( int degreeOfNeighbourhood, int distanceOfRelation, vector <Scaffold *> &scaffoldListToFill){
    // degreeOfNeighbourhood is 0 or 1 : it is an index (we only look at first or second neighbour
    // but distance of Relation can be any integer : the second neighbor can be at distance 3 of the particle for example
    vector <Site *> :: const_iterator it;
//    cout << "Will build neighbors of distance "<<distanceOfRelation<<endl;
    
    //cout<<"A total of "<<scaffoldListToFill.size()<<" scaffolds existed at the beging of build_neighbour_relation function"<<endl;
    // go across the list of all neighbours and connect it to its neighbours according to coordinates
    for (it = siteList.begin(); it != siteList.end() ; ++it){
        //**it is the site for which we will set the neighbours
        vector <int> real_position = (**it).get_real_space_position(NhashPosition);
        
        for (int edge_ref = 0 ; edge_ref<(typeOfThisSystem.get_bonds()).size() ; edge_ref++){
            int iNeighbour=uModulus(real_position[0]+distanceOfRelation*(typeOfThisSystem.get_bonds())[edge_ref][0],Lx);
            int jNeighbour=uModulus(real_position[1]+distanceOfRelation*(typeOfThisSystem.get_bonds())[edge_ref][1],Ly);
            int kNeighbour=uModulus(real_position[2]+distanceOfRelation*(typeOfThisSystem.get_bonds())[edge_ref][2],Lz);
            // this is the direction of the bond
            //cout << "target neighbour position :"<<iNeighbour<<" "<<jNeighbour<<" "<<kNeighbour<<endl;
            int neighbour_position_siteList = RealPositionMappedToSiteList.at(hash_ints(iNeighbour,jNeighbour,kNeighbour,NhashPosition,NhashPosition*NhashPosition));//we find the path to the neighbour we want to bind
            
            if ((**it).setNeighbourRelation(*siteList[neighbour_position_siteList],degreeOfNeighbourhood)){ //this sets the neighbou relation if it didn't already exist
                /* if (degreeOfNeighbourhood==0){
                cout<<"edge "<<(typeOfThisSystem.get_bonds())[edge_ref][0]<<" "<<(typeOfThisSystem.get_bonds())[edge_ref][1]<<endl;
                cout<<"start"<<real_position[0]<<" "<<real_position[1]<<endl;
                cout<<" neighbour_position_siteList "<<neighbour_position_siteList<<endl;
                cout << "to "<<(siteList[neighbour_position_siteList])->get_real_space_position(NhashPosition)[0]<<" "<<(siteList[neighbour_position_siteList])->get_real_space_position(NhashPosition)[0]<<endl;
                cout<<"meaning "<<iNeighbour<<" "<<jNeighbour<<endl;
                
                //PRINT FOR DEBUG
                string positionString = to_string(real_position[0])+ " " + to_string(real_position[1])+ " "+ to_string(real_position[2]);
                cout<<"building a new scaffold between "<<positionString <<" and "<<iNeighbour<<" "<<jNeighbour<<" "<<kNeighbour<< endl;
                } */
                //if this relation didn't exist, we also create a 2-particles-scaffold with these two sites, that we add to the list. 
                Scaffold* currentScaffold = new Scaffold(edge_ref+1, typeOfThisSystem);
                
//                cout <<" don "<< degreeOfNeighbourhood<<"->";
                (*it)->setScaffoldAppartenance(currentScaffold, 0, degreeOfNeighbourhood); //the considered site is the 1st member of the new scaffold
                siteList[neighbour_position_siteList]->setScaffoldAppartenance(currentScaffold, 1, degreeOfNeighbourhood); //its neighbour is the 2nd member of the same scaffold
                scaffoldListToFill.push_back(currentScaffold);
            }
        }
    
    }
    //cout<<"A total of "<<scaffoldListToFill.size()<<" scaffolds were build at the end of build_neighbour_relation function"<<endl;
    
}


void System::reorder_all_neighbors(){
    vector <Site *> :: const_iterator it;

    vector <vector<int>>  order_defaults ; 
    vector< vector<int >> order_first0 = (*siteList.begin())->reorder_neighbors(0, NhashPosition, Lx, Ly, Lz, order_defaults);
    vector< vector<int >> order_first1 = (*siteList.begin())->reorder_neighbors(1, NhashPosition, Lx, Ly, Lz, order_defaults);
     
    for (it = siteList.begin(); it != siteList.end() ; ++it){
        (*it)->reorder_neighbors(0, NhashPosition,Lx, Ly, Lz, order_first0);
        (*it)->reorder_neighbors(1, NhashPosition, Lx, Ly, Lz, order_first1);
    }

}





void System::fill_the_system(){
    cout << "The density of particles is "<<((float)NparticlesTot/(float)Nsites)<< " -> "<< NparticlesTot <<" particles" <<endl;
    //LOG_DURATION("System is initialized with particles");
    
    for (int particle_type=0 ; particle_type<NparticleTypes ; particle_type++){
        for (int i=0 ; i<Nparticles[particle_type] ; i++){
            //draw a random site in siteList
            //remove its adress from emptysites and put it in fullsites
            //give it a random orientation
            //actualize the corresponding scaffold
            if (emptySites.size()<0){
                cout << "Trying to fill the system with too much sites"<<endl;
            }
            
            int index = rand()%emptySites.size() ;
            Site* SiteToFill = remove_at(emptySites, index);
            int site_state = typeOfThisSystem.random_state(particle_type);
            //cout<<" New site filled with orientation "<<site_state<<endl;
            (*SiteToFill).setFullSite(site_state); //this function also actualizes the scaffold value
            fullSites.push_back(SiteToFill);
        }
    }

    /* cout<<" IN FILL THE SYSTEM "<<endl;
    cout <<endl<<"empty "<<emptySites.size()<<" ";
    for (int index_empty1=0; index_empty1<emptySites.size(); index_empty1++){
        cout<<emptySites[index_empty1]<<" ";
    }
    cout<<endl<<"full "<<fullSites.size()<<" ";
    for (int index_full1=0; index_full1<fullSites.size(); index_full1++){
        cout<<fullSites[index_full1]<<" ";
    }
    cout<<endl; */

    
    Etot=compute_Etot();
    print_energy_per_particle();
    
}

void System::print_energy_per_particle(){
    cout << "Energy per particle : "<<round(Etot*100/NparticlesTot)/100<<endl;
}



double System::compute_Etot(){
    double E = 0;
    vector <Site*> :: iterator it;
    for (it = fullSites.begin() ;it!=fullSites.end() ; it++){
        E+=(*it)->getSiteEnergy();
        
    }
    for (it = emptySites.begin() ;it!=emptySites.end() ; it++){
        E+=(*it)->getSiteEnergy();
    }
    E = E /2.0; //each link was counted twice
    //cout << E<< " ";
    return E;
}



vector <int> System::get_scaffold_list_2sites(Site* s1, Site* s2, int degreeOfNeighborhood){
//    Collects all the sturctures involved in two sites
    
    vector <int> Nvector (N_Structures, 0) ;
    set <Scaffold* > changingScaffolds ;
    
    for ( Scaffold* scaffAdress : s1->getBelongingScaffolds(degreeOfNeighborhood))
    {changingScaffolds.insert(scaffAdress);
    }
    for ( Scaffold* scaffAdress : s2->getBelongingScaffolds(degreeOfNeighborhood))
    {changingScaffolds.insert(scaffAdress);
    }
    
    for ( Scaffold* scaffAdress : changingScaffolds){
        Nvector[scaffAdress->getStructure()]++;
    }
    return Nvector ;
}


vector <int> System::get_scaffold_list_1site(Site* s1, int degreeOfNeighborhood){
//    Collects all the sturctures involved in one sites
    
    vector <int> Nvector (N_Structures, 0) ;
    
    for ( Scaffold* scaffAdress : s1->getBelongingScaffolds(degreeOfNeighborhood)){
        Nvector[scaffAdress->getStructure()]++;
    }
    
    return Nvector ;
}


pair<bool, float> System::MC_step_flip_plaquette(float Beta, bool saveStat){
    pair<bool, float> successAndDeltaE;
    
    ///Choose the radius of the plaquette
    int plaquette_radius=1;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    for (int i=0; i<cluster_probabilities.size(); i++){
        if (r<= cluster_probabilities[i]){
            plaquette_radius = i+1;
            break;  
        }
    }


    //cout<<"radius "<<plaquette_radius<<" ";

    int indexInPlaquettes = rand()%Nsites;
    Plaquette* thePlaquette = plaquettes_list[plaquette_radius-1].at(indexInPlaquettes);
    // We make sure that the center of the plaquette is full 
    //while (thePlaquette->centerIsEmpty()){
    while(false){
        indexInPlaquettes++;
        indexInPlaquettes = indexInPlaquettes%Nsites;
        thePlaquette = plaquettes_list[plaquette_radius-1].at(indexInPlaquettes);
    }
    
    /* if (plaquette_radius<=4){
        plaquette_radius=10;
    } */

    //cout<<endl;
    //cout <<"orientation plaquette center "<<theFullSite->getOrientationRef();
    //cout<<"nhash"<<NhashPosition<<endl;
    //cout<< ", position "<<theFullSite->get_real_space_position(NhashPosition)[0]<< " "<< theFullSite->get_real_space_position(NhashPosition)[1]<<endl;

    vector <Site *> full_to_empty_sites ; 
    vector <Site *> empty_to_full_sites ; 
    int rotation = rand()%5 + 1;
    //cout<<"rotation "<<rotation<<endl;
    successAndDeltaE = thePlaquette->attemptFlipPlaquette(rotation, Beta, index_correspondance_plaquette[plaquette_radius-1][rotation-1], full_to_empty_sites, empty_to_full_sites, N_Orientation);



    if (successAndDeltaE.first){
        //cout<<"radius "<<plaquette_radius<<endl;

        if (saveStat){
            cVector = scaffold_stats(scaffoldList);
            dVector = scaffold_stats(extendedScaffoldList);
        }

        //cout<<"check full sites"<<endl;

        //cout<<endl<<"treat empty to full "; 
        for (int index_empty=0; index_empty<empty_to_full_sites.size(); index_empty++){
            fullSites.push_back(empty_to_full_sites[index_empty]);
            emptySites.erase(std::remove(emptySites.begin(), emptySites.end(), (empty_to_full_sites[index_empty])), emptySites.end());
            //cout<<empty_to_full_sites[index_empty]<<" ";
        } 

        //cout<<endl<<"treat full to empty "; 
        for (int index_full=0; index_full<full_to_empty_sites.size(); index_full++){
            emptySites.push_back(full_to_empty_sites[index_full]);
            fullSites.erase(std::remove(fullSites.begin(), fullSites.end(), (full_to_empty_sites[index_full])), fullSites.end());
            //cout<<full_to_empty_sites[index_full]<<" ";
        }

    }

    return successAndDeltaE;

}



pair<bool, float> System::MC_step_swap_plaquette(float Beta, bool saveStat){
   //vector <int> Ci, Cf (N_Structures, 0) ;
    //vector <int> Di, Df (N_Structures, 0) ;
    pair<bool, float> successAndDeltaE;
    

    //Choose the radius of the plaquette
    int plaquette_radius=1;
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    for (int i=0; i<cluster_probabilities.size(); i++){
        if (r<= cluster_probabilities[i]){
            plaquette_radius = i+1;
            break;  
        }
    }
    
    
    //choose which sites to exchange
    int indexInPlaquettes = rand()%compatible_plaquettes[plaquette_radius-1].size();
    
    Plaquette* plaquette1 = compatible_plaquettes[plaquette_radius-1][indexInPlaquettes].first;
    Plaquette* plaquette2 = compatible_plaquettes[plaquette_radius-1][indexInPlaquettes].second;
    //while (plaquette1->centerIsEmpty() && plaquette2->centerIsEmpty()){
    while (false){
        indexInPlaquettes++;
        indexInPlaquettes = indexInPlaquettes%compatible_plaquettes[plaquette_radius-1].size();
        plaquette1 = compatible_plaquettes[plaquette_radius-1][indexInPlaquettes].first;
        plaquette2 = compatible_plaquettes[plaquette_radius-1][indexInPlaquettes].second;
    }
    

    vector <Site *> full_to_empty_sites ; 
    vector <Site *> empty_to_full_sites ; 
    successAndDeltaE = plaquette1->attemptSwapPlaquetteTo(plaquette2, Beta, full_to_empty_sites, empty_to_full_sites);
    //attemptSwapPlaquetteTo(Site &destinationSite, float Temperature, int new_state, int plaquette_radius)
    

    
    //cout<< "plaquette of size "<<plaquette_radius<<" -> success "<<successAndDeltaE.first<<" ,de="<<successAndDeltaE.second<<" "<<endl;
    if (successAndDeltaE.first){
        // if (plaquette_radius>=3){
            //cout<< "plaquette of size "<<plaquette_radius<<" ";//-> success  ,de="<<successAndDeltaE.second<<" "<<endl;
        
        //} */
        //cout<<"radius "<<plaquette_radius<<endl;

        if (saveStat){
            cVector = scaffold_stats(scaffoldList);
            dVector = scaffold_stats(extendedScaffoldList);
        }

        //cout<<"check full sites"<<endl;

       // cout<<endl<<"treat empty to full "; 
        for (int index_empty=0; index_empty<empty_to_full_sites.size(); index_empty++){
            fullSites.push_back(empty_to_full_sites[index_empty]);
            emptySites.erase(std::remove(emptySites.begin(), emptySites.end(), (empty_to_full_sites[index_empty])), emptySites.end());
            //cout<<empty_to_full_sites[index_empty]<<" ";
        } 

        //cout<<endl<<"treat full to empty "; 
        for (int index_full=0; index_full<full_to_empty_sites.size(); index_full++){
            emptySites.push_back(full_to_empty_sites[index_full]);
            fullSites.erase(std::remove(fullSites.begin(), fullSites.end(), (full_to_empty_sites[index_full])), fullSites.end());
           // cout<<full_to_empty_sites[index_full]<<" ";
        } 


        /* for (std::vector<Site*>::iterator it2 = full_to_empty_sites.begin(); it2 != full_to_empty_sites.end(); ++it2) {
            emptySites.push_back(*it2);
            fullSites.erase(std::remove(fullSites.begin(), fullSites.end(), (*it2)), fullSites.end());
        }
        for (std::vector<Site*>::iterator it1 = empty_to_full_sites.begin(); it1 != empty_to_full_sites.end(); ++it1) {
            fullSites.push_back(*it1);
            emptySites.erase(std::remove(emptySites.begin(), emptySites.end(), (*it1)), emptySites.end());
        } */

    }
    
    //cout<<"empty-full sites updated ";
    return successAndDeltaE;

}




pair <bool, float> System::MC_step_swap(float Beta, bool saveStat, bool isFlipSwap){
    vector <int> Ci, Cf (N_Structures, 0) ;
    vector <int> Di, Df (N_Structures, 0) ;
    pair <bool, float> successAndDeltaE;
    
    //choose which sites to exchange
    int indexInFull = rand()%NparticlesTot;
    int indexInEmpty = rand()%(Nsites-NparticlesTot);
    Site* theFullSite = fullSites.at(indexInFull);
    Site* theEmptySite = emptySites.at(indexInEmpty);
    
    if (saveStat){
    Ci = get_scaffold_list_2sites(theFullSite, theEmptySite, 0);
    Di = get_scaffold_list_2sites(theFullSite, theEmptySite, 1);
    }
            
    // attempt to and compute energy
    
    //Get the state of the full Site
    int current_site_state = theFullSite->getSiteStateTemporary() ;
    
    if (not isFlipSwap){
        successAndDeltaE = theFullSite->attemptSwapTo(*theEmptySite, Beta, current_site_state);
    }
    else {
        // Find a new state for the full Site (but this is not a mutation)
        int new_state = typeOfThisSystem.possible_new_state(current_site_state, false);
        successAndDeltaE = theFullSite->attemptSwapTo(*theEmptySite, Beta, new_state);
    }
    
    if (successAndDeltaE.first){
        //the swap was a success
//            Then we will look at the scaffolds that have changed and actualize cVector and dVector
    if (saveStat){

    LOG_DURATION("saveStat in MC_step_swap duration");
    Cf = get_scaffold_list_2sites(theFullSite, theEmptySite, 0);
    Df = get_scaffold_list_2sites(theFullSite, theEmptySite, 1);
    for (int s=0; s<N_Structures;s++){
        cVector[s] += Cf[s];
        cVector[s] -= Ci[s];
        dVector[s] += Df[s];
        dVector[s] -= Di[s];
        if (cVector[s]<0){
            cout<< "We found a negative cVector for structure s="<<s<<" in swap move, isFlipSwap="<<isFlipSwap<<endl;
        }
    }
    }
    
    emptySites.push_back(remove_at(fullSites, indexInFull));
    fullSites.push_back(remove_at(emptySites, indexInEmpty));
   // cout << "Swap done from "<<theFullSite<<" to "<<theEmptySite<<" de="<<successAndDeltaE.second<<endl;
    }
    //else{
    //cout<<"Swap not done"<<theFullSite<<" remains and de="<<successAndDeltaE.second;
    //}
    return successAndDeltaE;
}

pair <bool, float> System::MC_step_flip(float Beta, bool saveStat, bool isMutation){
    vector <int> Ci, Cf (N_Structures, 0) ;
    vector <int> Di, Df (N_Structures, 0) ;
    pair <bool, float> successAndDeltaE;
    
    int indexInFull = rand()%NparticlesTot; // randomly choose which site to flip
    Site* theFlippingSite = fullSites.at(indexInFull);
    
    if (saveStat){
    Ci = get_scaffold_list_1site(theFlippingSite, 0);
    Di = get_scaffold_list_1site(theFlippingSite, 1);
    }
    
    //attempt the flip
    int current_site_state = theFlippingSite->getSiteStateTemporary() ;
    int new_state = typeOfThisSystem.possible_new_state(current_site_state, isMutation);
//    cout<<"old state was "<<current_site_state<<" and new state would be "<<new_state<<endl;
    successAndDeltaE = theFlippingSite->attemptFlip(new_state, Beta);
    
    if (successAndDeltaE.first){
//          It was a success, then we will look at the scaffolds that have changed and actualize cVector and dVector
        if (saveStat){
        Cf = get_scaffold_list_1site(theFlippingSite, 0);
        Df = get_scaffold_list_1site(theFlippingSite, 1);
        for (int s=0; s<N_Structures;s++){
            cVector[s] += Cf[s];
            cVector[s] -= Ci[s];
            dVector[s] += Df[s];
            dVector[s] -= Di[s];
            if (cVector[s]<0){
                cout<< "We found a negative cVector for structure s="<<s<<" in flip move, isMutation="<<isMutation<<endl;
            }
        }
        }
//        cout << "Flip done from "<<current_site_state<<" to "<<new_state<<" de="<<successAndDeltaE.second<<endl;
    }
    else{
//    cout<<"Flip not done "<<current_site_state<<" remains and de="<<successAndDeltaE.second;
    }
    
    return successAndDeltaE;
}




pair <bool, float> System::MC_step(float Beta, bool saveStat){

    //vector <int> Ci, Cf (N_Structures, 0) ;
    //vector <int> Di, Df (N_Structures, 0) ;
    
    //float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    static thread_local std::mt19937 rng(std::random_device{}());
    static thread_local std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r = dist(rng);

    pair <bool, float> successAndDeltaE;
    
    /* T1 = SwapRate;
    T2 = T1 + FlipRate;
    T3 = T2 + FlipSwapRate;
    T4 = T3 + MutationRate;
    T5 = T4 + PlaquetteRate / 2.0f; */

    if (r<=T1){ 
        //LOG_DURATION("Swap duration");
        successAndDeltaE = MC_step_swap(Beta, saveStat, false);
    } else if  (r <=T2) {
        //LOG_DURATION("Flip duration");
        successAndDeltaE = MC_step_flip(Beta, saveStat, false);
    } else if (r <= T3){ 
        //LOG_DURATION("Flip-swap duration");
        successAndDeltaE = MC_step_swap(Beta, saveStat, true);
    } else if  (r <= T4) {
        successAndDeltaE = MC_step_flip(Beta, saveStat, true);
    }  else if (r <= T5) {
        successAndDeltaE = MC_step_swap_plaquette(Beta, saveStat);        
    } else{
        successAndDeltaE = MC_step_flip_plaquette(Beta, saveStat);
    }

    if (successAndDeltaE.first){
        Etot+=successAndDeltaE.second;
    }
    return successAndDeltaE;
}

void System::create_system(int firstNeighborDistance, int secondNeighborDistance){
    //initialization
    build_lattice(Lx, Ly, Lz);
    build_neighbour_relations(0, firstNeighborDistance, scaffoldList);
    build_neighbour_relations(1, secondNeighborDistance, extendedScaffoldList);
    reorder_all_neighbors();
    fill_the_system();
    //build_plaquettes();
    //create_rotations_map();
    LOG_DURATION(directoryForResults+" was created");
    

}


void System::create_rotations_map(){
    // CAREFUL THIS FUNCTION IS ONLY VALID FOR 2D SYSTEM 

    
    vector <int> sizes = {Lx, Ly, Lz};
    /* 
    vector <vector <int>> rotations ;
    for (int sign=-1 ; sign<=1 ; sign+=2){
        for (int edge_ref = 0 ; edge_ref<(typeOfThisSystem.get_bonds()).size() ; edge_ref++){
            
            vector <int> edge_coordinate = typeOfThisSystem.get_bonds()[edge_ref];
            for (int d=0; d<=3 ; d++){
                edge_coordinate[d] = (sign*edge_coordinate[d] + sizes[d])%sizes[d];
            }
            rotations.push_back(edge_coordinate);
        }
    }

    cout<<"rotations"<<endl;
    for (int j=0 ; j<rotations.size(); j++){
        cout<<rotations[j][0]<<" "<<rotations[j][1]<<endl;
    } */

    vector< vector <vector <int>>> rotation_matrices = {};

    if (latticeTypeName=="Triangular"){
    // the 5 rotation matrices in triangular lattice coordinates
    cout<<"triagular lattice rotation matrices"<<endl;
    vector <vector <int>> R = {{0,-1},{1, 1}};
    vector <vector <int>> R2 = {{-1,-1},{1, 0}};
    vector <vector <int>> R3 = {{-1,0},{0, -1}};
    vector <vector <int>> R4 = {{0,1},{-1, -1}};
    vector <vector <int>> R5 = {{1,1},{-1, 0}};
    rotation_matrices = {R, R2, R3, R4, R5};
    }

    else if(latticeTypeName=="Square"){
    vector <vector <int>> R = {{0,-1},{1, 0}};
    vector <vector <int>> R2 = {{-1,0},{0, -1}};
    vector <vector <int>> R3 = {{0,1},{-1, 0}};
    rotation_matrices = {R, R2, R3};
    }

    cout<<rotation_matrices.size()<<" rotation matrices built"<<endl;
    

    for (int level_max=1 ; level_max < cluster_probabilities.size()+1 ; level_max++){
        // Get the list of positions of the neighbor relative to center for a plaquette of size level_max
        //set<Site*> plaquette_sites  ;  
        vector<vector<int>> relative_positions;
        vector <int> center_position = fullSites[0]->get_real_space_position(NhashPosition);
        //siteList[0]->get_plaquette_relative_positions(0, level_max, NhashPosition, plaquette_sites, relative_positions, center_position, Lx, Ly, Lz);
        fullSites[0]->get_plaquette_relative_positions2(level_max, NhashPosition, relative_positions, center_position, Lx, Ly, Lz);


        int n_sites = relative_positions.size();
        //cout<<"plaquette radius "<<level_max<<" "<<n_sites<<" sites in plaquette"<<endl;

        // Find the correspondance between indices for each rotation
        vector <vector <int>> this_size_index_correspondance ; 

        for (int rotation_id=0; rotation_id< rotation_matrices.size() ; rotation_id ++){
            vector<int> this_size_this_rotation_index_correspondance(n_sites, 0); 
            for (int k=0 ; k<n_sites; k++){
                vector <int> pos_initiale = relative_positions[k];
                vector <int> new_pos = matrix_product(rotation_matrices[rotation_id], pos_initiale) ;
                
                // Implement periodic boundary condition 
                for (int d=0 ; d<3 ; d++){
                    new_pos[d] = (new_pos[d]+2*sizes[d])%sizes[d];
                }
                /* if (level_max==2){
                cout << "pos init "<<pos_initiale[0] <<" "<<pos_initiale[1];
                cout << "new pos "<<new_pos[0] <<" "<<new_pos[1]<<endl;
                } */

                // find the position after rotation
                bool found_corresp = false;
                for (int kprime=0 ; kprime< n_sites; kprime++){
                    // find the indices to which it corresponds 
                    if (relative_positions[kprime][0] == new_pos[0] && relative_positions[kprime][1] == new_pos[1] && relative_positions[kprime][2] == new_pos[2]){
                        this_size_this_rotation_index_correspondance[kprime] = k;
                        found_corresp = true;
                        //cout <<"k "<<k<<" kprime "<<kprime<<endl;
                    }
                    
                }
                if (!found_corresp){
                        cout<<"NOT FOUND ";
                        cout << "pos init "<<pos_initiale[0] <<" "<<pos_initiale[1];
                        cout << "new pos "<<new_pos[0] <<" "<<new_pos[1]<<endl;
                    } 
            }
            /* if (level_max==2){
            cout<<"this_size_this_rotation_index_correspondance"<<endl;
            for (int kprime2=0; kprime2<n_sites; kprime2++){
                cout<<this_size_this_rotation_index_correspondance[kprime2]<<" ";
            }
            cout<<endl;
            } */
            this_size_index_correspondance.push_back(this_size_this_rotation_index_correspondance);
        }
        index_correspondance_plaquette.push_back(this_size_index_correspondance);
    }   



    

}


void System::build_plaquettes(){

    for (int radius=1; radius<max_cluster_radius ; radius++){
        plaquettes_list.push_back(vector<Plaquette*>());
        compatible_plaquettes.push_back(vector<pair<Plaquette*, Plaquette*>>());
        cout<<"radius "<< radius;
        for (int i=0 ; i<Nsites ; i++){
            Plaquette* current_plaquette = new Plaquette(siteList[i], radius);
            
            plaquettes_list[radius-1].push_back(current_plaquette);
        }
        cout<<" ok ";
        for (int ip1=0 ; ip1<plaquettes_list[radius-1].size() ;  ip1++){
            for (int ip2=0 ; ip2<ip1 ;  ip2++){
                if (plaquettes_list[radius-1][ip1]->intersectionEmpty(plaquettes_list[radius-1][ip2])){
                    compatible_plaquettes[radius-1].push_back({plaquettes_list[radius-1][ip1], plaquettes_list[radius-1][ip2]});
                }
            }
        }
        cout<<compatible_plaquettes[radius-1].size()<<" pairs of compatible plaquettes"<<endl;
    }
    cout<<"We built ";
    for (int radius=1; radius<max_cluster_radius ; radius++){
        cout<<plaquettes_list[radius-1].size()<<" ";
    }
    cout<<"plaquettes" <<endl;
}

void System::anneal_system(std::vector<float> Betas){
    
    vector <float> :: iterator Beta ;
    pair <bool, float> successAndDeltaE;

    int indexOfTemperature = 0 ;
    //int Time = 0;
    
    int NTemperatures = Betas.size();
    
    int nextIMes = 0 ;
//    vector <int> allImes ;
    int freqOfMeasure =0 ;

     if (Nmes_per_T!=0){
        freqOfMeasure = (Nsteps_per_T)/(Nmes_per_T);
        
         if (freqOfMeasure==1){
             nextIMes=0;
         }
         else{
        nextIMes = Nsteps_per_T-1 - (Nmes_per_T-1)*(freqOfMeasure);
         }
     }
    float averageEtot = Etot;
    int MCsuccess = 0;
    int n_step_saveE= 1;
    int n_step_saveS= 1;

    {
    LOG_DURATION("Annealing");
    for (Beta =Betas.begin(); Beta!=Betas.end(); Beta++){
        //cout<<"total energy at T="<<(*T)<<" is "<<Etot<<endl;
        LOG_DURATION("Beta="+to_string(*Beta)+ ",Etot"+to_string(Etot));
//        cout <<"T"<<*T<<" do "<<Nsteps_per_T<<" steps"<<endl;
        

        for (int indexAtT = 0 ; indexAtT<Nsteps_per_T ;  indexAtT++){
            
            //cout<<"i "<<indexAtT<<" begins ";
            successAndDeltaE = MC_step(*Beta, false);
            //Time++;
            

            if (saveEnergies){
                averageEtot+=compute_Etot();
                n_step_saveE++;
                if (indexAtT%int(Nsteps_per_T/1000)==0){
                    //cout<<"averageEtot "<<averageEtot<<" n "<<n_step_saveE<<endl;
                    AllEnergies.push_back(averageEtot/n_step_saveE);
                    n_step_saveE=1;
                    averageEtot=Etot;
                }
            }
            
            if (saveMCsuccess){
                MCsuccess += int(successAndDeltaE.first);
                n_step_saveS++;
                if (indexAtT%int(Nsteps_per_T/100)==0){
                    AllMCsuccess.push_back(float(MCsuccess)/n_step_saveS);
                    //cout<<"success "<<float(MCsuccess)/n_step_saveS<<endl;
                    n_step_saveS=1;
                    MCsuccess = int(successAndDeltaE.first);;
                }
            }
            
            if (saveSnapshot){
                if (*Beta < Betas[NTemperatures-2] and indexAtT > std::max( int(Nsteps_per_T*0.9), Nsteps_per_T-20) ){
                    try_save_ground_state();
                }
            }
            
            if (Nmes_per_T !=0) {
                if (indexAtT == nextIMes || indexAtT == Nsteps_per_T-1){
                    save_system_state(indexOfTemperature, indexAtT);
//                    if (indexOfTemperature==0){
//                        allImes.push_back(nextIMes);}
                    nextIMes += freqOfMeasure ;
                }}
            //cout <<" i"<<indexAtT<<" finished"<<endl;
        }
        
        indexOfTemperature ++;
        if (Nmes_per_T !=0) {
            if (freqOfMeasure==1){
                nextIMes=0;
            }
            else{
            nextIMes = Nsteps_per_T-1 - (Nmes_per_T-1)*(freqOfMeasure);;
            }
        }
    }
    }
    print_energy_per_particle();
    if (saveTemperatures){
    vectorToFile(Betas, directoryForResults+"_Temperatures.txt");
    }
}

void System::quench_system(float Temperature_quenching, int Nsteps_quenching){
    LOG_DURATION("Quenching");
    pair <bool, float> successAndDeltaE;
    int nextIMes = 0 ;
    int freqOfMeasure =0 ;
     if (Nmes_per_T!=0){
        freqOfMeasure = (Nsteps_per_T)/(Nmes_per_T);
        nextIMes = freqOfMeasure;
     }
    
    for (int index=0 ; index<Nsteps_quenching ; index++){
        successAndDeltaE = MC_step(Temperature_quenching, false);
        if (saveSnapshot){
            try_save_ground_state();
        }
        
         if (Nmes_per_T !=0) {
            if (index == nextIMes || index == Nsteps_per_T-1){
                save_system_state(Ntemperatures+1, index);
                nextIMes += freqOfMeasure ;
                }}
        
    }
}

void System::save_equilibration(){
    
    if (saveSnapshot){
        if (not oneGroundStateSaved || Etot<Emin){
            // we stave the ground state state if it was not done before, or if the energy reached now is the lowest measured
            save_system_state(-1, -1);
            Emin = Etot ;
        }
//    vectorToFile(allImes, directoryForResults+"_MeasurementIndex.txt");
    }
    if (saveEnergies){
    vectorToFile(AllEnergies, directoryForResults+"_Energies.txt");
    }
    if (saveMCsuccess){
    vectorToFile(AllMCsuccess, directoryForResults+"_MCsuccess.txt");
    }
}



double System::get_Etot(){
    return Etot ;
}

double System::get_Emin(){
    return Emin ;
}

int System::get_MC_steps_per_Temperature(){
    
    return Nsteps_per_T ;
}

int System::get_Nsites(){
    return Nsites ;
}

//float System::get_final_temperature(){
//    return Temperatures[NTemperatures-1];
//}

bool System::get_saveAnnealing(){
    return saveSnapshot ;
}


void System::try_save_ground_state(){
    if (Etot<Emin){
        save_system_state(-1, -1);
        oneGroundStateSaved = true ;
        Emin = Etot ;
    }
}

void System::save_system_state(int indexOfTemperature, int current_step){
    
   
    string  directory1, directory2;
    if (indexOfTemperature!=-1 and current_step!=-1){
        // We save the system state at a given temperature and time step
        string indexT = to_string(indexOfTemperature);
        string indexS = to_string(current_step);
        directory1 = directoryForResults+"adressFulls_T" +indexT + "_step"+indexS+".txt";
        directory2 = directoryForResults+"orientationFulls_T" +indexT + "_step"+indexS+".txt";
    }
    
    else {
        // We save this system state because it is the lowest measured energy so far -> the ground state
        directory1 = directoryForResults+"adressFulls_groundState.txt";
        directory2 = directoryForResults+"orientationFulls_groundState.txt";
    }
    
    ofstream myStream1 (directory1.c_str());
    ofstream myStream2 (directory2.c_str());
    vector <Site*> :: iterator i;
    for (i=fullSites.begin(); i!=fullSites.end() ; i++){
        myStream1<<(*i)<<" ";
        myStream2<<(*i)->getOrientationRef()<<" ";
    }
    myStream1<<"#END";
    myStream1.flush();
    myStream1.close();
    myStream1.clear();
    myStream2<<"#END";
    myStream2.flush();
    myStream2.close();
    myStream2.clear();
}


vector <int>  System::scaffold_stats(vector <Scaffold*> aScaffoldList){
    vector <Scaffold*> :: iterator it;
    vector <int> Nvector = vector <int> (N_Structures,0);
    for (it=aScaffoldList.begin(); it!=aScaffoldList.end(); it++){
        Nvector[(*it)->getStructure()]++;

      }
    return Nvector ;
}






void System::do_statistics(float Temperature_statistics, int Naverage, bool saveCovariance, bool saveExtendedScaffold){
    //Creates covariance matrix and statistics vectors c and d, average over Naverage steps of MC simulation
    //the covariance matrix is symetric, and multiplied by Nsites to make it independent of the system size
    
    
    
    cVector = scaffold_stats(scaffoldList);
    dVector = scaffold_stats(extendedScaffoldList);
    
    
    LOG_DURATION("Statistics during "+to_string(Naverage)+" MC steps at T="+to_string(Temperature_statistics)+"kT");
 
    vector <long int> currentNalphas (N_Structures, 0), Nalphas (N_Structures, 0), currentNalphasPrime (N_Structures, 0), NalphasPrime (N_Structures, 0);
    
    vector<vector<long int>> secondMoment ;
    if (saveCovariance){
    //Initialize the secondMoment matrix
    for (int alpha = 0 ; alpha <N_Structures ; alpha++){
        secondMoment.push_back(vector<long int>(N_Structures, 0));
    }
    }
    
    
    //extract the values and compute their sum
    for (int i=0; i<Naverage ; i++){
        MC_step(Temperature_statistics, true);
        
        vector<int> scaffoldStats = cVector; //should create a copy without changing cVector
        std::transform(scaffoldStats.begin(), scaffoldStats.end(), currentNalphas.begin(), [](int x) { return (long int)x;});
        
        vector<int> scaffoldStatsExt = dVector;
        std::transform(scaffoldStatsExt.begin(), scaffoldStatsExt.end(), currentNalphasPrime.begin(), [](int x) { return (long int)x;});
        
        

        for (int alpha = 0 ; alpha <N_Structures ; alpha++){
            Nalphas[alpha]+= currentNalphas[alpha];
            NalphasPrime[alpha]+= currentNalphasPrime[alpha];
            if (saveCovariance){
            for (int beta = alpha ; beta <N_Structures ; beta++){
                secondMoment[alpha][beta] += currentNalphas[alpha]*currentNalphas[beta];
            }
            }
        }
    }
    
    if (saveCovariance){
    //substract the square of the means to the mean of the squares
    for (int alpha = 0 ; alpha <N_Structures ; alpha++){
        for (int beta = alpha ; beta <N_Structures ; beta++){
            if (alpha>N_Orientation){
                long int temp = secondMoment[alpha][beta];
            secondMoment[alpha][beta] = temp*(long int)Naverage - Nalphas[alpha]*Nalphas[beta];

            }
        }
    }
    }
        
    
    //Normalize
    vector <double> cVector (N_Structures, 0);
    vector <double> dVector (N_Structures, 0);
    
    if (saveCovariance){
    vector <vector<double>> Cmatrix ;
    }
        
    int N1particleScaffold = Nsites ;
    int N2particlesScaffold = (int) (scaffoldList.size() - N1particleScaffold) ; // number of particles * number of bonds
        
    for (int alpha = 0 ; alpha <N_Structures ; alpha++){
        if (alpha <= N_Orientation){
            cVector[alpha] = (double)Nalphas[alpha] / (double)N1particleScaffold / (double) Naverage;
            dVector[alpha] = (double)NalphasPrime[alpha] / (double)N1particleScaffold / (double) Naverage;
            
        }
        else{
            cVector[alpha] = (double)Nalphas[alpha] / (double)N2particlesScaffold / (double) Naverage;
            dVector[alpha] = (double)NalphasPrime[alpha] / (double)N2particlesScaffold / (double) Naverage;
        }
    }
    
        
    vectorToFile(cVector, directoryForResults+"_c.txt");
    if (saveExtendedScaffold){
    vectorToFile(dVector, directoryForResults+"_d.txt");
    }
    if (saveCovariance){
    matrixToFile(secondMoment, directoryForResults+"_Cab.txt");
    }
    vectorToFile(typeOfThisSystem.get_energies(), directoryForResults+"_LELchosen.txt");
    print_energy_per_particle();
       
}











void System::test_site_swap(){
    this->build_lattice(Lx, Ly, Lz);
    this->build_neighbour_relations(0, 1,scaffoldList);
    this->fill_the_system();
    
    int temp_n;
    temp_n= (int) fullSites.size();
    cout <<endl<<"We have "<<temp_n<<" sites in fullSites"<<endl;
    
    Site* aFullSite = fullSites[temp_n-1];
    int pos1 = RealPositionMappedToSiteList.at(aFullSite->get_hashedPosition());
    vector <int> RSpos = aFullSite->get_real_space_position(NhashPosition);
    cout<< "Site 1 : "<<endl;
    (*siteList[pos1]).printAttributes(NhashPosition);
    cout << "Structure of the belonging scaffolds :"<<endl;
    printVectorOfAdress((*siteList[pos1]).getBelongingScaffolds(0));
    
    int temp_n2;
    temp_n2= (int) fullSites.size();
    Site* anEmptySite = emptySites[temp_n2-1];
    int pos2 = RealPositionMappedToSiteList.at(anEmptySite->get_hashedPosition());
     vector <int> RSpos2 = aFullSite->get_real_space_position(NhashPosition);
    cout<< "Site 2 : "<<endl;
     (*siteList[pos2]).printAttributes(NhashPosition);
     cout << "Structure of the belonging scaffolds :"<<endl;
    
    pair <bool,float> SuccessAndDeltaE = aFullSite->attemptSwapTo(*anEmptySite, 1, false);
    
    cout << "Swap ? "<< (SuccessAndDeltaE.first)<<endl;
    cout<< "Site 1 : "<<endl;
    (*siteList[pos1]).printAttributes(NhashPosition);
    cout <<"Site 2 : "<<endl;
    (*siteList[pos2]).printAttributes(NhashPosition);
    
   
    cout <<"Common neighbours of this two sites"<<endl;
    printVectorOfAdress((*siteList[pos1]).commonNeigbhours(*siteList[pos2],0));
    cout <<"Common scaffold of this two sites"<<endl;
    vector <Scaffold*> commonScaf = (*siteList[pos1]).commonScaffold(*siteList[pos2],0);
    printVectorOfAdress(commonScaf);
    if (commonScaf.size()>0){
        Scaffold* theCommonScaf = commonScaf[0];
        theCommonScaf->printAttributes();
    }
}


void System::test_site_flip(){
    build_lattice(Lx, Ly, Lz);
    this->build_neighbour_relations(0, 1,scaffoldList);
    this->fill_the_system();

    int temp_n;
    temp_n= (int) fullSites.size();
    cout <<endl<<"We have "<<temp_n<<" sites in fullSites"<<endl;

    Site* aFullSite = fullSites[temp_n-1];
    int pos1 = RealPositionMappedToSiteList.at(aFullSite->get_hashedPosition());
    vector <int> RSpos = aFullSite->get_real_space_position(NhashPosition);
    cout<< "Site 1 : "<<endl;
    (*siteList[pos1]).printAttributes(NhashPosition);
    cout << "Structure of the belonging scaffolds :"<<endl;
    printVectorOfAdress((*siteList[pos1]).getBelongingScaffolds(0));
    
    pair <bool, float> successAndDeltaE = aFullSite->attemptFlip(N_Orientation, 100);
    cout << "Delta Energy ="<<successAndDeltaE.second<<endl;
    
    (*siteList[pos1]).printAttributes(NhashPosition);
    cout << "Structure of the belonging scaffolds :"<<endl;
    printVectorOfAdress((*siteList[pos1]).getBelongingScaffolds(0));
    
}
    

void extract_value(string* token, string const expectation){
    size_t found = (*token).find(expectation);
    if (found!=std::string::npos){
        (*token).erase(0, expectation.size());
        //(*token).pop_back();
    }
}

int get_int_in_line(string token, string  const expectation){
    int result ;
    extract_value(&token, expectation);
    std::stringstream stream(token);
    if (!(stream >> result)){
        cout << "Problem of conversion reading " <<expectation <<endl;
    }
    
    return result ;
    }
        


int * build_2particles_scaffold(string scaffold_filename){
    ifstream file(scaffold_filename.c_str());
    if (file){
        string token;
        getline(file, token);
        int N_Orientations (get_int_in_line(token, "N_Orientations="));
        
        string token2;
        getline(file, token2);
        int N_Edges (get_int_in_line(token2, "N_Edges="));
        
        string token3;
        getline(file, token3);
        int N_Structures (get_int_in_line(token3, "N_Structures="));
        
        cout << N_Orientations << N_Edges << N_Structures <<endl; 
}
    else {
    cout << "Couldn't read file"<<endl;
    }
    return 0;
}

