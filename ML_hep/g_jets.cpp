#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<ctime>

#include "Pythia8Plugins/FastJet3.h"
#include "Pythia8/Pythia.h"

// Declaration of namespaces
using namespace Pythia8;
using namespace std;
using namespace fastjet;


// to compile the program: g++ -I /home/kamil/Desktop/pythia8310/include pp_gluon_event.cpp -o pp_gluon_event -lpythia8 -L/home/kamil/Desktop/pythia8310/lib
// to access the library path: export LD_LIBRARY_PATH=~/Desktop/pythia8310/lib:$LD_LIBRARY_PATH
// also for fastjet library: export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

ofstream fdata;

int main(){
    // Here we declare all the variables that we are going to use
    int nevents = 10000; // number of events

    // Parameters used in the fastjet reconstruction
    double R = 0.4; // jet size
    double pTMin = 5.0; // minimum pT of the jet
    double pTMax = 200; // maximum ..
    double particle_eta, particle_phi, particle_pt, particle_id;
    double center_phi, center_eta;
    valarray<double> fourvec(4);
        
    // Declaration of pythia object and shorthand for event
    Pythia pythia;
    Event& event = pythia.event;

    // Read random seed:
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 72649");

    // characteristichs of each event:
    //processes, CM energy, pT energy trashold
    pythia.readString("HardQCD:gg2gg = on");
    pythia.readString("HardQCD:qqbar2gg = on");

    pythia.readString("Beams:eCM = 13000"); 
    pythia.readString("PhaseSpace:pTHatMin = 200.");

    // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    
    // Initialize pythia (?)
    pythia.init();

    // Set up the FastJet jet finder
    JetAlgorithm algorithm;
    algorithm = antikt_algorithm;
    JetDefinition jetDef(algorithm, R);
    vector<PseudoJet> fjImputs;

    int jet_count = 0;
    fdata.open("jets_g.dat");

    // Loop to generate nevents, skip if error occures
    for(int iEvent = 0; iEvent < nevents; iEvent++){
        if (!pythia.next()) continue; // In case of errors
        cout << " event " << iEvent << " : " << endl;
        // Begin FastJet analisys
        fjImputs.resize(0);
        for (int i = 0; i < 4; i++) fourvec[i] = 0; // Reset 
        Vec4    pTemp;
        int     nAnalyze = 0;

        // Loop over each particle in the event, 
        //if the particle is a final state particle add it to PseudoJet
        for(int i = 0; i < event.size(); i++) if (event[i].isFinal()) {

            //if(event[i].isNeutral()) continue; // Select charged particle
            if(!event[i].isVisible()) continue;// Select visible particle
            if(abs(event[i].eta()) > 2.5) continue; // Select |eta|<2.5
           
            for (int i = 0; i < 4; i++) fourvec[i] = 0; // Reset px, py, pz, e
            
            // Create a pseudojet for the pythia particle
            fourvec[0] = event[i].px();
            fourvec[1] = event[i].py();
            fourvec[2] = event[i].pz();
            fourvec[3] = event[i].e();
            PseudoJet particleTemp(fourvec);

            // Store acceptable particle as imput to fastjet
            fjImputs.push_back(particleTemp);
            nAnalyze++;
        }

        // Run FastJet algorithm and sort jets in pT order
        vector<PseudoJet> inclusiveJets, sortedJets;
        ClusterSequence clustSeq(fjImputs, jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMin);
        sortedJets    = sorted_by_pt(inclusiveJets);
        
        // loop through each jet and extract information about eta and phi of single particles
        for(unsigned int i = 0; i < sortedJets.size(); i++){
            // Print out the selected jet caracteristichs
            cout << "pT = " << sortedJets[i].pt() << " eta = " << sortedJets[i].eta();
            cout << " phi = " << sortedJets[i].phi() << endl;
            // If pt in 100, 150 Gev, write the constituents into a file
            if(sortedJets[i].pt() >= 200.0 && sortedJets[i].pt() <= 220.0){
                vector<PseudoJet> constituents = sortedJets[i].constituents();
                cout << "Jet selected for imaging"<< endl;
                
                jet_count++;
                // Centroid of the histogram
                double sum_phi = 0; // weighted sum with pt
                double sum_eta = 0;
                double sum_pt  = 0;
                for(unsigned int k = 0; k < constituents.size(); k++){
                    sum_phi += constituents[k].phi() * constituents[k].pt();
                    sum_eta += constituents[k].eta() * constituents[k].pt();
                    sum_pt  += constituents[k].pt();
                }

                double centroid_eta = sum_eta / sum_pt;
                double centroid_phi = sum_phi / sum_pt;

                // Write file for histogram, centering image centroid in (eta, phi) = (0, 0)
                for( unsigned int j = 0; j < constituents.size(); j++){
                    particle_eta = constituents[j].eta();
                    particle_phi = constituents[j].phi();
                    particle_pt  = constituents[j].pt();

                    // Eliminate particles too distant from center of the jet
                    if (fabs(particle_eta - sortedJets[i].eta()) > 0.5) continue;
                    if (fabs(particle_phi - sortedJets[i].phi()) > 0.5) continue;
                    
                    // Centering image in the centroid
                    particle_eta -= centroid_eta;
                    particle_phi -= centroid_phi;

                    // Print eta, phi, pt of each particle of the jet in a file
                    fdata << jet_count << " " << particle_eta << " " << particle_phi;
                    fdata << " " << particle_pt << endl;
                }

                // Separate one jet from the other
                fdata << endl;            
            }
        }
    }
    fdata.close();
    return 0;
}