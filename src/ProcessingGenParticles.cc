#include "UserCode/ttbar-leptons-80X/interface/ProcessingGenParticles.h"

using namespace std;

/*
 * const reco::Candidate* find_W_decay(const reco::Candidate * W) {
 *
 * returns the final W, daughters of which are decay products
 * or returns NULL if something weird happend
 *
 * Usage:
 * const reco::Candidate * W = p.daughter( W_num );
 * const reco::Candidate * W_final = find_W_decay(W);
 * W_final->daughter(0)->pdgId();
 * W_final->daughter(1)->pdgId();
*/
const reco::Candidate* find_W_decay(const reco::Candidate * W) {
	const reco::Candidate * p = W; // eeh.. how C passes const arguments? are they local to function
	int d0_id, d1_id; // ids of decay daughters
	//bool found_decay=false;
	// follow W, whatever happens to it (who knows! defensive programming here)
	// until leptonic/quarkonic decay is found
	// 
	// assume there can be only W->W transitions inbetween (TODO: to check actually)
	//while (!found_decay) {
	while (true) {
		int n = p->numberOfDaughters();
		switch(n) {
		case 0: return NULL;
		case 1: // it should be another W->W transition
			p = p->daughter(0);
			break; // there should not be no infinite loop here!
		case 2: // so, it should be the decay
			//found_decay = true;
			return p;
			/* all this is derived from p
			d0_id = fabs(p->daughter(0)->pdgId());
			d1_id = fabs(p->daughter(1)->pdgId());
			if (d0_id == 15 || d1_id == 15 ) return string("tau");
			if (d0_id == 11 || d1_id == 11 ) return string("el");
			if (d0_id == 13 || d1_id == 13 ) return string("mu");
			return string("q"); // FiXME: quite dangerous control-flow!
			*/
		default: // and this is just crazy
			return NULL;
		}
	}
}

string parse_W_decay(const reco::Candidate * W) {
	int d0_id = fabs(W->daughter(0)->pdgId());
	int d1_id = fabs(W->daughter(1)->pdgId());
	if (d0_id == 15 || d1_id == 15 ) return string("tau");
	if (d0_id == 11 || d1_id == 11 ) return string("el");
	if (d0_id == 13 || d1_id == 13 ) return string("mu");
	return string("q"); // FiXME: quite dangerous control-flow!
}
