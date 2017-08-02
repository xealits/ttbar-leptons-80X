#ifndef PROCESSINGGENPARTICLES_H
#define PROCESSINGGENPARTICLES_H

#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <string>

const reco::Candidate* find_W_decay(const reco::Candidate * W);
std::string parse_W_decay(const reco::Candidate * W);
unsigned int simple_tau_decay_id(const reco::Candidate * tau);
double top_pT_SF(double x);


#endif /* PROCESSINGGENPARTICLES_H */

