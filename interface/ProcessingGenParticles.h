#ifndef PROCESSINGGENPARTICLES_H
#define PROCESSINGGENPARTICLES_H

#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <string>

const reco::Candidate* find_W_decay(const reco::Candidate * W);
std::string parse_W_decay(const reco::Candidate * W);
int simple_tau_decay_id(const reco::Candidate * tau);


#endif /* PROCESSINGGENPARTICLES_H */

