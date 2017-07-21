#ifndef PROCESSINGGENPARTICLES_H
#define PROCESSINGGENPARTICLES_H

#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <string>

const reco::Candidate* find_W_decay(const reco::Candidate * W);
std::string parse_W_decay(const reco::Candidate * W);


#endif /* PROCESSINGGENPARTICLES_H */

