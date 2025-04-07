#ifndef COMBINESIMILARITYENERGY
#define COMBINESIMILARITYENERGY

#include "foldalign.hxx"

inline scoreType combineSimilarityEnergy(const scoreType& similarityScore,
										 const scoreType& energyScore) {
	if (similarityScore == big_neg || energyScore == big_neg) {
		return big_neg;
	}

	return similarityScore + energyScore;
}

#endif /* COMBINESIMILARITYENERGY */
