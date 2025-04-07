#ifndef PVALUE
#define PVALUE

#include "nohit.hxx"
#include "datalist.cxx"
#include "alignmentCells.cxx"
#include <math.h>

class pValue {

public:
	pValue(const lengthType length1, const lengthType length2, 
		const double gcContent1, const double gcContent2)
		:	len1(length1), len2(length2),
			gc1(gcContent1), gc2(gcContent2),
			lambda(setLambda()), k(setK()) {}
	
	double calcPvalue(const scoreType score) const {
		double p = 1 - exp(-1*k*len1*len2*exp(-1*lambda*double(score)));
		return p;
	};
	
	double getLambda() const {return lambda;}
	double getK() const {return k;}
	
	// Assumes that all alignments in the hitlist are random
	void estimateParameters(datalist<alignmentcells>& hitlist, scoreType scoreCutOff);
private:
	const lengthType len1;
	const lengthType len2;
	const double gc1;
	const double gc2;

	double lambda;
	double k;
	
	pValue();
	double setLambda();
	double setK();
	void calcNewLambda(long long sumOfScores, scoreType scoreCutoff, long numberOfScores);
	void calcNewK(long numberOfScores, scoreType scoreCutoff);
};

inline void pValue::estimateParameters(datalist<alignmentcells>& hitlist, scoreType scoreCutOff) {

	scoreType sa;
	alignmentcells a;
	hitlist.initCheckNext();

	hitlist.checkNext(sa, a);
	scoreType badScore = scoreCutOff;
	if (badScore < bad) {
		badScore = bad;
	}
	
	long numberOfScores = 0;
	long long sumOfScores = 0;
	while (sa > badScore) {
		sumOfScores += sa;
		numberOfScores++;
		
		hitlist.checkNext(sa, a);
	}
	calcNewLambda(sumOfScores, scoreCutOff, numberOfScores);
	calcNewK(numberOfScores, scoreCutOff);
}

inline void pValue::calcNewLambda(long long sumOfScores, scoreType scoreCutoff, long numberOfScores) {

	double averageAboveCut = double(sumOfScores)/double(numberOfScores) - double(scoreCutoff);
	
	if (averageAboveCut != 0) {
		lambda = log(1 + 1/averageAboveCut);
	}
	else {
		std::cerr << "The lambda value could not be estimated" << std::endl;
		lambda = 1000;
	}
}

inline void pValue::calcNewK(long numberOfScores, scoreType scoreCutoff) {

	if (lambda == 1000) {
		k = 1000;
	}
	
	k = numberOfScores*exp(lambda*scoreCutoff)/(len1*len2);
}
	

inline double pValue::setLambda() {

//return 0.0167201347394947;

	const double gcMax = gc1 > gc2 ? gc1 : gc2;
	const double gcMin = gc1 > gc2 ? gc2 : gc1;

        if (gcMax < 0.20) {
                return 0.0147;
        } else if (gcMax < 0.25) {
                if (gcMin < 0.20) {
                        return 0.0171;
                } else {
                        return 0.01939;
                }
        } else if (gcMax < 0.30) {
                if (gcMin < 0.20) {
                        return 0.0193;
                } else if (gcMin < 0.25) {
                        return 0.0211;
                } else {
                        return 0.0226;
                }
        } else if (gcMax < 0.35) {
                if (gcMin < 0.20) {
                        return 0.0208;
                } else if (gcMin < 0.25) {
                        return 0.0224;
                } else if (gcMin < 0.30) {
                        return 0.0235;
                } else {
                        return 0.0236;
                }
        } else if (gcMax < 0.40) {
                if (gcMin < 0.20) {
                        return 0.0220;
                } else if (gcMin < 0.25) {
                        return 0.0232;
                } else if (gcMin < 0.30) {
                        return 0.0235;
                } else if (gcMin < 0.35) {
                        return 0.0232;
                } else {
                        return 0.0223;
                }
        } else if (gcMax < 0.45) {
                if (gcMin < 0.20) {
                        return 0.0226;
                } else if (gcMin < 0.25) {
                        return 0.0233;
                } else if (gcMin < 0.30) {
                        return 0.0230;
                } else if (gcMin < 0.35) {
                        return 0.0223;
                } else if (gcMin < 0.40) {
                        return 0.0210;
                } else {
                        return 0.019398;
                }
        } else if (gcMax < 0.50) {
                if (gcMin < 0.20) {
                        return 0.0222;
                } else if (gcMin < 0.25) {
                        return 0.0226;
                } else if (gcMin < 0.30) {
                        return 0.0220;
                } else if (gcMin < 0.35) {
                        return 0.0209;
                } else if (gcMin < 0.40) {
                        return 0.0192;
                } else if (gcMin < 0.45) {
                        return 0.0175;
                } else {
                        return 0.0159;
                }
        } else if (gcMax < 0.55) {
                if (gcMin < 0.20) {
                        return 0.0217;
                } else if (gcMin < 0.25) {
                        return 0.0213;
                } else if (gcMin < 0.30) {
                        return 0.0205;
                } else if (gcMin < 0.35) {
                        return 0.0191;
                } else if (gcMin < 0.40) {
                        return 0.0174;
                } else if (gcMin < 0.45) {
                        return 0.0157;
                } else if (gcMin < 0.50) {
                        return 0.0137;
                } else {
                        return 0.010;
                }
        } else if (gcMax < 0.60) {
                if (gcMin < 0.20) {
                        return 0.0201;
                } else if (gcMin < 0.25) {
                        return 0.0197;
                } else if (gcMin < 0.30) {
                        return 0.0186;
                } else if (gcMin < 0.35) {
                        return 0.0171;
                } else if (gcMin < 0.40) {
                        return 0.0154;
                } else if (gcMin < 0.45) {
                        return 0.01340;
                } else if (gcMin < 0.50) {
                        return 0.0105;
                } else if (gcMin < 0.55) {
                        return 0.006;
                } else {
                        return 0.0044;
                }
        } else if (gcMax < 0.65) {
                if (gcMin < 0.20) {
                        return 0.0185;
                } else if (gcMin < 0.25) {
                        return 0.0180;
                } else if (gcMin < 0.30) {
                        return 0.0166;
                } else if (gcMin < 0.35) {
                        return 0.0152;
                } else if (gcMin < 0.40) {
                        return 0.01316;
                } else if (gcMin < 0.45) {
                        return 0.010;
                } else if (gcMin < 0.50) {
                        return 0.006;
                } else if (gcMin < 0.55) {
                        return 0.0041;
                } else if (gcMin < 0.60) {
                        return 0.0031;
                } else {
                        return 0.0020;
                }
        } else {
                if (gcMin < 0.20) {
                        return 0.0175;
                } else if (gcMin < 0.25) {
                        return 0.0166;
                } else if (gcMin < 0.30) {
                        return 0.0156;
                } else if (gcMin < 0.35) {
                        return 0.0141;
                } else if (gcMin < 0.40) {
                        return 0.01169;
                } else if (gcMin < 0.45) {
                        return 0.008;
                } else if (gcMin < 0.50) {
                        return 0.0044;
                } else if (gcMin < 0.55) {
                        return 0.0040;
                } else if (gcMin < 0.60) {
                        return 0.0023;
                } else if (gcMin < 0.65) {
                        return 0.00168;
                } else {
                        return 0.0018674140;
                }
        }

}
	

inline double pValue::setK() {

//return 0.000198574695044439;

	const double gcMax = gc1 > gc2 ? gc1 : gc2;
	const double gcMin = gc1 > gc2 ? gc2 : gc1;
	
        if (gcMax < 0.20) {
                return 0.000333;
        } else if (gcMax < 0.25) {
                if (gcMin < 0.20) {
                        return 0.00059;
                } else {
                        return 0.00095;
                }
        } else if (gcMax < 0.30) {
                if (gcMin < 0.20) {
                        return 0.00092;
                } else if (gcMin < 0.25) {
                        return 0.00121;
                } else {
                        return 0.001447;
                }
        } else if (gcMax < 0.35) {
                if (gcMin < 0.20) {
                        return 0.00114;
                } else if (gcMin < 0.25) {
                        return 0.001397;
                } else if (gcMin < 0.30) {
                        return 0.00155;
                } else {
                        return 0.00159;
                }
        } else if (gcMax < 0.40) {
                if (gcMin < 0.20) {
                        return 0.001286;
                } else if (gcMin < 0.25) {
                        return 0.00149;
                } else if (gcMin < 0.30) {
                        return 0.001553;
                } else if (gcMin < 0.35) {
                        return 0.00156;
                } else {
                        return 0.00146;
                }
        } else if (gcMax < 0.45) {
                if (gcMin < 0.20) {
                        return 0.001330;
                } else if (gcMin < 0.25) {
                        return 0.001458;
                } else if (gcMin < 0.30) {
                        return 0.00147837;
                } else if (gcMin < 0.35) {
                        return 0.00145;
                } else if (gcMin < 0.40) {
                        return 0.001330;
                } else {
                        return 0.00114;
                }
        } else if (gcMax < 0.50) {
                if (gcMin < 0.20) {
                        return 0.00126079;
                } else if (gcMin < 0.25) {
                        return 0.001367;
                } else if (gcMin < 0.30) {
                        return 0.001350;
                } else if (gcMin < 0.35) {
                        return 0.001268;
                } else if (gcMin < 0.40) {
                        return 0.0011006;
                } else if (gcMin < 0.45) {
                        return 0.00088;
                } else {
                        return 0.00065;
                }
        } else if (gcMax < 0.55) {
                if (gcMin < 0.20) {
                        return 0.001166;
                } else if (gcMin < 0.25) {
                        return 0.00118;
                } else if (gcMin < 0.30) {
                        return 0.00115;
                } else if (gcMin < 0.35) {
                        return 0.00103;
                } else if (gcMin < 0.40) {
                        return 0.00084;
                } else if (gcMin < 0.45) {
                        return 0.00062;
                } else if (gcMin < 0.50) {
                        return 0.00036;
                } else {
                        return 0.000111;
                }
        } else if (gcMax < 0.60) {
                if (gcMin < 0.20) {
                        return 0.00097;
                } else if (gcMin < 0.25) {
                        return 0.00098;
                } else if (gcMin < 0.30) {
                        return 0.000923;
                } else if (gcMin < 0.35) {
                        return 0.000762;
                } else if (gcMin < 0.40) {
                        return 0.0005426;
                } else if (gcMin < 0.45) {
                        return 0.00033;
                } else if (gcMin < 0.50) {
                        return 0.00015;
                } else if (gcMin < 0.55) {
                        return 0.00002;
                } else {
                        return 0.0000018;
                }
        } else if (gcMax < 0.65) {
                if (gcMin < 0.20) {
                        return 0.00076862;
                } else if (gcMin < 0.25) {
                        return 0.0007633;
                } else if (gcMin < 0.30) {
                        return 0.000638;
                } else if (gcMin < 0.35) {
                        return 0.000502;
                } else if (gcMin < 0.40) {
                        return 0.00030;
                } else if (gcMin < 0.45) {
                        return 0.00014;
                } else if (gcMin < 0.50) {
                        return 0.00002;
                } else if (gcMin < 0.55) {
                        return 0.000002;
                } else if (gcMin < 0.60) {
                        return 0.0000009;
                } else {
                        return 0.0000005;
                }
        } else {
                if (gcMin < 0.20) {
                        return 0.00063;
                } else if (gcMin < 0.25) {
                        return 0.00058;
                } else if (gcMin < 0.30) {
                        return 0.00050;
                } else if (gcMin < 0.35) {
                        return 0.00036;
                } else if (gcMin < 0.40) {
                        return 0.00020;
                } else if (gcMin < 0.45) {
                        return 0.00006;
                } else if (gcMin < 0.50) {
                        return 0.000003;
                } else if (gcMin < 0.55) {
                        return 0.0000013;
                } else if (gcMin < 0.60) {
                        return 0.0000006;
                } else if (gcMin < 0.65) {
                        return 0.0000004;
                } else {
                        return 0.0000005020;
                }
        }
}
		

#endif /* PVALUE */
