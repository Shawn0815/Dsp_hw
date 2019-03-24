//
// Created by yang on 2019-03-21.
//

#ifndef NUMOFTEST
#	define NUMOFTEST	2500
#endif

#ifndef TIME_LENGTH
#	define TIME_LENGTH	50
#endif

#ifndef STATE
#	define STATE	6
#endif

#ifndef MODEL
#	define MODEL	5
#endif


#ifndef DSP_HW1_VITERBI_H
#define DSP_HW1_VITERBI_H

#include "hmm.h"

typedef struct {
    short sequence[NUMOFTEST][TIME_LENGTH];
    double delta[NUMOFTEST][TIME_LENGTH][STATE][MODEL];
    short psi[NUMOFTEST][TIME_LENGTH][STATE][MODEL];
    double maxProbTable[NUMOFTEST][MODEL];
    int modelResult[NUMOFTEST];
}ViterbiBlock;

static void LoadTestSeq(ViterbiBlock *viterbiBlock, char *seqFile){
    FILE *fp = fopen(seqFile, "r");
    char token[NUMOFTEST] = "";
    int i=0;
    while( fscanf( fp, "%s", token ) > 0 )
    {
        if( token[0] == '\0' || token[0] == '\n' ) continue;
        for (int t=0; t<TIME_LENGTH; t++){
            viterbiBlock->sequence[i][t] = (short)token[t] - (short)'A';
        }
        i++;
    }
}

static void CalculateDelta(ViterbiBlock *viterbiBlock, HMM *hmm, int num_model){

    double maxPreState;
    double tempStatePro;
    int tempState;

    for (int u=0; u<NUMOFTEST; u++){
        for (int m=0; m<num_model; m++){
            for (int t=0; t<TIME_LENGTH; t++){
                for (int s=0; s<STATE; s++) {
                    maxPreState = 0;
                    tempStatePro = 0;
                    tempState = -1;
                    if (t == 0) {
                        viterbiBlock->delta[u][t][s][m] =
                                hmm[m].initial[s] * hmm[m].observation[viterbiBlock->sequence[u][t]][s];
                        //printf("%f ", viterbiBlock->delta[u][t][s]);
                    } else {
                        for (int ss = 0; ss < STATE; ss++) {
                            tempStatePro = viterbiBlock->delta[u][t-1][ss][m] * hmm[m].transition[ss][s];
                            if (maxPreState < tempStatePro) {
                                maxPreState = tempStatePro;
                                tempState = ss;
                            }
                        }
                        viterbiBlock->delta[u][t][s][m] =
                                maxPreState * hmm[m].observation[viterbiBlock->sequence[u][t]][s];
                        viterbiBlock->psi[u][t][s][m] = tempState;
                    }
                }

            }
        }



        //printf("\n");

    }

}

static void FindTheMaxProb(ViterbiBlock *viterbiBlock, int num_model){

    double maxProb;
    for (int u=0; u<NUMOFTEST; u++){
        for (int m=0; m<num_model; m++){
            maxProb=0;
            for (int s=0; s<STATE; s++){
                if (maxProb<viterbiBlock->delta[u][TIME_LENGTH-1][s][m]){
                    maxProb=viterbiBlock->delta[u][TIME_LENGTH-1][s][m];
                }
                viterbiBlock->maxProbTable[u][m]=maxProb;
            }
        }
    }

}

static void MatchTheModel(ViterbiBlock *viterbiBlock, int num_model){
    int modelNum=-1;
    double modelPro=0;
    for (int u=0; u<NUMOFTEST; u++){
        modelNum=0;
        modelPro=0;
        for (int m=0; m<num_model; m++){
            if (viterbiBlock->maxProbTable[u][m]>modelPro){
                modelPro = viterbiBlock->maxProbTable[u][m];
                modelNum=m+1;
            }
        }
        viterbiBlock->modelResult[u]=modelNum;
    }
}
#endif //DSP_HW1_VITERBI_H
