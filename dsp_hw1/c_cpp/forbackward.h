//
// Created by yang on 2019-03-19.
//

#ifndef DSP_HW1_FORWARD_H
#define DSP_HW1_FORWARD_H

#ifndef TIME_LENGTH
#	define TIME_LENGTH	50
#endif

#ifndef STATE
#	define STATE	6
#endif

#ifndef OBSERVE
#	define OBSERVE	6
#endif

#ifndef NUMOFTEST
#	define NUMOFTEST	10000
#endif


#include "hmm.h"
#include "math.h"
#include "float.h"
#include "stdio.h"

typedef struct {

    int sequence[NUMOFTEST][TIME_LENGTH+1];
    double alphaN[NUMOFTEST][TIME_LENGTH][STATE];
    double betaN[NUMOFTEST][TIME_LENGTH][STATE];
    double gammaN[NUMOFTEST][TIME_LENGTH][STATE];
    double epsilonN[NUMOFTEST][TIME_LENGTH][STATE][STATE];
    double newInitial[STATE];

    double alpha[TIME_LENGTH][STATE];
    double beta[TIME_LENGTH][STATE];
    double gamma[TIME_LENGTH][STATE];
    double epsilon[TIME_LENGTH][STATE][STATE];
    double newTransition[STATE][STATE];
    double newObservation[STATE][OBSERVE];

}ForBackwardBlock;

// Load the test sequence into numbers, and store into sequence[][]
static void LoadSequence(ForBackwardBlock *forwardBlock, char *seqFile){

    FILE *fp = fopen(seqFile, "r");

    char *temp;
    temp = (char *)malloc(52* sizeof(char));

    // Make temp array to be all ""
    for (int st=0; st<52; st++){
            temp[st]="";
    }

    for (int i=0; i<10000; i=i+1)
    {
        fgets(temp, 52, fp);
        for (int j=0; j<50; j++){
            forwardBlock->sequence[i][j] = (int)temp[j] - (int)'A';
        }
        forwardBlock->sequence[i][50]='\0';

    }

    free(temp);
    //printf("%s", sequence[1]);
    fclose(fp);
}

//static void CalculateAlpha(ForBackwardBlock *forwardBlock, HMM *hmm){
//
//    // Calculate alpha_0(i)
//    int multiple=0;
//    double init_alpha=0;
//    // Make alpha to be all -1
//        for (int tt = 0; tt < TIME_LENGTH; tt++) {
//            for (int st = 0; st < STATE; st++) {
//                forwardBlock->alpha[tt][st] = -1;
//            }
//        }
//
//        printf("Initial Alpha[0].\n");
//        for (int j = 0; j < STATE; j++) {
//
//            int obse = forwardBlock->sequence[0][0];
//            init_alpha = hmm->initial[j] * hmm->observation[j][obse];
//            multiple++;
//            forwardBlock->alpha[0][j] = init_alpha;
//            //printf("%f ",hmm->initial[j] * hmm->observation[j][o]);
//            //printf("\n");
//        }
//
//        // Calculate alpha_t(i)
//        printf("Calculate and update Alpha[t] (t>1).\n");
//        for (int t = 1; t < TIME_LENGTH; t++) {
//
//            //Calculate initial matrix for time=t>0
//            double temp_initial[STATE];
//            for (int st = 0; st < STATE; st++) {
//                temp_initial[st] = 0;
//            }
//
//            // Sum the alpha_i * transition_ij
//            for (int s1 = 0; s1 < STATE; s1++) {
//                for (int s2 = 0; s2 < STATE; s2++) {
//                    temp_initial[s2] += forwardBlock->alpha[t - 1][s1] * hmm->transition[s2][s1];
//                }
//            }
//
//            //Calculate alpha_t(i)= sum(alpha_i * transition_ij) * observation_ij
//            int obse = forwardBlock->sequence[0][t];
//            for (int j = 0; j < STATE; j++) {
//                double init_alpha = temp_initial[j] * hmm->observation[j][obse];
//                multiple++;
//                forwardBlock->alpha[t][j] = init_alpha;
//                //printf("%f\n",hmm->initial[j] * hmm->observation[j][o]);
//            }
//        }
//
//        //printf("Update Alpha for sequence[0] completed! For multiple %d times.\n\n", multiple);
//
//}

static void CalculateAlphaN(ForBackwardBlock *forwardBlock, HMM *hmm){

    int multiple=0; // Calculate the multiple times

    // Make alpha array  to be all -1
    for (int u=0; u<NUMOFTEST; u++){
        for (int tt=0; tt<TIME_LENGTH; tt++){
            for (int st=0; st<STATE; st++){
                forwardBlock->alphaN[u][tt][st]=-1;
            }
        }
    }

    //printf("Start to calculate Alpha[][][].\n");
    for (int u=0; u<NUMOFTEST; u++){

        double alpha_initial=0; // Initial array in alpha value

        // Calculate alpha_0(i)
        for (int j=0; j<STATE; j++){

            int obse = forwardBlock->sequence[u][0];
            alpha_initial = hmm->initial[j] * hmm->observation[obse][j];
            multiple++;
            forwardBlock->alphaN[u][0][j] = alpha_initial;
            //printf("%f ",hmm->initial[j] * hmm->observation[j][o]);
            //printf("\n");
        }

        // Calculate alpha_t(i)
        //printf("Calculate and update Alpha[][t] (t>1).\n");
        for (int t=1; t<TIME_LENGTH; t++){

            //Calculate initial matrix for time=t>0
            double *temp_initial;
            temp_initial = (double *)malloc(STATE* sizeof(double));

            for (int st=0; st<STATE; st++){
                temp_initial[st]=0;
            }

            // Sum the alpha_i * transition_ij
            for (int s1=0; s1<STATE; s1++){
                for (int s2=0; s2<STATE; s2++){
                    temp_initial[s2] += forwardBlock->alphaN[u][t-1][s1] * hmm->transition[s1][s2];
                }
            }

            //Calculate alpha_t(i)= sum(alpha_i * transition_ij) * observation_ij
            int obse = forwardBlock->sequence[u][t];
            for (int j=0; j<STATE; j++){
                double init_alpha = temp_initial[j] * hmm->observation[obse][j];
                multiple++;
                forwardBlock->alphaN[u][t][j] = init_alpha;
                //printf("%f\n",hmm->initial[j] * hmm->observation[j][o]);
            }
            free(temp_initial);
        }
    }

    //printf("Update all Alpha for sequence[] completed! For multiple %d times.\n\n", multiple);

}

static void CalculateBeta(ForBackwardBlock *forwardBlock, HMM *hmm){
    // Calculate beta
    int multiple=0; // multiple times
    double init_beta=1.0;

    // Make alpha to be all -1
    for (int tt=0; tt<TIME_LENGTH; tt++){
        for (int st=0; st<STATE; st++){
            forwardBlock->beta[tt][st]=0;
        }
    }

    printf("Initial Beta[T]=1.\n");
    for (int j=0; j<STATE; j++){
        forwardBlock->beta[TIME_LENGTH-1][j] = init_beta;
        //printf("Initial Beta[TIME_LENGTH-1][%d] = %f\n", j, init_beta);
    }

    // Calculate alpha_t(i)
    printf("Calculate and update Beta[t] t<T).\n");
    for (int t=TIME_LENGTH-2; t>=0; t--){

        //Calculate initial matrix for time=t>0
        double temp_initial[STATE];
        for (int st=0; st<STATE; st++){
            temp_initial[st]=0;
        }

        //Multiple of transition_ij * beta_j, the state probability at t
        int obse = forwardBlock->sequence[0][t+1];
        for (int j=0; j<STATE; j++){
            temp_initial[j] = forwardBlock->beta[t+1][j] * hmm->observation[obse][j];
        }

        // The Sum of observation_io * transition_ij * beta_j
        for (int s1=0; s1<STATE; s1++){
            for (int s2=0; s2<STATE; s2++){
                forwardBlock->beta[t][s1] += temp_initial[s2] * hmm->transition[s1][s2];
                multiple++;
            }
            //printf("%f\n",hmm->initial[j] * hmm->observation[j][o]);
        }

    }
    printf("Update Beta for sequence[0] completed! For multiple %d times.\n\n", multiple);

}

static void CalculateBetaN(ForBackwardBlock *forwardBlock, HMM *hmm){
    // Calculate beta
    int multiple=0; // multiple times
    double init_beta=1.0;

    // Make alpha to be all -1
    for (int u=0; u<NUMOFTEST; u++){
        for (int tt=0; tt<TIME_LENGTH; tt++){
            for (int st=0; st<STATE; st++){
                forwardBlock->betaN[u][tt][st]=0;
            }
        }
    }

    //printf("Initial Beta[][T]=1.\n");

    for (int u=0; u<NUMOFTEST; u++){

        for (int j=0; j<STATE; j++){
            forwardBlock->betaN[u][TIME_LENGTH-1][j] = init_beta;
            //printf("Initial Beta[TIME_LENGTH-1][%d] = %f\n", j, init_beta);
        }

        // Calculate alpha_t(i)
        //printf("Calculate and update Beta[t] t<T).\n");
        for (int t=TIME_LENGTH-2; t>=0; t--){

            //Calculate initial matrix for time=t>0
            double *temp_initial;
            temp_initial = (double *)malloc(STATE* sizeof(double));

            for (int st=0; st<STATE; st++){
                temp_initial[st]=0;
            }

            //Multiple of transition_ij * beta_j, the state probability at t
            int obse = forwardBlock->sequence[u][t+1];
            for (int j=0; j<STATE; j++){
                temp_initial[j] = forwardBlock->betaN[u][t+1][j] * hmm->observation[obse][j];
            }

            // The Sum of observation_io * transition_ij * beta_j
            for (int s1=0; s1<STATE; s1++){
                for (int s2=0; s2<STATE; s2++){
                    forwardBlock->betaN[u][t][s1] += temp_initial[s2] * hmm->transition[s1][s2];
                    multiple++;
                }
                //printf("%f\n",hmm->initial[j] * hmm->observation[j][o]);
            }
            free(temp_initial);

        }
    }


    //printf("Update Beta for sequence[] completed! For multiple %d times.\n\n", multiple);

}

static void CalculateGammaN(ForBackwardBlock *forwardBlock){

    //printf("Start to calculate GammaN for sequence[].\n");

    for (int u=0; u<NUMOFTEST; u++){

        double *state_i;
        state_i = (double *)malloc(STATE* sizeof(double));

        double sumofState_i=0;
        for (int t=0; t<TIME_LENGTH; t++){

            sumofState_i=0;
            for (int i=0; i<STATE; i++){
                state_i[STATE]=-1;
            }

            for (int i=0; i<STATE; i++){
                state_i[i]=forwardBlock->alphaN[u][t][i] * forwardBlock->betaN[u][t][i];
                sumofState_i += state_i[i];
            }

            for (int i=0; i<STATE; i++){
                forwardBlock->gammaN[u][t][i]=state_i[i] / sumofState_i;
            }

        }

        free(state_i);
    }

    //printf("Calculating GammaN for sequence[] completed!\n\n");
}

static void CalculateGamma(ForBackwardBlock *forwardBlock){
    for (int t=0; t<TIME_LENGTH; t++){
        double *state_i;
        state_i = (double *)malloc(STATE* sizeof(double));
        //double state_i[STATE];
        double sumofState_i=0;
        for (int i=0; i<STATE; i++){
            state_i[STATE]=-1;
        }

        for (int i=0; i<STATE; i++){
            state_i[i]=forwardBlock->alpha[t][i] * forwardBlock->beta[t][i];
            sumofState_i += state_i[i];
        }

        for (int i=0; i<STATE; i++){
            forwardBlock->gamma[t][i]=state_i[i] / sumofState_i;
        }
        free(state_i);

    }
}

static void CalculateEpsilon(ForBackwardBlock *forwardBlock, HMM *hmm){

    double state_ij[STATE][STATE];
    double sumofState_ij;
    for (int t=0; t<TIME_LENGTH-1; t++){
        // Calculate Epsilon
        sumofState_ij=0;
        for (int i=0; i<STATE; i++){
            for (int j=0; j<STATE; j++){
                state_ij[i][j]=-1;
            }
        }

        for (int i=0; i<STATE; i++){
            for (int j=0; j<STATE; j++){
                double alpha_temp = forwardBlock->alpha[t][i];
                double beta_temp = forwardBlock->beta[t+1][j];
                double trans_temp = hmm->transition[i][j];
                double obser_temp = hmm->observation[forwardBlock->sequence[0][t+1]][j];
                state_ij[i][j] = alpha_temp * beta_temp * trans_temp * obser_temp;
                sumofState_ij = sumofState_ij + state_ij[i][j];
            }
        }

        for (int i=0; i<STATE; i++){
            for (int j=0; j<STATE; j++){
                forwardBlock->epsilon[t][i][j]= state_ij[i][j] / sumofState_ij;
            }
        }

    }
}

static void CalculateEpsilonN(ForBackwardBlock *forwardBlock, HMM *hmm){

    for (int u=0; u<NUMOFTEST; u++){

        double state_ij[STATE][STATE];
        double sumofState_ij;

        for (int t=0; t<TIME_LENGTH-1; t++){
            // Calculate Epsilon
            sumofState_ij=0;
            for (int i=0; i<STATE; i++){
                for (int j=0; j<STATE; j++){
                    state_ij[i][j]=-1;
                }
            }

            for (int i=0; i<STATE; i++){
                for (int j=0; j<STATE; j++){
                    double alpha_temp = forwardBlock->alphaN[u][t][i];
                    double beta_temp = forwardBlock->betaN[u][t+1][j];
                    double trans_temp = hmm->transition[i][j];
                    double obser_temp = hmm->observation[forwardBlock->sequence[u][t+1]][j];
                    state_ij[i][j] = alpha_temp * beta_temp * trans_temp * obser_temp;
                    sumofState_ij = sumofState_ij + state_ij[i][j];
                }
            }

            for (int i=0; i<STATE; i++){
                for (int j=0; j<STATE; j++){
                    forwardBlock->epsilonN[u][t][i][j]= state_ij[i][j] / sumofState_ij;
                }
            }

        }
    }

}

static void CalculateNewInitial(ForBackwardBlock *forwardBlock){

    //printf("Start to calculate the new Initial matrix!\n");

    double sumofInitial;
    for (int i=0; i<STATE; i++){
        forwardBlock->newInitial[i]=0;
    }
    for (int i=0; i<STATE; i++){
        sumofInitial=0;
        for (int u=0; u<NUMOFTEST; u++){
            sumofInitial += forwardBlock->gammaN[u][0][i];
        }
        forwardBlock->newInitial[i] = sumofInitial / NUMOFTEST;
    }
    //printf("Calculating the new initial matrix completed!\n\n");

}

static void CalculateNewTransition(ForBackwardBlock *forwardBlock){
    printf("Start to calculate the new transition matrix!\n");
    double sumofEpsilon=0;
    double sumofGamma=0;
    //double sumofState_i=0;
    for (int i=0; i<STATE; i++){
        //sumofState_i=0;
        for (int j=0; j<STATE; j++){
            sumofEpsilon=0;
            sumofGamma=0;
            for (int t=0; t<TIME_LENGTH-1; t++){
                sumofEpsilon += forwardBlock->epsilon[t][i][j];
                sumofGamma += forwardBlock->gamma[t][i];
            }
            forwardBlock->newTransition[i][j]=sumofEpsilon/sumofGamma;
            //sumofState_i += forwardBlock->newTransition[i][j];
        }
        //printf("The sum of row %d ", i);
        //printf(" is %f \n", sumofState_i);
    }
    printf("Calculating the new transition matrix completed!\n\n");
}

static void CalculateNewTransitionN(ForBackwardBlock *forwardBlock){
    //printf("Start to calculate the new transition matrix!\n");
    double sumofEpsilon=0;
    double sumofGamma=0;

    for (int i=0; i<STATE; i++){
        for (int j=0; j<OBSERVE; j++){
            forwardBlock->newTransition[j][i]=0;
        }
    }

    //double sumofState_i=0;
    for (int i=0; i<STATE; i++){
        //sumofState_i=0;
        for (int j=0; j<STATE; j++){
            sumofEpsilon=0;
            sumofGamma=0;
            for (int u=0; u<NUMOFTEST; u++){
                for (int t=0; t<TIME_LENGTH-1; t++){
                sumofEpsilon += forwardBlock->epsilonN[u][t][i][j];
                sumofGamma += forwardBlock->gammaN[u][t][i];
                }
            }

            forwardBlock->newTransition[i][j]=sumofEpsilon/sumofGamma;
            //sumofState_i += forwardBlock->newTransition[i][j];
        }
        //printf("The sum of row %d ", i);
        //printf(" is %f \n", sumofState_i);
    }
    //printf("Calculating the new transition matrix completed!\n\n");
}

static void CalculateNewObservation(ForBackwardBlock *forwardBlock){
    printf("Start to calculate the new observation matrix!\n");
    double sumofObse=0;
    double sumofAll=0;
    //double sumofState_j=0;
    for (int j=0; j<STATE; j++){
        //sumofState_j=0;
        for (int k=0; k<OBSERVE; k++){
            sumofObse=0;
            sumofAll=0;
            for (int t=0; t<TIME_LENGTH; t++){
                if (k==forwardBlock->sequence[0][t]){
                    sumofObse = sumofObse + forwardBlock->gamma[t][j];
                }
                sumofAll = sumofAll + forwardBlock->gamma[t][j];
            }
            forwardBlock->newObservation[k][j] = sumofObse/sumofAll;

            //sumofState_j += forwardBlock->newObservation[j][k];
        }
        //printf("The sum of row %d ", j);
        //printf(" is %f \n", sumofState_j);
    }
    printf("Calculating the new observation matrix completed!\n\n");
}

static void CalculateNewObservationN(ForBackwardBlock *forwardBlock){
    //printf("Start to calculate the new observation matrix!\n");
    double sumofObse=0;
    double sumofAll=0;

    for (int i=0; i<STATE; i++){
        for (int j=0; j<OBSERVE; j++){
            forwardBlock->newObservation[i][j]=0;
        }
    }
    //double sumofState_j=0;
    for (int j=0; j<STATE; j++){
        //sumofState_j=0;
        for (int k=0; k<OBSERVE; k++){

            sumofObse=0;
            sumofAll=0;
            for (int u=0; u<NUMOFTEST; u++){
                for (int t=0; t<TIME_LENGTH; t++){
                    if (k==forwardBlock->sequence[u][t]){
                        sumofObse = sumofObse + forwardBlock->gammaN[u][t][j];
                    }
                    sumofAll = sumofAll + forwardBlock->gammaN[u][t][j];
                }
            }
//            if (sumofObse==0){
//                printf("The observation of state is %d observation is %d\n", j, k);
//                exit(1);
//            }

            forwardBlock->newObservation[k][j] = sumofObse/sumofAll;
            //sumofState_j += forwardBlock->newObservation[j][k];
        }
        //printf("The sum of row %d ", j);
        //printf(" is %f \n", sumofState_j);
    }
    //printf("Calculating the new observation matrix completed!\n\n");
}

static void CheckResultMatrix(ForBackwardBlock *forwardBlock, int iteration){

    int transition=0;
    int observation=0;
    int initial=0;
    for (int i=0; i<STATE; i++){
        for (int j=0; j<OBSERVE; j++){
            if (forwardBlock->newTransition[i][j] >= 0 &&  forwardBlock->newTransition[i][j] < 1) continue;
            else{
                transition=1;
                break;
            }
        }
    }

    for (int i=0; i<STATE; i++){
        for (int j=0; j<OBSERVE; j++){
            if (forwardBlock->newObservation[i][j] >= 0 &&  forwardBlock->newObservation[i][j] < 1) continue;
            else{
                observation=1;
                break;
            }
        }
    }

    for (int j=0; j<STATE; j++){
        if (forwardBlock->newInitial[j] >= 0 && forwardBlock->newInitial[j]<1) continue;
        else{
            initial=1;
            break;
        }
    }
    if (transition){
        printf("Transition matrix computation error!\n");
    }
    if (observation){
        printf("Observation matrix computation error!\n");
    }
    if (transition){
        printf("Initial matrix computation error!\n");
    }

    if (transition+observation+initial!=0){
        printf("Zero components show.");
        //exit(1);
    }
    else if( iteration%10==0 ){
        printf("Check completed without error for iteration #%d!\n", iteration);
    }


}
#endif //DSP_HW1_FORWARD_H
