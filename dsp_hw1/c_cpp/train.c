//
// Created by yang on 2019-02-28.
//

#ifndef TIME_LENGTH
#	define TIME_LENGTH	50
#endif

#ifndef STATE
#	define STATE	6
#endif

#ifndef OBSERVE
#	define OBSERVE	6
#endif



#include "hmm.h"
#include "math.h"
#include "float.h"
#include "stdio.h"
#include "forbackward.h"

HMM hmm_initial;
ForBackwardBlock forwardBlock;

int main(int argc, char *argv[]){
    if (argc != 5){
        printf("The number of arguments is not correct");
        exit(1);
    }

    //HMM hmm_initial;

    //get the number of iterations
    int num_Iteration = atoi(argv[1]);
    //printf("%d\n",num_Iteration);

    //init the model
    char* initModel = argv[2];
    loadHMM(&hmm_initial, initModel);
    initModel="model_init_temp.txt";
	//dumpHMM( stderr, &hmm_initial );

    //load the training sequences
    char *seqFile = argv[3];


    // Training
    LoadSequence(&forwardBlock, seqFile);
    short iteration=1;

    do {

        CalculateAlphaN(&forwardBlock, &hmm_initial);
        CalculateBetaN(&forwardBlock, &hmm_initial);
        CalculateGammaN(&forwardBlock);
        CalculateEpsilonN(&forwardBlock, &hmm_initial);
        CalculateNewTransitionN(&forwardBlock);
        CalculateNewObservationN(&forwardBlock);
        CalculateNewInitial(&forwardBlock);

        FILE *fp = NULL;

        fp = fopen("model_init_temp.txt", "w+");

        fprintf(fp, "initial: %d \n", STATE);
        for (int tt=0; tt<STATE; tt++){
            fprintf(fp, "%.15f ", forwardBlock.newInitial[tt]);
        }

        fprintf(fp, "\n\ntransition: %d \n", STATE);
        for (int tt=0; tt<STATE; tt++){
            for (int st=0; st<STATE; st++){
                fprintf(fp, "%.15f ", forwardBlock.newTransition[tt][st]);
            }
            fprintf(fp, "\n");
        }

        fprintf(fp, "\n\nobservation: %d \n", STATE);
        for (int tt=0; tt<STATE; tt++){
            for (int st=0; st<STATE; st++){
                fprintf(fp, "%.15f ", forwardBlock.newObservation[tt][st]);
            }
            fprintf(fp, "\n");
        }

        fclose(fp);

        loadHMM(&hmm_initial, initModel);

        if (iteration%10==0){
            printf("Iteration %d compelted! \n", iteration);
        }


        iteration++;

    }while (iteration<num_Iteration+1);





//    CalculateBeta(&forwardBlock, &hmm_initial);
//    CalculateEpsilon(&forwardBlock, &hmm_initial);
//    CalculateGamma(&forwardBlock);
//    CalculateNewTransition(&forwardBlock);
//
//    CalculateNewObservation(&forwardBlock);


    //printf("%d", alpha[0][0]);



    // Calculate Gamma and Epsilon

    // Calculate new Transaction matrix and Observe Probability

    // Get new Initial probability

    //Iteration

    //write out the model file

}




