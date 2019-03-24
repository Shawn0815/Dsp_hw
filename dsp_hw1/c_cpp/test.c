//
// Created by yang on 2019-02-28.
//
#include "hmm.h"
#include "math.h"
#include "float.h"
#include "viterbi.h"
#include "stdio.h"

ViterbiBlock viterbiBlock;
HMM trainedHMM[5];

int main(int argc, char *argv[]){

    if (argc != 4){
        printf("The number of arguments is not correct");
        exit(1);
    }

    char* modellist = argv[1];
    char* testing_data = argv[2];
    char* result_file = argv[3];
    //HMM hmms[5];
	//loadHMM(&trainedHMM, "model_01.txt");
    load_models(modellist, trainedHMM, 5);

    LoadTestSeq(&viterbiBlock, testing_data);
    CalculateDelta(&viterbiBlock, &trainedHMM, 5);
    FindTheMaxProb(&viterbiBlock, 5);
    MatchTheModel(&viterbiBlock, 5);


    FILE *fp = NULL;
    fp = fopen(result_file, "w+");

    //fprintf(fp, "initial: %d \n", STATE);
    for (int u=0; u<NUMOFTEST; u++){
        fprintf(fp, "model_0%d.txt     %g\n", viterbiBlock.modelResult[u], viterbiBlock.maxProbTable[u][viterbiBlock.modelResult[u]-1]);
    }

    fclose(fp);
    char *answer_file = "testing_answer.txt";
    char answer[NUMOFTEST];
    short correct=0;
    fp = fopen(answer_file, "r");
    int u=0;
    while( fscanf( fp, "%s", answer) > 0 ){

        if ((int)(answer[7]-'0')==viterbiBlock.modelResult[u]){
            correct++;
        }
        u++;
    }

    printf("Correct one number is %d, Accuracy is %f %%", correct,  ((double)correct / NUMOFTEST) );
    exit(1);
    //*init_model[];

}
