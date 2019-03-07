#include <stdio.h>
#include "cgp.h"
#include <math.h>
double functionName(const int numInputs, const double *inputs, const double *connectionWeights)
{
    if (numInputs == 1)
    {   
        return 0.512384;
    }else
    {
        return 0;
    }
}

double meanSquareError(struct parameters *params, struct chromosome *chromo, struct dataSet *data){

    int i,j;
    double squareError = 0;

    if(getNumChromosomeInputs(chromo) !=getNumDataSetInputs(data)){
        printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
        printf("Terminating.\n");
        exit(0);
    }

    if(getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)){
        printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
        printf("Terminating.\n");
        exit(0);
    }

    
    
    for(i=0; i<getNumDataSetSamples(data); i++){

        executeChromosome(chromo, getDataSetSampleInputs(data, i));

        for(j=0; j<getNumChromosomeOutputs(chromo); j++){

            squareError += pow(getDataSetSampleOutput(data,i,j) - getChromosomeOutput(chromo,j), 2);
        }
    }

    

    return squareError / (getNumDataSetSamples(data) * getNumDataSetOutputs(data));
}



int main(void){

    struct parameters *params = NULL;
    struct dataSet *trainingData = NULL;
    struct chromosome *chromo = NULL;

    int numInputs = 4;
    int numNodes = 20;
    int numOutputs = 1;
    int nodeArity = 4;

    int numGens = 100000;
    int updateFrequency = 500;
    double targetFitness = 0;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,mul,abs,exp,sin,cos,tan,logic,and,nand,or,nor,xor,xnor,not,neuron,sig,gauss,step,softsign,tanh,rand,pi,1,0,wire");
//    addCustomNodeFunction(params, functionName , "hypt", -1);

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);
    setMutationRate(params,0.01);
//    setMutationType(params, "onlyActive");
    printParameters(params);

    double recurrentConnectionProbability = 0.2;

        
    trainingData = initialiseDataSetFromFile("./dataSets/symbolic.data");
    setRecurrentConnectionProbability(params, recurrentConnectionProbability);

    setNumThreads(params, 3); 
    
    setCustomFitnessFunction(params, meanSquareError, "MSE");

    chromo = runCGP(params, trainingData, numGens);

    printChromosome(chromo, 0);
    saveChromosomeDot(chromo, 0, "chromo.dot");
    saveChromosomeLatex(chromo, 0, "chromo.tex");
    
    double a1[] = {5.1,3.5,1.4,0.2}; //0 
    double a2[] = {5.1,3.5,1.4,0.2}; // 1
    double a3[] = {5.1,3.5,1.4,0.2}; //2
    double o1 ; 
    double o2 ; 
    double o3 ; 
    executeChromosome(chromo, a1);
    o1 = getChromosomeOutput(chromo, 0);
    executeChromosome(chromo, a2);
    o2 = getChromosomeOutput(chromo, 0);
    executeChromosome(chromo, a3);
    o3 = getChromosomeOutput(chromo, 0);
//      printf("ASCII value = %f\n", o1 );

    printf("actual is 0 and predicted is %f \n" , o1);
    printf("actual is 1 and predicted is %f \n" ,o2);
    printf("actual is 2 and predicted is %f \n" ,o3);
    
            
    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);

    return 0;
}