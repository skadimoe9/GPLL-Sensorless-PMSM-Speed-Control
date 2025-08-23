#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME ClarkeTF
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 


static void mdlInitializeSizes(SimStruct *S){ 
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 3); 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 2); 
    ssSetNumSampleTimes(S, 1); 

    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); } 

static void mdlInitializeSampleTimes(SimStruct *S) { 
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
    ssSetOffsetTime(S, 0, 0.0); } 

static void mdlOutputs(SimStruct *S, int_T tid) { 
    real_T *Y = ssGetOutputPortRealSignal(S,0); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    real_T t = ssGetT(S);

    real_T vi_a, vi_b, vi_c, vi_alfa, vi_beta;
    real_T K = 0.8164966;
    real_T L = 0.866025;

    vi_a = U(0);
    vi_b = U(1);
    vi_c = U(2);

    vi_alfa = K*(vi_a-0.5*vi_b-0.5*vi_c);  // Transformasi Clarke
    vi_beta = K*L*(vi_b-vi_c);

    Y[0] = vi_alfa; // Output alfa
    Y[1] = vi_beta; // Output beta
}


static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif