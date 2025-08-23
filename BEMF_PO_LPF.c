#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME BEMF_PO_LPF
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumDiscStates(S, 4); // i_hat_alpha, i_hat_beta, e_hat_alpha, e_hat_beta
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 4); // u_alpha, u_beta, i_alpha, i_beta
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 2); // e_alpha, e_beta
    ssSetNumSampleTimes(S, 1); 

    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
    | SS_OPTION_DISCRETE_VALUED_OUTPUT));
} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
    ssSetSampleTime(S, 0, 1e-4); 
    ssSetOffsetTime(S, 0, 0.0);
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S){ 
    real_T *X0 = ssGetRealDiscStates(S); 
    int_T nXStates = ssGetNumDiscStates(S); 
    for (int_T i = 0; i < nXStates; i++) { 
        X0[i] = 0.0; 
    } 
} 

static void mdlOutputs(SimStruct *S, int_T tid){ 
    real_T *Y = ssGetOutputPortRealSignal(S,0); 
    real_T *X = ssGetRealDiscStates(S); 

    // Output estimated BEMF
    Y[0] = X[2]; // e_hat_alpha
    Y[1] = X[3]; // e_hat_beta
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) { 
    real_T *X = ssGetRealDiscStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

    real_T dt = 1e-4;

    // Inputs
    real_T u_alpha = U(0);
    real_T u_beta  = U(1);
    real_T i_alpha = U(2);
    real_T i_beta  = U(3);
    real_T theta; // Not used in this implementation
    // Not used in this implementation, but can be used for further enhancements

    // Constants (from paper)
    real_T Rs    = 0.15;         // Ohm
    real_T Ls    = 0.193e-3;     // H
    real_T K     = 1.0;          // SMO gain
    real_T eps   = 5.6;   // wider transition zone
    real_T alpha = 0.055; // slower LPF

    real_T i_hat_alpha = X[0];
    real_T i_hat_beta  = X[1];
    real_T e_hat_alpha = X[2];
    real_T e_hat_beta  = X[3];

    // Current errors -> Jadi Sliding Surface
    real_T err_alpha = i_alpha - i_hat_alpha;
    real_T err_beta  = i_beta  - i_hat_beta;

    // Sliding mode gains 
    real_T G_alpha = K * tanh(err_alpha / eps);
    real_T G_beta  = K * tanh(err_beta / eps);

    // BEMF Estimasi pake LPF
    e_hat_alpha += alpha * (G_alpha - e_hat_alpha);
    e_hat_beta  += alpha * (G_beta  - e_hat_beta);

    // Update arus estimasi (Euler integration) -> Karena Discrete
    i_hat_alpha += dt * (1.0 / Ls) * (u_alpha - Rs * i_hat_alpha - e_hat_alpha + G_alpha);
    i_hat_beta  += dt * (1.0 / Ls) * (u_beta  - Rs * i_hat_beta  - e_hat_beta  + G_beta);

    X[0] = i_hat_alpha;
    X[1] = i_hat_beta;
    X[2] = e_hat_alpha;
    X[3] = e_hat_beta;
} 

static void mdlTerminate(SimStruct *S) 
{ } 

#ifdef MATLAB_MEX_FILE 
#include "simulink.c" 
#else 
#include "cg_sfun.h" 
#endif
