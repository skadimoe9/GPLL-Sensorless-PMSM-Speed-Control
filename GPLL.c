#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME GPLL
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumDiscStates(S, 6); // [theta_hat, omega_pll, omega_ffilt, integrator, polarity]
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 2); // [e_alpha, e_beta] 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 4); // [theta_hat, omega_total, omega_pll]
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

    real_T theta_hat   = X[0];
    real_T omega_pll   = X[1];
    real_T omega_ffilt = X[2];
    real_T omega_all    = X[5];

    real_T theta_norm = fmod(theta_hat, 2 * M_PI);
    if (theta_norm < 0) theta_norm += 2 * M_PI;

    Y[0] = theta_hat;
    Y[1] = theta_norm;
    Y[2] = omega_all; // total estimated speed
    Y[3] = omega_pll;               // PLL-only component

} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) { 
    real_T *X = ssGetRealDiscStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    real_T t = ssGetT(S);
    // Sample time
    real_T dt = 1e-4;

    // Inputs
    real_T e_alpha = U(0);
    real_T e_beta  = U(1);

    // Constants (from paper)
    const real_T psi_f  = 0.0116;
    const real_T wc     = 80.0 * 2.0 * M_PI;
    const real_T wpll   = 40.0 * 2.0 * M_PI; // 40 Hz in rad/s
    
    //const real_T kp     = (sqrt(3.0)*(wpll))/2.0; // 40 Hz in rad/s, divided by 2 for PLL
    const real_T kp     = 50.5;
    //const real_T ki     = (wpll*wpll)/2.0; // 40 Hz in rad/s, divided by 2 for PLL
    const real_T ki     = 30000.0; 

    const real_T alpha_ff = dt / (dt + (1.0 / wc)); // LPF smoothing factor
    const real_T Pcor   = 2.3;
    const real_T EPS    = 1e-9; // smal3 value to prevent divide-by-zero
    real_T alpha_pll    = 1.0; // Scale contribution of PLL
    real_T N = 5.0;
   
    // State variables
    real_T theta_hat = X[0];
    real_T omega_pll = X[1];
    real_T omega_ffilt = X[2];
    real_T x_int = X[3];
    real_T polarity = X[4];
    real_T omega_all = X[5];

    real_T bemf_mag = sqrt(e_alpha * e_alpha + e_beta * e_beta);
    real_T omega_ff = bemf_mag / psi_f;
    omega_ffilt += alpha_ff * (omega_ff - omega_ffilt);
    
    real_T omega_total = omega_ffilt + omega_pll;

    real_T cos2t = cos(2.0 * theta_hat);
    real_T sin2t = sin(2.0 * theta_hat);
    real_T e_alpha2 = e_alpha * e_alpha;
    real_T e_beta2  = e_beta * e_beta;
    real_T denom = e_alpha2 + e_beta2 + EPS;
    real_T num = (-e_alpha * e_beta * cos2t) + (0.5 * (e_alpha2 - e_beta2) * sin2t);
    real_T error = num / denom;

    real_T fcos = (-e_alpha * sin(theta_hat) + e_beta * cos(theta_hat)) * omega_total;
    if (t > 0.0) { 
        polarity = Pcor * ((fcos < 0 ?1: 1.0 - 2.0 * Pcor) - 1.0) + 1.0;
    } else {
        polarity = 1.0; 
    }

    real_T error_c = polarity * error;

    x_int += dt * error_c;
    omega_pll = kp * error_c + ki * x_int;
    
    omega_all = (omega_pll + omega_ffilt)/N;
    theta_hat += dt * (omega_pll + omega_ffilt);

    // Store updated states
    X[0] = theta_hat;
    X[1] = omega_pll;
    X[2] = omega_ffilt;
    X[3] = x_int;
    X[4] = polarity;
    X[5] = omega_all; // Store total estimated speed
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
#include "simulink.c" 
#else 
#include "cg_sfun.h" 
#endif 
