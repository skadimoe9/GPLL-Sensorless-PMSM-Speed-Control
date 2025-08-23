#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME rfocPMSM 
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
ssSetNumDiscStates(S, 8); 
if (!ssSetNumInputPorts(S, 1)) return; 
ssSetInputPortWidth(S, 0, 6); 
ssSetInputPortDirectFeedThrough(S, 0, 1); 
ssSetInputPortOverWritable(S, 0, 1); 
if (!ssSetNumOutputPorts(S, 1)) return; 
ssSetOutputPortWidth(S, 0, 5); 
ssSetNumSampleTimes(S, 1); 

ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
| SS_OPTION_DISCRETE_VALUED_OUTPUT));} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
ssSetSampleTime(S, 0, 1e-4); 
ssSetOffsetTime(S, 0, 0.0);} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S){ 
real_T *X0 = ssGetRealDiscStates(S); 
int_T nXStates = ssGetNumDiscStates(S); 
InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
int_T i; 

/* initialize the states to 0.0 */ 
for (i=0; i < nXStates; i++) { 
X0[i] = 0.0; } 

} 

static void mdlOutputs(SimStruct *S, int_T tid){ 
real_T *Y = ssGetOutputPortRealSignal(S,0); 
real_T *X = ssGetRealDiscStates(S); 

real_T va, vb, vc;

va  = X[5];
vb  = X[6];
vc  = X[7];

 Y[0] = va;
 Y[1] = vb;
 Y[2] = vc;
 Y[3] = X[2];
 Y[4] = X[3];
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) { 
   real_T *X = ssGetRealDiscStates(S); 
   InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

   real_T dt = 1e-4; 

   real_T ia, ib, ic, theta, theta_prv;
   real_T xid, xid_prv, xiq, xiq_prv, isd_1, isd_1_prv, isq_1, isq_1_prv;
   real_T isd, isq, ialfa, ibeta, isdref, isqref;
   real_T vsd, vsq, valfa, vbeta, va, vb, vc;
   real_T Upi_d, Upi_q, wr;

   real_T K	= 0.8164966;
   real_T L	= 0.866025;

   real_T N = 5.0;
   real_T psi = 0.0116;
   real_T Lsd = 0.193e-3;
   real_T Lsq = 0.193e-3;
   real_T Rs = 0.15;
   real_T Td = 0.01;
   
   real_T Kp_id = Lsd/Td;
   real_T Ki_id = Rs/Td;
   real_T Kp_iq = Lsq/Td;
   real_T Ki_iq = Rs/Td;        
            

   ia = U(0);
   ib = U(1);
   ic = U(2);
   theta = U(5);
   isdref = U(3);
   isqref = U(4);


    xid_prv = X[0];
    xiq_prv =  X[1];
    isd_1_prv = X[2];
    isq_1_prv = X[3];
    theta_prv = X[4];
    
wr = ((theta - theta_prv)/dt)/N;

   ialfa = (ia-0.5*ib-0.5*ic)*K;
   ibeta = (L*(ib-ic))*K;

  isd = ialfa*cos(theta) + ibeta*sin(theta);
  isq = -ialfa*sin(theta) + ibeta*cos(theta);

/* d-axis current control */
xid = xid_prv + dt*(isdref-isd);
Upi_d = Kp_id*(isdref-isd)+Ki_id*xid;

/* q-axis current control */
xiq = xiq_prv + dt*(isqref-isq);
Upi_q = Kp_iq*(isqref-isq)+Ki_iq*xiq;

real_T vsmax = 5;
/* Anti-Windup 1
if (Upi_d >= vsmax) { 
    Upi_d = vsmax; xid = xid_prv; 
} 
if (Upi_d <= -vsmax) { 
    Upi_d = -vsmax; xid = xid_prv;
}

if (Upi_q >= vsmax) { 
    Upi_q = vsmax; xiq = xiq_prv; 
} 
if (Upi_q <= -vsmax) { 
    Upi_q = -vsmax; xiq = xiq_prv;
}*/

/* decoupling equation */
isd_1 = isd_1_prv + dt*(isdref-isd_1_prv)/Td;
isq_1 = isq_1_prv + dt*(isqref-isq_1_prv)/Td;

vsd = -N*wr*Lsq*isq_1 + Upi_d;
vsq = N*wr*Lsd*isd_1 + N*wr*psi +Upi_q;

   /* Anti-Windup 2
if (vsd >= vsmax) { 
    vsd = vsmax; isd_1 = isd_1_prv;
}
if (vsd <= -vsmax) { 
    vsd = -vsmax; isd_1 = isd_1_prv;
}
if (vsq >= vsmax) { 
    vsq = vsmax; isq_1 = isq_1_prv;
}
if (vsq <= -vsmax) { 
    vsq = -vsmax; isq_1 = isq_1_prv;
}*/

   valfa =vsd*cos(theta)-vsq*sin(theta);
   vbeta =vsd*sin(theta)+vsq*cos(theta);

   va = valfa*K;
   vb = (-0.5*valfa+L*vbeta)*K;
   vc = (-0.5*valfa-L*vbeta)*K;

   X[0] = xid;
   X[1] = xiq;
   X[2] = isd_1;
   X[3] = isq_1;
   X[4] = theta;
   X[5] = va;
   X[6] = vb;
   X[7] = vc;
} 

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 