#define S_FUNCTION_LEVEL 2                                                       
#define S_FUNCTION_NAME PMSM                                             
#include "simstruc.h"                                                                           
#include <math.h>                                                                                

#define U(element) (*uPtrs0[element])	/*Pointer to Input Port0*/			

static void mdlInitializeSizes(SimStruct *S)
{                                   
  ssSetNumContStates(S, 4);                                                               

  if (!ssSetNumInputPorts(S, 1)) return;                                            
  ssSetInputPortWidth(S, 0, 4);                                                          
  ssSetInputPortDirectFeedThrough(S, 0, 1);                                   
  ssSetInputPortOverWritable(S, 0, 1);                                            

  if (!ssSetNumOutputPorts(S, 1)) return;                                        
  ssSetOutputPortWidth(S, 0, 7);                                                         

  ssSetNumSampleTimes(S, 1);                                                            
                                                                                                              
  ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); 
}    
                                                                                                              
static void mdlInitializeSampleTimes(SimStruct *S) 
{                   
  ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);         

  ssSetOffsetTime(S, 0, 0.0); 
}                                                              
                                                                                                             
#define MDL_INITIALIZE_CONDITIONS                                  
static void mdlInitializeConditions(SimStruct *S) 
{                                                                                                                  
  real_T 	*X0 = ssGetContStates(S);                                       
  int_T 	nStates = ssGetNumContStates(S);                         
  int_T 	i;                                                                                  
                                                                                                            
  /* initialize the states to 0.0 */                                                            
  for (i=0; i < nStates; i++) 
	{
	X0[i] = 0.0;
	} 
//  X0[3] = 100.0;
}                                       
                                                                                                              
static void mdlOutputs(SimStruct *S, int_T tid) 
{                          
  real_T            *Y = ssGetOutputPortRealSignal(S,0);                  
  real_T            *X = ssGetContStates(S);                                         
  InputRealPtrsType uPtrs0 = ssGetInputPortRealSignalPtrs(S,0);  

  real_T ia, ib, ic, wr, theta, ialfa, ibeta, isd, isq, te, theta_out;
  
  /* transformation constant */
  real_T K	= 0.8164966;
  real_T L	= 0.866025;
  
  wr = X[0];
  isd = X[1];
  isq = X[2];
  theta = X[3];
  
    /* dq to alfa beta transformation */
 ialfa = isd*cos(theta) - isq*sin(theta);
 ibeta = isd*sin(theta) + isq*cos(theta);

 /* alfa beta to abc transformation */
 ia = K*ialfa; 
 ib = K*( -0.5*ialfa + L*ibeta ); 
 ic = K*( -0.5*ialfa - L*ibeta ); 

/* PMSM parameters */
real_T N = 5.0;
real_T psi = 0.0116;
real_T Lsd = 0.193e-3;
real_T Lsq = 0.193e-3;
real_T Rs = 0.15;
real_T J = 0.015;//0.01;

te = N*(psi+(Lsd-Lsq)*isd)*isq;
  
 /* Normalize theta to 0 to 2π for output */
    theta_out = fmod(theta, 2 * M_PI);
    if (theta_out < 0) {
        theta_out += 2 * M_PI; // Ensure theta is positive
    }

Y[0]= ia; 
Y[1]= ib;
Y[2]= ic;
Y[3]= wr; 
Y[4]= theta;
Y[5]= te; 
Y[6]= theta_out; // Output the normalized theta value

}                                                                                           



#define MDL_DERIVATIVES                                                          
static void mdlDerivatives(SimStruct *S) 
{                                      
  real_T  *dX = ssGetdX(S);                                                               
  real_T  *X = ssGetContStates(S);                                                    
  InputRealPtrsType uPtrs0 = ssGetInputPortRealSignalPtrs(S,0);  
  
  
  real_T va, vb, vc;
  real_T valfa, vbeta, vsd, vsq, theta;
  real_T isd_dot, isq_dot, isd, isq;
  real_T wr_dot, te, tl;
  real_T ia, ib, ic, wr, ialfa, ibeta, theta_dot;
  
  /* transformation constant */
  real_T K	= 0.8164966;
  real_T L	= 0.866025;
  
 /* PMSM parameters */
  real_T N = 5.0;
  real_T psi = 0.0116;
  real_T Lsd = 0.193e-3;
  real_T Lsq = 0.193e-3;
  real_T Rs = 0.15;
  real_T J = 0.015;//0.01;
  
  va = U(0);
  vb = U(1);
  vc = U(2);
  
  tl = U(3);
  
  wr = X[0];
  isd = X[1];
  isq = X[2];
  theta = X[3];
  
  valfa   = K*(va-0.5*vb-0.5*vc);  //transformasi clarke
  vbeta  = K*L*(vb-vc);

  vsd	=   valfa*cos(theta)+vbeta*sin(theta); //transformasi park
  vsq	=  -valfa*sin(theta)+vbeta*cos(theta);

  /* PMSM model */
  isd_dot = -(Rs/Lsd)*isd +N*wr*(Lsq/Lsd)*isq +vsd/Lsd;
  isq_dot = -(Rs/Lsq)*isq -N*wr*((Lsd/Lsq)*isd +psi/Lsq) +vsq/Lsq;
  
  /* torque motor */
  te = N*(psi+(Lsd-Lsq)*isd)*isq;
  
  /* Mechanical model */
  wr_dot = (te - tl)/(J);
 
  theta_dot = N*wr;
  
  dX[0] = wr_dot;
  dX[1] = isd_dot;
  dX[2] = isq_dot;
  dX[3] = theta_dot;
}

   
static void mdlTerminate(SimStruct *S)                                         
{} /*Keep this function empty since no memory is allocated*/    
                                                                                                            
#ifdef  MATLAB_MEX_FILE                                                        
/* Is this file being compiled as a MEX-file? */                             
#include "simulink.c"   /* MEX-file interface mechanism */      
#else                                                                                                   
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif      
                                                                                           