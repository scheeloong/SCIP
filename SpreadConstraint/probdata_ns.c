#include <string.h>
 
#include "probdata_ns.h"
/** sets up the problem data using CIP */


//--------------------------------------------------------------------------------------------------------
SCIP_RETCODE SCIPprobNSdataCreateCIP
//--------------------------------------------------------------------------------------------------------
( // These are the parameters to this function 
        SCIP*                 	scip,               /**< SCIP data structure */
	const char*           	probname,           /**< problem name */
  	int			N, 
	int			stdDevLB,
	int			stdDevUB,  
	int 			meanLB, 
	int			meanUB,
	int**			values
)
{
   int i, j, k;  
   SCIP_Real value; //this variable is used to pass integer variables to functions like createconslinear, addcoeflinear.. 
	            // DEBUG: for some reason using integer variables doesn't work

   // Create an empty problem which:
   // i) Frees old problem 
   // ii) switch SCIP to SCIP_STAGE_PROBLEM
   // iii) Create solution pool for original solution candidates
   SCIP_CALL(SCIPcreateProb(scip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL)); // defined in <scip.h>
   // set objective sense 
   // Precondition: In SCIP_STAGE_PROBLEM stage
   // Sets the objective to minimize not maximize, it's either one of the two. 	
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   // Limit the amount of memory and time used for solving. 
   SCIP_CALL(SCIPsetRealParam(scip, "limits/memory", 2500));
   SCIP_CALL(SCIPsetRealParam(scip, "limits/time", 3600));
   // Names of all the linear constraints that will be created 
   char name[100];
double a, b; // to fix weird bug 
// 1. Create variable X
   SCIP_VAR* X[N]; // To store the 10 possible X values, each are integer variables   
   for (i = 0; i < N; i++) 
   {
      a = values[i][0]; 
      b = values[i][1]; 
      printf("Values %d is [%f, %f] \n", i, a, b); 
      SCIP_VAR *varX; 
      snprintf(name, sizeof name, "X_%d", i); 
	// Set both lower and upper bound the same value, values[i] 
     SCIP_CALL(SCIPcreateVar(scip, & varX, name, a,b, 1.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE,NULL, NULL, NULL, NULL, NULL)); 
     SCIP_CALL(SCIPaddVar(scip, varX)); 
printf("Inside Create Variable X,  %d is %s with [%f, %f] \n", i, SCIPvarGetName(varX), SCIPvarGetLbLocal(varX), SCIPvarGetUbLocal(varX)); 
     X[i] = varX; 
   }


// 2. Create domain variable stdDev
// Create Y variable, stdDev = [stdDevlb, stDevub], stdDev points to var which contains both the lower and upper bound. 
   SCIP_VAR* stdDev;//integer;
   SCIP_VAR *varStdDev;       
   snprintf(name, sizeof name, "stdDev");		// INTEGER => {lb, lb+1, ... , ub} 
   SCIP_CALL(SCIPcreateVar(scip, &varStdDev, name, stdDevLB, stdDevUB, 1.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE,NULL, NULL, NULL, NULL, NULL));
   SCIP_CALL(SCIPaddVar(scip, varStdDev)); 
   stdDev = varStdDev;


// 3. Create domain variable mean
   SCIP_VAR* mean; 
   SCIP_VAR *varMean; 
   snprintf(name, sizeof name, "mean"); 
   // printf("mean's lowerbound is %d\n", meanLB); 
   SCIP_CALL(SCIPcreateVar(scip, &varMean, name, meanLB, meanUB, 1.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL)); 
   SCIP_CALL(SCIPaddVar(scip, varMean)); 
   mean =  varMean;  

// 4. Create Spread Constraint
   SCIP_CONS* cons_spread; 
   snprintf(name, sizeof name, "Spread"); 
   SCIP_VAR* allVars[N+2]; // To store the N possible X values, stdDev, and Mean, each are integer variables   
   for (i = 0; i < N; i++) 
   {
	allVars[i] = X[i];
   }
   allVars[N] = stdDev;
   allVars[N+1] = mean;
   SCIP_CALL (SCIPcreateConsBasicSpread(scip, allVars, &cons_spread, name, N+2)); 
   SCIP_CALL( SCIPaddCons(scip, cons_spread));
   return SCIP_OKAY;
}

