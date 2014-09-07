#include <string.h>

// NOTE: THIS FILE HAS BEEN CHANGED FOR ELEMENT CONSTRAINT TESTING, A BACKUP OF ORIGINAL PROBDATA_NS HAS BEEN CREATED
// ONLY SETS UP FOR CIP FOR NOW, THE MILP AND CP ARE NOT USED 
 
#include <stdio.h> 
#include "probdata_ns.h"

// TODO: DO I HAVE TO INCLUDE CONS_SPREAD.H? IN THIS CASE, CONS_ELEMENT.H? 


// Various Variable Types
// 1. SCIP_VARTYPE_BINARY 	
	//binary variable: x in {0,1}
// 2. SCIP_VARTYPE_INTEGER 	
	// integer variable: x in {lb, ..., ub}
// 3. SCIP_VARTYPE_IMPLINT 	
	// implicit integer variable: continuous variable, that is always integral
// 4. SCIP_VARTYPE_CONTINUOUS 	
	// continuous variable: x in [lb,ub]
 

/** sets up the problem data using CIP */
//--------------------------------------------------------------------------------------------------------
SCIP_RETCODE SCIPprobNSdataCreateCIP
//--------------------------------------------------------------------------------------------------------
( // These are the parameters to this function 
    SCIP*               scip,               /**< SCIP data structure */
    const char*         probname,           /**< problem name */
    int			N,
    int			Ylb,
    int			Yub, 
    int			Zlb, 
    int			Zub, 
    int**		values, // X[0] -> X[N]
    int 		method // the method used to solve this problem (CIP , MIQP, or CP) 
    )
{
   int i;  // for forloops

   SCIP_CALL(SCIPcreateProb(scip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
   // set objective sense 
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   SCIP_CALL(SCIPsetRealParam(scip, "limits/memory", 2500));
   SCIP_CALL(SCIPsetRealParam(scip, "limits/time", 3600));
// Note: For SCIPgetSolVal(scip, sol, vars[i]) will get you a fixed value from X, Y , and Z 

// 1. Create variable X
   // TODO: Check if can just use one name and not all 4 different names like what past student did 
   char name[100];
double a, b; // to fix weird bug 
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

// Note: Y and Z will be stored in vars in struct SCIP_ConsData in cons_element.cpp 

// 2. Create domain variable Y 
// Create Y variable, Y = [Ylb, Yub], Y points to var which contains both the lower and upper bound. 
   SCIP_VAR* Y;//integer;
   SCIP_VAR *varY;       
   snprintf(name, sizeof name, "Y");		// INTEGER => {lb, lb+1, ... , ub} 
   SCIP_CALL(SCIPcreateVar(scip, &varY, name, Ylb, Yub, 1.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE,NULL, NULL, NULL, NULL, NULL));
   SCIP_CALL(SCIPaddVar(scip, varY)); 
   Y = varY;

// 3. Create domain variable Z
   SCIP_VAR* Z; 
   SCIP_VAR *varZ; 
   snprintf(name, sizeof name, "Z"); 
printf("Z's lowerbound is %d\n", Zlb); 
   SCIP_CALL(SCIPcreateVar(scip, &varZ, name, Zlb, Zub, 1.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL)); 
   SCIP_CALL(SCIPaddVar(scip, varZ)); 
   Z =  varZ;  
// 4. Create Element Constraint
   SCIP_CONS* cons_element; 
   snprintf(name, sizeof name, "Element"); 
// TODO: Call the function once you implemented it 

 SCIP_VAR* allVars[N+2]; // To store the 10 possible X values, each are integer variables   
 for (i = 0; i < N; i++) 
   {
	allVars[i] = X[i];
   }
   allVars[N] = Y;
   allVars[N+1] = Z;
   SCIP_CALL (SCIPcreateConsBasicElement(scip, allVars, &cons_element, name, N+2));  // FIXME: Segmentation Fault Occurs here 
   SCIP_CALL( SCIPaddCons(scip, cons_element));
   return SCIP_OKAY;
}
