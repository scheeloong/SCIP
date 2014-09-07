/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_Element.c
 * @brief  constraint handler for Element constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "cons_element.h"
#include "scip/cons_linking.h" // DONT NEED THIS
#include "scip/cons_knapsack.h" // DONT NEED THIS
#include "scip/scipdefplugins.h"
//headers used by bound consistency algorithm
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h> // also used for debugging
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>


//fundamental constraint handler properties 
#define CONSHDLR_NAME          "Element"
#define CONSHDLR_DESC          "Element constraint handler"
#define CONSHDLR_SEPAPRIORITY   2100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -3000000 //< priority of the constraint handler for constraint enforcing 
#define CONSHDLR_CHECKPRIORITY -3000000 //< priority of the constraint handler for checking feasibility 
#define CONSHDLR_EAGERFREQ          100 //< frequency for using all instead of only the useful constraints in separation,                                                 propagation and enforcement, -1 for no eager evaluations, 0 for first only 
// Do you need this constraint handler? if no constraint handler available ?
#define CONSHDLR_NEEDSCONS         TRUE //< should the constraint handler be skipped, if no constraints are available? 

// optional constraint handler properties 
// TODO: remove properties which are never used because the corresponding routines are not supported 
#define CONSHDLR_SEPAFREQ            1 //< frequency for separating cuts; zero means to separate only in the root node 
#define CONSHDLR_DELAYSEPA        FALSE //< should separation method be delayed, if other separators found cuts? 
#define CONSHDLR_PROPFREQ            1 //< frequency for propagating domains; zero means only preprocessing propagation 
#define CONSHDLR_DELAYPROP        FALSE //< should propagation method be delayed, if other propagators found reductions? 
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP//< propagation timing mask of the constraint handler
#define CONSHDLR_MAXPREROUNDS        -1 //< maximal number of presolving rounds the constraint handler participates in (-1: no limit) 
#define CONSHDLR_DELAYPRESOL      FALSE //< should presolving method be delayed, if other presolvers found reductions? 

// (optional) enable linear or nonlinear constraint upgrading 
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"

/* default parameter values */

/* separation */
#define DEFAULT_USEBINVARS             FALSE //< should the binary representation be used? 
#define DEFAULT_LOCALCUTS              FALSE //< should cuts be added only locally? 
#define DEFAULT_USECOVERCUTS            TRUE //< should covering cuts be added? 
#define DEFAULT_CUTSASCONSS             TRUE //< should the cuts be created as knapsack constraints? 
#define DEFAULT_SEPAOLD                 TRUE //< shall old sepa algo be applied? 

// propagation 
#define DEFAULT_CORETIMES               TRUE //< should core-times be propagated (time tabling)? 
#define DEFAULT_OVERLOAD                TRUE //< should edge finding be used to detect an overload? 
#define DEFAULT_EDGEFINDING             TRUE //< should edge-finding be executed? 
#define DEFAULT_USEADJUSTEDJOBS        FALSE //< should during edge-finding jobs be adusted which run on the border of the effective time horizon? 

// presolving 
#define DEFAULT_DUALPRESOLVE            TRUE //< should dual presolving be applied? 
#define DEFAULT_COEFTIGHTENING         FALSE //< should coefficient tightening be applied? 
#define DEFAULT_NORMALIZE               TRUE //< should demands and capacity be normalized? 
#define DEFAULT_MAXNODES             10000LL //< number of branch-and-bound nodes to solve an independent Element constraint  (-1: no limit) 

// enforcement 
#define DEFAULT_FILLBRANCHCANDS        FALSE //< should branching candidates be added to storage? 

// conflict analysis 
#define DEFAULT_USEBDWIDENING           TRUE //< should bound widening be used during conflict analysis? 

#define EVENTHDLR_NAME         "Element"
#define EVENTHDLR_DESC         "bound change event handler for Element constraints"

const double numerical_gap = 0.000001;

// Data structure used in the Element constraint

// struct ... 

// data for Element constraints 
// Element(X, Y, Z) 
struct SCIP_ConsData
{
   // All your variables, its mean and the allowed maximum deviation 
   SCIP_VAR**     	 vars;         		//< array of variables for X. Y and Z
   int		 	 nvars;			//number of variables for X , Y and Z
};

// Data for constraint handler during Propagation for Element Constraint
// constraint handler data
// NOTE: IGNORE ALL OF THESE
// NOTE: DID NOT CHANGE AT ALL 
struct SCIP_ConshdlrData
{
// To handle bound change events 
   SCIP_EVENTHDLR*       eventhdlr;          //< event handler for bound change events 
// Do you use binary variables? 
   SCIP_Bool             usebinvars;         //< should the binary variables be used? 
   SCIP_Bool             edgefinding;        //< should edge-finding be executed? 
   SCIP_Bool             localcuts;          //< should cuts be added only locally? 
   SCIP_Bool             usecovercuts;       //< should covering cuts be added? 
   SCIP_Bool             sepaold;            //< shall old sepa algo be applied? 
   SCIP_Bool             fillbranchcands;    //< should branching candidates be added to storage? 
   SCIP_Bool             dualpresolve;       //< should dual presolving be applied? 
   SCIP_Bool             coeftightening;     //< should coeffisient tightening be applied? 
   SCIP_Bool             normalize;          //< should demands and capacity be normalized? */
   SCIP_Bool             usebdwidening;      //< should bound widening be used during conflict analysis? */
   SCIP_Longint          maxnodes;           //< number of branch-and-bound nodes to solve an independent disjunctive constraint  (-1: no limit) */
   SCIP_Bool             cutsasconss;        //< should the disjunctive constraint create cuts as knapsack constraints? */
   SCIP_Bool             coretimes;          //< should core-times be propagated (time tabling)? */
   SCIP_Bool             overload;           //< should edge finding be used to detect an overload? */
   SCIP_Bool             useadjustedjobs;    //< should during edge-finding jobs be adusted which run on the border of the effective time horizon? */
   SCIP_NODE*            lastenfolpnode;     /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenfolprounds;      /**< counter on number of enforcement rounds for the current node */
};

// Element Constraint methods 
/** Local methods */

// This function creates a constraint data  (sort of like the constructor but in C, not C++)
// used by (SCIPcreateConsElement() below  ) 
static SCIP_RETCODE consdataCreate (
// Note: These are the parameters to this function 
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to consdata for element constraint*/
   SCIP_VAR**            vars,               /**< array of integer variables */
   int 			 nvars
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(vars != NULL);	
   assert(nvars > 0); 

   // create constraint data
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   
   // Assign parameter values to actual data 
   (*consdata)->nvars = nvars; 
   (*consdata)->vars = vars;
printf("Nvars is  %d \n", nvars); 
int i = 0; 
for (i = 0; i < nvars; i++)
// Note: To print the values, always print with %f and not %d 
	printf("vars %d is %s with [%f, %f] \n", i, SCIPvarGetName((*consdata)->vars[i]), SCIPvarGetLbLocal((*consdata)->vars[i]), SCIPvarGetUbLocal((*consdata)->vars[i])); 
// DEBUGGING FUNCTIONS!!!!!   

   if( nvars > 0 )
   {
      assert(vars != NULL); /* for flexelint */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );

      // You need this (W.Y.) 
      /* transform variables, if they are not yet transformed */
      if( SCIPisTransformed(scip) )
      {
         SCIPdebugMessage("get tranformed variables and constraints\n");
         /* get transformed variables and do NOT capture these */
         SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
      }
   }
   else
   {
printf("Consdata is null\n"); 
      (*consdata)->vars = NULL;
   }
   return SCIP_OKAY;
}

// Note: Removed relaxation for mean, not sure if need relaxation for Element Constraint
//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
// Constructor for Element Constraint handler
/** creates constaint handler data for Element constraint handler */
static SCIP_RETCODE conshdlrdataCreate
(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
)
{
   /* create precedence constraint handler data */
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* set event handler for checking if bounds of start time variables are tighten */
   (*conshdlrdata)->eventhdlr = eventhdlr;

#ifdef SCIP_STATISTIC
   (*conshdlrdata)->nirrelevantjobs = 0;
   (*conshdlrdata)->nalwaysruns = 0;
   (*conshdlrdata)->ndualfixs = 0;
   (*conshdlrdata)->nremovedlocks = 0;
   (*conshdlrdata)->ndecomps = 0;
   (*conshdlrdata)->nallconsdualfixs = 0;
#endif
   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
// This function is used to catch events
static SCIP_RETCODE consdataCatchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< Element constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   int v;
   /* catch event for every single variable */
   for( v = 0; v < consdata->nvars; ++v )
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[v],
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );

 /*  SCIP_CALL( SCIPcatchVarEvent(scip, consdata->Ylb,
            SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, (SCIP_EVENTDATA*)consdata, NULL) );*/
   return SCIP_OKAY;
}
//---------------------------------------------------------------------------------------------------------------------
// This function is used to create the Element Constraint (used by prob_data_ns.c to create Element Constraint)
SCIP_RETCODE SCIPcreateConsElement(

   SCIP* 		 scip, /**< SCIP data structure */
   SCIP_VAR** 		 vars,                /**< variables */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
)
{ 
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(vars != NULL);
   assert(cons != NULL);
   assert(nvars > 0);

   // find the Element constraint handler
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage(""CONSHDLR_NAME" constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   printf("create Element constraint <%s> with %d variables\n", name, nvars);
   SCIPdebugMessage("create Element constraint <%s> with %d variables\n", name, nvars);

   int i;	
   for(i=0;i<nvars;i++)
   {
      printf("Element constraint var %s with [%f, %f] \n", SCIPvarGetName(vars[i]), SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]));
   }

   // create constraint data , defined above		
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nvars)); 

   //  create constraint (defined in <scip.h>)
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata,
         initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) ); 
// Note: Did not create relaxation for mean here 

   // Check if in correct stage, if it is, do catch events. 
   // Note: NOTE: DID NOT CHANGE AT ALL 
   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( consdataCatchEvents(scip, consdata, conshdlrdata->eventhdlr) );
   }
   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyElement)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrElement(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopyElement NULL
#endif
//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
/** frees constraint handler data for logic or constraint handler */
static void conshdlrdataFree
(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
)
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);
}

//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static SCIP_DECL_CONSFREE(consFreeElement)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   conshdlrdataFree(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);
   return SCIP_OKAY;
}

// NOTE: DID NOT CHANGE AT ALL 
/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   printf("\n\nCONSInitPre");
   return SCIP_OKAY;
}
#else
#define consInitpreElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreElement)
{  
   printf("\n\nCONSExitPre");

   return SCIP_OKAY;
}
#else
#define consExitpreElement NULL
#endif

//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
//#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolElement)
{  /*lint --e{715}*/
/*   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); */ /*lint --e{527}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->lastenfolpnode = NULL;
   conshdlrdata->nenfolprounds = 0;
   return SCIP_OKAY;
}
/*#else
#define consInitsolElement NULL
#endif*/

// NOTE: DID NOT CHANGE AT ALL 
/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolElement)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* free rows */
   //   SCIP_CALL( consdataFreeRows(scip, &consdata) );
   }

   return SCIP_OKAY;
}
#else
#define consExitsolElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpElement NULL
#endif

//---------------------------------------------------------------------------------------------------------------------
// NOTE: DID NOT CHANGE AT ALL 
// separation method of constraint handler for LP solutions 
static
SCIP_DECL_CONSSEPALP(consSepalpElement)
{  
  return SCIP_OKAY;
}

// NOTE: DID NOT CHANGE AT ALL 
/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolElement NULL
#endif
//---------------------------------------------------------------------------------------------------------------------
// This function checks that the solution found by other constraints does not violate the element constraint 
// Note: used in SCIP_DECL_CONSCHECK() below 
// (called by CheckConsElement ()) 
// It checks if: 
// 1. 
// 2. 
// returns SCIP_OKAY if current solution did not violate Element Constraint (defined in <scip.h>) 
// otherwise, returns a suitable error (refer to <scip.h>), note: past student just return SCIP_OKAY always but changes
// violated = 1 to mean the constraint was violated and violated = 0 to mean that the constraint was not violated 
static SCIP_RETCODE CheckCons
(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   int&			 violated	            /**< pointer to store whether the constraint is violated */
)
{ 
   assert( consdata != NULL );
   violated = 0; // Initialize as not violated 

   SCIP_VAR** vars;
   vars = consdata->vars; // get the X values 
   int nvars = consdata->nvars; 
   assert( consdata->nvars != NULL );	

// Get the solution values 
   int Zsol = SCIPgetSolVal(scip,sol,vars[nvars-1]); 
   int Ysol = SCIPgetSolVal(scip,sol,vars[nvars-2]);
   int Xsol[consdata->nvars-2]; 

   for(int i = 0; i < (int)consdata->nvars-2; i++ )
   	Xsol[i] = SCIPgetSolVal(scip, sol, vars[i]); // SCIPgetSolVal defined in <scip.h>

   // Note: SCIPgetSolVal gets value of the current solution 

// Currently assumes Zsol is the right value. 
// Initialize as violated
   violated = 1; 

printf("Zsol: %d  Ysol:%d Xsol[Ysol]: %d\n", Zsol, Ysol, Xsol[Ysol]);  
     
   // for every current Y value in the solution bound
// if Z is one of these values, return SCIP_OKAY with no violation
   if (Zsol == Xsol[Ysol]) // Note: Assumes Y starts from 0
   {
      violated = 0; // not violated anymore
      return SCIP_OKAY; 
   }
   // If here, means it is violated 
   return SCIP_OKAY;
}


// NOTE: DID NOT CHANGE AT ALL 
// constraint enforcing method of constraint handler for LP solutions 
#if 0
static
SCIP_DECL_CONSENFOLP(consEnfolpElement)
{  


   return SCIP_OKAY;

}
#else
#define consEnfolpElement NULL
#endif

// NOTE: DID NOT CHANGE AT ALL 
// constraint enforcing method of constraint handler for pseudo solutions 
#if 0
static
SCIP_DECL_CONSENFOPS(consEnfopsElement)
{  


   return SCIP_OKAY;
}
#else
#define consEnfopsElement NULL
#endif

// feasibility check method of constraint handler for solutions 
static SCIP_DECL_CONSCHECK(consCheckElement)
{  
   int violated = 0; // Initialize with no violation 
   int c;
   SCIPdebugMessage("\nConsCheck");   

   assert( conss != NULL );
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   for (c = 0; c < nconss; ++c)
   {
	 SCIP_CONSDATA* consdata;
         consdata = SCIPconsGetData(conss[c]);
         SCIP_CALL ( CheckCons(scip, consdata, sol, violated) ); // Check if violated by passing by reference. 
	 if (violated != 0)  break;
   }
   if( violated != 0 )
   {
      *result = SCIP_INFEASIBLE;
      if (violated == 1) //violated something 
        SCIPinfoMessage(scip, NULL, "Violated Element constraint\nRefer to static SCIP_DECL_CONSCHECK(consCheckElement) function \n");
   }
   else
      *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


//-------------------------------------------------------------------------------------------------------------------------------
// Propagation function for element constraint 
/** domain propagation method of constraint handler */
// nconss -> prolly global variable, just  use it. 
// conss  -> prolly global variable, just use it. 

// Functions to use:
// SCIPconsGetData() 
// SCIPvarGetLbLocal() // <pub_var.h> 
// SCIPvarGetUbLocal() // <pub_var.h> 
// 	   *result = SCIP_CUTOFF;
// 	   *result = SCIP_REDUCEDDOM;
// SCIPtightenVarUb() // <scip.h> 
//SCIP_RETCODE SCIPtightenVarUb (SCIP* scip, SCIP_VAR* var, SCIP_Real newbound, SCIP_Bool force, SCIP_Bool* infeasible, SCIP_Bool* tightened ) 	
// SCIP_CALL ( SCIPtightenVarUb	(scip, vars[k], iv[k].max, TRUE, &infeasible, NULL ) );
// SCIP_CALL ( SCIPtightenVarLb	(scip, vars[k], iv[k].min, TRUE, &infeasible, NULL ) );
// 

// double globalD, global variable goes into constraint data. 
static SCIP_DECL_CONSPROP(consPropElement)
{
printf(" INSIDE PROPAGATOR!!!\n"); 
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
  
/* // TO RETURN WITHOUT CHANGING THE BOUNDS 
   *result = SCIP_CUTOFF;
   return SCIP_OKAY;
*/

   SCIP_Bool infeasible; // to check if solution is feasible 

   (*result) = SCIP_DIDNOTFIND; // Initialize result as no result found 
// loop through all Element constraints 
// Note: nconss will be 1 since only have 1 Element Constraint for now. 
// May have multiple element constraints for the same problem, so have to go through each of them. 
   for (int c = 0; c < nconss; ++c) 
   {
     	SCIP_CONSDATA* consdata;
	SCIP_CONS* cons;
 	SCIP_VAR** vars;
     	cons = conss[c]; // Get the current constraint  in the loop 
     	assert( cons != NULL );	 // Make sure it is not NULL 
     	consdata = SCIPconsGetData(cons); // Get it's data 
	vars = consdata->vars; 
	// Note: Need to distinguish the usual X constraint with Y and Z constraint by not looping through the final two constraints. 
	// STEPS TODO: 
	// GET X1, ..., Xn, Ylb, Yub, Zlb, Zub 
	double Ylb, Yub, Zlb, Zub; 	
	Zlb = SCIPvarGetLbLocal(consdata->vars[consdata->nvars-1]); 
	Zub = SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1]); 
	Ylb = SCIPvarGetLbLocal(consdata->vars[consdata->nvars-2]);
	Yub = SCIPvarGetUbLocal(consdata->vars[consdata->nvars-2]);
		// printf("CHECKING Propagation is: %f \n", 	SCIPvarGetLbLocal(consdata->vars[consdata->nvars-1])); 
 
	// Keep booleans of whether or not Ylb, and Yub or Zlb, Zub changes. 		
	bool yLB, yUB, zLB, zUB; // to know if these changes 
	yLB = true; // yLB did change
	yUB = true; 
	zLB = true; 
	zUB = true; // assume all values did change before this
	double Xmin, Xmax; // to prune out Z values 
	double n; // to know number of X variables to prune out Y values, but Xn may be infeasible totally, so need update n. 		
	n = consdata->nvars -2 ; 
	double **X; // malloc enough memory for all X[i][0] for Xilb, and X[i][1] for Xiub 
	int i; // for for loop 	
int a,b,c,d; // for debugging , remove later, temporary
a = b = c = d = 0;  
	while (yLB || yUB || zLB || zUB)
	{
// TEMP DEBUGGING START
printf("in while Loop\n");
for (i = 0; i < n; i++)
	printf("X[%d] Lb = %f  Ub = %f \n", i,  SCIPvarGetLbLocal(consdata->vars[i]),SCIPvarGetUbLocal(consdata->vars[i])); 
printf("Ylb: %f Yub: %f Zlb: %f Zub: %f \n", Ylb, Yub, Zlb, Zub); 
if(yLB == true) a = 1 ;
if(yUB == true) b = 1 ;
if(zLB == true) c = 1 ;
if(zUB == true) d = 1 ;
printf("yLB: %d, yUB: %d, zLB: %d, zUB: %d\n", a, b, c, d); 
a = b = c = d = 0;  
// TEMP DEBUGGING END 
	    yLB = false; yUB = false; zLB = false; zUB = false; // initialize all to nothing changed 
	    i = Ylb; 
	    if ( i < n && i >= 0) // If i is one of the values available on n
	    {
printf("Enter  normal case\n"); 
		    // Initialize minimum and maximum of all the available X values to prune Z
		    Xmin = SCIPvarGetLbLocal(consdata->vars[i]); //X[i][0]; 
		    Xmax = SCIPvarGetUbLocal(consdata->vars[i]); //X[i][1]; 
	    }
	    else 
	    {
		printf("Enter  special case\n"); 
// NOWSS // for quick reference to this file 
 		  *result = SCIP_CUTOFF;
 		  return SCIP_OKAY;
	    } 
	    if (Ylb == Yub) // for case where Ylb == Yub which means can change value of X 
	    {
		int j = Ylb; 
		if (SCIPvarGetLbLocal(consdata->vars[j]) < Zlb) // X[i][0] < Zlb
		{
		   if (Zlb <= SCIPvarGetUbLocal(consdata->vars[j]))
		   {
		      SCIP_CALL(SCIPtightenVarLb(scip, consdata->vars[j], Zlb, TRUE, &infeasible, NULL )); //  X[i][0] = Zlb; 
		   }
		   else
		   {	
			printf("Special case 2\n"); // Debug
	  		 *result = SCIP_CUTOFF;
	 		  return SCIP_OKAY;
		   }
		}
		if (SCIPvarGetUbLocal(consdata->vars[j])< Zub)  // X[i][1] < Zub
		{
		   if (Zub >= SCIPvarGetLbLocal(consdata->vars[j]))
		   {
		     SCIP_CALL(SCIPtightenVarLb(scip, consdata->vars[j], Zub, TRUE, &infeasible, NULL )); // X[i][1] = Zub; 
		   }
		   else
		   {	
			printf("Special case 3\n"); // Debug
     		        *result = SCIP_CUTOFF;
	  		return SCIP_OKAY;
		   }
		}
	    }
	    for ( i = Ylb; i <= Yub; i++)
	    {
		// If no more X values to loop through 			// NOTE: if Y = 0 => Z = X0, => X starts from 0!!!!! 
		if (i >= n)	// Yub = 10 and X = 10 , n = number of variables in X -> X0,X1,...,X9 => n = 10
		{
		   // 
		   Yub = n-1; 
		   SCIP_CALL ( SCIPtightenVarUb	(scip, consdata->vars[consdata->nvars-2], Yub, TRUE, &infeasible, NULL )); 
		   yUB = true;
		   break; 
		}


		if (SCIPvarGetLbLocal(consdata->vars[i]) < Xmin) // X[i][0] < Xmin
		{
		    Xmin = SCIPvarGetLbLocal(consdata->vars[i]); 
		}
		if (SCIPvarGetUbLocal(consdata->vars[i]) > Xmax) // X[i][1] < Xmax
		{
		    Xmax = SCIPvarGetUbLocal(consdata->vars[i]); 
		}
	    }	 // End of for loop around all the X's 
	    if (Xmin > Zlb)
	    {
		Zlb = Xmin; 
		SCIP_CALL(SCIPtightenVarLb(scip, consdata->vars[consdata->nvars-1], Xmin, TRUE, &infeasible, NULL )); 
		zLB = true; 
	    }
	    if (Xmax < Zub)
	    {
		Zub = Xmax; 
		SCIP_CALL (SCIPtightenVarUb	(scip, consdata->vars[consdata->nvars-1], Xmax, TRUE, &infeasible, NULL )); 
		zUB = true; 
	    }
	}	// End of while loop
// Print final values
printf("The final assigned values are: \n"); 
for (i = 0; i < n; i++)
	printf("X[%d] Lb = %f  Ub = %f \n", i,  SCIPvarGetLbLocal(consdata->vars[i]),SCIPvarGetUbLocal(consdata->vars[i])); 
printf("Ylb: %f Yub: %f Zlb: %f Zub: %f \n", Ylb, Yub, Zlb, Zub); 
	// Here, finished pruning the first Element constraint, may loop to the 2nd element constraint if it exists. 
   } // End of for loop 
      // Exit all the Element Constraints  that exist
      return SCIP_OKAY; 
} // End of propagation function 





//-------------------------------------------------------------------------------------------------------------------------------
/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolElement NULL
#endif


/** propagation conflict resolving method of constraint handler */
//#if 0
static
SCIP_DECL_CONSRESPROP(consRespropElement)
{  /*lint --e{715}*/
   SCIPdebugMessage("\nConsresprop");
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
//#else
//#define consRespropElement NULL
//#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockElement)
{  
    SCIP_CONSDATA* consdata;
    SCIP_VAR** vars;

    SCIPdebugMessage("lock Element constraint <%s> with nlockspos = %d, nlocksneg = %d\n", SCIPconsGetName(cons), nlockspos, nlocksneg);

    assert(scip != NULL);
    assert(cons != NULL);

    SCIPdebugMessage("\nConsLock");


    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    int i;

    vars = consdata->vars;
    assert(vars != NULL);

    /* In a general way, a variable can make the solution infeasible if rounded both up and down. Here, we tell scip not to round variables in both ways. */
    for( i = 0; i < consdata->nvars; i++)
    {
        SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
    }
	// TODO: FIXME: No idea about this yet. (W.Y and Soon) 
    //SCIP_CALL( SCIPaddVarLocks(scip, consdata->max_deviation, nlockspos, nlocksneg) ); // Question: No idea if this is right
    return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   return SCIP_OKAY;
}


/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   return SCIP_OKAY;
}

/** variable deletion of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   return SCIP_OKAY;
}

/**/
static
SCIP_RETCODE consdataPrint (
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
   FILE*                 file               /**< output file (or NULL for standard output) */
   )
{
   assert(consdata != NULL);

/* print coefficients */
   SCIPinfoMessage( scip, file, "Element({");
   for( int v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars[v] != NULL);
      if( v > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(consdata->vars[v]) );
   }
   SCIPinfoMessage( scip, file, "Z with Bound [ %f , %f ]",SCIPvarGetLbLocal(consdata->vars[consdata->nvars-1]),SCIPvarGetUbLocal(consdata->vars[consdata->nvars-1])); 
   return SCIP_OKAY;

}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintElement)
{ 
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   consdataPrint(scip, SCIPconsGetData(cons), file);
   return SCIP_OKAY;

}


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyElement)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyElement NULL
#endif

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseElement)
{  /*lint --e{715}*/
   SCIPdebugMessage("\nConsParse");
   SCIPerrorMessage("method of Element constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsElement)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }
   return SCIP_OKAY;}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsElement)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   (*nvars) = consdata->nvars;
   (*success) = TRUE;
   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecElement)
{  //lint --e{715}
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);
   return SCIP_OKAY;
}


// Note: This is called in cmain.c so SCIP knows i have this constraint handler
/** creates the handler for Element constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrElement(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;
   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecElement, NULL) );
   /* create Element constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpElement, consEnfopsElement, consCheckElement, consLockElement,
         conshdlrdata) );


   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyElement, consCopyElement) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteElement) );
//#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreElement) );
//#endif
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolElement) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeElement) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsElement) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsElement) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreElement) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpElement) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseElement) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolElement, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintElement) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropElement, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropElement) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpElement, consSepasolElement, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransElement) );

   /* add Element constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/coretimes", "should core-times be propagated (time tabling)?",
         &conshdlrdata->coretimes, FALSE, DEFAULT_CORETIMES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/overload", "should edge finding be used to detect an overload?",
         &conshdlrdata->overload, FALSE, DEFAULT_OVERLOAD, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/edgefinding", "should edge-finding be executed?",
         &conshdlrdata->edgefinding, FALSE, DEFAULT_EDGEFINDING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/useadjustedjobs", "should edge-finding be executed?",
         &conshdlrdata->useadjustedjobs, FALSE, DEFAULT_USEADJUSTEDJOBS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usebinvars", "should the binary representation be used?",
         &conshdlrdata->usebinvars, FALSE, DEFAULT_USEBINVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/localcuts", "should cuts be added only locally?",
         &conshdlrdata->localcuts, FALSE, DEFAULT_LOCALCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usecovercuts", "should covering cuts be added every node?",
         &conshdlrdata->usecovercuts, FALSE, DEFAULT_USECOVERCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/cutsasconss",
         "should the Element constraint create cuts as knapsack constraints?",
         &conshdlrdata->cutsasconss, FALSE, DEFAULT_CUTSASCONSS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/sepaold",
         "shall old sepa algo be applied?",
         &conshdlrdata->sepaold, FALSE, DEFAULT_SEPAOLD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/fillbranchcands", "should branching candidates be added to storage?",
         &conshdlrdata->fillbranchcands, FALSE, DEFAULT_FILLBRANCHCANDS, NULL, NULL) );

   /* presolving parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/dualpresolve", "should dual presolving be applied?",
         &conshdlrdata->dualpresolve, FALSE, DEFAULT_DUALPRESOLVE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/coeftightening", "should coefficient tightening be applied?",
         &conshdlrdata->coeftightening, FALSE, DEFAULT_COEFTIGHTENING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/normalize", "should demands and capacity be normalized?",
         &conshdlrdata->normalize, FALSE, DEFAULT_NORMALIZE, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/"CONSHDLR_NAME"/maxnodes",
         "number of branch-and-bound nodes to solve an independent Element constraint (-1: no limit)?",
         &conshdlrdata->maxnodes, FALSE, DEFAULT_MAXNODES, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   /* conflict analysis parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/usebdwidening", "should bound widening be used during the conflict analysis?",
         &conshdlrdata->usebdwidening, FALSE, DEFAULT_USEBDWIDENING, NULL, NULL) );

   return SCIP_OKAY;
}

//---------------------------------------------------------------------------------------------------------------------
 
/** creates and captures an Element constraint in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsElement(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsElement() for information about the basic constraint flag configuration*/
SCIP_RETCODE SCIPcreateConsBasicElement(
   SCIP*                 scip,               //< SCIP data structure 
   SCIP_VAR**            vars,               //< array with variables of constraint entries 
   SCIP_CONS**           cons,               //< pointer to hold the created constraint 
   const char*           name,               //< name of constraint 
   int                   nvars              //< number of variables in the constraint 

   )
{
   SCIP_CALL( SCIPcreateConsElement(scip, vars, cons, name, nvars, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   return SCIP_OKAY;
}


