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

/**@file   cons_Element.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for Element constraints
 * @author Tobias Achterberg
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_Element_H__
#define __SCIP_CONS_Element_H__
#include "scip/scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for Spread constraints and includes it in SCIP */
EXTERN SCIP_RETCODE SCIPincludeConshdlrElement(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a Spread constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
// This function creates the Spread Constaint 
EXTERN SCIP_RETCODE SCIPcreateConsElement
(
// First, create scip data structure, with the variables, and hold created constraint, name of constraint to identify it, number of variables
// max deviation allowed, 
   SCIP* 		 scip, /**< SCIP data structure */
   SCIP_VAR** 		 vars,                /**< variables */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables  including X, Y and Z */
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
   );

/** creates and captures an Element constraint in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsElement(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsElement() for information about the basic constraint flag configuration*/

// Creates a very basic element constraint 
EXTERN SCIP_RETCODE SCIPcreateConsBasicElement(
   SCIP*                 scip,               //< SCIP data structure 
   SCIP_VAR**            vars,               //< array with variables of constraint entries 
   SCIP_CONS**           cons,               //< pointer to hold the created constraint 
   const char*           name,               //< name of constraint 
   int                   nvars              //< number of variables in the constraint 
   );


#ifdef __cplusplus
}
#endif

#endif
