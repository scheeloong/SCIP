/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// This file basically initializes the constaints needed and starts running the SCIP shell 


/**@file   cmain.c
 * @brief  Main file for disjunctive constraint handler example
 * @author Stefan Heinz
 *
 *  This the file contains the \ref main() main function of the projects. This includes all the default plugins of
 *  \SCIP and the once which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */
#include <stdio.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "cons_spread.h"
#include "reader_myreader.h"
#include "reader_ruler.h"
#include "reader_csp.h"
#include "reader_ns.h"

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   /* include JSP reader */

   /* include GR reader */
//   SCIPincludeReaderGR(scip);

   /* include CSP reader */
//   SCIPincludeReaderCSP(scip);

   /* include NS reader */
   SCIPincludeReaderNS(scip);

// Used to ensure no nurse is in the same two zones
  /* include binpacking pricer  */
   //SCIP_CALL( SCIPincludeConshdlrDisjunctive(scip) );

// Used to ensure 
   /* include AllDiff Constraint Handler */
//   SCIP_CALL( SCIPincludeConshdlrAllDiff(scip) );

// Used to limit each nurse with the maximum number of patients to handle 
   /* include Global Cardinality Constraint Handler */
//   SCIP_CALL( SCIPincludeConshdlrGlobalCardinality(scip) );


// Used to get
   /* include Spread Constraint Handler */
   SCIP_CALL( SCIPincludeConshdlrSpread(scip) );

// Get the 13 available linear constraints 
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();
   return SCIP_OKAY;
}

int
main(int argc,char** argv)
{
   SCIP_RETCODE retcode;
   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }
   return 0;
}
