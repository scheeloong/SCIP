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

// For quick opening all of these in gedit 
// gedit cmain.c  cons_element.cpp cons_element.h probdata_ns.c probdata_ns.h reader_ns.h reader_ns.cpp &

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
#include "cons_element.h" 
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
 //  SCIPincludeReaderXyz(scip);
   /* include GR reader */
  // SCIPincludeReaderGR(scip);

   /* include CSP reader */
 //  SCIPincludeReaderCSP(scip);

   /* include NS reader */
   SCIPincludeReaderNS(scip);

   SCIP_CALL( SCIPincludeConshdlrElement(scip) ); // defined in cons_element.cpp 

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
