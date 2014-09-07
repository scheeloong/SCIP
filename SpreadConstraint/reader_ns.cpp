/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving CoNStraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer InformatioNStechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic LiceNSe.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic LiceNSe              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_xyz.c
 * @brief  XYZ file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "reader_ns.h"
#include "scip/reader_xyz.h"
#include "probdata_ns.h" // to create the problem 

// ns => Nurse Scheduling 

#define READER_NAME             "NSreader"
#define READER_DESC             "Nurse Scheduling Problem file reader"
#define READER_EXTENSION        "NS"

/*
method: 	   	0 -> MILP model
        	   	1 -> MIQP model
			2 -> MIQP model (binary obj func)
        	   	3 -> CIP model
        	   	4 -> CIP model  (binary obj func)
*/

/** copy method for reader plugiNS (called when SCIP copies plugiNS) */
#define readerCopyNS NULL

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeNS NULL

// QUESTION: What is this function called? How is this function called? 

// This function is used to read problem
// It first reads a file to understand the 5 constraints at the first line and 2 constraints at the last line
// Then it reads a 2nd file to understand which method is needed to solve it
// Then it creates the problem and deallocates all memories needed to create the problem. 
// and returns. 
/** problem reading method of reader */
static SCIP_DECL_READERREAD(readerReadNS)
{  

int numX, stdDevLB, stdDevUB, meanLB, meanUB; 

    FILE* file; // To read the file 
    int status;
    *result = SCIP_DIDNOTRUN; // Initialize result to SCIP did not even run 
    int i; 
    /* open problem file */
    file = fopen(filename, "r");
    if( file == NULL )
    {
       SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
       return SCIP_NOFILE;
    }

   // Scan the first line of data, it has 3 numbers, nZones, nNurses, nPatients
   //first line of data
   status = fscanf(file, "%d %d %d %d %d", &numX, &stdDevLB, &stdDevUB, &meanLB , &meanUB);

   int **values;
   values = (int **) malloc(numX * sizeof(int *));
   if(values == NULL)
   {
      fprintf(stderr, "out of memory\n");
   }
   for(i = 0; i < numX; i++)
   {
      values[i] = (int *) malloc(2 * sizeof(int)); // 2 columns each int 
      if(values[i] == NULL)
      {
         fprintf(stderr, "out of memory\n");
      }
   }

   // O(n) 
   // Continue scanning number of patients lines of the file (Guaranteed to have this as assumption) 
   // and scans the zone and acuity for this patient. 
   for (i =0; i < numX; i++)
   {
      status = fscanf(file, "%d %d", &values[i][0], &values[i][1]);
      if ( ! status )
      {
         SCIPerrorMessage("Reading failed.\n");
        return SCIP_READERROR;
      }
   }
  // close the file  for understanding the 5 parameters 
   fclose( file );
   SCIPprobNSdataCreateCIP(scip, filename, numX, stdDevLB, stdDevUB, meanLB, meanUB, values); 
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

// TO include this reader into SCIP 
#define readerWriteNS NULL
/** includes the xyz file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderNS(
   SCIP*                 scip                
   )
{
   SCIP_READERDATA* readerdata;

   /* create xyz reader data */
   readerdata = NULL;
   // TODO: (optional) create reader specific data here 
   
   // include xyz reader 
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyNS, readerFreeNS, readerReadNS, readerWriteNS, readerdata) );
   return SCIP_OKAY;
}

