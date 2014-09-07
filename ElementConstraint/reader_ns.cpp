
// NOTE: THIS FILE HAS BEEN CHANGED FOR ELEMENTTEST AND NO LONGER USED FOR NURSE SCHEDULING USING SPREAD 

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
#include <stdio.h> 
#include "reader_ns.h"
#include "scip/reader_xyz.h"
#include "probdata_ns.h" // to create the problem 

// ns => Nurse Scheduling 

#define READER_NAME             "NSreader"
#define READER_DESC             "Nurse Scheduling Problem file reader"
#define READER_EXTENSION        "NS"

// CHANGED: ALWAYS USES THE CIP MODEL 

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
// For now, it only reads 1 file, and assumes method 0 is always used. 
/** problem reading method of reader */
static SCIP_DECL_READERREAD(readerReadNS)
{  
    int i;
    int N ; // Number of values in the set of constants or set of variables X 
    int Ylb ; // the rank of the ordered set of N values that Z has to take (lower bound) 
    int Yub ; // rank of the ordered set of N values that Z has to take (upper bound) 
    int Zlb, Zub ; // Z = aY or Z= xY , a -> const. , x => variable 
    

    FILE* file; // To read the file 
    int status;
    *result = SCIP_DIDNOTRUN; // Initialize result to SCIP did not even run 

    /* open problem file */
    file = fopen(filename, "r");
    if( file == NULL )
    {
       SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
       return SCIP_NOFILE;
    }

   // Scan the first line of data, it has 3 numbers, n ylb yub 
   //first line of data
   status = fscanf(file, "%d %d %d %d %d", &N, &Ylb, &Yub, &Zlb, &Zub);
printf("values gotten are N:%d Ylb: %d, Yub: %d, Zlb: %d, Zub: %d \n", N, Ylb, Yub, Zlb, Zub); 
   // Create integer arrays to store the array of values 
//   int values[N][2];

   int **values;
   values = (int **) malloc(N * sizeof(int *));
   if(values == NULL)
   {
      fprintf(stderr, "out of memory\n");
   }
   for(i = 0; i < N; i++)
   {
      values[i] = (int *) malloc(2 * sizeof(int)); // 2 columns each int 
      if(values[i] == NULL)
      {
         fprintf(stderr, "out of memory\n");
      }
   }


   for (i =0; i < N; i++)
   {
      status = fscanf(file, "%d %d", &values[i][0], &values[i][1]);// X 
      if ( ! status )
      {
         SCIPerrorMessage("Reading failed.\n");
        return SCIP_READERROR;
      }
   }
  // close the file  for understanding the 2 parameters and all the values in the set 
   fclose( file );

  int method = 0; // Instead of opening a second file, just initialize method to 0 here. 
  // Regardless of what method is, create a cip problem 

// The bottom global functions are defined in probdata_ns.h (note: Included above at declaration) 
// If it is pure MILP method, create an MILP problem 
// If it is pure MIQP method, create a MIQP problem 
// If it is a CIP method, create a CIP problem 
   SCIPprobNSdataCreateCIP(scip, filename, N, Ylb, Yub, Zlb, Zub, values, method); 
// Done creating problem, deallocate memory made for values
//   delete [] values;
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

// END OF FUNCTION 


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

