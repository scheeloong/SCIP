
#ifndef __SCIP_PROBDATA_NS__
#define __SCIP_PROBDATA_NS__

#include "scip/scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */

// Creates the CIP prolem dat 
EXTERN SCIP_RETCODE SCIPprobNSdataCreateCIP(
        SCIP*                 	scip,               /**< SCIP data structure */
	const char*           	probname,           /**< problem name */
  	int			N, 
	int			stdDevLB,
	int			stdDevUB,  
	int 			meanLB, 
	int			meanUB,
	int**			values
    );

#ifdef __cplusplus
}
#endif
#endif
