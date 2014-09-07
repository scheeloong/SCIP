
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
	int			Ylb,
	int			Yub,  
	int 			Zlb, 
	int			Zub,
	int**			values,
	int			method
    );

#ifdef __cplusplus
}
#endif
#endif

