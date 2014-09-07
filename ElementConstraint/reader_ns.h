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

/**@file   reader_xyz.h
 * @ingroup FILEREADERS
 * @brief  XYZ file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_NS_H__
#define __SCIP_READER_NS_H__

#include "scip/scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** includes the NS file reader into SCIP */
EXTERN
SCIP_RETCODE SCIPincludeReaderNS(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
