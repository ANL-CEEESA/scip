/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_sol.c
 * @brief  public methods for solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Robert Lion Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
#include "scip/exprinterpret.h"
#include "scip/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branchexact.h"
#include "scip/bounding_exact.h"
#include "scip/branch_nodereopt.h"
#include "scip/certificate.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cutpool.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/def.h"
#include "scip/dialog.h"
#include "scip/dialog_default.h"
#include "scip/disp.h"
#include "scip/event.h"
#include "scip/heur.h"
#include "scip/heur_ofins.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heuristics.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/nodesel.h"
#include "scip/paramset.h"
#include "scip/presol.h"
#include "scip/presolve.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/retcode.h"
#include "scip/sepastoreexact.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scipgithash.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"

#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/pub_lpexact.h"
#include "scip/struct_certificate.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** returns whether the solution process should be probably correct
 *
 *  @note This feature is not supported yet!
 *
 *  @return Returns TRUE if \SCIP is exact solving mode, otherwise FALSE
 */
SCIP_Bool SCIPisExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->exact_enabled);
}

/** returns whether aggreagtion is allowed to use negative slack */
SCIP_Bool SCIPallowNegSlack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (!SCIPisExactSolve(scip)) || (scip->set->exact_allownegslack);
}

/** returns which method is used for computing truely valid dual bounds at the nodes ('n'eumaier and shcherbina,
 *  'v'erify LP basis, 'r'epair LP basis, 'p'roject and scale, 'e'xact LP,'i'nterval neumaier and shcherbina,
 *  e'x'act neumaier and shcherbina, 'a'utomatic); only relevant for solving the problem provably correct
 */
char SCIPdualBoundMethod(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->exact_safedbmethod);
}

/** returns whether the certificate output is activated? */
SCIP_Bool SCIPisCertificateActive(
   SCIP*                 scip                /**< certificate information */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   return (scip->stat->certificate != NULL && scip->stat->certificate->transfile != NULL);
}

/** returns certificate data structure
 *
 *  @return certificate data structure
 */
SCIP_CERTIFICATE* SCIPgetCertificate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   return scip->stat->certificate;
}

/** free information that is possibly still stored about this row in the certifacte structure */
SCIP_RETCODE SCIPfreeRowCertInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< a SCIP row */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_AGGREGATIONINFO* aggrinfo;
   SCIP_MIRINFO* mirinfo;

   if( !SCIPisExactSolve(scip) || !SCIPisCertificateActive(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   /* only do something if row does not already exist*/
   if( SCIPhashmapExists(certificate->rowdatahash, (void*) SCIProwGetRowExact(row)) )
      return SCIP_OKAY;

   if( certificate->workingaggrinfo || certificate->workingmirinfo )
      return SCIP_OKAY;

   SCIPdebugMessage("Removing information stored in certificate for row \n");

   if( (certificate->aggrinfohash != NULL) && SCIPhashmapExists(certificate->aggrinfohash, (void*) row) )
   {
      aggrinfo = (SCIP_AGGREGATIONINFO*) SCIPhashmapGetImage(certificate->aggrinfohash, (void*) row);
      SCIP_CALL( SCIPcertificateFreeAggrInfo(scip->set, certificate, scip->lp, aggrinfo, row) );
   }

   if( (certificate->mirinfohash != NULL) && SCIPhashmapExists(certificate->mirinfohash, (void*) row) )
   {
      mirinfo = (SCIP_MIRINFO*) SCIPhashmapGetImage(certificate->mirinfohash, (void*) row);
      SCIP_CALL( SCIPcertificateFreeMirInfo(scip->set, certificate, scip->lp, mirinfo, row) );
   }

   return SCIP_OKAY;
}

/** agg aggregation information to certificate for one row */
SCIP_RETCODE SCIPaddCertificateAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL( SCIPcertificateNewAggrInfo(scip, aggrrow, aggrrows, weights, naggrrows, negslackrows, negslackweights, nnegslackrows) );

   return SCIP_OKAY;
}

/** agg aggregation information to certificate for one row */
SCIP_RETCODE SCIPaddCertificateMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL( SCIPcertificateNewMirInfo(scip) );

   return SCIP_OKAY;
}

/** compute a safe bound for the current lp solution */
SCIP_RETCODE SCIPcomputeSafeBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             proveinfeas,        /**< should infeasibility be proven */
   SCIP_Real*            safebound           /**< safe the resulting bound */
   )
{
   SCIP_Bool lperror = false;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->lp != NULL && scip->lpexact != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcomputeSafeBound", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPlpExactComputeSafeBound(scip->lp, scip->lpexact, scip->set, scip->messagehdlr, SCIPblkmem(scip),
         scip->stat, scip->eventqueue, scip->eventfilter, scip->transprob, scip->lpexact->lpiitlim,
         &lperror, proveinfeas, safebound, NULL, NULL) );

   if( lperror )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/** force the next lp to be solved exactly */
SCIP_RETCODE SCIPforceExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->lpexact != NULL);

   scip->lpexact->forceexactsolve = TRUE;

   return SCIP_OKAY;
}

/** check exact integrality of lp solution */
SCIP_RETCODE SCIPcheckIntegralityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   assert(scip != NULL);
   assert(scip->lp != NULL && scip->lpexact != NULL);

   SCIP_CALL( SCIPlpExactcheckIntegralityExact(scip->lp, scip->lpexact,
         scip->set, scip->stat, result) );

   return SCIP_OKAY;
}

/** branches on an LP solution exactly; does not call branching rules, since fractionalities are assumed to small;
 *  if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchLPexact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchLPexact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecLPexact(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->sepastore, scip->branchcand, scip->eventqueue, scip->primal->cutoffbound,
         TRUE, result) );

   return SCIP_OKAY;
}

/** adds row to exact separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        rowexact            /**< exact row to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsepastoreExactAddCut(scip->sepastoreexact, SCIPblkmem(scip), scip->set, scip->stat, scip->eventqueue,
            scip->lpexact, rowexact) );

   return SCIP_OKAY;
}
