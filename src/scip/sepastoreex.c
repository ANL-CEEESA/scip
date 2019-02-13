/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepastoreex.c
 * @brief  internal methods for storing separated exact cuts
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/sepastoreex.h"
#include "scip/event.h"
#include "scip/sepa.h"
#include "scip/cons.h"
#include "scip/debug.h"
#include "scip/scip.h"
#include "scip/cuts.h"
#include "scip/struct_sepastore.h"
#include "scip/misc.h"
#include "scip/lpex.h"
#include "scip/rational.h"
#include "scip/pub_lpex.h"

/** resizes cuts and score arrays to be able to store at least num entries */
static
SCIP_RETCODE sepastoreexEnsureCutsMem(
   SCIP_SEPASTOREEX*     sepastoreex,        /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(sepastoreex != NULL);
   assert(set != NULL);

   if( num > sepastoreex->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastoreex->cuts, newsize) );
      sepastoreex->cutssize = newsize;
   }
   assert(num <= sepastoreex->cutssize);

   return SCIP_OKAY;
}

/** creates separation storage */
SCIP_RETCODE SCIPsepastoreexCreate(
   SCIP_SEPASTOREEX**    sepastoreex,          /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepastoreex != NULL);

   SCIP_ALLOC( BMSallocMemory(sepastoreex) );

   (*sepastoreex)->cuts = NULL;
   (*sepastoreex)->cutssize = 0;
   (*sepastoreex)->ncuts = 0;
   (*sepastoreex)->ncutsfound = 0;
   (*sepastoreex)->ncutsfoundround = 0;
   (*sepastoreex)->ncutsapplied = 0;

   (*sepastoreex)->initiallp = FALSE;

   return SCIP_OKAY;
}

/** frees separation storage */
extern
SCIP_RETCODE SCIPsepastoreexFree(
   SCIP_SEPASTOREEX**    sepastoreex,          /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(sepastoreex != NULL);
   assert(*sepastoreex != NULL);
   assert((*sepastoreex)->ncuts == 0);

   BMSfreeMemoryArrayNull(&(*sepastoreex)->cuts);
   BMSfreeMemory(sepastoreex);

   return SCIP_OKAY;
}

/** informs separation storage that the setup of the initial LP starts now */
void SCIPsepastoreexStartInitialLP(
   SCIP_SEPASTOREEX*     sepastoreex         /**< separation storage */
   )
{
   assert(sepastoreex != NULL);
   assert(!sepastoreex->initiallp);
   assert(sepastoreex->ncuts == 0);

   sepastoreex->initiallp = TRUE;
}

/** informs separation storage that the setup of the initial LP is now finished */
void SCIPsepastoreexEndInitialLP(
   SCIP_SEPASTOREEX*     sepastoreex         /**< separation storage */
   )
{
   assert(sepastoreex != NULL);
   assert(sepastoreex->initiallp);
   assert(sepastoreex->ncuts == 0);

   sepastoreex->initiallp = FALSE;
}

/** adds cut to separation storage and captures it */
extern
SCIP_RETCODE SCIPsepastoreexAddCut(
   SCIP_SEPASTOREEX*     sepastoreex,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEX*            lp,                 /**< LP data */
   SCIP_ROWEX*           cut,                /**< separated cut */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut is infeasible */
   )
{
   SCIP_Bool redundant;
   int pos;

   assert(sepastoreex != NULL);
   assert(set != NULL);
   assert(cut != NULL);
   assert(!RisNegInfinity(SCIProwexGetLhs(cut)) || !RisInfinity(SCIProwexGetRhs(cut)));
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);

   /* debug: check cut for feasibility */
   SCIP_CALL( SCIPdebugCheckRow(set, cut) ); /*lint !e506 !e774*/

   /* update statistics of total number of found cuts */
   if( !sepastoreex->initiallp )
   {
      sepastoreex->ncutsfound++;
      sepastoreex->ncutsfoundround++;
   }

   /* check cut for redundancy or infeasibility */
   // todo exip: we will need that probably :redundant = sepastoreIsCutRedundantOrInfeasible(sepastore, set, stat, cut, infeasible);
   /* Note that we add infeasible rows in any case, since we cannot be sure whether the return values are handled
    * correctly. In this way, the LP becomes infeasible. */


   /* get enough memory to store the cut */
   SCIP_CALL( sepastoreexEnsureCutsMem(sepastoreex, set, sepastoreex->ncuts+1) );
   assert(sepastoreex->ncuts < sepastoreex->cutssize);

   SCIPsetDebugMsg(set, "adding cut <%s> to exact separation storage of size %d (forcecut=%u, len=%d)\n",
      SCIProwGetName(cut->fprow), sepastoreex->ncuts, SCIProwGetNNonz(cut->fprow));
   /*SCIP_CALL( SCIPprintRow(set->scip, cut, NULL) );*/

   /* capture the cut */
   // todo exip: might need this SCIProwCapture(cut);

   pos = sepastoreex->ncuts;

   sepastoreex->cuts[pos] = cut;
   sepastoreex->ncuts++;

   /* check, if the row addition to separation storage events are tracked if so, issue ROWADDEDSEPA event */
   /* if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWADDEDSEPA) != 0 )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowAddedSepa(&event, blkmem, cut) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   } */

   /* If the duals need to be collected, then the infeasible flag is set to FALSE. This ensures that the LP is solved */
   if( set->lp_alwaysgetduals && sepastoreex->initiallp )
      (*infeasible) = FALSE;

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
extern
SCIP_RETCODE SCIPsepastoreexApplyCuts(
   SCIP_SEPASTOREEX*     sepastoreex,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEX*            lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_LP* fplp;
   SCIP_ROW** fprows;
   SCIP_ROWEX* rowex;
   int nrowsfp;
   int nrowsex;
   int nreleases;
   int nadded;
   int i;

   fplp = lp->fplp;
   nreleases = 0;
   nadded = 0;

   assert(fplp != NULL);

   fprows = SCIPlpGetRows(fplp);
   nrowsfp = SCIPlpGetNRows(fplp);
   nrowsex = SCIPlexGetNRows(lp);

   assert(fprows != NULL);

   /** this method should sync the fp-lp withe the exact lp */

   /** remove all rows from exact lp that are not in the floating point lp */
   for( i = 0; i < nrowsex; ++i )
   {
      SCIP_ROW* fprow =lp->rows[i]->fprow;
      if( !SCIProwIsInLP(fprow) )
      {
         SCIPlpexReleaseRow(lp->rows[i]);
         nreleases++;
      }
   }

   for( i = 0; i < nrowsfp; ++i )
   {
      rowex = SCIProwGetExRow(lp, lp->rows[i]);
      if( rowex != NULL )
      {
         /* corresponding exact row has already been created, add it to lp */
         assert(...);         
      }
      else
      {
         /* corresponding excact row still needs to be created */
         if( fprows[i]->origintype != SCIP_ROWORIGINTYPE_CONS )
         {
            SCIPerrorMessage("under construction \n");
            SCIPABORT();
         }
         
         SCIProwGetOriginCons(fprow[i]);
      }
   }


   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
extern
SCIP_RETCODE SCIPsepastoreexClearCuts(
   SCIP_SEPASTOREEX*     sepastoreex,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEX*            lp                  /**< LP data */
   )
{
   int c;

   assert(sepastoreex != NULL);

   SCIPsetDebugMsg(set, "clearing %d cuts\n", sepastoreex->ncuts);

   /* release cuts */
   for( c = 0; c < sepastoreex->ncuts; ++c )
   {
      /* check, if the row deletions from separation storage events are tracked if so, issue ROWDELETEDSEPA event */
      /* if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
      {
         SCIP_EVENT* event;

         SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[c]) );
         SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
      } */

      SCIP_CALL( SCIProwexRelease(&sepastoreex->cuts[c], blkmem, set, lp) );
   }

   /* reset counters */
   sepastoreex->ncuts = 0;
   sepastoreex->ncutsfoundround = 0;

   /* if we have just finished the initial LP construction, free the (potentially large) cuts array */
   if( sepastoreex->initiallp )
   {
      BMSfreeMemoryArrayNull(&sepastoreex->cuts);
      sepastoreex->cutssize = 0;
   }

   return SCIP_OKAY;
}


/** indicates whether a cut is applicable
 *
 *  A cut is applicable if it is modifiable, not a bound change, or a bound change that changes bounds by at least epsilon.
 */
extern
SCIP_Bool SCIPsepastoreexIsCutApplicable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             cut                 /**< cut to check */
   )
{
   return TRUE; // todo: exip still do this SCIProwexIsModifiable(cut) || SCIProwexGetNNonz(cut) != 1 || sepastoreIsBdchgApplicable(set, cut);
}

/** get cuts in the separation storage */
extern
SCIP_ROWEX** SCIPsepastoreexGetCuts(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   )
{
   assert(sepastoreex != NULL);

   return sepastoreex->cuts;
}

/** get number of cuts in the separation storage */
extern
int SCIPsepastoreexGetNCuts(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   )
{
   assert(sepastoreex != NULL);

   return sepastoreex->ncuts;
}

/** get total number of cuts found so far */
extern
int SCIPsepastoreexGetNCutsFound(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   )
{
   assert(sepastoreex != NULL);

   return sepastoreex->ncutsfound;
}

/** get number of cuts found so far in current separation round */
extern
int SCIPsepastoreexGetNCutsFoundRound(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   )
{
   assert(sepastoreex != NULL);

   return sepastoreex->ncutsfoundround;
}

/** get total number of cuts applied to the LPs */
extern
int SCIPsepastoreexGetNCutsApplied(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   )
{
   assert(sepastoreex != NULL);

   return sepastoreex->ncutsapplied;
}
