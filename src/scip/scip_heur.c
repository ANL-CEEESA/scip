/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_heur.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for primal heuristic plugins and divesets
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
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

#include "scip/debug.h"
#include "scip/heur.h"
#include "scip/pub_message.h"
#include "scip/scip_heur.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a primal heuristic and includes it in SCIP.
 *
 *  @note method has all heuristic callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeHeurBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPincludeHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   SCIP_HEURTIMING       timingmask,         /**< positions in the node solving loop where heuristic should be executed;
                                              *   see definition of SCIP_HEURTIMING for possible values */
   SCIP_Bool             usessubscip,        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_DECL_HEURCOPY    ((*heurcopy)),      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol)),   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol)),   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   SCIP_HEUR* heur;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeHeur", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether heuristic is already present */
   if( SCIPfindHeur(scip, name) != NULL )
   {
      SCIPerrorMessage("heuristic <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPheurCreate(&heur, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, dispchar, priority, freq, freqofs, maxdepth, timingmask, usessubscip,
         heurcopy, heurfree, heurinit, heurexit, heurinitsol, heurexitsol, heurexec, heurdata) );

   SCIP_CALL( SCIPsetIncludeHeur(scip->set, heur) );

   return SCIP_OKAY;
}

/** creates a primal heuristic and includes it in SCIP with its most fundamental callbacks.
 *  All non-fundamental (or optional) callbacks
 *  as, e. g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetHeurCopy(), SCIPsetHeurFree(),
 *  SCIPsetHeurInit(), SCIPsetHeurExit(), SCIPsetHeurInitsol(), and SCIPsetHeurExitsol()
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeHeur() instead
 *
 */
SCIP_RETCODE SCIPincludeHeurBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   SCIP_HEURTIMING       timingmask,         /**< positions in the node solving loop where heuristic should be executed;
                                              *   see definition of SCIP_HEURTIMING for possible values */
   SCIP_Bool             usessubscip,        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   SCIP_HEUR* heurptr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeHeurBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether heuristic is already present */
   if( SCIPfindHeur(scip, name) != NULL )
   {
      SCIPerrorMessage("heuristic <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPheurCreate(&heurptr, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, dispchar, priority, freq, freqofs, maxdepth, timingmask, usessubscip,
         NULL, NULL, NULL, NULL, NULL, NULL, heurexec, heurdata) );

   assert(heurptr != NULL);

   SCIP_CALL( SCIPsetIncludeHeur(scip->set, heurptr) );

   if( heur != NULL )
      *heur = heurptr;

   return SCIP_OKAY;
}

/* new callback/method setter methods */

/** sets copy method of primal heuristic */
SCIP_RETCODE SCIPsetHeurCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURCOPY    ((*heurcopy))       /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetHeurCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(heur != NULL);

   SCIPheurSetCopy(heur, heurcopy);

   return SCIP_OKAY;
}

/** sets destructor method of primal heuristic */
SCIP_RETCODE SCIPsetHeurFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURFREE    ((*heurfree))       /**< destructor of primal heuristic */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetHeurFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(heur != NULL);

   SCIPheurSetFree(heur, heurfree);

   return SCIP_OKAY;
}

/** sets initialization method of primal heuristic */
SCIP_RETCODE SCIPsetHeurInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit))       /**< initialize primal heuristic */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetHeurInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(heur != NULL);

   SCIPheurSetInit(heur, heurinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of primal heuristic */
SCIP_RETCODE SCIPsetHeurExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit))       /**< deinitialize primal heuristic */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetHeurExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(heur != NULL);

   SCIPheurSetExit(heur, heurexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of primal heuristic */
SCIP_RETCODE SCIPsetHeurInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol))    /**< solving process initialization method of primal heuristic */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetHeurInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(heur != NULL);

   SCIPheurSetInitsol(heur, heurinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of primal heuristic */
SCIP_RETCODE SCIPsetHeurExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol))    /**< solving process deinitialization method of primal heuristic */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetHeurExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(heur != NULL);

   SCIPheurSetExitsol(heur, heurexitsol);

   return SCIP_OKAY;
}

/** returns the primal heuristic of the given name, or NULL if not existing */
SCIP_HEUR* SCIPfindHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of primal heuristic */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindHeur(scip->set, name);
}

/** returns the array of currently available primal heuristics */
SCIP_HEUR** SCIPgetHeurs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->heurs;
}

/** returns the number of currently available primal heuristics */
int SCIPgetNHeurs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nheurs;
}

/** sets the priority of a primal heuristic */
SCIP_RETCODE SCIPsetHeurPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   priority            /**< new priority of the primal heuristic */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPheurSetPriority(heur, scip->set, priority);

   return SCIP_OKAY;
}

/** create a diving set associated with a primal heuristic. The primal heuristic needs to be included
 *  before this method can be called. The diveset is installed in the array of divesets of the heuristic
 *  and can be retrieved later by accessing SCIPheurGetDivesets()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPcreateDiveset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET**        diveset,            /**< pointer to created diving heuristic settings, or NULL if not needed */
   SCIP_HEUR*            heur,               /**< primal heuristic to which the diveset belongs */
   const char*           name,               /**< name for the diveset, or NULL if the name of the heuristic should be used */
   SCIP_Real             minreldepth,        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth,        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot,      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   SCIP_Real             maxdiveubquot,      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot,     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol, /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol,/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             lpresolvedomchgquot,/**< percentage of immediate domain changes during probing to trigger LP resolve */
   int                   lpsolvefreq,        /**< LP solve frequency for (0: only if enough domain reductions are found by propagation)*/
   int                   maxlpiterofs,       /**< additional number of allowed LP iterations */
   unsigned int          initialseed,        /**< initial seed for random number generation */
   SCIP_Bool             backtrack,          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Bool             onlylpbranchcands,  /**< should only LP branching candidates be considered instead of the slower but
                                              *   more general constraint handler diving variable selection? */
   SCIP_Bool             ispublic,           /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */
   SCIP_Bool             specificsos1score,  /**< should SOS1 variables be scored by the diving heuristics specific score function;
                                              *   otherwise use the score function of the SOS1 constraint handler */
   SCIP_DECL_DIVESETGETSCORE((*divesetgetscore)), /**< method for candidate score and rounding direction */
   SCIP_DECL_DIVESETAVAILABLE((*divesetavailable)) /**< callback to check availability of dive set at the current stage, or NULL if always available */
   )
{
   SCIP_DIVESET* divesetptr = NULL;
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateDiveset", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* create the diveset (this will add diving specific parameters for this heuristic) */
   SCIP_CALL( SCIPdivesetCreate(&divesetptr, heur, name, scip->set, scip->messagehdlr, scip->mem->setmem,
         minreldepth, maxreldepth, maxlpiterquot, maxdiveubquot, maxdiveavgquot, maxdiveubquotnosol,
         maxdiveavgquotnosol, lpresolvedomchgquot, lpsolvefreq, maxlpiterofs, initialseed, backtrack,
         onlylpbranchcands, ispublic, specificsos1score, divesetgetscore, divesetavailable) );

   assert(divesetptr != NULL);
   if( diveset != NULL )
      *diveset = divesetptr;

   return SCIP_OKAY;
}

/** check specific preconditions for diving, e.g., if an incumbent solution is available */
SCIP_RETCODE SCIPisDivesetAvailable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diving heuristic settings */
   SCIP_Bool*            available           /**< pointer to store if the diving can run at the current solving stage */
   )
{
   assert(scip != NULL);
   assert(diveset != NULL);
   assert(available != NULL);

   SCIP_CALL( SCIPdivesetIsAvailable(diveset, scip->set, available) );

   return SCIP_OKAY;
}
