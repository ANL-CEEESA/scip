/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   branch_vanillafullstrong.c
 * @ingroup DEFPLUGINS_BRANCH
 * @brief  vanilla full strong LP branching rule
 * @author Tobias Achterberg
 * @author Maxime Gasse
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/branch_vanillafullstrong.h"
#include "scip/pub_branch.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include <string.h>


#define BRANCHRULE_NAME            "vanillafullstrong"
#define BRANCHRULE_DESC            "vanilla full strong branching"
#define BRANCHRULE_PRIORITY        -2000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define DEFAULT_INTEGRALCANDS      FALSE   /**< should integral variables in the current LP solution be considered as
                                            *   branching candidates ? */
#define DEFAULT_SCOREALL           FALSE   /**< should strong branching scores be computed for all candidates, or can
                                            *   we early stop when a variable has infinite score ? */
#define DEFAULT_IDEMPOTENT         FALSE   /**< should strong branching side-effects be prevented (e.g., domain
                                            *   changes, stat updates etc.) ? */
#define DEFAULT_COLLECTSCORES      FALSE   /**< should strong branching scores be collected ? */
#define DEFAULT_DONOTBRANCH        FALSE   /**< should branching be done ? */
#define DEFAULT_UPDATEPARENTBOUND  TRUE    /**< should parent lower bound be updated based on strong branching? */
#define DEFAULT_UPDATECHILDBOUND   TRUE    /**< should child lower bound be updated based on strong branching? */

#define DEFAULT_MAXDEPTHBINS       10      /**< maximum number of bins for dividing the branching depths */
#define DEFAULT_BRANCHSTATSFILENAME "-"   /**< name of the text output file for branching stats, or "-" if no text
                                            *   output should be created */


/** variable branching stats */
struct VarBranchStats
{
   int nbranchings;
   int mindepth;
   int maxdepth;
   int* depthhist;
};
typedef struct VarBranchStats VARBRANCHSTATS;

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             integralcands;         /**< should integral variables in the current LP solution be considered
                                                 *   as branching candidates ? */
   SCIP_Bool             scoreall;              /**< should strong branching scores be computed for all candidates, or
                                                 *   can we early stop when a node is detected infeasible ? */
   SCIP_Bool             idempotent;            /**< should strong branching side-effects be prevented (e.g., domain
                                                 *   changes, stat updates etc.) ? */
   SCIP_Bool             collectscores;         /**< should strong branching scores be collected ? */
   SCIP_Bool             donotbranch;           /**< should branching be done ? */
   SCIP_Bool             updateparentbound;     /**< should parent lower bound be updated based on strong branching? */
   SCIP_Bool             updatechildbound;      /**< should child lower bound be updated based on strong branching? */
   SCIP_VAR**            cands;                 /**< candidate variables */
   SCIP_Real*            candscores;            /**< candidate scores */
   int                   ncands;                /**< number of candidates */
   int                   npriocands;            /**< number of priority candidates */
   int                   bestcand;              /**< best branching candidate */
   int                   candcapacity;          /**< capacity of candidate arrays */
   SCIP_HASHMAP*         varstatsmap;
   int                   depthbins[DEFAULT_MAXDEPTHBINS];
   int                   ndepthbins;
   SCIP_Bool             printtoconsole;
   char*                 branchstatsfilename;
};


/*
 * local methods
 */


/** selects a variable from a set of candidates by strong branching */
static
SCIP_RETCODE runVanillaStrongBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< branching candidates */
   int                   ncands,             /**< number of branching candidates */
   int                   npriocands,         /**< number of branching candidates with highest priority */
   SCIP_Bool             scoreall,           /**< should strong branching scores be computed for all candidates, or can
                                              *   we early stop when a node is detected infeasible ? */
   SCIP_Bool             idempotent,         /**< should strong branching side-effects be prevented (e.g., domain
                                              *   changes, stat updates etc.) ? */
   SCIP_Real*            scores,             /**< candidate scores */
   int*                  bestcand,           /**< best candidate for branching */
   SCIP_Real*            bestdown,           /**< objective value of the down branch for bestcand */
   SCIP_Real*            bestup,             /**< objective value of the up branch for bestcand */
   SCIP_Real*            bestscore,          /**< score for bestcand */
   SCIP_Bool*            bestdownvalid,      /**< is bestdown a valid dual bound for the down branch? */
   SCIP_Bool*            bestupvalid,        /**< is bestup a valid dual bound for the up branch? */
   SCIP_Real*            provedbound         /**< proved dual bound for current subtree */
   )
{  /*lint --e{715}*/
   SCIP_Real lpobjval;
   int nsbcalls;
   int c;

   assert(scip != NULL);
   assert(cands != NULL);
   assert(bestcand != NULL);
   assert(bestdown != NULL);
   assert(bestup != NULL);
   assert(bestscore != NULL);
   assert(bestdownvalid != NULL);
   assert(bestupvalid != NULL);
   assert(provedbound != NULL);
   assert(ncands > 0);

   /* get current LP objective bound of the local sub problem and global cutoff bound */
   lpobjval = SCIPgetLPObjval(scip);
   *provedbound = lpobjval;

   *bestcand = 0;
   *bestdown = lpobjval;
   *bestup = lpobjval;
   *bestdownvalid = FALSE;
   *bestupvalid = FALSE;
   *bestscore = -SCIPinfinity(scip);

   if( scores != NULL )
      for( c = 0; c < ncands; ++c )
         scores[c] = -SCIPinfinity(scip);

   /* if only one candidate exists, choose this one without applying strong branching; also, when SCIP is about to be
    * stopped, all strongbranching evaluations will be aborted anyway, thus we can return immediately
    */
   if( (!scoreall && ncands == 1) || SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* this assert may not hold if SCIP is stopped, thus we only check it here */
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* initialize strong branching without propagation */
   SCIP_CALL( SCIPstartStrongbranch(scip, FALSE) );

   /* compute strong branching scores */
   nsbcalls = 0;
   for( c = 0; c < ncands ; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Bool integral;
      SCIP_Real down, up;
      SCIP_Real downgain, upgain;
      SCIP_Bool downvalid, upvalid;
      SCIP_Bool downinf, upinf;
      SCIP_Bool downconflict, upconflict;
      SCIP_Bool lperror;
      SCIP_Real gains[3];
      SCIP_Real score;

      var = cands[c];
      assert(var != NULL);

      val = SCIPvarGetLPSol(var);
      integral = SCIPisFeasIntegral(scip, val);

      up = -SCIPinfinity(scip);
      down = -SCIPinfinity(scip);

      SCIPdebugMsg(scip, "applying vanilla strong branching on variable <%s> with solution %g\n",
         SCIPvarGetName(var), val);

      /* apply strong branching */
      if( integral )
      {
         SCIP_CALL( SCIPgetVarStrongbranchInt(scip, cands[c], INT_MAX, idempotent,
               &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );
      }
      else
      {
         SCIP_CALL( SCIPgetVarStrongbranchFrac(scip, cands[c], INT_MAX, idempotent,
               &down, &up, &downvalid, &upvalid, &downinf, &upinf, &downconflict, &upconflict, &lperror) );
      }
      nsbcalls++;

      /* check for an error in strong branching */
      if( lperror )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "(node %" SCIP_LONGINT_FORMAT ") error in strong branching call for variable <%s> with solution %g\n",
            SCIPgetNNodes(scip), SCIPvarGetName(var), val);
         break;
      }

      /* evaluate strong branching */
      down = MAX(down, lpobjval);
      up = MAX(up, lpobjval);
      downgain = down - lpobjval;
      upgain = up - lpobjval;

      assert(!SCIPallColsInLP(scip) || SCIPisExactSolve(scip) || !downvalid || downinf == SCIPisGE(scip, down, SCIPgetCutoffbound(scip)));
      assert(!SCIPallColsInLP(scip) || SCIPisExactSolve(scip) || !upvalid || upinf == SCIPisGE(scip, up, SCIPgetCutoffbound(scip)));
      assert(downinf || !downconflict);
      assert(upinf || !upconflict);

      if( !idempotent )
      {
         /* display node information line */
         if( SCIPgetDepth(scip) == 0 && nsbcalls % 100 == 0 )
         {
            SCIP_CALL( SCIPprintDisplayLine(scip, NULL, SCIP_VERBLEVEL_HIGH, TRUE) );
         }
         /* update variable pseudo cost values */
         if( !downinf && downvalid )
         {
            SCIP_CALL( SCIPupdateVarPseudocost(scip, var, integral ? -1.0 : 0.0 - SCIPfrac(scip, val), downgain, 1.0) );
         }
         if( !upinf && upvalid )
         {
            SCIP_CALL( SCIPupdateVarPseudocost(scip, var, integral ? +1.0 : 1.0 - SCIPfrac(scip, val), upgain, 1.0) );
         }
      }

      /* compute strong branching score */
      gains[0] = downgain;
      gains[1] = upgain;
      gains[2] = 0.0;
      score = SCIPgetBranchScoreMultiple(scip, var, integral ? 3 : 2, gains);

      /* collect scores if requested */
      if( scores != NULL )
         scores[c] = score;

      /* check for a better score */
      if( score > *bestscore )
      {
         *bestcand = c;
         *bestdown = down;
         *bestup = up;
         *bestdownvalid = downvalid;
         *bestupvalid = upvalid;
         *bestscore = score;
      }

      SCIPdebugMsg(scip, " -> cand %d/%d (prio:%d) var <%s> (solval=%g, downgain=%g, upgain=%g, score=%g) -- best: <%s> (%g)\n",
         c, ncands, npriocands, SCIPvarGetName(var), val, downgain, upgain, score,
         SCIPvarGetName(cands[*bestcand]), *bestscore);

      /* node is infeasible -> early stopping (highest score) */
      if( !integral && !scoreall && downinf && upinf )
      {
         /* we should only detect infeasibility if the LP is a valid relaxation */
         assert(SCIPallColsInLP(scip));
         assert(!SCIPisExactSolve(scip));
         assert(*bestcand == c);

         SCIPdebugMsg(scip, " -> variable <%s> is infeasible in both directions\n", SCIPvarGetName(var));
         break;
      }
   }

   /* end strong branching */
   SCIP_CALL( SCIPendStrongbranch(scip) );

   /* update proved bound */
   if( *bestdownvalid && *bestupvalid && !SCIPisFeasIntegral(scip, SCIPvarGetLPSol(cands[*bestcand])) )
   {
      SCIP_Real minbound = MIN(*bestdown, *bestup);

      *provedbound = MAX(*provedbound, minbound);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyVanillafullstrong)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleVanillafullstrong(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeVanillafullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeBlockMemoryNull(scip, &branchruledata);

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitVanillafullstrong)
{  /*lint --e{715}*/
#ifndef NDEBUG
   SCIP_BRANCHRULEDATA* branchruledata;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
#endif
   assert(branchruledata != NULL);
   assert(branchruledata->candscores == NULL);
   assert(branchruledata->cands == NULL);

   SCIP_CALL( SCIPhashmapCreate(&(branchruledata->varstatsmap), SCIPblkmem(scip), SCIPgetNVars(scip)) );

   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
static
SCIP_DECL_BRANCHEXIT(branchExitVanillafullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   VARBRANCHSTATS* stats;
   SCIP_HASHMAPENTRY* entry;
   FILE* file;
   int nvars;
   int i;
   int j;

   /* initialize branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* write the branch stats file as needed */
   file = fopen(branchruledata->branchstatsfilename, "w");
   if( file != NULL )
   {
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);

      fprintf(file, "varname nbranched mindepth maxdepth");
      for( i = 0; i <= branchruledata->ndepthbins; i++ )
         fprintf(file, " depth[%d,%d]", (i == 0 ? 0 : branchruledata->depthbins[i-1] + 1), (i < branchruledata->ndepthbins ? branchruledata->depthbins[i] : INT_MAX));
      fprintf(file, "\n");

      for( i = 0; i < nvars; i++ )
      {
         var = vars[i];
         stats = (VARBRANCHSTATS*) SCIPhashmapGetImage(branchruledata->varstatsmap, (void*)var);

         if( stats == NULL )
            continue;

         if( stats->nbranchings > 0 )
         {
            fprintf(file, "%s %d %d %d", SCIPvarGetName(var), stats->nbranchings, stats->mindepth, stats->maxdepth);
            for( j = 0; j <= branchruledata->ndepthbins; j++ )
               fprintf(file, " %d", stats->depthhist[j]);
            fprintf(file, "\n");

            if( branchruledata->printtoconsole )
            {
               SCIPinfoMessage(scip, NULL, "%s branched %d times [depths %d-%d]:", SCIPvarGetName(var), stats->nbranchings, stats->mindepth, stats->maxdepth);
               for( j = 0; j <= branchruledata->ndepthbins; j++ )
                  SCIPinfoMessage(scip, NULL, " depth%d=%d", j, stats->depthhist[j]);
               SCIPinfoMessage(scip, NULL, "\n");
            }
         }
      }
   }
   fclose(file);

   /* free branch stats hashmap and arrays if any */
   for( i = 0; i < SCIPhashmapGetNEntries(branchruledata->varstatsmap); i++ )
   {
      entry = SCIPhashmapGetEntry(branchruledata->varstatsmap, i);

      if( entry == NULL )
         continue;

      stats = (VARBRANCHSTATS*) SCIPhashmapEntryGetImage(entry);

      /* if stats has been added to the hashmap, it can't be empty */
      assert(stats->depthhist != NULL);

      SCIPfreeBlockMemoryArray(scip, &stats->depthhist, (branchruledata->ndepthbins + 1));
      SCIPfreeBlockMemory(scip, &stats);
   }
   SCIPhashmapFree(&branchruledata->varstatsmap);

   /* free candidate arrays if any */
   if( branchruledata->candscores != NULL )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->candscores, branchruledata->candcapacity);
   }
   if( branchruledata->cands != NULL )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->cands, branchruledata->candcapacity);
   }

   branchruledata->candcapacity = -1;
   branchruledata->ncands = -1;
   branchruledata->npriocands = -1;
   branchruledata->bestcand = -1;

   return SCIP_OKAY;
}

/** branching execution method */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpVanillafullstrong)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Real bestdown;
   SCIP_Real bestup;
   SCIP_Real bestscore;
   SCIP_Real provedbound;
   SCIP_Bool bestdownvalid;
   SCIP_Bool bestupvalid;
   SCIP_VAR** cands;
   int ncands;
   int npriocands;
   int i;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Execlp method of vanilla fullstrong branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates, either all non-fixed variables or only the
    * fractional ones */
   if( branchruledata->integralcands )
   {
      SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, &ncands, &npriocands) );
   }
   else
   {
      SCIP_CALL( SCIPgetLPBranchCands(scip, &cands, NULL, NULL, &ncands, &npriocands, NULL) );
   }

   assert(ncands > 0);
   assert(npriocands > 0);

   /* increase candidate arrays capacity if needed */
   if( ncands > branchruledata->candcapacity )
   {
      /* free previously allocated arrays if any */
      if( branchruledata->candscores != NULL)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->candscores, branchruledata->candcapacity);
         branchruledata->candscores = NULL;
      }
      if( branchruledata->cands != NULL)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &branchruledata->cands, branchruledata->candcapacity);
         branchruledata->cands = NULL;
      }

      /* update capacity */
      branchruledata->candcapacity = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);
   }
   assert(branchruledata->candcapacity >= ncands);

   /* allocate new candidate arrays if needed */
   if( branchruledata->cands == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->cands, branchruledata->candcapacity) );
   }
   if( branchruledata->candscores == NULL && branchruledata->collectscores )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &branchruledata->candscores, branchruledata->candcapacity) );
   }

   /* copy candidates */
   branchruledata->ncands = ncands;
   branchruledata->npriocands = npriocands;

   for( i = 0; i < ncands; i++ )
      branchruledata->cands[i] = cands[i];

   SCIP_CALL( runVanillaStrongBranching(scip, branchruledata->cands, branchruledata->ncands, branchruledata->npriocands,
         branchruledata->scoreall, branchruledata->idempotent, branchruledata->candscores,
         &branchruledata->bestcand, &bestdown, &bestup, &bestscore, &bestdownvalid,
         &bestupvalid, &provedbound) );

   if( !branchruledata->donotbranch )
   {
      assert(0 <= branchruledata->bestcand && branchruledata->bestcand < branchruledata->ncands);

      if( SCIPisGE(scip, provedbound, SCIPgetCutoffbound(scip)) )
      {
         SCIPdebugMsg(scip, " -> variable <%s> is infeasible in both directions\n",
               SCIPvarGetName(branchruledata->cands[branchruledata->bestcand]));

         *result = SCIP_CUTOFF;
      }
      else
      {
         SCIP_VAR* var;
         SCIP_Real val;
         SCIP_NODE* downchild;
         SCIP_NODE* eqchild;
         SCIP_NODE* upchild;
         SCIP_Bool allcolsinlp;
         SCIP_Bool exactsolve;
         VARBRANCHSTATS* stats;
         int depth;
         int bin;
         int i;

         /* check, if we want to solve the problem exactly, meaning that strong branching information is not useful
          * for cutting off sub problems and improving lower bounds of children
          */
         exactsolve = SCIPisExactSolve(scip);

         /* check, if all existing columns are in LP, and thus the strong branching results give lower bounds */
         allcolsinlp = SCIPallColsInLP(scip);

         if( branchruledata->updateparentbound && allcolsinlp && !exactsolve )
         {
            SCIP_CALL( SCIPupdateLocalLowerbound(scip, provedbound) );
         }

         assert(SCIPisLT(scip, provedbound, SCIPgetCutoffbound(scip)));

         var = branchruledata->cands[branchruledata->bestcand];
         val = SCIPvarGetLPSol(var);

         /* perform the branching */
         SCIPdebugMsg(scip, " -> %d candidates, selected candidate %d: variable <%s>[%g,%g] (solval=%g, down=%g, up=%g, score=%g)\n",
               branchruledata->ncands, branchruledata->bestcand, SCIPvarGetName(var), SCIPvarGetLbLocal(var),
               SCIPvarGetUbLocal(var), val, bestdown, bestup, bestscore);
         SCIP_CALL( SCIPbranchVarVal(scip, var, val, &downchild, &eqchild, &upchild) );

         /* update the lower bounds in the children */
         if( branchruledata->updatechildbound && allcolsinlp && !exactsolve )
         {
            if( downchild != NULL && bestdownvalid )
            {
               SCIP_CALL( SCIPupdateNodeLowerbound(scip, downchild, bestdown) );
               SCIPdebugMsg(scip, " -> down child's lowerbound: %g\n", SCIPnodeGetLowerbound(downchild));
            }
            if( upchild != NULL && bestupvalid )
            {
               SCIP_CALL( SCIPupdateNodeLowerbound(scip, upchild, bestup) );
               SCIPdebugMsg(scip, " -> up child's lowerbound:   %g\n", SCIPnodeGetLowerbound(upchild));
            }
         }

         /* init stats */
         stats = (VARBRANCHSTATS*) SCIPhashmapGetImage(branchruledata->varstatsmap, (void*)var);
         if( stats == NULL )
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &stats) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &stats->depthhist, branchruledata->ndepthbins + 1) );
            stats->nbranchings = 0;
            stats->mindepth = INT_MAX;
            stats->maxdepth = -1;
            for( i = 0; i <= branchruledata->ndepthbins; i++ )
               stats->depthhist[i] = 0;

            SCIP_CALL( SCIPhashmapInsert(branchruledata->varstatsmap, (void*)var, (void*)stats) );
         }
         assert(stats != NULL);
         /* update stats */
         depth = SCIPgetDepth(scip);
         stats->nbranchings++;
         stats->mindepth = MIN(stats->mindepth, depth);
         stats->maxdepth = MAX(stats->maxdepth, depth);
         bin = 0;
         while( bin < branchruledata->ndepthbins && depth > branchruledata->depthbins[bin] )
            bin++;
         stats->depthhist[bin]++;

         *result = SCIP_BRANCHED;
      }
   }

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the vanilla full strong LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleVanillafullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;
   int i;

   /* create fullstrong branching rule data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );
   branchruledata->cands = NULL;
   branchruledata->candscores = NULL;
   branchruledata->candcapacity = -1;
   branchruledata->ncands = -1;
   branchruledata->npriocands = -1;
   branchruledata->bestcand = -1;
   branchruledata->varstatsmap = NULL;
   branchruledata->branchstatsfilename = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non-fundamental callbacks via specific setter functions*/
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyVanillafullstrong) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeVanillafullstrong) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitVanillafullstrong) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitVanillafullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpVanillafullstrong) );

   /* fullstrong branching rule parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "branching/vanillafullstrong/ndepthbins",
         "number of depth bins while tracking the variable branching statistics (max 10)",
         &branchruledata->ndepthbins, TRUE, 5, 1, DEFAULT_MAXDEPTHBINS, NULL, NULL) );
   for( i = 0; i < DEFAULT_MAXDEPTHBINS; i++ )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "branching/vanillafullstrong/depthbin%d", i+1);
      SCIP_CALL( SCIPaddIntParam(scip,
               paramname,
               "depth upper bound for the current bin",
               &branchruledata->depthbins[i], TRUE, (i + 1) * 10, 0, INT_MAX, NULL, NULL) );
   }
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/integralcands",
         "should integral variables in the current LP solution be considered as branching candidates?",
         &branchruledata->integralcands, FALSE, DEFAULT_INTEGRALCANDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/idempotent",
         "should strong branching side-effects be prevented (e.g., domain changes, stat updates etc.)?",
         &branchruledata->idempotent, FALSE, DEFAULT_IDEMPOTENT, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/scoreall",
         "should strong branching scores be computed for all candidates, or can we early stop when a variable has infinite score?",
         &branchruledata->scoreall, TRUE, DEFAULT_SCOREALL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/collectscores",
         "should strong branching scores be collected?",
         &branchruledata->collectscores, TRUE, DEFAULT_COLLECTSCORES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/donotbranch",
         "should candidates only be scored, but no branching be performed?",
         &branchruledata->donotbranch, TRUE, DEFAULT_DONOTBRANCH, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/updateparentbound",
         "should parent lower bound be updated based on strong branching?",
         &branchruledata->updateparentbound, TRUE, DEFAULT_UPDATEPARENTBOUND, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/updatechildbound",
         "should child lower bound be updated based on strong branching?",
         &branchruledata->updatechildbound, TRUE, DEFAULT_UPDATECHILDBOUND, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "branching/vanillafullstrong/printtoconsole",
         "should variable branching stats be printed to the console at the end?",
         &branchruledata->printtoconsole, TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddStringParam(scip,
         "branching/vanillafullstrong/branchstatsfilename",
         "name of the text output file for branching stats, or - if no text output should be created",
         &branchruledata->branchstatsfilename, TRUE, DEFAULT_BRANCHSTATSFILENAME,
         NULL, NULL) );

   return SCIP_OKAY;
}


/** recovers candidate variables and their scores from last vanilla full strong branching call */
SCIP_RETCODE SCIPgetVanillafullstrongData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           cands,              /**< pointer to store candidate variables; or NULL */
   SCIP_Real**           candscores,         /**< pointer to store candidate scores; or NULL */
   int*                  ncands,             /**< pointer to store number of candidates; or NULL */
   int*                  npriocands,         /**< pointer to store number of priority candidates; or NULL */
   int*                  bestcand            /**< pointer to store best branching candidate; or NULL */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert( branchrule != NULL );
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert( branchruledata != NULL );

   if( cands )
   {
      *cands = branchruledata->cands;
   }
   if( candscores && branchruledata->collectscores )
   {
      *candscores = branchruledata->candscores;
   }
   if( ncands )
   {
      *ncands = branchruledata->ncands;
   }
   if( npriocands )
   {
      *npriocands = branchruledata->npriocands;
   }
   if( bestcand )
   {
      *bestcand = branchruledata->bestcand;
   }

   return SCIP_OKAY;
}
