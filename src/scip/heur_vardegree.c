/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define SCIP_DEBUG
#define SCIP_STATISTIC
/**@file   heur_vardegree.c
 * @brief  LNS heuristic that tries to delimit the search region to a neighborhood in the constraint graph
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_vardegree.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "vardegree"
#define HEUR_DESC             "vardegree works on k-neighborhood in a variable-constraint graph"
#define HEUR_DISPCHAR         'K'
#define HEUR_PRIORITY         -1103000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          8
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      500           /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_MAXNODES      5000          /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINIMPROVE    0.01          /**< factor by which Vardegree should at least improve the incumbent */
#define DEFAULT_MINNODES      500           /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_MINFIXINGRATE 0.66          /**< minimum percentage of integer variables that have to be fixed */
#define DEFAULT_NODESQUOT     0.15          /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_NWAITINGNODES 20            /**< number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS     FALSE         /**< should subproblem be created out of the rows in the LP rows,
                                             **< otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE          /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                             **< cutpool of the original scip be copied to constraints of the subscip */
#define DEFAULT_BESTSOLLIMIT   3            /**< limit on number of improving incumbent solutions in sub-CIP */
#define DEFAULT_FIXCONTVARS FALSE           /**< should continuous variables outside the neighborhoods get fixed? */
#define DEFAULT_POTENTIAL      'r'          /**< the reference point to compute the neighborhood potential: (r)oot or (p)seudo solution */
#define DEFAULT_MAXDISTANCE     3           /**< maximum distance to selected variable to enter the subproblem, or -1 to
                                             *   select the distance that best approximates the minimum fixing rate from below */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which Vardegree should at least improve the incumbent */
   SCIP_Longint          usednodes;          /**< nodes already used by Vardegree in earlier calls */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   unsigned int          randseed;           /**< seed value for random number generator */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows? */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem? */
   SCIP_Bool             fixcontvars;        /**< should continuous variables outside the neighborhoods get fixed? */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP */
   int                   maxdistance;        /**< maximum distance to selected variable to enter the subproblem, or -1 to
                                              *   select the distance that best approximates the minimum fixing rate from below */
   int                   sumneighborhoodvars;/**< neighborhood variables sum over all seen neighborhoods */
   int                   sumdiscneighborhoodvars; /**< neighborhood discrete variables sum over all seen neighboorhoods */
   int                   nneighborhoods;     /**< number of calculated neighborhoods */
   char                  potential;          /**< the reference point to compute the neighborhood potential: (r)oot or (p)seudo solution */
};

/** variable graph data structure;
 *
 *  a bipartite graph with variables and global constraints as nodes, and edges if a variable is part of the constraint */
struct VariableGraph
{
   SCIP_CONS***          varconss;           /**< constraints of each variable */
   int*                  nvarconss;          /**< number of constraints for each variable */
   int*                  varconssize;        /**< size array for every varconss entry */
};
typedef struct VariableGraph VARIABLEGRAPH;

/*
 * Local methods
 */

/* breadth first search on the variable constraint graph; uses a special kind of queue data structure that holds two queue levels
 * at the same time: the variables at the current distance and the ones at the next distance
 */
static
SCIP_RETCODE variablegraphBreadthFirst(
   SCIP*                scip,                /**< SCIP data structure */
   VARIABLEGRAPH*       vargraph,            /**< pointer to the variable graph */
   SCIP_VAR*            startvar,            /**< variable to calculate distance from */
   int*                 distances,           /**< array to keep distance in vargraph from startvar for every variable */
   int                  maxdistance          /**< maximum distance from start variable */
   )
{
   SCIP_VAR** vars;

   int* queue;
   int  nvars;
   int i;
   int currlvlidx;
   int nextlvlidx;
   int increment = 1;
   int currentdistance;
   SCIP_VAR** varbuffer;

   assert(scip != NULL);
   assert(vargraph != NULL);
   assert(startvar != NULL);
   assert(distances != NULL);

   /* get variable data */
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   SCIP_CALL( SCIPallocClearBufferArray(scip, &queue, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &varbuffer, nvars) );

   /* initialize distances to -1 */
   for( i = 0; i < nvars; ++i )
   {
      queue[i] = -1;
      distances[i] = -1;
   }

   /* start variable has a distance of 0 */
   distances[SCIPvarGetProbindex(startvar)] = 0;


   queue[0] = SCIPvarGetProbindex(startvar);
   currlvlidx = 0;
   nextlvlidx = nvars - 1;

   /* loop over the queue and pop the next variable, starting with start variable */
   do
   {
      SCIP_VAR* currvar;
      int c;
      int varpos;

      assert(distances[queue[currlvlidx]] >= 0);
      currvar = vars[queue[currlvlidx]];
      varpos = SCIPvarGetProbindex(currvar);
      currentdistance = distances[currlvlidx];

      /* check for termination because maximum distance has been reached */
      if( currentdistance > maxdistance )
         break;

      assert(varpos >= 0);

      /* loop over variable constraints and enqueue variables that were not visited yet */
      for( c = 0; c < vargraph->nvarconss[varpos]; ++c )
      {
         int nconsvars;
         int v;
         SCIP_Bool success;
         SCIP_CONS* cons = vargraph->varconss[varpos][c];

         /* request number of variables */
         SCIP_CALL( SCIPgetConsNVars(scip, cons, &nconsvars, &success) );

         if( !success )
            continue;

         /* request constraint to write its variables into this buffer here */
         SCIP_CALL( SCIPgetConsVars(scip, cons, varbuffer, nvars, &success) );

         if( !success )
            continue;

         /* loop over variables of this constraint */
         for( v = 0; v < nconsvars; ++v )
         {
            SCIP_VAR* nextvar = varbuffer[v];
            int nextvarpos;
            assert(nextvar != NULL);
            if( !SCIPvarIsActive(nextvar) )
               continue;

            nextvarpos = SCIPvarGetProbindex(nextvar);
            assert(nextvarpos >= 0);

            /* insert variables that were not considered yet into the next level queue */
            if( distances[nextvarpos] == -1 )
            {
               distances[nextvarpos] = distances[varpos] + 1;
               queue[nextlvlidx] = nextvarpos;
               nextlvlidx -= increment;
            }
         } /* end constraint variables loop */
      } /* end constraint loop */

      queue[currlvlidx] = -1;
      currlvlidx += increment;

      /* check if we need to swap current and next level index and reverse the increment */
      if( currlvlidx == nvars || currlvlidx == 0 || queue[currlvlidx] == -1 || currlvlidx == nextlvlidx )
      {
         /* increment knows whether we are currently looping upwards (all variables with odd distance) or downwards the queue */
         if( increment == +1 )
         {
            currlvlidx = nvars - 1;
            nextlvlidx = 0;
            increment = -1;
         }
         else
         {
            currlvlidx = 0;
            nextlvlidx = nvars - 1;
            increment = +1;
         }
      }
   }
   while( queue[currlvlidx] != -1 && distances[queue[currlvlidx]] >= currentdistance );



   SCIPfreeBufferArray(scip, &varbuffer);
   SCIPfreeBufferArray(scip, &queue);

   return SCIP_OKAY;
}

/* fills variable graph data structure
 *
 * loops over global problem constraints and updates a mapping from the variables to their respective constraints
 */
static
SCIP_RETCODE fillVariableGraph(
   SCIP*                scip,                /**< SCIP data structure */
   VARIABLEGRAPH*       vargraph             /**< pointer to the variable graph */
   )
{
   SCIP_CONS** conss;
   int nconss;
   int nvars;
   int c;
   SCIP_VAR** varbuffer;

   assert(scip != NULL);
   assert(vargraph != NULL);

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuffer, nvars) );

   for( c = 0; c < nconss; ++c )
   {
      int nconsvars;
      int v;
      SCIP_Bool success;
      SCIP_CONS* cons = conss[c];

      /* we only consider constraints that are checkable */
      if( !SCIPconsIsChecked(cons) || !SCIPconsIsInitial(cons) )
         continue;

      /* request number of variables */
      SCIP_CALL( SCIPgetConsNVars(scip, cons, &nconsvars, &success) );

      if( !success )
         continue;

      /*request constraint to write its variables into this buffer here */
      SCIP_CALL( SCIPgetConsVars(scip, cons, varbuffer, nvars, &success) );

      if( !success )
         continue;

      /* loop over constraint variables and add this constraint to them if they are active */
      for( v = 0; v < nconsvars; ++v )
      {
         int varpos = SCIPvarGetProbindex(varbuffer[v]);

         /* skip inactive variables */
         if( varpos == -1 )
            continue;

         /* ensure array size */
         if( vargraph->varconssize[varpos] == vargraph->nvarconss[varpos]  )
         {
            int newmem = SCIPcalcMemGrowSize(scip, vargraph->nvarconss[varpos] + 1);

            assert(newmem > vargraph->varconssize[varpos]);

            if( vargraph->varconss[varpos] == NULL )
            {
               SCIP_CALL( SCIPallocMemoryArray(scip, &vargraph->varconss[varpos], newmem) );
            }
            else
            {
               SCIP_CALL( SCIPreallocMemoryArray(scip, &vargraph->varconss[varpos], newmem) );
            }
            vargraph->varconssize[varpos] = newmem;
         }

         assert(vargraph->nvarconss[varpos] < vargraph->varconssize[varpos]);

         /* add constraint to constraint array for this variable */
         vargraph->varconss[varpos][vargraph->nvarconss[varpos]] = cons;
         vargraph->nvarconss[varpos] += 1;
      }

   }

   /* free the buffer */
   SCIPfreeBufferArray(scip, &varbuffer);

   return SCIP_OKAY;
}

/** initialization method of variable graph data structure */
static
SCIP_RETCODE variableGraphCreate(
   SCIP*                scip,                /**< SCIP data structure */
   VARIABLEGRAPH**      vargraph             /**< pointer to the variable graph */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(vargraph != NULL);

   nvars = SCIPgetNVars(scip);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocMemory(scip, vargraph) );

   /* allocate and clear memory */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(*vargraph)->varconss, nvars) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(*vargraph)->nvarconss, nvars) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &(*vargraph)->varconssize, nvars) );

   SCIP_CALL( fillVariableGraph(scip, *vargraph) );

   return SCIP_OKAY;
}

/** deinitialization method of variable graph data structure */
static
SCIP_RETCODE variableGraphFree(
   SCIP*                scip,                /**< SCIP data structure */
   VARIABLEGRAPH**      vargraph             /**< pointer to the variable graph */
   )
{
   int nvars;
   int v;
   assert(scip != NULL);
   assert(vargraph != NULL);

   nvars = SCIPgetNVars(scip);

   for( v = nvars - 1; v >= 0; --v )
   {
      SCIPfreeMemoryArrayNull(scip, &(*vargraph)->varconss[v]);
   }

   /* allocate and clear memory */
   SCIPfreeMemoryArray(scip, &(*vargraph)->varconss);
   SCIPfreeMemoryArray(scip, &(*vargraph)->nvarconss);
   SCIPfreeMemoryArray(scip, &(*vargraph)->varconssize);

   SCIPfreeMemory(scip, vargraph);

   return SCIP_OKAY;
}
/** get the potential of a subset of variables (distance to a reference point such as the pseudo-solution or root LP solution) */
static
SCIP_Real getPotential(
   SCIP*                scip,                /**< SCIP data structure */
   SCIP_HEURDATA*       heurdata,            /**< heuristic data */
   SCIP_SOL*            sol,                 /**< solution */
   SCIP_VAR**           vars,                /**< variable array */
   int                  nvars                /**< length of variable array */
   )
{
   SCIP_Real potential;
   int i;
   assert(scip != NULL);
   assert(vars != NULL);
   assert(sol != NULL);


   if( nvars == 0 )
      return 0.0;

   potential = 0.0;

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real objdelta;
      SCIP_VAR* var;
      SCIP_Real referencepoint;
      SCIP_Real varobj;

      var = vars[i];
      assert(var != NULL);
      varobj = SCIPvarGetObj(var);

      if( SCIPisZero(scip, varobj) )
         continue;

      /* determine the reference point for potential computation */
      switch( heurdata->potential )
      {
         /* use difference to pseudo solution using the bound in the objective direction */
         case 'p':
            referencepoint = SCIPvarGetObj(var) > 0.0 ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
            break;
         /* use root LP solution difference */
         case 'r':
            referencepoint = SCIPvarGetRootSol(var);
            break;
         default:
            SCIPerrorMessage("Unknown potential computation %c specified\n", heurdata->potential);
            SCIPABORT();
            referencepoint = SCIPgetSolVal(scip, sol, var);
            break;
      }
      /* calculate the delta to the variables best bound */

      if( SCIPisInfinity(scip, REALABS(referencepoint)) )
         continue;

      objdelta = (SCIPgetSolVal(scip, sol, var) - referencepoint) * SCIPvarGetObj(var);
      potential += objdelta;
   }

   return potential;
}

/** gets the average neighborhood size of all selected variables */
static
SCIP_Real heurdataAvgNeighborhoodSize(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   return heurdata->sumneighborhoodvars / (MAX(1.0, (SCIP_Real)heurdata->nneighborhoods));
}

/** gets the average size of a discrete neighborhood over all variables tested */
static
SCIP_Real heurdataAvgDiscreteNeighborhoodSize(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   return heurdata->sumdiscneighborhoodvars / (MAX(1.0, (SCIP_Real)heurdata->nneighborhoods));
}

/** is the variable in the current neighborhood which is given by the breadth-first distances from a central variable? */
static
SCIP_Bool isVariableInNeighborhood(
   SCIP_VAR*             var,                /**< problem variable */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   int                   maxdistance         /**< maximum distance (inclusive) to be considered for neighborhoods */
   )
{
   assert(var != NULL);
   assert(distances != NULL);
   assert(maxdistance >= 0);
   assert(SCIPvarGetProbindex(var) >= 0);

   return (distances[SCIPvarGetProbindex(var)] != -1 && distances[SCIPvarGetProbindex(var)] <= maxdistance);
}

/** creates subproblem out of LP rows */
static
SCIP_RETCODE createConstraintsFromRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars             /**< the variables of the subproblem */
   )
{
   SCIP_ROW** rows;   /* original scip rows */
   int nrows;
   int i;
   int j;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      SCIP_COL** cols;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real* vals;
      int nnonz;

      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert( lhs <= rhs );

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** fixes variables in subproblem based on long breadth-first distances in variable graph */
static
SCIP_RETCODE fixNonNeighborhoodVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_SOL*             sol,                /**< solution in main SCIP for fixing values */
   SCIP_VAR**            vars,               /**< variables in the main SCIP */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   int                   maxdistance,        /**< maximum distance (inclusive) to be considered for neighborhoods */
   int*                  nfixings            /**< pointer to store number of fixed variables */
   )
{
   int i;
   int nbinvars;
   int nintvars;
   int nvars;
   int nvarstofix;

   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nvarstofix = heurdata->fixcontvars ? nvars : nbinvars + nintvars;
   *nfixings = 0;
   /* change bounds of variables of the subproblem */
   for( i = 0; i < nvarstofix; i++ )
   {
      /* fix all variables that are too far away from this variable according to maxdistance */
      if( distances[i] == -1 || distances[i] > maxdistance )
      {
         SCIP_Real solval;
         SCIP_Real lb;
         SCIP_Real ub;

         solval = SCIPgetSolVal(scip, sol, vars[i]);
         lb = SCIPvarGetLbGlobal(subvars[i]);
         ub = SCIPvarGetUbGlobal(subvars[i]);
         assert(SCIPisLE(scip, lb, ub));

         /* due to dual reductions, it may happen that the solution value is not in
            the variable's domain anymore */
         if( SCIPisLT(scip, solval, lb) )
            solval = lb;
         else if( SCIPisGT(scip, solval, ub) )
            solval = ub;

         /* perform the bound change */
         if( !SCIPisInfinity(scip, solval) && !SCIPisInfinity(scip, -solval) )
         {
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );
            ++(*nfixings);
         }
      }
   }

   return SCIP_OKAY;
}

/** determine the maximum allowed distance to stay within the restriction to fix at least minfixingrate many non neighborhood variables */
static
SCIP_RETCODE determineMaxDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   int*                  choosevardistance   /**< pointer to store the computed maximum distance */
   )
{
   int* distancestmp;
   int nrelevantdistances;
   int criticalidx;
   int zeropos;
   int nvars;
   int nbinvars;
   int nintvars;

   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nrelevantdistances = (heurdata->fixcontvars ?  nvars : (nbinvars + nintvars));

   /* copy the relevant distances of either the discrete or all problem variables and sort them */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &distancestmp, distances, nrelevantdistances) );
   SCIPsortInt(distancestmp, nrelevantdistances);

   /* distances can be infinite in the case of multiple connected components; therefore, search for the index of the
    * zero entry, which is the unique representative of the chosen variable in the sorted distances */
   zeropos = -1;
   SCIPsortedvecFindInt(distancestmp, 0, nrelevantdistances, &zeropos);
   assert(zeropos >= 0);

   /* determine the critical index to look for an appropriate neighborhood distance, starting from the zero position */
   criticalidx = zeropos + (int)((1.0 - heurdata->minfixingrate) * nrelevantdistances);
   criticalidx = MIN(criticalidx, nrelevantdistances - 1);

   /* determine the maximum breadth-first distance such that the neighborhood size stays small enough (below 1-minfixingrate)*/
   *choosevardistance = distancestmp[criticalidx];

   /* we set the distance to exactly the distance at the critical index. If the distance at critical index is not the
    * last one before the distance increases, we decrease the choosevardistance such that the entire neighborhood
    * fits into the minfixingrate restriction
    */
   if( criticalidx != nrelevantdistances - 1 && distancestmp[criticalidx + 1] == *choosevardistance )
      (*choosevardistance)--;

   SCIPfreeBufferArray(scip, &distancestmp);

   return SCIP_OKAY;
}

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed */
   unsigned int*         randseed,           /**< a seed value for the random number generator */
   SCIP_Bool             uselprows,          /**< should subproblem be created out of the rows in the LP rows? */
   SCIP_Bool*            success             /**< used to store whether the creation of the subproblem worked */
   )
{
   SCIP_VAR** vars;                          /* original scip variables */
   SCIP_VAR** varscopy;
   SCIP_SOL* sol;                            /* pool of solutions */
   int* distances;
   VARIABLEGRAPH* vargraph;
   SCIP_VAR* selvar;
   SCIP_Real maxpotential;
   int nintegralvarsleft; /* the number of remaining discrete variables after deletion of previous neighborhood from graph */
   int nvars;
   int nfixings;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int nintegralvarsbound;
   int nsearched;
   int searchlimit;
   int fixthreshold;
   int v;
   int maxdistance;
   int selvarmaxdistance;


   *success = TRUE;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, NULL) );
   sol = SCIPgetBestSol(scip);
   assert(sol != NULL);

   /* create variable graph */
   SCIPdebugMessage("Creating variable constraint graph\n");
   SCIP_CALL( variableGraphCreate(scip, &vargraph) );

   /* allocate buffer memory to hold distances */
   SCIP_CALL( SCIPallocBufferArray(scip, &distances, nvars) );



   /* copy SCIP variables */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &varscopy, vars, nbinvars + nintvars + nimplvars) );
   nsearched = 0;
   maxpotential = SCIP_REAL_MIN;

   /* determine upper bound on neighborhood size */
   nintegralvarsbound = (int)((1.0 - minfixingrate) * (nbinvars + nintvars));

   nintegralvarsleft  = nbinvars + nintvars + nimplvars;
   selvar = NULL;

   /* sort inactive variables to the end of the array */
   for( v = nintegralvarsleft - 1; v >= 0; --v )
   {
      if( ! SCIPvarIsActive(varscopy[v]) )
      {
         varscopy[v] = varscopy[nintegralvarsleft - 1];
         --nintegralvarsleft;
      }
   }

   /* maximum distance from selected variable for breadth-first search (if set to -1, we compute an exhaustive breadth-first search and sort the variables by their distance) */
   maxdistance = (heurdata->maxdistance == - 1 ? INT_MAX : heurdata->maxdistance);

   /* adjust the search limit */
   searchlimit = heurdata->nneighborhoods < 10 ? 5 : (int)(nintegralvarsleft / heurdataAvgDiscreteNeighborhoodSize(heurdata));
   searchlimit = MIN(searchlimit, 10);
   /* multi variable potential: choose different disjoint neighborhoods, compare their potential */
   while( nsearched < searchlimit && nintegralvarsleft > 0 )
   {
      SCIP_VAR** neighborhood;
      SCIP_VAR* choosevar;
      int neighborhoodsize;
      int ndiscvarsneighborhood;
      int choosevardistance;

      ++nsearched;

      /* select a variable to start with randomly, but make sure it is active */
      choosevar = NULL;
      do
      {
         int index = SCIPgetRandomInt(0, nintegralvarsleft - 1, randseed);
         choosevar = varscopy[index];
         /* sort inactive variables to the end */
         if( SCIPvarGetProbindex(choosevar) < 0 )
         {
            varscopy[index] = varscopy[nintegralvarsleft - 1];
            --nintegralvarsleft;
         }
      }
      while( choosevar != NULL && SCIPvarGetProbindex(choosevar) < 0 && nintegralvarsleft > 0);

      /* if there was no variable chosen, there are no active variables left */
      if( choosevar == NULL || SCIPvarGetProbindex(choosevar) < 0 )
      {
         SCIPdebugMessage("No active variable left to perform breadth first search\n");
         break;
      }

      assert(SCIPvarIsIntegral(choosevar));

      /* get neighborhood storage */
      SCIP_CALL( SCIPallocBufferArray(scip, &neighborhood, nvars) );

      /* determine breadth-first distances from the chosen variable */
      variablegraphBreadthFirst(scip, vargraph, choosevar, distances, maxdistance);

      /* use either automatic or user-defined max-distance for neighborhood in variable constraint graph */
      if( heurdata->maxdistance != -1 )
      {
         choosevardistance = heurdata->maxdistance;
      }
      else
      {
         SCIP_CALL( determineMaxDistance(scip, heurdata, distances, &choosevardistance) );
      }

      ndiscvarsneighborhood = 0;
      neighborhoodsize = 0;

      /* loop over variables and determine neighborhood */
      for( v = nvars - 1; v >= 0; --v )
      {
         SCIP_VAR* currvar;
         currvar = vars[v];

         /* put variable in the neighborhood */
         if( isVariableInNeighborhood(currvar, distances, choosevardistance) )
         {
            neighborhood[neighborhoodsize++] = currvar;

            /* increase discrete variables counter */
            if( SCIPvarIsIntegral(currvar) )
               ++ndiscvarsneighborhood;
         }
      }

      /* check if neighborhood contains too many integer variables in order to satisfy the minimum fixing rate */
      if( ndiscvarsneighborhood >= nintegralvarsbound || ndiscvarsneighborhood <= 1 )
      {
         SCIPdebugMessage("Too many or too few discrete variables in neighboorhood: %d (%d)\n", ndiscvarsneighborhood, nbinvars + nintvars);
      }
      else
      {
         /* compare the neighborhood potential to the previous one */
         SCIP_Real potential;

         /* compare the neighborhood potential to the best potential found so far */
         potential = getPotential(scip, heurdata, sol, neighborhood, neighborhoodsize);

         /* big potential, take this variable */
         if( potential > maxpotential )
         {
            maxpotential = potential;
            selvar = choosevar;
            selvarmaxdistance = choosevardistance;
         }
      }

      /* sort neighborhood variables to the end so that this neighborhood is not considered in further samples */
      for( v = nintegralvarsleft - 1; v >= 0; --v )
      {
         SCIP_VAR* currvar;
         currvar = vars[v];

         if( isVariableInNeighborhood(currvar, distances, choosevardistance) )
         {
            varscopy[v] = varscopy[nintegralvarsleft - 1];
            --nintegralvarsleft;
         }
      }

      heurdata->sumdiscneighborhoodvars += ndiscvarsneighborhood;
      heurdata->sumneighborhoodvars += neighborhoodsize;
      ++heurdata->nneighborhoods;

      /* free current neighborhood */
      SCIPfreeBufferArray(scip, &neighborhood);
   }

   SCIPfreeBufferArray(scip, &varscopy);

   /* maybe no variable has a positive delta */
   if( !SCIPisPositive(scip, maxpotential) || selvar == NULL )
   {
      SCIPdebugMessage("Stopping with maxpotential %15.9f and selected variable %s\n",
         maxpotential, selvar != NULL ? SCIPvarGetName(selvar) : "none");
      *success = FALSE;
      goto TERMINATE;
   }

   assert(selvar != NULL);
   SCIPdebugMessage("Selected variable <%s> as central variable for a <%d>-neighborhood\n", SCIPvarGetName(selvar), selvarmaxdistance);

   /* collect distances in the variable graph of all variables to the selected variable */
   SCIP_CALL( variablegraphBreadthFirst(scip, vargraph, selvar, distances, selvarmaxdistance) );

   /* fix variables that are not in the neighborhood around the selected variable */
   SCIP_CALL( fixNonNeighborhoodVariables(scip, subscip, heurdata, sol, vars, subvars, distances, selvarmaxdistance, &nfixings) );

   fixthreshold = (int)(minfixingrate * (heurdata->fixcontvars ? nvars : (nbinvars + nintvars)));

   /* compare actual number of fixings to limit; if we fixed not enough variables we terminate here; we also terminate if no discrete variables are left */
   if( nfixings < fixthreshold )
   {
      SCIPdebugMessage("Fixed %d/%d variables in vardegree heuristic, stopping", nfixings, fixthreshold);
      *success = FALSE;
      goto TERMINATE;

   }

   /* create problem constraints from rows. If uselprows is FALSE, the entire problem was copied already, including constraints */
   if( uselprows )
   {
      SCIP_CALL( createConstraintsFromRows(scip, subscip, subvars) );
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &distances);
   SCIP_CALL( variableGraphFree(scip, &vargraph) );
   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< vardegree heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
)
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyVardegree)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurVardegree(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeVardegree)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitVardegree)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->randseed = 0;
   heurdata->sumdiscneighborhoodvars = heurdata->sumneighborhoodvars = 0;
   heurdata->nneighborhoods = 0;

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEUREXIT(heurExitVardegree)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPstatistic(
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Vardegree: Avg Neighborhood size: %.1f Avg. discrete neighboorhood vars: %.1f\n", heurdataAvgNeighborhoodSize(heurdata), heurdataAvgDiscreteNeighborhoodSize(heurdata));
      )

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecVardegree)
{  /*lint --e{715}*/
   SCIP_Longint maxnnodes;
   SCIP_Longint nsubnodes;                   /* node limit for the subproblem */

   SCIP_HEURDATA* heurdata;                  /* heuristic's data */
   SCIP* subscip;                            /* the subproblem created by vardegree */
   SCIP_VAR** vars;                          /* original problem's variables */
   SCIP_VAR** subvars;                       /* subproblem's variables */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */

   SCIP_Real cutoff;                         /* objective cutoff for the subproblem */
   SCIP_Real maxnnodesr;
   SCIP_Real memorylimit;
   SCIP_Real timelimit;                      /* timelimit for the subproblem */
   SCIP_Real upperbound;

   int nvars;                                /* number of original problem's variables */
   int i;

   SCIP_Bool success;

   SCIP_RETCODE retcode;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* only call heuristic, if feasible solution is available */
   if( SCIPgetNSols(scip) <= 0 )
      return SCIP_OKAY;

   /* only call heuristic, if the best solution comes from transformed problem */
   assert( SCIPgetBestSol(scip) != NULL );
   if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip,SCIPgetBestSol(scip))  < heurdata->nwaitingnodes)
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if discrete variables are present */
   if( SCIPgetNBinVars(scip) == 0 && SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodesr = heurdata->nodesquot * SCIPgetNNodes(scip);

   /* reward vardegree if it succeeded often, count the setup costs for the sub-MIP as 100 nodes */
   maxnnodesr *= 1.0 + 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0);
   maxnnodes = (SCIP_Longint) maxnnodesr - 100 * SCIPheurGetNCalls(heur);
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nsubnodes < heurdata->minnodes )
       return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_vardegreesub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_vardegreesub", SCIPgetProbName(scip));

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
   }
   else
   {
      SCIP_Bool valid;
      valid = FALSE;

      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "vardegree", TRUE, FALSE, TRUE, &valid) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }

      SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");
   }

   for( i = 0; i < nvars; i++ )
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* create a new problem, by fixing all variables except for a small neighborhood */
   SCIP_CALL( createSubproblem(scip, subscip, subvars, heurdata, heurdata->minfixingrate, &heurdata->randseed,
         heurdata->uselprows, &success) );

   /* terminate if it was not possible to create the subproblem */
   if( !success )
   {
      SCIPdebugMessage("Could not create the subproblem -> skip call\n");
      goto TERMINATE;
   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

  /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      goto TERMINATE;

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   if( !SCIPisParamFixed(subscip, "conflict/useprop") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useinflp") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/usesb") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/usepseudo") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );
   }

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* add an objective cutoff */
   cutoff = SCIPinfinity(scip);
   assert( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) );

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
   {
      cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip)
            + heurdata->minimprove * SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound(scip) >= 0 )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
      else
         cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
   }
   cutoff = MIN(upperbound, cutoff);
   SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));

   /* solve the subproblem */
   SCIPdebugMessage("Solve Vardegree subMIP\n");
   retcode = SCIPsolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in Vardegree heuristic; sub-SCIP terminated with code <%d>\n",retcode);
   }
   else
   {
      /* transfer variable statistics from sub-SCIP */
      SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );
   }

   heurdata->usednodes += SCIPgetNNodes(subscip);

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      }
      if( success )
         *result = SCIP_FOUNDSOL;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the vardegree primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurVardegree(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Vardegree primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecVardegree, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyVardegree) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeVardegree) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitVardegree) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitVardegree) );

   /* add vardegree primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, SCIPsumepsilon(scip), 1.0-SCIPsumepsilon(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/fixcontvars",
         "should continuous variables outside the neighborhoods be fixed?",
         &heurdata->fixcontvars, TRUE, DEFAULT_FIXCONTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxdistance",
         "maximum distance to selected variable to enter the subproblem, or -1 to select the distance "
         "that best approximates the minimum fixing rate from below",
         &heurdata->maxdistance, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/potential",
         "the reference point to compute the neighborhood potential: (r)oot or (p)seudo solution",
         &heurdata->potential, TRUE, DEFAULT_POTENTIAL, "pr", NULL, NULL) );

   return SCIP_OKAY;
}
