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

/**@file   presol_donotmultagg.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  donotmultaggr presolver: mark donotmultaggr=TRUE for certain variables
 * @author Suresh Bolusani
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/presol_donotmultaggr.h"
#include "scip/pub_message.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_presol.h"
#include "scip/scip_prob.h"
#include "scip/scip_var.h"
#include <string.h>

#define PRESOL_NAME            "donotmultaggr"
#define PRESOL_DESC            "mark donotmultaggr=TRUE for certain variables"
#define PRESOL_PRIORITY       -10000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDonotmultaggr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDonotmultaggr(scip) );

   return SCIP_OKAY;
}


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecDonotmultaggr)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* get the problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* scan the variables for marking donotmultaggr=TRUE as needed
    * (loop backwards, since the variables that need to be marked are present at the end)
    */
   for( v = nvars-1; v >= 0; --v )
   {
      /* is variable integral with its name starting with "t_zdisj" or "zdisj"? */
      if( (SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS) && ((strncmp(SCIPvarGetName(vars[v]), "t_zdisj", 7) == 0)
               || (strncmp(SCIPvarGetName(vars[v]), "zdisj", 5) == 0)) )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, vars[v]) );
         SCIPdebugMsg(scip, "Marked donotmultaggr=TRUE for the variable <%s>\n", SCIPvarGetName(vars[v]));
      }
   }

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the donotmultaggr presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDonotmultaggr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presolptr;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING, presolExecDonotmultaggr, NULL) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyDonotmultaggr) );

   return SCIP_OKAY;
}
