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

/**@file   nodesel_vanillabfs.c
 * @ingroup DEFPLUGINS_NODESEL
 * @brief  node selector for vanilla best first search
 * @author Suresh Bolusani
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nodesel_vanillabfs.h"
#include "scip/pub_message.h"
#include "scip/pub_nodesel.h"
#include "scip/pub_tree.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/type_misc.h"
#include <string.h>

#define NODESEL_NAME             "vanillabfs"
#define NODESEL_DESC             "vanilla best first search"
#define NODESEL_STDPRIORITY      -20000
#define NODESEL_MEMSAVEPRIORITY       0

/*
 * Callback methods
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyVanillabfs)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeNodeselVanillabfs(scip) );

   return SCIP_OKAY;
}

/** destructor of node selector to free user data (called when SCIP is exiting) */
/**! [SnippetNodeselFreeVanillabfs] */
static
SCIP_DECL_NODESELFREE(nodeselFreeVanillabfs)
{  /*lint --e{715}*/
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   return SCIP_OKAY;
}
/**! [SnippetNodeselFreeVanillabfs] */


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectVanillabfs)
{  /*lint --e{715}*/
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   *selnode = NULL;

   /* get the best node */
   *selnode = SCIPgetBestNode(scip);
   SCIPdebugMsg(scip, "  -> best node   : lower=%g\n", *selnode != NULL ? SCIPnodeGetLowerbound(*selnode) : SCIPinfinity(scip));

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompVanillabfs)
{  /*lint --e{715}*/
   SCIP_Real lowerbound1;
   SCIP_Real lowerbound2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   lowerbound1 = SCIPnodeGetLowerbound(node1);
   lowerbound2 = SCIPnodeGetLowerbound(node2);
   if( SCIPisLT(scip, lowerbound1, lowerbound2) )
      return -1;
   else if( SCIPisGT(scip, lowerbound1, lowerbound2) )
      return +1;
   else
   {
      SCIP_Real estimate1;
      SCIP_Real estimate2;

      estimate1 = SCIPnodeGetEstimate(node1);
      estimate2 = SCIPnodeGetEstimate(node2);
      if( (SCIPisInfinity(scip,  estimate1) && SCIPisInfinity(scip,  estimate2)) ||
          (SCIPisInfinity(scip, -estimate1) && SCIPisInfinity(scip, -estimate2)) ||
          SCIPisEQ(scip, estimate1, estimate2) )
      {
         SCIP_NODETYPE nodetype1;
         SCIP_NODETYPE nodetype2;

         nodetype1 = SCIPnodeGetType(node1);
         nodetype2 = SCIPnodeGetType(node2);
         if( nodetype1 == SCIP_NODETYPE_CHILD && nodetype2 != SCIP_NODETYPE_CHILD )
            return -1;
         else if( nodetype1 != SCIP_NODETYPE_CHILD && nodetype2 == SCIP_NODETYPE_CHILD )
            return +1;
         else if( nodetype1 == SCIP_NODETYPE_SIBLING && nodetype2 != SCIP_NODETYPE_SIBLING )
            return -1;
         else if( nodetype1 != SCIP_NODETYPE_SIBLING && nodetype2 == SCIP_NODETYPE_SIBLING )
            return +1;
         else
         {
            int depth1;
            int depth2;

            depth1 = SCIPnodeGetDepth(node1);
            depth2 = SCIPnodeGetDepth(node2);
            if( depth1 < depth2 )
               return -1;
            else if( depth1 > depth2 )
               return +1;
            else
               return 0;
         }
      }

      if( SCIPisLT(scip, estimate1, estimate2) )
         return -1;

      assert(SCIPisGT(scip, estimate1, estimate2));
      return +1;
   }
}


/*
 * vanillabfs specific interface methods
 */

/** creates the node selector for vanilla best first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselVanillabfs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESEL* nodesel;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectVanillabfs, nodeselCompVanillabfs, NULL) );

   assert(nodesel != NULL);

   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyVanillabfs) );
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreeVanillabfs) );

   return SCIP_OKAY;
}
