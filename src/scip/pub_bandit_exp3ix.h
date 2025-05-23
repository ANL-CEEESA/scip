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

/**@file   pub_bandit_exp3ix.h
 * @ingroup PublicBanditMethods
 * @brief  public methods for Exp.3-IX
 * @author Antonia Chmiela
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_PUB_BANDIT_EXP3IX_H_
#define SRC_SCIP_PUB_BANDIT_EXP3IX_H_

#include "scip/def.h"
#include "scip/type_bandit.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBanditMethods
 *
 * ## Exp.3-IX TODO: update
 *
 * Exp.3 is a randomized selection method for the multi-armed bandit problem
 *
 * Exp3 maintains a probability distribution
 * according to which an action is drawn
 * in every iteration.
 * The probability distribution is a mixture between
 * a uniform distribution and a softmax distribution
 * based on the cumulative rewards of the actions.
 * The weight of the uniform distribution in the mixture
 * is controlled by the parameter \f$ \gamma \f$, ie.,
 * setting \f$ \gamma = 1\f$ uses a uniform distribution
 * in every selection step.
 * The cumulative reward for the actions can be
 * fine-tuned by adding a general bias for all actions.
 * The bias is given by the parameter \f$ \beta \f$.
 *
 * @{
 */

/** creates and resets an Exp.3-IX bandit algorithm using \p scip pointer */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateBanditExp3IX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         exp3ix,             /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial seed for random number generation */
   );

/** returns probability to play an action */
SCIP_EXPORT
SCIP_Real SCIPgetProbabilityExp3IX(
   SCIP_BANDIT*          exp3ix,             /**< bandit algorithm */
   int                   action              /**< index of the requested action */
   );

/** @}*/

#ifdef __cplusplus
}
#endif

#endif
