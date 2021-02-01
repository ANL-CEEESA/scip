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

/**@file   lpexact.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for exact LP management
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LPEXACT_H__
#define __SCIP_LPEXACT_H__


#include <stdio.h>

#include "scip/def.h"
#include "scip/rational.h"
#include "blockmemshell/memory.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_misc.h"
#include "scip/type_lpexact.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/pub_lp.h"

#include "scip/struct_lpexact.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Column methods
 */

/** checks if the exact column and its fpcol are consistent */
SCIP_Bool colExactInSync(
   SCIP_COLEXACT*        colexact,           /**< exact column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   );

/** creates an LP column */
SCIP_RETCODE SCIPcolExactCreate(
   SCIP_COLEXACT**       col,                /**< pointer to column data */
   SCIP_COL*             fpcol,              /**< the corresponding fp col */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROWEXACT**       rows,               /**< array with rows of column entries */
   SCIP_Rational**       vals,               /**< array with coefficients of column entries */
   SCIP_Bool             removable           /**< should the column be removed from the LP due to aging or cleanup? */
   );

/** frees an LP column */
SCIP_RETCODE SCIPcolExactFree(
   SCIP_COLEXACT**       col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   );

/** output column to file stream */
void SCIPcolExactPrint(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** adds a previously non existing coefficient to an LP column */
SCIP_RETCODE SCIPcolExactAddCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   );

/** deletes coefficient from column */
SCIP_RETCODE SCIPcolExactDelCoef(
   SCIP_COLEXACT*        col,                /**< column to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row                 /**< coefficient to be deleted */
   );

/** changes or adds a coefficient to an LP column */
SCIP_RETCODE SCIPcolExactChgCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        val                 /**< value of coefficient */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
SCIP_RETCODE SCIPcolExactIncCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        incval              /**< value to add to the coefficient */
   );

/** changes objective value of column */
SCIP_RETCODE SCIPcolExactChgObj(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_Rational*        newobj              /**< new objective value */
   );

/** changes lower bound of column */
SCIP_RETCODE SCIPcolExactChgLb(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_Rational*        newlb               /**< new lower bound value */
   );

/** changes upper bound of column */
SCIP_RETCODE SCIPcolExactChgUb(
   SCIP_COLEXACT*        col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_Rational*        newub               /**< new upper bound value */
   );

/*
 * Row methods
 */

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwExactSort(
   SCIP_ROWEXACT*        row                 /**< row to be sorted */
   );

/** increases usage counter of LP row */
void SCIProwExactCapture(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** output column to file stream */
void SCIProwExactPrint(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** get the index of an exact row */
int SCIProwExactGetIndex(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** get the length of a row */
int SCIProwExactGetNNonz(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** returns TRUE iff row is member of current LP */
SCIP_Bool SCIProwExactIsInLP(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** return TRUE iff row is modifiable */
SCIP_Bool SCIProwExactIsModifiable(
   SCIP_ROWEXACT*        row                 /**< LP row */
   );

/** returns true, if an exact row for this fprow was already created */
SCIP_Bool SCIProwHasExRow(
   SCIP_LPEXACT*         lpexact,            /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   );

/** returns exact row corresponding to fprow, if it exists. Otherwise returns NULL */
SCIP_ROWEXACT* SCIProwGetExRow(
   SCIP_LPEXACT*         lpexact,            /**< exact lp data structure */
   SCIP_ROW*             row                 /**< SCIP row */
   );

/** changes left hand side of exact LP row */
SCIP_RETCODE SCIProwExactChgLhs(
   SCIP_ROWEXACT*        row,                /**< exact LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_Rational*        lhs                 /**< new left hand side */
   );

/** changes right hand side of exact LP row */
SCIP_RETCODE SCIProwExactChgRhs(
   SCIP_ROWEXACT*        row,                /**< exact LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_Rational*        rhs                 /**< new right hand side */
   );

/** returns exact col corresponding to fpcol, if it exists. Otherwise returns NULL */
SCIP_COLEXACT* SCIPcolGetColExact(
   SCIP_COL*             col                 /**< SCIP col */
   );

/** calculates the Farkas coefficient or reduced cost of a column i using the given dual Farkas vector y */
void SCIPcolExactCalcFarkasRedcostCoef(
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_SET*             set,                /**< SCIP settings pointer */
   SCIP_Rational*        result,             /**< rational to store the result */
   SCIP_Rational**       dual,               /**< dense dual Farkas vector, NULL to use internal row-values */
   SCIP_Bool             usefarkas           /**< should the farkas coefficient be computed ? */
   );

/** creates and captures an LP row */
SCIP_RETCODE SCIProwCreateExact(
   SCIP_ROWEXACT**       row,                /**< pointer to LP row data */
   SCIP_ROW*             fprow,              /**< corresponding fp row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COLEXACT**       cols,               /**< array with columns of row entries */
   SCIP_Rational**       vals,               /**< array with coefficients of row entries */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_ROWORIGINTYPE    origintype,         /**< type of origin of row */
   void*                 origin              /**< pointer to constraint handler or separator who created the row (NULL if unkown) */
   );

/** applies all cached changes to the LP solver */
SCIP_RETCODE SCIPlpExactFlush(
   SCIP_LPEXACT*         lp,                 /**< current exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/*
 * lp methods
 */

/** returns whether it is possible to use neumair-shcherbina bounding method */
SCIP_Bool SCIPlpExactBSpossible(
   SCIP_LPEXACT*         lp                  /**< pointer to LP data object */
   );

/** returns whether it is possible to use project and shift bounding method */
SCIP_Bool SCIPlpExactPSpossible(
   SCIP_LPEXACT*         lp                  /**< pointer to LP data object */
   );

/** checks that lp and fplp are properly synced */
SCIP_Bool SCIPlpExactIsSynced(
   SCIP_LPEXACT*         lp,                 /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     msg                 /**< message handler */
   );

/** creates empty LP data object */
SCIP_RETCODE SCIPlpExactCreate(
   SCIP_LPEXACT**        lpexact,            /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              fplp,               /**< the normal floating point lp */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   );

/** frees LP data object */
SCIP_RETCODE SCIPlpExactFree(
   SCIP_LPEXACT**        lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** adds a column to the LP and captures the variable */
SCIP_RETCODE SCIPlpExactAddCol(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COLEXACT*        col,                /**< LP column */
   int                   depth               /**< depth in the tree where the column addition is performed */
   );

/** adds a row to the LP and captures it */
SCIP_RETCODE SCIPlpExactAddRow(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_ROWEXACT*        rowexact,           /**< LP row */
   int                   depth               /**< depth in the tree where the row addition is performed */
   );

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpExactShrinkRows(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newnrows            /**< new number of rows in the LP */
   );

/** returns the feasibility of a row for the given solution */
void SCIProwExactGetSolFeasibility(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   );

/** returns the activity of a row for a given solution */
void SCIProwExactGetSolActivity(
   SCIP_ROWEXACT*        rowexact,           /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             useexact,           /**< should an exact solution be used */
   SCIP_Rational*        result              /**< resulting activity */
   );

/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIProwExactRelease(
   SCIP_ROWEXACT**       row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   );

/** frees an LP row */
SCIP_RETCODE SCIProwExactFree(
   SCIP_ROWEXACT**       row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   );

/** ensuresr, that column array of row can store at least num entries */
SCIP_RETCODE SCIProwExactEnsureSize(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** add constant value to a row */
SCIP_RETCODE SCIProwExactAddConstant(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_Rational*        addval              /**< constant value to add to the row */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
SCIP_RETCODE SCIProwExactIncCoef(
   SCIP_ROWEXACT*        row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_COLEXACT*        col,                /**< LP column */
   SCIP_Rational*        incval              /**< valpelue to add to the coefficient */
   );

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
void SCIProwExactGetLPFeasibility(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_Rational*        result              /**< rational pointer to store the result */
   );

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
void SCIProwExactGetPseudoFeasibility(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Rational*        result              /**< rational pointer to store the result */
   );

/** returns the activity of a row in the current LP solution */
SCIP_Rational* SCIProwExactGetLPActivity(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lp                  /**< current LP data */
   );

/** returns the pseudo activity of a row in the current pseudo solution */
SCIP_Rational* SCIProwExactGetPseudoActivity(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** enables delaying of row sorting */
void SCIProwExactDelaySort(
   SCIP_ROWEXACT*        rowexact            /**< LP rowexact */
   );

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
void SCIProwExactForceSort(
   SCIP_ROWEXACT*        rowexact,           /**< LP rowexact */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** recalculates the current activity of a row */
void SCIProwExactRecalcLPActivity(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

 /** calculates the current pseudo activity of a row */
void SCIProwExactRecalcPseudoActivity(
   SCIP_ROWEXACT*        row,                /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** gets objective value of column */
SCIP_Rational* SCIPcolExactGetObj(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/** gets lower bound of column */
SCIP_Rational* SCIPcolExactGetLb(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/** gets upper bound of column */
SCIP_Rational* SCIPcolExactGetUb(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/** gets best bound of column with respect to the objective function */
SCIP_Rational* SCIPcolExactGetBestBound(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/** gets the primal LP solution of a column */
SCIP_Rational* SCIPcolExactGetPrimsol(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/** gets the minimal LP solution value, this column ever assumed */
SCIP_Rational* SCIPcolExactGetMinPrimsol(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/** gets the maximal LP solution value, this column ever assumed */
SCIP_Rational* SCIPcolExactGetMaxPrimsol(
   SCIP_COLEXACT*        col                 /**< LP column */
   );

/*
 * lp update methods
 */

/** updates current pseudo and loose objective value for a change in a variable's objective value */
SCIP_RETCODE SCIPlpExactUpdateVarObj(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldobj,             /**< old objective value of variable */
   SCIP_Rational*        newobj              /**< new objective value of variable */
   );

/** updates current root pseudo objective value for a global change in a variable's lower bound */
SCIP_RETCODE SCIPlpExactUpdateVarLbGlobal(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb               /**< new lower bound of variable */
   );

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
SCIP_RETCODE SCIPlpExactUpdateVarLb(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldlb,              /**< old lower bound of variable */
   SCIP_Rational*        newlb               /**< new lower bound of variable */
   );

/** updates current root pseudo objective value for a global change in a variable's upper bound */
SCIP_RETCODE SCIPlpExactUpdateVarUbGlobal(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub               /**< new upper bound of variable */
   );

/** updates current pseudo objective value for a change in a variable's upper bound */
SCIP_RETCODE SCIPlpExactUpdateVarUb(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Rational*        oldub,              /**< old upper bound of variable */
   SCIP_Rational*        newub               /**< new upper bound of variable */
   );

/** informs LP, that given variable was added to the problem */
SCIP_RETCODE SCIPlpExactUpdateAddVar(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   );

/** informs LP, that given variable is to be deleted from the problem */
SCIP_RETCODE SCIPlpExactUpdateDelVar(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   );

/** informs LP, that given formerly loose problem variable is now a column variable */
SCIP_RETCODE SCIPlpExactUpdateVarColumn(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   );

/** informs LP, that given formerly column problem variable is now again a loose variable */
SCIP_RETCODE SCIPlpExactUpdateVarLoose(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   );

/** decrease the number of loose variables by one */
void SCIPlpExactDecNLoosevars(
   SCIP_LPEXACT*         lp                  /**< current LP data */
   );

SCIP_RETCODE SCIPlexGetNRows(
   SCIP_LPEXACT*         lp                  /**< current LP data */
   );

/** stores the LP solution in the columns and rows */
SCIP_RETCODE SCIPlpExactGetSol(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible,       /**< pointer to store whether the solution is dual feasible, or NULL */
   SCIP_Bool             overwritefplp       /**< should the floating point values be overwritten, e.g. if fp lp was infeasible */
   );

/** stores LP solution with infinite objective value in the columns and rows */
SCIP_RETCODE SCIPlpExactGetUnboundedSol(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            rayfeasible         /**< pointer to store whether the primal ray is a feasible unboundedness proof, or NULL */
   );

/** returns primal ray proving the unboundedness of the current LP */
SCIP_RETCODE SCIPlpExactGetPrimalRay(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Rational**       ray                 /**< array for storing primal ray values, they are stored w.r.t. the problem index of the variables,
                                              *   so the size of this array should be at least number of active variables
                                              *   (all entries have to be initialized to 0 before) */
   );

/** stores the dual Farkas multipliers for infeasibility proof in rows. besides, the proof is checked for validity if
 *  lp/checkfarkas = TRUE.
 *
 *  @note the check will not be performed if @p valid is NULL.
 */
SCIP_RETCODE SCIPlpExactGetDualfarkas(
   SCIP_LPEXACT*         lpexact,            /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            valid,              /**< pointer to store whether the Farkas proof is valid  or NULL */
   SCIP_Bool             overwritefplp       /**< should the floating point values be overwritten, e.g. if fp lp was infeasible */
   );

/** get number of iterations used in last LP solve */
SCIP_RETCODE SCIPlpExactGetIterations(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   int*                  iterations          /**< pointer to store the iteration count */
   );

/** gets objective value of current LP
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status is
 *        SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 */
void SCIPlpExactGetObjval(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational*        res                 /**< result pointer to store rational */
   );

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 */
void SCIPlpExactGetPseudoObjval(
   SCIP_LPEXACT*         lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Rational*        res                 /**< result pointer to store rational */
   );

/** removes all columns after the given number of cols from the LP */
SCIP_RETCODE SCIPlpExactshrinkCols(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newncols            /**< new number of columns in the LP */
   );

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpExactshrinkRows(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   newnrows            /**< new number of rows in the LP */
   );

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpExactReset(
   SCIP_LPEXACT*         lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** removes all columns and rows from LP, releases all rows */
SCIP_RETCODE SCIPlpExactClear(
   SCIP_LPEXACT*         lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** checks whether primal solution satisfies all integrality restrictions exactly.
 * This checks either the fp solution exactly or checks the exact solution, if one exists.
 */
SCIP_RETCODE SCIPlpExactcheckIntegralityExact(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_RESULT*          result              /**< result pointer */
   );

/** forces an exact lp to be solved in the next exact bound computation */
void SCIPlpExactForceExactSolve(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   );


/** gets solution status of current exact LP */
SCIP_LPSOLSTAT SCIPlpExactGetSolstat(
   SCIP_LPEXACT*         lpexact              /**< current LP data */
   );

/** sets the upper objective limit of the exact LP solver */
SCIP_RETCODE SCIPlpExactSetCutoffbound(
   SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             cutoffbound         /**< new upper objective limit */
   );

/** solves the LP with simplex algorithm, and copy the solution into the column's data */
SCIP_RETCODE SCIPlpExactSolveAndEval(
   SCIP_LPEXACT*         lpexact,            /**< LP data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool             usefarkas           /**< are we aiming to prove infeasibility? */
   );

/** stores exact LP state (like basis information) into LP state object */
SCIP_RETCODE SCIPlpExactGetState(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** loads exact LP state (like basis information) into solver */
SCIP_RETCODE SCIPlpExactSetState(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPISTATE*        lpistate,           /**< LP state information (like basis information) */
   SCIP_Bool             wasprimfeas,        /**< primal feasibility when LP state information was stored */
   SCIP_Bool             wasprimchecked,     /**< true if the LP solution has passed the primal feasibility check */
   SCIP_Bool             wasdualfeas,        /**< dual feasibility when LP state information was stored */
   SCIP_Bool             wasdualchecked      /**< true if the LP solution has passed the dual feasibility check */
   );

/** frees exact LP state information */
SCIP_RETCODE SCIPlpExactFreeState(
   SCIP_LPEXACT*         lpexact,            /**< exact LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** gets solution status of current exact LP */
SCIP_LPSOLSTAT SCIPlpExactGetSolstat(
   SCIP_LPEXACT*         lpexact             /**< current LP exact data */
    );

/** starts exact LP diving and saves bounds and objective values of columns to the current nodes's values */
SCIP_RETCODE SCIPlpExactStartDive(
    SCIP_LPEXACT*         lpexact,            /**< current exact LP data */
    BMS_BLKMEM*           blkmem,             /**< block memory */
    SCIP_SET*             set,                /**< global SCIP settings */
    SCIP_STAT*            stat                /**< problem statistics */
    );

/** quits exact LP diving and resets bounds and objective values of columns to the current node's values */
SCIP_RETCODE SCIPlpExactEndDive(
    SCIP_LPEXACT*         lpexact,            /**< current LP data */
    BMS_BLKMEM*           blkmem,             /**< block memory */
    SCIP_SET*             set,                /**< global SCIP settings */
    SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
    SCIP_STAT*            stat,               /**< problem statistics */
    SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
    SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
    SCIP_PROB*            prob,               /**< problem data */
    SCIP_VAR**            vars,               /**< array with all active variables */
    int                   nvars               /**< number of active variables */
    );

/** returns whether the exact LP is in diving mode */
SCIP_EXPORT
SCIP_Bool SCIPlpExactDiving(
   SCIP_LPEXACT*         lpexact             /**< current exact LP data */
   );

#ifdef __cplusplus
}
#endif

#endif
