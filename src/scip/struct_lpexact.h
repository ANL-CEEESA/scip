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

/**@file   struct_lpexact.h
 * @ingroup INTERNALAPI
 * @brief  data structures for exact LP management
 * @author Tobias Achterberg
 *
 *  In SCIP, the LP is defined as follows:
 *
 *     min       obj * x
 *        lhs <=   A * x + const <= rhs
 *        lb  <=       x         <= ub
 *
 *  The row activities are defined as activity = A * x + const and must
 *  therefore be in the range of [lhs,rhs].
 *
 *  Mathematically, each range constraint would account for two dual
 *  variables, one for each inequality. Since in an optimal solution (at
 *  least) one of them may be chosen to be zero, we may define one dual
 *  multiplier for each row as the difference of those two.
 *
 *  Let y be the vector of dual multipliers for the rows, then the reduced
 *  costs are defined as
 *
 *     redcost = obj - A^T * y.
 *
 *  In an optimal solution, y must be
 *
 *     - nonnegative, if the corresponding row activity is not tight at its rhs
 *     - nonpositive, if the corresponding row activity is not tight at its lhs
 *     - zero, if the corresponding row activity is not at any of its sides
 *
 *  and the reduced costs must be
 *
 *     - nonnegative, if the corresponding variable is not tight at its ub
 *     - nonpositive, if the corresponding variable is not tight at its lb
 *     - zero, if the corresponding variable is not at any of its bounds.
 *
 *  The main datastructures for storing an LP are the rows and the columns.
 *  A row can live on its own (if it was created by a separator), or as SCIP_LP
 *  relaxation of a constraint. Thus, it has a nuses-counter, and is
 *  deleted, if not needed any more.
 *  A column cannot live on its own. It is always connected to a problem
 *  variable. Because pricing is always problem specific, it cannot create
 *  LP columns without introducing new variables. Thus, each column is
 *  connected to exactly one variable, and is deleted, if the variable
 *  is deleted.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_LPEXACT_H__
#define __SCIP_STRUCT_LPEXACT_H__


#include "scip/def.h"
#include "scip/scip_exact.h"
#include "scip/type_lp.h"
#include "scip/type_lpexact.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"
#include "scip/type_event.h"
#include "lpi/type_lpi.h"
#include "lpiexact/type_lpiexact.h"
#include "scip/rational.h"
#include "rectlu/rectlu.h"

#ifdef __cplusplus
extern "C" {
#endif

/** LP column;
 *  The row vector of the LP column is partitioned into two parts: The first col->nlprows rows in the rows array
 *  are the ones that belong to the current LP (col->rows[j]->lppos >= 0) and that are linked to the column
 *  (col->linkpos[j] >= 0). The remaining col->len - col->nlprows rows in the rows array are the ones that
 *  don't belong to the current LP (col->rows[j]->lppos == -1) or that are not linked to the column
 *  (col->linkpos[j] == -1).
 */
struct SCIP_ColExact
{
   SCIP_COL*             fpcol;              /**< the floating point column corresponding to this exact column */
   SCIP_ROWEXACT**       rows;               /**< the exact rows corresponding the this exact column */
   SCIP_Rational*        obj;                /**< current objective value of column in LP (might be changed in diving or probing) */
   SCIP_Rational*        lb;                 /**< current lower bound of column in LP */
   SCIP_Rational*        ub;                 /**< current upper bound of column in LP */
   SCIP_Rational*        flushedobj;         /**< objective value of column already flushed to the LP solver */
   SCIP_Rational*        flushedlb;          /**< lower bound of column already flushed to the LP solver */
   SCIP_Rational*        flushedub;          /**< upper bound of column already flushed to the LP solver */
   SCIP_Rational*        primsol;            /**< primal solution value in LP, is 0 if col is not in LP */
   SCIP_Rational*        redcost;            /**< reduced cost value in LP, or SCIP_INVALID if not yet calculated */
   SCIP_Rational*        farkascoef;         /**< coefficient in dual Farkas infeasibility proof (== dualfarkas^T A_c) */
   SCIP_VAR*             var;                /**< variable, this column represents; there cannot be a column without variable */
   SCIP_Rational**       vals;               /**< coefficients of column entries */
   SCIP_Longint          validredcostlp;     /**< LP number for which reduced cost value is valid */
   SCIP_Longint          validfarkaslp;      /**< LP number for which Farkas coefficient is valid */
   int*                  linkpos;            /**< position of col in col vector of the row, or -1 if not yet linked */
   int                   index;              /**< consecutively numbered column identifier */
   int                   size;               /**< size of the row- and val-arrays */
   int                   len;                /**< number of nonzeros in column */
   int                   nlprows;            /**< number of linked rows in column, that belong to the current LP */
   int                   nunlinked;          /**< number of column entries, where the rows don't know about the column */
   int                   lppos;              /**< column position number in current LP, or -1 if not in current LP */
   int                   lpipos;             /**< column position number in LP solver, or -1 if not in LP solver */
   unsigned int          basisstatus:2;      /**< basis status of column in last LP solution, invalid for non-LP columns */
   unsigned int          lprowssorted:1;     /**< are the linked LP rows in the rows array sorted by non-decreasing index? */
   unsigned int          nonlprowssorted:1;  /**< are the non-LP/not linked rows sorted by non-decreasing index? */
   unsigned int          objchanged:1;       /**< has objective value changed, and has data of LP solver to be updated? */
   unsigned int          lbchanged:1;        /**< has lower bound changed, and has data of LP solver to be updated? */
   unsigned int          ubchanged:1;        /**< has upper bound changed, and has data of LP solver to be updated? */
   unsigned int          coefchanged:1;      /**< has the coefficient vector changed, and has LP solver to be updated? */
};

/** LP row
 *  The column vector of the LP row is partitioned into two parts: The first row->nlpcols columns in the cols array
 *  are the ones that belong to the current LP (row->cols[j]->lppos >= 0) and that are linked to the row
 *  (row->linkpos[j] >= 0). The remaining row->len - row->nlpcols columns in the cols array are the ones that
 *  don't belong to the current LP (row->cols[j]->lppos == -1) or that are not linked to the row
 *  (row->linkpos[j] == -1).
 */
struct SCIP_RowExact
{
   SCIP_ROW*             fprow;              /**< pointer to the corresponding row in the fp lp */
   SCIP_Rational*        constant;           /**< constant shift c in row lhs <= ax + c <= rhs */
   SCIP_Rational*        lhs;                /**< left hand side of row */
   SCIP_Rational*        rhs;                /**< right hand side of row */
   SCIP_Rational*        flushedlhs;         /**< left hand side minus constant of row already flushed to the LP solver */
   SCIP_Rational*        flushedrhs;         /**< right hand side minus constant of row already flushed to the LP solver */
   SCIP_Rational*        objprod;            /**< scalar product of row vector with objective function */
   SCIP_Rational*        dualsol;            /**< dual solution value in LP, is 0 if row is not in LP */
   SCIP_Rational*        activity;           /**< row activity value in LP, or SCIP_INVALID if not yet calculated */
   SCIP_Rational*        dualfarkas;         /**< multiplier value in dual Farkas infeasibility proof */
   SCIP_Rational*        pseudoactivity;     /**< row activity value in pseudo solution, or SCIP_INVALID if not yet calculated */
   SCIP_Rational**       vals;               /**< coefficients of row entries */
   SCIP_INTERVAL*        valsinterval;       /**< interval-array of coefficients rounded up and down, respectively */
   SCIP_COLEXACT**       cols;               /**< columns of row entries, that may have a nonzero primal solution value */
   int*                  cols_index;         /**< copy of cols[i]->index for avoiding expensive dereferencing */
   int*                  linkpos;            /**< position of row in row vector of the column, or -1 if not yet linked */
   SCIP_Longint          validactivitylp;    /**< LP number for which activity value is valid */
   int                   index;              /**< consecutively numbered row identifier */
   int                   size;               /**< size of the col- and val-arrays */
   int                   len;                /**< number of nonzeros in row */
   int                   nlpcols;            /**< number of linked columns in row, that belong to the current LP */
   int                   nunlinked;          /**< number of row entries, where the columns don't know about the row */
   int                   nuses;              /**< number of times, this row is referenced */
   int                   lppos;              /**< row position number in current LP, or -1 if not in current LP */
   int                   lpipos;             /**< row position number in LP solver, or -1 if not in LP solver */
   int                   lpdepth;            /**< depth level at which row entered the LP, or -1 if not in current LP */
   unsigned int          basisstatus:2;      /**< basis status of row in last LP solution, invalid for non-LP rows */
   unsigned int          lpcolssorted:1;     /**< are the linked LP columns in the cols array sorted by non-decreasing index? */
   unsigned int          nonlpcolssorted:1;  /**< are the non-LP/not linked columns sorted by non-decreasing index? */
   unsigned int          delaysort:1;        /**< should the row sorting be delayed and done in a lazy fashion? */
   unsigned int          lhschanged:1;       /**< was left hand side or constant changed, and has LP solver to be updated? */
   unsigned int          rhschanged:1;       /**< was right hand side or constant changed, and has LP solver to be updated? */
   unsigned int          coefchanged:1;      /**< was the coefficient vector changed, and has LP solver to be updated? */
   unsigned int          integral:1;         /**< is activity (without constant) of row always integral in feasible solution? */
   unsigned int          nlocks:15;          /**< number of sealed locks of an unmodifiable row */
};

struct SCIP_ProjShiftData
{
   SCIP_Rational**       interiorpoint;      /**< stores S-interior point for root node dual problem */
   SCIP_Rational**       interiorray;        /**< stores S-interior ray for root node dual problem */
   SCIP_Rational**       violation;          /**< needed on every iteration, so only construct once and possibly resize */
   SCIP_Rational**       correction;         /**< needed on every iteration, so only construct once and possibly resize */
   int*                  includedrows;       /**< 1 if constraints dual variable is included in original S-interior point/ray */
   int*                  projshiftbasis;     /**< mapping for basis used in factorization */
#ifdef SCIP_WITH_GMP
   qsnum_factor_work*    rectfactor;         /**< stores factorized matrix for project-and-shift */
#endif
   SCIP_Rational*        commonslack;        /**< slack by which S-interior point/ray satisfies inequalities */
   int                   projshiftbasisdim;  /**< length of projshiftbasis */
   int                   nextendedrows;      /**< dimension of S-interior point/ray = 2*(ncols+nrows) */
   int                   violationsize;      /**< size of violation array */
   unsigned int          projshiftdatacon:1; /**< was project-and-shift data structure constructed? */
   unsigned int          projshiftdatafail:1;/**< did the construction of the project-and-shift root node data fail? */
   unsigned int          projshifthaspoint:1;/**< has an S-interior point successfully been constructed? */
   unsigned int          projshifthasray:1;  /**< has an S-interior ray successfully been constructed? */
   unsigned int          projshiftobjweight:1;/**< weight of the original objective function in lp to compute interior point */
   unsigned int          scaleobj:1;         /**< should the objective be scaled to be integral if possible? */
   unsigned int          projshiftuseintpoint:1;/**< should correction shift use an interior pt? (otherwise use interior ray of recession cone) */
   unsigned int          psdualcolselection:2;/**< strategy to select which dual columns to use for lp to compute interior point
                                              *   ('n'o sel, 'a'ctive rows of exact primal LP, 'A'ctive rows of inexact primal LP) */
   unsigned int          psintpointselection:3;/**< method to select interior point ('a'rbitrary interior point, 'o'ptimized interior point
                                              *   'A'rbitrary interior point in dual form, 't'wo stage optimized interior point) */
};


/** current LP data */
struct SCIP_LpExact
{
   SCIP_LP*              fplp;               /**< pointer to the floating point lp */
   SCIP_PROJSHIFTDATA*   projshiftdata;      /**< data stored for usage in project+shift, NULL if ps not used */
   SCIP_Rational*        lpobjval;           /**< objective value of LP without loose variables, or SCIP_INVALID */
   SCIP_Rational*        looseobjval;        /**< current solution value of all loose variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Rational*        glbpseudoobjval;    /**< global pseudo solution value with all variables set to their best global bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Rational*        pseudoobjval;       /**< current pseudo solution value with all variables set to their best bounds,
                                              *   ignoring variables, with infinite best bound */
   SCIP_Real             cutoffbound;        /**< upper objective limit of LP (copy of primal->cutoffbound) */
   SCIP_Real             lpiobjlim;          /**< current objective limit in LPI */
   SCIP_LPIEXACT*        lpiexact;           /**< exact LP solver interface */
   SCIP_COLEXACT**       lpicols;            /**< array with columns currently stored in the LP solver */
   SCIP_ROWEXACT**       lpirows;            /**< array with rows currently stored in the LP solver */
   SCIP_COLEXACT**       chgcols;            /**< array of changed columns not yet applied to the LP solver */
   SCIP_ROWEXACT**       chgrows;            /**< array of changed rows not yet applied to the LP solver */
   SCIP_COLEXACT**       cols;               /**< array with current LP columns in correct order */
   SCIP_ROWEXACT**       rows;               /**< array with current LP rows in correct order */
   int                   lpicolssize;        /**< available slots in lpicols vector */
   int                   nlpicols;           /**< number of columns in the LP solver */
   int                   lpifirstchgcol;     /**< first column of the LP which differs from the column in the LP solver */
   int                   lpirowssize;        /**< available slots in lpirows vector */
   int                   nlpirows;           /**< number of rows in the LP solver */
   int                   lpifirstchgrow;     /**< first row of the LP which differs from the row in the LP solver */
   int                   chgcolssize;        /**< available slots in chgcols vector */
   int                   nchgcols;           /**< current number of chgcols (number of used slots in chgcols vector) */
   int                   chgrowssize;        /**< available slots in chgrows vector */
   int                   nchgrows;           /**< current number of chgrows (number of used slots in chgrows vector) */
   int                   colssize;           /**< available slots in cols vector */
   int                   ncols;              /**< current number of LP columns (number of used slots in cols vector) */
   int                   nremovablecols;     /**< number of removable columns in the LP */
   int                   firstnewcol;        /**< first column added at the current node */
   int                   rowssize;           /**< available slots in rows vector */
   int                   nrows;              /**< current number of LP rows (number of used slots in rows vector) */
   int                   nremovablerows;     /**< number of removable rows in the LP */
   int                   firstnewrow;        /**< first row added at the current node */
   int                   looseobjvalinf;     /**< number of loose variables with infinite best bound in current solution */
   int                   nloosevars;         /**< number of loose variables in LP */
   int                   glbpseudoobjvalinf; /**< number of variables with infinite best bound in global pseudo solution */
   int                   pseudoobjvalinf;    /**< number of variables with infinite best bound in current pseudo solution */
   int                   lpiitlim;           /**< current iteration limit setting in LPI */
   int                   lpitiming;          /**< current timing type in LPI */
   int                   lpirandomseed;      /**< current initial random seed in LPI */
   int                   lpiscaling;         /**< current SCALING setting in LPI */
   int                   lpirefactorinterval;/**< current refactorization interval */
   int                   interleavedbfreq;   /**< frequency at which dual bounding strategy is interleaved (-1: never, 0: if prommising, x: xth node) */
   SCIP_PRICING          lpipricing;         /**< current pricing setting in LPI */
   SCIP_LPSOLSTAT        lpsolstat;          /**< solution status of last LP solution */
   SCIP_LPALGO           lastlpalgo;         /**< algorithm used for last LP solve */
   SCIP_Bool             lpisolutionpolishing;/**< LP solution polishing method (0: disabled, 1: enabled) */
   SCIP_Bool             flushdeletedcols;   /**< have LPI-columns been deleted in the last lpFlush() call? */
   SCIP_Bool             flushaddedcols;     /**< have LPI-columns been added in the last lpFlush() call? */
   SCIP_Bool             flushdeletedrows;   /**< have LPI-rows been deleted in the last lpFlush() call? */
   SCIP_Bool             flushaddedrows;     /**< have LPI-rows been added in the last lpFlush() call? */
   SCIP_Bool             updateintegrality;  /**< does integrality information need to be updated? */
   SCIP_Bool             flushed;            /**< are all cached changes applied to the LP solver? */
   SCIP_Bool             solved;             /**< is current LP solved? */
   SCIP_Bool             primalfeasible;     /**< is current LP solution (rather LPI state) primal feasible? */
   SCIP_Bool             primalchecked;      /**< was current LP solution checked for primal feasibility?? */
   SCIP_Bool             dualfeasible;       /**< is current LP solution (rather LPI state) dual feasible? */
   SCIP_Bool             dualchecked;        /**< was current LP solution checked for primal feasibility?? */
   SCIP_Bool             solisbasic;         /**< is current LP solution a basic solution? */
   SCIP_Bool             resolvelperror;     /**< an error occured during resolving the LP after diving or probing */
   SCIP_Bool             lpihasscaling;      /**< does the LPI support the SCALING parameter? */
   SCIP_Bool             lpihaspresolving;   /**< does the LPI support the PRESOLVING parameter? */
   SCIP_Bool             projshiftpossible;  /**< can a safe bound be computed with project-and-shift? */
   SCIP_Bool             boundshiftviable;   /**< is bound-shift viable? set to FALSE if success rate too low */
   SCIP_Bool             forceexactsolve;    /**< should the next safe bounding step be forced to solve the lp exactly? */
};

#ifdef __cplusplus
}
#endif

#endif
