/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_matrix.h
 * @ingroup INTERNALAPI
 * @brief  data structure for MIP matrix
 * @author Dieter Weninger
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_MATRIX_H__
#define __SCIP_STRUCT_MATRIX_H__

#include "scip/def.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

struct SCIP_MatrixValsExact
{
   SCIP_Rational**       lbexact;            /**< exact lower bound per variable */
   SCIP_Rational**       ubexact;            /**< exact upper bound per variable */
   SCIP_Rational**       lhsexact;            /**< exact lhs per constraint */
   SCIP_Rational**       rhsexact;            /**< exact rhs per constraint */
   SCIP_Rational**       colmatvalexact;     /**< exact coefficients in column major format */
   SCIP_Rational**       rowmatvalexact;     /**< exact coefficients in column major format */
   SCIP_Rational**       minacivityexact;    /**< exact min activity per row */
   SCIP_Rational**       maxacivityexact;    /**< exact max activity per row */
   int                   buffersize;         /**< necessary because rational buffer arrays need to be freed with a size */
   int                   buffersizenconss;   /**< necessary because rational buffer arrays need to be freed with a size */
};

/** constraint matrix data structure in column and row major format */
struct SCIP_Matrix
{
   SCIP_Real*            colmatval;          /**< coefficients in column major format */
   int*                  colmatind;          /**< row indexes in column major format */
   int*                  colmatbeg;          /**< column storage offset */
   int*                  colmatcnt;          /**< number of row entries per column */
   int                   ncols;              /**< complete number of columns */
   SCIP_Real*            lb;                 /**< lower bound per variable */
   SCIP_Real*            ub;                 /**< upper bound per variable */
   int*                  nuplocks;           /**< number of up locks per variable */
   int*                  ndownlocks;         /**< number of down locks per variable */

   SCIP_VAR**            vars;               /**< variables pointer */

   SCIP_Real*            rowmatval;          /**< coefficients in row major format */
   int*                  rowmatind;          /**< column indexed in row major format */
   int*                  rowmatbeg;          /**< row storage offset */
   int*                  rowmatcnt;          /**< number of column entries per row */

   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */

   SCIP_CONS**           cons;               /**< constraints pointer */

   SCIP_Bool*            isrhsinfinite;      /**< is right hand side infinity */
   int                   nnonzs;             /**< sparsity counter */
   SCIP_Real*            minactivity;        /**< min activity per row */
   SCIP_Real*            maxactivity;        /**< max activity per row */
   int*                  minactivityneginf;  /**< min activity negative infinity counter */
   int*                  minactivityposinf;  /**< min activity positive infinity counter */
   int*                  maxactivityneginf;  /**< max activity negative infinity counter */
   int*                  maxactivityposinf;  /**< max activity positive infinity counter */
   SCIP_MATRIXVALSEXACT* matrixvalsexact;    /**< exact matrix data, or NULL if matrix is not exact */
};

#ifdef __cplusplus
}
#endif

#endif
