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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_certificate.h
 * @brief  data structures for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CERTIFICATE_H__
#define __SCIP_STRUCT_CERTIFICATE_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_certificate.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for hashing bounds of variables in a certificate file */
struct SCIP_CertificateBound
{
   int                   varindex;           /**< index of the variable */
   SCIP_Rational*        boundval;           /**< value of the bound */
   SCIP_BOUNDTYPE        boundtype;          /**< is it the upper bound? */
   SCIP_Bool             isbound;            /**< is the last printed index a bound? if it is not, the other information is not useful */
};

struct SCIP_AggregationInfo
{
   SCIP_AGGRROW*         aggrrow;            /**< aggregation row to be saved */
   SCIP_ROW**            aggrrows;           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights;            /**< array of weights */
   int                   naggrrows;          /**< length of the arrays */
   SCIP_Longint          fileindex;          /**< index of the aggregated row in the certificate file */
   SCIP_Longint          arpos;              /**< position in the aggrinfo array, so we can access it from the hashmap */
};

struct SCIP_MirInfo
{
   SCIP_Real*            splitcoefs;         /**< coefficients in the split, saved in the complemented variable space */
   SCIP_Real*            contcoefs;          /**< coefficients of continous part, saved in the complemented variable space */
   int*                  splitvarinds;       /**< indices of variables in split */
   int*                  contvarinds;        /**< indices of variables in split */
   SCIP_Bool*            splitupperused;     /**< TRUE if ub was used to complemented variable, FALSE if lb was used */
   SCIP_Bool*            contupperused;      /**< TRUE if ub was used to complemented variable, FALSE if lb was used */
   int                   nsplitvars;         /**< number of variables in the split */
   int                   ncontvars;          /**< number of variables in the split */
   SCIP_Real             rhs;                /**< rhs of the split disjunction */
   SCIP_Longint          arpos;              /**< position in the mirinfo array, so we can access it from the hashmap */
   SCIP_Bool             global;             /**< is the cut valid globally? (relevant for which bounds to use in rounding)*/
};

struct SCIP_Certnodedata
{
   SCIP_Longint          assumptionindex_self;/**< Line Index where assumption is printed */
   SCIP_Longint          assumptionindex_left;/**< Line Index of left branch assumption */
   SCIP_Longint          derindex_left;      /**< Line Index of derivation assuming assumption left */
   SCIP_Rational*        derbound_left;      /**< Bound of left derivation */
   SCIP_Longint          assumptionindex_right;/**< Line Index of right branch assumption */
   SCIP_Longint          derindex_right;     /**< Line Index of derivation assuming assumption right */
   SCIP_Rational*        derbound_right;     /**< Bound of right derivation */
   SCIP_Longint          derindex_inherit;   /**< Line index of bound inherited from parent */
   SCIP_Rational*        derbound_inherit;   /**< inherited bound */
   unsigned int          leftfilled:1;       /**< Is the left node filled ? */
   unsigned int          leftinfeas:1;       /**< Is the left node infeasible ? */
   unsigned int          rightfilled:1;      /**< Is the node right filled ? */
   unsigned int          rightinfeas:1;      /**< Is the node right infeasible ? */
   unsigned int          inheritedbound:1;   /**< did the node inherit its bound from its parent node? */
};

/** certificate data structure */
struct SCIP_Certificate
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler to use */
   SCIP_HASHMAP*         nodedatahash;       /**< Hashmap storing pointer to data of each node */
   SCIP_HASHMAP*         aggrinfohash;       /**< Hashmap storing aggregation information of rows */
   SCIP_HASHMAP*         mirinfohash;      /**< Hashmap storing split disjunctions */
   SCIP_AGGREGATIONINFO** aggrinfo;          /**< array to store the aggregation info to avoid memory leaks */
   SCIP_MIRINFO**        mirinfo;          /**< array to store the split info to avoid memory leaks */
   SCIP_Longint          aggrinfosize;       /**< size of aggrinfo array */
   SCIP_Longint          naggrinfos;         /**< number of elements in aggrinfo array */
   SCIP_Longint          mirinfosize;      /**< size of mirinfo array */
   SCIP_Longint          nmirinfos;        /**< number of elements in mirinfo array */
   SCIP_CERTIFICATEBOUND* lastinfo;          /**< information on last printed certificate index */
   BMS_BLKMEM*           blkmem;             /**< SCIP block memory */
   SCIP_Longint          indexcounter;       /**< counter for line indices in file */
   SCIP_Longint          indexcounter_ori;   /**< counter for line indices in origial problem vipr file */
   SCIP_Longint          conscounter;        /**< counter for line indices in constraint section */
   SCIP_FILE*            origfile;           /**< file to store original problem definition */
   SCIP_FILE*            transfile;          /**< file to store transformed problem (after presolving) */
   SCIP_FILE*            derivationfile;     /**< file to store derivations temporarily */
   char*                 derivationfilename; /**< name of the derivation file */
   char*                 origfilename;       /**< name of the original problem file */
   SCIP_Real             filesize;           /**< size of derivation file in MB */
   SCIP_HASHMAP*         rowdatahash;        /**< Hashmap storing mapping between rows and file index */
   SCIP_Rational*        rootbound;          /**< the bound for the root node */
   SCIP_Rational*        finalbound;         /**< the final dual bound value */
   SCIP_Longint          derindex_root;      /**< index of root bound in certificate */
   SCIP_Bool             rootinfeas;         /**< is the root node infeasible */
   SCIP_Bool             objintegral;        /**< is the objective always integral? copy this so we don't need the prob everywhere */
   SCIP_Bool             workingmirinfo;     /**< true if mirinfo is under construction and not sparsely stored, false otherwise */
   SCIP_Bool             workingaggrinfo;     /**< true if aggrinfo is under construction (last entry not in hashmap), false otherwise */
   SCIP_Rational**       vals;               /**< we maintain an array for solvals so we don't have to reallocate at every bounding call */
   int                   valssize;           /**< the size of the vals array */
};

#ifdef __cplusplus
}
#endif

#endif
