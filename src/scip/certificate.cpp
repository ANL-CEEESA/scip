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

/**@file   certificate.c
 * @brief  methods for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <map>

#include "blockmemshell/memory.h"
#include "scip/cuts.h"
#include "scip/cons_exactlp.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/pub_lpexact.h"
#include "lpiexact/lpiexact.h"
#include "scip/pub_misc.h"
#include "scip/prob.h"
#include "scip/certificate.h"
#include "scip/struct_certificate.h"
#include "scip/sol.h"
#include "scip/struct_scip.h"
#include "scip/var.h"

#define SCIP_HASHSIZE_CERTIFICATE    500 /**< size of hash map for certificate -> nodesdata mapping used for certificate output */

/** updates file size */
static
void updateFilesize(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Real             nchars              /**< number of characters printed */
   )
{
   certificate->filesize += nchars/1048576.0;
}

/** checks whether node is a left node or not */
static
SCIP_Bool certificateIsLeftNode(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node                /**< node from branch and bound tree */
   )
{
   SCIP_CERTNODEDATA* nodedata;
   SCIP_CERTNODEDATA* nodedataparent;

   assert(node != NULL);
   assert(SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE);

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( SCIPnodeGetParent(node) == NULL )
      return FALSE;

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

   assert(nodedata->assumptionindex_self != -1);
   if( nodedataparent->assumptionindex_left == nodedata->assumptionindex_self )
      return TRUE;
   else
      return FALSE;
}

/** print variable bound assumption into certificate,
 *  return index of this bound in the certificate file
 */
static
SCIP_Longint printBoundAssumption(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_VAR*             var,                /**< variable to print assumption for */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_BOUNDTYPE        boundtype           /**< is it the upper bound? */
   )
{
   /* check whether output should be created */
   if ( certificate->origfile == NULL )
      return SCIP_OKAY;

#ifndef NDEBUG
   certificate->lastinfo->isbound = TRUE;
   certificate->lastinfo->boundtype = boundtype;
   certificate->lastinfo->varindex = SCIPvarGetCertificateIndex(var);
   RatSet(certificate->lastinfo->boundval, boundval);
#endif

   /** @todo: it could be better to seperate the printing from insertion of variable bound */
   SCIPcertificatePrintProofMessage(certificate, "A%lld %c ", certificate->indexcounter, (boundtype == SCIP_BOUNDTYPE_LOWER) ? 'G' : 'L');

   SCIPcertificatePrintProofRational(certificate, boundval, 10);
   SCIPcertificatePrintProofMessage(certificate, " 1 %d 1 { asm } -1\n", SCIPvarGetCertificateIndex(var));
   certificate->indexcounter++;

   return certificate->indexcounter - 1;
}

/** free nodedata of corresponding node */
static
SCIP_RETCODE certificateFreeNodeData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node                /**< focus node */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   assert(node != NULL);
   assert(certificate != NULL);

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);
   RatFreeBlock(certificate->blkmem, &nodedata->derbound_left);
   RatFreeBlock(certificate->blkmem, &nodedata->derbound_right);
   RatFreeBlock(certificate->blkmem, &nodedata->derbound_inherit);
   BMSfreeBlockMemory(certificate->blkmem, &nodedata);
   SCIP_CALL( SCIPhashmapRemove(certificate->nodedatahash, node) );

   return SCIP_OKAY;
}

/** print the best solution found */
static
SCIP_RETCODE SCIPcertificatePrintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_SOL*             sol                 /**< solution to be printed */
   )
{
   SCIP_VAR** vars;
   SCIP_Rational** vals;
   int nvars;
   int nnonz;
   int i;

   /* check if certificate output should be created */
   if( certificate->origfile == NULL )
      return SCIP_OKAY;

   assert(scip != NULL);

   if( sol == NULL )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "SOL 0\n");
      return SCIP_OKAY;
   }
   else if( !SCIPsolIsExact(sol) )
   {
      SCIP_CALL( SCIPmakeSolExact(scip, sol) );
   }

   /* get variables and number of the transformed problem */
   if( isorigfile )
   {
      SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   }

   SCIP_CALL( RatCreateBufferArray(SCIPbuffer(scip), &vals, nvars) );

   /* get number of non-zero coefficient in the solution */
   nnonz = 0;
   for( i = 0; i < nvars; i++)
   {
      SCIPsolGetValExact(vals[i], sol, scip->set, scip->stat, vars[i]);
      if( !RatIsZero(vals[i]) )
         nnonz++;
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "SOL 1\nbest %d", nnonz);

   for( i = 0; i < nvars; i++ )
   {
      if( !RatIsZero(vals[i]) )
      {
         /* print the solution into certificate */
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d ", SCIPvarGetCertificateIndex(vars[i]));
         SCIPcertificatePrintProblemRational(certificate, isorigfile, vals[i], 10);
      }
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   RatFreeBufferArray(SCIPbuffer(scip), &vals, nvars);

   return SCIP_OKAY;
}

/** set the node to have its own bound proof */
SCIP_RETCODE SCIPcertificateSetInheritanceData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_Rational*        newbound            /**< the inherited bound */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   /* do nothing if the newbound is worse than the inherited bound */
   if( RatIsLT(newbound, nodedata->derbound_inherit) )
      return SCIP_OKAY;

   nodedata->inheritedbound = FALSE;
   nodedata->derindex_inherit = fileindex;
   RatSet(nodedata->derbound_inherit, newbound);

   return SCIP_OKAY;
}

/** creates certificate data structure */
SCIP_RETCODE SCIPcertificateCreate(
   SCIP_CERTIFICATE**    certificate,        /**< pointer to store the certificate information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_ALLOC( BMSallocMemory(certificate) );

   (*certificate)->messagehdlr = messagehdlr;
   (*certificate)->lastinfo = NULL;
   (*certificate)->blkmem = NULL;
   (*certificate)->indexcounter = 0;
   (*certificate)->indexcounter_ori = 0;
   (*certificate)->conscounter = 0;
   (*certificate)->origfile = NULL;
   (*certificate)->transfile = NULL;
   (*certificate)->origfilename = NULL;
   (*certificate)->derivationfile = NULL;
   (*certificate)->derivationfilename = NULL;
   (*certificate)->filesize = 0.0;
   (*certificate)->rowdatahash = NULL;
   (*certificate)->naggrinfos = 0;
   (*certificate)->nmirinfos = 0;
   (*certificate)->aggrinfosize = 0;
   (*certificate)->mirinfosize = 0;
   (*certificate)->nodedatahash = NULL;
   (*certificate)->rootbound = NULL;
   (*certificate)->finalbound = NULL;
   (*certificate)->derindex_root = -1;
   (*certificate)->rootinfeas = FALSE;
   (*certificate)->objintegral = FALSE;
   (*certificate)->workingmirinfo = FALSE;
   (*certificate)->workingaggrinfo = FALSE;
   (*certificate)->vals = NULL;
   (*certificate)->valssize = 0;
   (*certificate)->aggrinfo = NULL;
   (*certificate)->mirinfo = NULL;

   return SCIP_OKAY;
}

/** frees certificate data structure */
void SCIPcertificateFree(
   SCIP_CERTIFICATE**    certificate         /**< pointer to store the certificate information */
   )
{
   assert(certificate != NULL);
   assert(*certificate != NULL);
   assert((*certificate)->origfile == NULL);
   assert((*certificate)->transfile == NULL);
   assert((*certificate)->derivationfile == NULL);

   BMSfreeMemory(certificate);
}

/** initializes certificate information and creates files for certificate output */
SCIP_RETCODE SCIPcertificateInit(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   int filenamelen;
   int bufferlen;
   int nvars;
   int nintvars;
   int nbinvars;
   int nboundconss;
   int ncertcons;
   int j;
   char* name = NULL;
   char* compression = NULL;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_Rational* lb;
   SCIP_Rational* ub;

   assert(certificate != NULL);
   assert(set != NULL);
   assert(set->certificate_filename != NULL);
   assert(certificate->derivationfile == NULL);
   assert(certificate->nodedatahash == NULL);
   assert(certificate->rowdatahash == NULL);

   if( !(set->exact_enabled) || (set->certificate_filename[0] == '-' && set->certificate_filename[1] == '\0') )
      return SCIP_OKAY;

   filenamelen = strlen(set->certificate_filename);
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, filenamelen + 1) );
   BMScopyMemoryArray(name, set->certificate_filename, filenamelen);
   name[filenamelen] = '\0';

   /** @todo exip: this currently strips the .vipr from the filename if there is no compression...
    * the whole compression code currently makes no sense, since the certificate has to be unpacked to be read by vipr
    * again anyway, so it is just disabled.
    */
   //SCIPsplitFilename(name, NULL, NULL, NULL, &compression);

   SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "storing certificate information in file <%s>\n", set->certificate_filename);

   if( NULL != compression && 0 == strncmp(compression, "gz", 2) )
      certificate->transfile = SCIPfopen(set->certificate_filename, "wb");
   else
      certificate->transfile = SCIPfopen(set->certificate_filename, "wT");

   bufferlen = strlen(name);
   SCIP_ALLOC( BMSallocMemoryArray(&certificate->derivationfilename, filenamelen+5) );
   SCIP_ALLOC( BMSallocMemoryArray(&certificate->origfilename, filenamelen+5) );
   BMScopyMemoryArray(certificate->derivationfilename, name, bufferlen);
   BMScopyMemoryArray(certificate->origfilename, name, bufferlen);
   certificate->derivationfilename[bufferlen] = '_';
   certificate->derivationfilename[bufferlen+1] = 'd';
   certificate->derivationfilename[bufferlen+2] = 'e';
   certificate->derivationfilename[bufferlen+3] = 'r';
   certificate->origfilename[bufferlen] = '_';
   certificate->origfilename[bufferlen+1] = 'o';
   certificate->origfilename[bufferlen+2] = 'r';
   certificate->origfilename[bufferlen+3] = 'i';
   if( NULL != compression && 0 == strncmp(compression, "gz", 2) )
   {
      certificate->derivationfilename[bufferlen+4] = '.';
      certificate->derivationfilename[bufferlen+5] = 'g';
      certificate->derivationfilename[bufferlen+6] = 'z';
      certificate->derivationfilename[bufferlen+7] = '\0';

      certificate->origfilename[bufferlen+4] = '.';
      certificate->origfilename[bufferlen+5] = 'g';
      certificate->origfilename[bufferlen+6] = 'z';
      certificate->origfilename[bufferlen+7] = '\0';

      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wb");
      certificate->origfile = SCIPfopen(certificate->origfilename, "wb");
   }
   else
   {
      certificate->derivationfilename[bufferlen+4] = '\0';
      certificate->origfilename[bufferlen+4] = '\0';
      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wT");
      certificate->origfile = SCIPfopen(certificate->origfilename, "wT");
   }

   if( certificate->origfile == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and derivation file\n", set->certificate_filename);
      SCIPprintSysError(set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   /* initialisation of hashmaps and hashtables */
   SCIP_CALL( SCIPhashmapCreate(&certificate->nodedatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->rowdatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->aggrinfohash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->mirinfohash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(certificate->aggrinfo), SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(certificate->mirinfo), SCIP_HASHSIZE_CERTIFICATE) );
   certificate->aggrinfosize = SCIP_HASHSIZE_CERTIFICATE;
   certificate->mirinfosize = SCIP_HASHSIZE_CERTIFICATE;

   certificate->blkmem = blkmem;
   SCIPsetFreeBufferArray(set, &name);

   SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nboundconss = 0;
   for ( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !RatIsAbsInfinity(lb) )
         nboundconss++;
      if( !RatIsAbsInfinity(ub) )
         nboundconss++;
   }

   /* print the Version Header into certificate */
   SCIPcertificatePrintVersionHeader(certificate, TRUE);

   /* print the Variable Header into certificate */
   SCIPcertificatePrintVarHeader(certificate, TRUE, nvars);
   for( j = 0; j < nvars; j++ )
   {
      const char* varname;

      varname = SCIPvarGetName(vars[j]);
      SCIPvarSetCertificateIndex(vars[j], j);
      SCIPvarSetCertificateIndex(SCIPvarGetTransVar(vars[j]), j);
      if( strstr(varname, " ") != NULL || strstr(varname, "\t") != NULL || strstr(varname, "\n") != NULL
         || strstr(varname, "\v") != NULL || strstr(varname, "\f") != NULL || strstr(varname, "\r") != NULL )
      {
         SCIPerrorMessage("Variable name <%s> cannot be printed to certificate file because it contains whitespace.\n",
            varname);
         return SCIP_ERROR;
      }

      SCIPcertificatePrintProblemMessage(certificate, TRUE, "%s\n", varname);
   }

   /* print the Integer Variable Header into certificate */
   SCIPcertificatePrintIntHeader(certificate, TRUE, nintvars + nbinvars);
   for( j = 0; j < nvars; j++ )
   {
      if( SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[j]) == SCIP_VARTYPE_INTEGER )
      {
         SCIPcertificatePrintProblemMessage(certificate, TRUE, "%d \n", j);
      }
   }

   {
      SCIP_Rational** objcoefs;
      SCIP_CALL( RatCreateBufferArray(SCIPbuffer(scip), &objcoefs, nvars) );

      for( j = 0; j < nvars; j++)
         RatSet(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIPcertificateSetAndPrintObjective(certificate, TRUE, blkmem, objcoefs, nvars);

      RatFreeBufferArray(SCIPbuffer(scip), &objcoefs, nvars);
   }

   conss = SCIPgetOrigConss(scip);
   ncertcons = 0;
   for( j = 0; j < SCIPgetNOrigConss(scip); j++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSHDLR* conshdlr;

      cons = conss[j];
      conshdlr = SCIPconsGetHdlr(cons);

      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear-exact") == 0 )
      {
         lb = SCIPgetLhsExactLinear(scip, cons);
         ub = SCIPgetRhsExactLinear(scip, cons);

         if( !RatIsEqual(lb, ub) && !RatIsAbsInfinity(lb) && !RatIsAbsInfinity(ub) )
            ncertcons += 2;
         else
            ncertcons += 1;
      }
      else
      {
         SCIPerrorMessage("Cannot print certificate for non-exact constraints \n");
         SCIPABORT();
         return SCIP_ERROR;
      }
   }

   SCIPcertificatePrintConsHeader(certificate, TRUE, ncertcons, nboundconss);

   for( j = 0; j < nvars; j++ )
   {
      if( !RatIsAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, TRUE, NULL, vars[j], SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
      }
      if( !RatIsAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, TRUE, NULL, vars[j], SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
      }
   }

   SCIP_CALL( RatCreateBlock(blkmem, &certificate->rootbound) );
   SCIP_CALL( RatCreateBlock(blkmem, &certificate->finalbound) );
   certificate->valssize = SCIPgetNVars(scip) + SCIPgetNConss(scip);
   SCIP_CALL( RatCreateBlockArray(SCIPblkmem(scip), &(certificate->vals), certificate->valssize) );

   return SCIP_OKAY;
}

/** initializes certificate information and creates files for certificate output */
SCIP_RETCODE SCIPcertificateInitTransFile(
   SCIP*                 scip                /**< scip data structure */
   )
{
   int nvars;
   int nintvars;
   int nbinvars;
   int nboundconss;
   int ncertcons;
   int j;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_Rational* lb;
   SCIP_Rational* ub;
   SCIP_CERTIFICATE* certificate;
   BMS_BLKMEM* blkmem;

   assert(scip != NULL);
   assert(scip->set->certificate_filename != NULL);

   certificate = SCIPgetCertificate(scip);
   blkmem = SCIPblkmem(scip);

   assert(certificate != NULL);

   if( !(scip->set->exact_enabled) || (scip->set->certificate_filename[0] == '-' && scip->set->certificate_filename[1] == '\0') )
      return SCIP_OKAY;

   if( certificate->origfile == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and derivation file\n", scip->set->certificate_filename);
      SCIPprintSysError(scip->set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nboundconss = 0;
   for ( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !RatIsAbsInfinity(lb) )
         nboundconss++;
      if( !RatIsAbsInfinity(ub) )
         nboundconss++;
   }

   /* print the Version Header into certificate */
   SCIPcertificatePrintVersionHeader(certificate, FALSE);

   /* print the Variable Header into certificate */
   SCIPcertificatePrintVarHeader(certificate, FALSE, nvars);
   for( j = 0; j < nvars; j++ )
   {
      const char* varname;

      varname = SCIPvarGetName(vars[j]);
      SCIPvarSetCertificateIndex(vars[j], j);
      if( strstr(varname, " ") != NULL || strstr(varname, "\t") != NULL || strstr(varname, "\n") != NULL
         || strstr(varname, "\v") != NULL || strstr(varname, "\f") != NULL || strstr(varname, "\r") != NULL )
      {
         SCIPerrorMessage("Variable name <%s> cannot be printed to certificate file because it contains whitespace.\n",
            varname);
         return SCIP_ERROR;
      }

      SCIPcertificatePrintProblemMessage(certificate, FALSE, "%s\n", varname);
   }

   /* print the Integer Variable Header into certificate */
   SCIPcertificatePrintIntHeader(certificate, FALSE, nintvars + nbinvars);
   for( j = 0; j < nvars; j++ )
   {
      if( SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[j]) == SCIP_VARTYPE_INTEGER )
      {
         SCIPcertificatePrintProblemMessage(certificate, FALSE, "%d \n", SCIPvarGetCertificateIndex(vars[j]));
      }
   }

   {
      SCIP_Rational** objcoefs;
      SCIP_CALL( RatCreateBufferArray(SCIPbuffer(scip), &objcoefs, nvars) );

      for( j = 0; j < nvars; j++)
         RatSet(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIPcertificateSetAndPrintObjective(certificate, FALSE, blkmem, objcoefs, nvars);

      RatFreeBufferArray(SCIPbuffer(scip), &objcoefs, nvars);
   }

   conss = SCIPgetConss(scip);
   ncertcons = 0;
   for( j = 0; j < SCIPgetNConss(scip); j++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSHDLR* conshdlr;

      cons = conss[j];
      conshdlr = SCIPconsGetHdlr(cons);

      SCIPdebug(SCIPprintCons(scip, conss[j], NULL));

      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear-exact") == 0 )
      {
         lb = SCIPgetLhsExactLinear(scip, cons);
         ub = SCIPgetRhsExactLinear(scip, cons);

         if( !RatIsEqual(lb, ub) && !RatIsAbsInfinity(lb) && !RatIsAbsInfinity(ub) )
         {
            SCIPdebugMessage("constraint is a ranged constraint \n");
            ncertcons += 2;
         }
         else
         {
            SCIPdebugMessage("constraint only has one side \n");
            ncertcons += 1;
         }
      }
      else
      {
         SCIPerrorMessage("Cannot print certificate for non-exact constraints \n");
         SCIPABORT();
         return SCIP_ERROR;
      }
   }

   SCIPcertificatePrintConsHeader(certificate, FALSE, ncertcons, nboundconss);

   for( j = 0; j < nvars; j++ )
   {
      if( !RatIsAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, FALSE, NULL, vars[j], SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
      }
      if( !RatIsAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, FALSE, NULL, vars[j], SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
      }
   }

   return SCIP_OKAY;
}

/** Concatenate the certificate and the _der file and delete the _der file  */
static
void concatCert(
   SCIP_CERTIFICATE*     certificate,        /**< The certificate pointer */
   const char*           certificatefilename /**< The name of the certificate file */
   )
{
   SCIP_FILE* derivationfile;
   char buffer[SCIP_MAXSTRLEN];
   size_t size;

   derivationfile = SCIPfopen(certificate->derivationfilename, "r");

   /* append the derivation file to the problem file */
   while( 0 != (size = SCIPfread(buffer, sizeof(char), SCIP_MAXSTRLEN, derivationfile)) )
      SCIPfwrite(buffer, sizeof(char), size, certificate->transfile);

   SCIPfclose(derivationfile);
   /* delete the derivation file */
   remove(certificate->derivationfilename);
}

/** closes the certificate output files */
void SCIPcertificateExit(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   int i;

   assert(certificate != NULL);
   assert(set != NULL);

   if( certificate->origfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "closing CERTIFICATE information file (wrote %.1f MB)\n", certificate->filesize);

      if( certificate->derivationfile != NULL )
      {
         /* CERT TODO: DER line with counter and append two files */
         SCIPfclose(certificate->derivationfile);
         certificate->derivationfile = NULL;
         concatCert(certificate, set->certificate_filename);
         // if the file is empty (e.g. because we detected infeasibility in presolving) we delete it
         if( certificate->indexcounter == 0 )
         {
            SCIPdebugMessage("derivation file is empty; deleting it");
            remove(set->certificate_filename);
         }
      }
      SCIPfclose(certificate->origfile);
      SCIPfclose(certificate->transfile);
      certificate->origfile = NULL;
      certificate->transfile = NULL;

      BMSfreeMemoryArray(&certificate->derivationfilename);
      BMSfreeMemoryArray(&certificate->origfilename);

      if( certificate->lastinfo != NULL )
      {
         RatFreeBlock(certificate->blkmem, &certificate->lastinfo->boundval);
         BMSfreeBlockMemory(certificate->blkmem, &certificate->lastinfo);
      }

      if( certificate->rowdatahash)
         SCIPhashmapFree(&certificate->rowdatahash);

      if( certificate->nodedatahash )
      {
         assert(SCIPhashmapIsEmpty(certificate->nodedatahash));
         SCIPhashmapFree(&certificate->nodedatahash);
      }
      if( certificate->aggrinfohash )
      {
         SCIPcertificateClearAggrinfo(scip);
      }
      if( certificate->mirinfohash )
      {
         for( i = 0; i < certificate->nmirinfos; i++ )
         {
            RatFreeBlock(certificate->blkmem, &certificate->mirinfo[i]->frac);
            RatFreeBlock(certificate->blkmem, &certificate->mirinfo[i]->rhs);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->varinds), certificate->mirinfo[i]->nsplitvars);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->splitcoefficients), certificate->mirinfo[i]->nsplitvars);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->upperused), certificate->mirinfo[i]->nsplitvars);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->localbdused), certificate->mirinfo[i]->nsplitvars);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackrows), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slacksign), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackcoefficients), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackweight), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackroundeddown), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackscale), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackusedcoef), certificate->mirinfo[i]->nslacks);
            BMSfreeBlockMemory(certificate->blkmem, &certificate->mirinfo[i]);
         }
         BMSfreeBlockMemoryArray(certificate->blkmem, &certificate->mirinfo, certificate->mirinfosize);
         SCIPhashmapRemoveAll(certificate->mirinfohash);
         SCIPhashmapFree(&certificate->mirinfohash);
      }

      RatFreeBlock(certificate->blkmem, &certificate->rootbound);
      RatFreeBlock(certificate->blkmem, &certificate->finalbound);
      RatFreeBlockArray(certificate->blkmem, &certificate->vals, certificate->valssize);
   }
}

/** returns whether the certificate output is activated? */
SCIP_Bool SCIPcertificateIsActive(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   return (set->exact_enabled && certificate != NULL && (certificate->transfile != NULL || strcmp(set->certificate_filename, "-")) );
}

/** returns whether the certificate output is activated? */
SCIP_Bool SCIPsetCertificateEnabled(
   SCIP_SET*             set                 /**< SCIP settings */
   )
{
   if( !set->exact_enabled )
      return FALSE;
   if( (set->certificate_filename[0] == '-' && set->certificate_filename[1] == '\0') )
      return FALSE;

   return TRUE;
}

/** returns current certificate file size in MB */
SCIP_Real SCIPcertificateGetFilesize(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( certificate == NULL || certificate->transfile == NULL)
      return 0.0;
   else
      return certificate->filesize;
}

/** returns current certificate index (return -1 if certificate not active) */
SCIP_Longint SCIPcertificateGetCurrentIndex(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( certificate == NULL || certificate->transfile == NULL)
      return -1;
   else
      return certificate->indexcounter;
}

#ifndef NDEBUG
/** checks if information is consistent with printed certificate line */
SCIP_Bool SCIPcertificateEnsureLastBoundInfoConsistent(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_VAR*             var,                /**< variable that gets changed */
   SCIP_BOUNDTYPE        boundtype,          /**< lb or ub changed? */
   SCIP_Real             newbound            /**< new bound */
   )
{
   SCIP_Bool consistent;

   assert(certificate != NULL);

   consistent = TRUE;
   consistent = consistent && certificate->lastinfo->isbound;
   consistent = consistent && certificate->lastinfo->varindex == SCIPvarGetCertificateIndex(var);
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      consistent = consistent && certificate->lastinfo->boundtype == SCIP_BOUNDTYPE_LOWER;
      consistent = consistent && RatRoundReal(certificate->lastinfo->boundval, SCIP_R_ROUND_DOWNWARDS) == newbound;
   }
   else
   {
      consistent = consistent && certificate->lastinfo->boundtype == SCIP_BOUNDTYPE_UPPER;
      consistent = consistent && RatRoundReal(certificate->lastinfo->boundval, SCIP_R_ROUND_UPWARDS) == newbound;
   }

   return consistent;

}
#endif

/** sets the objective function used when printing dual bounds */
SCIP_RETCODE SCIPcertificateSetAndPrintObjective(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Rational**       coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   )
{
   char obj[SCIP_MAXSTRLEN - 2];
   int nnonz;
   int i;

   assert(coefs != NULL);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   /* create a hash table for the variable bounds (we work on the tranformed problem) */
   if( isorigfile )
   {
      /* create working memory for bound struct */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &certificate->lastinfo) );
      SCIP_CALL( RatCreateBlock(blkmem, &certificate->lastinfo->boundval) );
   }

   nnonz = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( !RatIsZero(coefs[i]) )
         nnonz++;
   }

   SCIPsnprintf(obj, SCIP_MAXSTRLEN - 2, "OBJ min\n %d ", nnonz);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s ", obj);

   for( i = 0; i < nvars; i++ )
   {
      if( !RatIsZero(coefs[i]) )
      {
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s%d ", (i > 0 ? " " : ""), i );
         SCIPcertificatePrintProblemRational(certificate, isorigfile, coefs[i], 10);
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   return SCIP_OKAY;
}

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificatePrintResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   SCIP_Rational* primalbound;
   SCIP_Rational* dualbound;
   SCIP_SOL* bestsol;

   assert(scip != NULL);

   if( certificate->transfile == NULL)
      return SCIP_OKAY;

   SCIP_CALL( RatCreateBuffer(set->buffer, &primalbound) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &dualbound) );

   /* check status first: OPTIMAL / INFEAS / NOT OPTIMAL */
   if( SCIPisInRestart(scip) )
      return SCIP_ERROR;

   if( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL )
   {
      bestsol = SCIPgetBestSol(scip);

      if( !SCIPisExactSol(scip, bestsol) )
      {
         SCIP_CALL( SCIPmakeSolExact(scip, bestsol) );
      }

      if( isorigfile )
      {
         SCIP_CALL( SCIPretransformSolExact(scip, bestsol) );
      }

      if( isorigfile )
         SCIPgetPrimalboundExact(scip, primalbound);
      else
      {
         SCIPgetUpperboundExact(scip, primalbound);
      }

      assert(!RatIsAbsInfinity(primalbound));

      /* for the orig file we only print the primal bound, since the derivation happens in the transformed problem */
      if( isorigfile )
      {
         RatSetString(dualbound, "-inf");
         /* print RTP range (same when optimal solution found) */
         SCIPcertificatePrintRtpRange(certificate, isorigfile, dualbound, primalbound);
      }
      else
         SCIPcertificatePrintRtpRange(certificate, isorigfile, primalbound, primalbound);

      /* print optimal solution into certificate */
      SCIP_CALL( SCIPcertificatePrintSol(scip, isorigfile, certificate, bestsol) );
   }
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPcertificatePrintRtpInfeas(certificate, isorigfile);
      SCIP_CALL(SCIPcertificatePrintSol(scip, isorigfile, certificate, NULL) );
   }
   else
   {
      /* two cases to distinguish: a primal bound has been found or not */
      if( SCIPisPrimalboundSol(scip) )
      {
         bestsol = SCIPgetBestSol(scip);
         if( isorigfile )
         {
            SCIP_CALL( SCIPretransformSolExact(scip, bestsol) );
         }

         if( isorigfile )
            SCIPgetPrimalboundExact(scip, primalbound);
         else
            SCIPgetUpperboundExact(scip, primalbound);
      }
      else
      {
         bestsol = NULL;
         RatSetString(primalbound, "inf");
      }

      if( isorigfile )
         RatSetString(dualbound, "-inf");
      else
         RatSet(dualbound, certificate->finalbound);

      SCIPcertificatePrintRtpRange(certificate, isorigfile, dualbound, primalbound);
      SCIP_CALL( SCIPcertificatePrintSol(scip, isorigfile, certificate, bestsol) );
   }

   SCIPcertificatePrintDerHeader(certificate, isorigfile);

   RatFreeBuffer(set->buffer, &dualbound);
   RatFreeBuffer(set->buffer, &primalbound);

   return SCIP_OKAY;
}

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificateSaveFinalbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   assert(scip != NULL);

   if( certificate->transfile == NULL)
      return SCIP_OKAY;

   SCIPgetLowerboundExact(scip, certificate->finalbound);

   return SCIP_OKAY;
}

/** prints a string to the problem section of the certificate file */
void SCIPcertificatePrintProblemMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   va_start(ap, formatstr);
   vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   if( isorigfile )
      SCIPfprintf(certificate->origfile, "%s", buffer);
   else
      SCIPfprintf(certificate->transfile, "%s", buffer);

   va_end(ap);
   updateFilesize(certificate, strlen(buffer));
}

/** prints a string to the proof section of the certificate file */
void SCIPcertificatePrintProofMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
      return;
   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   SCIPfprintf(certificate->derivationfile, "%s", buffer); // todo: is this correct?
   va_end(ap);
   updateFilesize(certificate, strlen(buffer));
}


/** prints a rational number to the problem section of the certificate file */
void SCIPcertificatePrintProblemRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   SCIP_Longint len = RatStrlen(val) + 1;
   char* formatstr;

   assert(len <= INT_MAX);

   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;

   BMSallocMemoryArray(&formatstr, len);
   RatToString(val, formatstr, len);
   if( isorigfile )
      SCIPfputs(formatstr, certificate->origfile);
   else
      SCIPfputs(formatstr, certificate->transfile);

   BMSfreeMemoryArray(&formatstr);
}


/** prints a rational number to the proof section of the certificate file */
void SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   SCIP_Longint len = RatStrlen(val) + 1;
   char* formatstr;

   assert(len <= INT_MAX);

   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;

   BMSallocMemoryArray(&formatstr, len);
   RatToString(val, formatstr, len);
   SCIPfputs(formatstr, certificate->derivationfile);
   BMSfreeMemoryArray(&formatstr);
}

/** prints a comment to the problem section of the certificate file */
void SCIPcertificatePrintProblemComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPfprintf(certificate->origfile, "# ");

   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   if( isorigfile )
      SCIPfprintf(certificate->origfile, "%s", formatstr); // todo: is this correct?
   else
      SCIPfprintf(certificate->transfile, "%s", formatstr); // todo: is this correct?

   va_end(ap);
   updateFilesize(certificate, 2 + strlen(formatstr));
}

/** prints a comment to the proof section of the certificate file */
void SCIPcertificatePrintProofComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
      return;

   SCIPfprintf(certificate->derivationfile, "# ");

   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   SCIPfprintf(certificate->derivationfile, "%s", formatstr); // todo: is this correct?
   va_end(ap);
   updateFilesize(certificate, 2 + strlen(formatstr));
}

/** prints version header */
void SCIPcertificatePrintVersionHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< shoud the line be printed to the origfile or the transfile */
   )
{
   assert(certificate != NULL);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "VER 1.0 \n");
}

/** prints variable section header */
void SCIPcertificatePrintVarHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   int                   nvars               /**< number of variables */
   )
{
   assert(certificate != NULL);
   assert(nvars >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "VAR %d \n", nvars);
}

/** prints integer section header */
void SCIPcertificatePrintIntHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   int                   nints               /**< number of integer variables */
   )
{
   assert(certificate != NULL);
   assert(nints >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "INT %d\n", nints);
}

/** prints constraint section header */
void SCIPcertificatePrintConsHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   int                   nconss,             /**< number of all constraints */
   int                   nboundconss         /**< number of bound constraints */
   )
{
   assert(certificate != NULL);
   assert(nconss >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "CON %d %d\n", nconss + nboundconss, nboundconss);
}

/** prints derivation section header */
void SCIPcertificatePrintDerHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< shoud the line be printed to the origfile or the transfile */
   )
{
   int nders;
   assert(certificate != NULL);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   if( !isorigfile )
      nders = certificate->indexcounter - certificate->conscounter;
   else
      nders = 0;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "DER %d\n", nders);
}

/** prints constraint */
void SCIPcertificatePrintCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   SCIP_Rational*        side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   SCIP_Rational**       val                 /**< coefficient array */
   )
{
   int i;
   SCIP_Longint index;

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   index = isorigfile ? certificate->indexcounter_ori : certificate->indexcounter;

   if( consname == NULL )
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "C%d %c ", index, sense);
   else
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", consname, sense);

   SCIPcertificatePrintProblemRational(certificate, isorigfile, side, 10);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d", len);

   for( i = 0; i < len; i++ )
   {
      /** @todo exip: perform line breaking before exceeding maximum line length */
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d ", ind[i]);
      SCIPcertificatePrintProblemRational(certificate, isorigfile, val[i], 10);
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   if( isorigfile )
      certificate->indexcounter_ori++;
   else
   {
      certificate->indexcounter++;
      certificate->lastinfo->isbound = FALSE;
   }

   if( !isorigfile )
   {
      certificate->conscounter++;
   }
}

/** print a line for an exact row to the certificate (without derivation), @param alternativerhs is used instead
 * of the real rhs of the row (infinity if real rhs should be used). This is necessary for integer cuts
 * where the rhs was rounded down from the original rhs */
static
SCIP_RETCODE SCIPcertificatePrintRow(
      SCIP_SET*          set,                /**< global SCIP settings */
      SCIP_CERTIFICATE*  certificate,        /**< certificate structure */
      SCIP_ROWEXACT*     rowexact,           /**< exact SCIP row */
      SCIP_Real          alternativerhs      /**< rhs to be used instead or rowexact->rhs (infinity to disable this) */
   )
{
   SCIP_ROW* row;
   SCIP_Rational* rhs;
   int i;

   assert(rowexact != NULL);

   row = SCIProwExactGetRow(rowexact);
   assert(SCIProwGetNNonz(row) == SCIProwExactGetNNonz(rowexact));

   SCIPcertificatePrintProofMessage(certificate, "L%d_%s %c ", certificate->indexcounter, row->name, 'L');

   /* if we rounded the rhs in cut tightening we need to first verify the cut without it */
   if( SCIPsetIsInfinity(set, alternativerhs) )
   {
      rhs = SCIProwExactGetRhs(rowexact);
      SCIPcertificatePrintProofRational(certificate, rhs, 10);
   }
   else
   {
      SCIP_CALL( RatCreateBuffer(set->buffer, &rhs) );
      RatSetReal(rhs, alternativerhs);
      SCIPcertificatePrintProofRational(certificate, rhs, 10);
      RatFreeBuffer(set->buffer, &rhs);
   }

   SCIPcertificatePrintProofMessage(certificate, " %d", SCIProwGetNNonz(row));

   for( i = 0; i < SCIProwGetNNonz(row); i++ )
   {
      SCIP_Rational* val;
      int varindex;
      /** @todo exip: perform line breaking before exceeding maximum line length */
      assert(rowexact->cols[i]->fpcol->var->index == rowexact->fprow->cols[i]->var->index);

      varindex = SCIPvarGetCertificateIndex(SCIPcolExactGetVar(SCIProwExactGetCols(rowexact)[i]));
      val = SCIProwExactGetVals(rowexact)[i];

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIPcertificatePrintProofRational(certificate, val, 10);
   }
   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   return SCIP_OKAY;
}

/** prints mir split for the specified aggrrow */
static
SCIP_RETCODE certificatePrintMirSplit(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_ROW*             row                 /**< row that split should be printed for */
   )
{
   SCIP_MIRINFO* mirinfo;
   SCIP_Rational* val;
   SCIP_VAR** vars;
   int i, j;
   SCIP_Real slackrhs;
   std::map<int, SCIP_Real> coefs;

   assert(SCIPhashmapExists(certificate->mirinfohash, (void*) row));

   mirinfo = (SCIP_MIRINFO*) SCIPhashmapGetImage(certificate->mirinfohash, (void*) row);

   vars = SCIPprobGetVars(prob);
   slackrhs = 0;

   RatDebugMessage("printing split disjunction <= %q \\/ >= %q+1  \n", mirinfo->rhs, mirinfo->rhs);

   for( i = 0; i < mirinfo->nsplitvars; ++i )
   {
      coefs[SCIPvarGetCertificateIndex(vars[mirinfo->varinds[i]])] += mirinfo->splitcoefficients[i];
      SCIPdebugMessage("+%g%s", mirinfo->splitcoefficients[i], SCIPvarGetName(vars[mirinfo->varinds[i]]));
   }
   for( i = 0; i < mirinfo->nslacks; ++i )
   {
      SCIP_ROW* slackrow = mirinfo->slackrows[i];
      SCIP_Real slackval = mirinfo->slackcoefficients[i];

      for( j = 0; j < SCIProwGetNNonz(slackrow); ++j)
      {
         SCIP_VAR* var = SCIPcolGetVar(SCIProwGetCols(slackrow)[j]);
         int varidx = SCIPvarGetCertificateIndex(var);
         SCIP_Real rowcoef = SCIProwGetVals(slackrow)[j];

         assert(SCIPisExactlyIntegral(set->scip, rowcoef * slackval));
         assert(SCIPvarIsBinary(var) || SCIPvarIsIntegral(var));
         coefs[varidx] += rowcoef  * slackval;
         assert(SCIPisExactlyIntegral(set->scip, coefs[varidx]));
      }
      if( mirinfo->slacksign[i] == 1 )
      {
         assert(SCIPisExactlyIntegral(set->scip, SCIProwGetRhs(slackrow) - SCIProwGetConstant(slackrow)));
         slackrhs += (SCIProwGetRhs(slackrow) - SCIProwGetConstant(slackrow)) * slackval;

         assert(SCIPisExactlyIntegral(set->scip, slackrhs));
         SCIPdebugMessage("+%gsrhs_%s", mirinfo->slackcoefficients[i], SCIProwGetName(mirinfo->slackrows[i]));
      }
      else
      {
         assert(SCIPisExactlyIntegral(set->scip, SCIProwGetLhs(slackrow) - SCIProwGetConstant(slackrow)));
         slackrhs += (SCIProwGetLhs(slackrow) - SCIProwGetConstant(slackrow)) * slackval;

         assert(SCIPisExactlyIntegral(set->scip, slackrhs));
         SCIPdebugMessage("+%gslhs_%s", mirinfo->slackcoefficients[i], SCIProwGetName(mirinfo->slackrows[i]));
      }
   }

   SCIP_CALL( RatCreateBuffer(set->buffer, &val) );

   /* transform the split back into original variable space -> undo the bound transformations */
   SCIPcertificatePrintProofMessage(certificate, "A%d_split %c ", certificate->indexcounter, 'L');

   /* add rhs change from integer slacks and print rhs */
   RatAddReal(val, mirinfo->rhs, slackrhs);
   SCIPcertificatePrintProofRational(certificate, val, 10);

   SCIPcertificatePrintProofMessage(certificate, " %d", coefs.size());

   for( auto it = coefs.begin(); it != coefs.end(); ++it )
   {
      /** @todo exip: perform line breaking before exceeding maximum line length */
      int varindex = it->first;
      RatSetReal(val, it->second);

      assert(RatIsIntegral(val));

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIPcertificatePrintProofRational(certificate, val, 10);
   }

   SCIPcertificatePrintProofMessage(certificate, " { asm } -1 \n");

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   SCIPcertificatePrintProofMessage(certificate, "A%d_split %c ", certificate->indexcounter, 'G');

   RatAddReal(val, mirinfo->rhs, slackrhs + 1.0);
   SCIPcertificatePrintProofRational(certificate, val, 10);

   SCIPcertificatePrintProofMessage(certificate, " %d", coefs.size());

   for( auto it = coefs.begin(); it != coefs.end(); ++it )
   {
      /** @todo exip: perform line breaking before exceeding maximum line length */
      int varindex = it->first;
      RatSetReal(val, it->second);

      assert(RatIsIntegral(val));

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIPcertificatePrintProofRational(certificate, val, 10);
   }

   SCIPcertificatePrintProofMessage(certificate, " { asm } -1 \n");

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   RatFreeBuffer(set->buffer, &val);

   return SCIP_OKAY;
}

/** prints all local bounds that differ from their global bounds as the bounds to take into account */
static
SCIP_RETCODE certificatePrintWeakDerStart(
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_Bool             local               /**< TRUE if the cut is only valid locally */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;
   int nboundentries;

   nboundentries = 0;

   vars = SCIPprobGetVars(prob);
   nvars = SCIPprobGetNVars(prob);

   /* count the number of needed entries */
   for( i = 0; i < nvars && local; i++ )
   {
      SCIP_VAR* var = vars[i];
      if( !RatIsEqual(var->exactdata->glbdom.lb, var->exactdata->locdom.lb) )
         nboundentries++;
      if( !RatIsEqual(var->exactdata->glbdom.ub, var->exactdata->locdom.ub) )
         nboundentries++;
   }

   SCIPcertificatePrintProofMessage(certificate, " { lin weak { %d", nboundentries);

   for( i = 0; i < nvars && nboundentries > 0; i++ )
   {
      SCIP_VAR* var = vars[i];

      if( !RatIsEqual(var->exactdata->glbdom.lb, var->exactdata->locdom.lb) )
      {
         SCIP_Longint index;
         SCIP_Rational* boundval;

         index = SCIPvarGetLbCertificateIndexLocal(var);
         boundval = SCIPvarGetLbLocalExact(var);
         SCIPcertificatePrintProofMessage(certificate, " L %d %d ", SCIPvarGetCertificateIndex(var), index);
         SCIPcertificatePrintProofRational(certificate, boundval, 10);
      }
      if( !RatIsEqual(var->exactdata->glbdom.ub, var->exactdata->locdom.ub) )
      {
         SCIP_Longint index;
         SCIP_Rational* boundval;

         index = SCIPvarGetUbCertificateIndexLocal(var);
         boundval = SCIPvarGetUbLocalExact(var);
         SCIPcertificatePrintProofMessage(certificate, " U %d %d ", SCIPvarGetCertificateIndex(var), index);
         SCIPcertificatePrintProofRational(certificate, boundval, 10);
      }
   }

   SCIPcertificatePrintProofMessage(certificate, " } ");

   return SCIP_OKAY;
}

/** prints verification of row as a MIR cut (viewed as a split cut) */
SCIP_RETCODE SCIPcertificatePrintMirCut(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_ROW*             row,                /**< the row to be printed */
   const char            sense               /**< sense of the constraint, i.e., G, L, or E */
   )
{
   SCIP_ROWEXACT* rowexact;
   SCIP_Rational* tmpval;
   SCIP_Rational* value;
   SCIP_Rational* oneminusf0;
   SCIP_AGGREGATIONINFO* aggrinfo;
   SCIP_Longint leftdisjunctionindex;
   SCIP_Longint rightdisjunctionindex;
   SCIP_MIRINFO* mirinfo;
   int arraypos;
   int i;

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   assert(row != NULL);
   assert(sense == 'L'); // only take care of this case for now
   assert(SCIPhashmapExists(certificate->aggrinfohash, (void*) row));

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmpval) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &oneminusf0) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &value) );

   rowexact = SCIProwGetRowExact(row);

   /* only do something if row does not already exist*/
   if( SCIPhashmapExists(certificate->rowdatahash, (void*) rowexact) )
      return SCIP_OKAY;

   SCIPdebugMessage("Printing MIR-certifiaction for row ");
   SCIPdebug(SCIProwExactPrint(rowexact,  set->scip->messagehdlr, NULL));

   /* get aggregation info and print aggregation row to certificate */
   aggrinfo = (SCIP_AGGREGATIONINFO*) SCIPhashmapGetImage(certificate->aggrinfohash, (void*) row);
   mirinfo = (SCIP_MIRINFO*) SCIPhashmapGetImage(certificate->mirinfohash, (void*) row);

   /* compute the correct split from the aggregation row, and print the two assumptions (\xi \le \lfloor \beta \rfloor), and  (\xi \ge \lfloor \beta + 1 \rfloor) */
   SCIP_CALL( certificatePrintMirSplit(set, prob, certificate, row) );

   leftdisjunctionindex = certificate->indexcounter - 2;
   rightdisjunctionindex = certificate->indexcounter - 1;

   /* if this aggregation depends on another no yet certified MIR cut, we need to print that first */
   for( i = 0; i < aggrinfo->naggrrows; i++ )
   {
      SCIP_ROWEXACT* aggrrow;
      aggrrow = SCIProwGetRowExact(aggrinfo->aggrrows[i]);
        assert(aggrrow != NULL);
      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) aggrrow) )
      {
         SCIP_CALL( SCIPcertificatePrintMirCut(set, lp, certificate, prob, aggrinfo->aggrrows[i], 'L') );
      }
   }

   /* print possibly missing derivations to the certifiacte */
   for( i = 0; i < aggrinfo->nnegslackrows; i++ )
   {
      SCIP_ROWEXACT* slackrow;
      slackrow = SCIProwGetRowExact(aggrinfo->negslackrows[i]);

      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow) )
      {
         SCIP_CALL( SCIPcertificatePrintMirCut(set, lp, certificate, prob, aggrinfo->negslackrows[i], 'L') );
      }
   }

   /* print possibly missing derivations to the certifiacte */
   for( i = 0; i < mirinfo->nslacks; i++ )
   {
      SCIP_ROWEXACT* slackrow;
      slackrow = SCIProwGetRowExact(mirinfo->slackrows[i]);

      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow) )
      {
         SCIP_CALL( SCIPcertificatePrintMirCut(set, lp, certificate, prob, mirinfo->slackrows[i], 'L') );
      }
   }

   /* print the mir cut with proof 1 * (\xi \le \lfloor \beta \rfloor) - (1/1-f)(\nu \ge 0) */
   /* we dont need the \nu \ge 0 part since it will be taken care of by the vipr completion part */
   assert(rowexact != NULL);

   SCIPcertificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs);

   /* calculate 1-f0 */
   RatSetReal(oneminusf0, 1.0);
   RatDiff(oneminusf0, oneminusf0, mirinfo->frac);

   certificatePrintWeakDerStart(certificate, prob, SCIProwIsLocal(row));

   /* 1 * (\xi \le \lfloor \beta \rfloor) we also have to add the correct multipliers for the negative slacks that were used here */
   SCIPcertificatePrintProofMessage(certificate, "%d %d ", 1 + aggrinfo->nnegslackrows + mirinfo->nslacks, leftdisjunctionindex);
   /* multiply with scaling parameter that was used during cut computation */
   RatSetReal(tmpval, mirinfo->scale);
   SCIPcertificatePrintProofRational(certificate, tmpval, 10);

   SCIPdebugMessage("Verifying left part of split disjunction, multipliers 1 and 1/%g \n", RatApproxReal(oneminusf0));

   SCIPdebugMessage("Correcting for negative continous slacks ( needed for v >= 0 part ) \n");
   for( i = 0; i < aggrinfo->nnegslackrows; i++ )
   {
       size_t key;
       SCIP_ROWEXACT* slackrow;
       slackrow = SCIProwGetRowExact(aggrinfo->negslackrows[i]);

       SCIPdebugMessage("adding (weight/(1-f0)) %g times row: ", aggrinfo->substfactor[i]);
       SCIPdebug(SCIProwExactPrint(slackrow, set->scip->messagehdlr, NULL));
       assert(slackrow != NULL);
       assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));
       RatSetReal(tmpval, aggrinfo->substfactor[i]);

       key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) slackrow);
       /* for ranged rows, the key always corresponds to the >= part of the row;
          therefore we need to increase it by one to get the correct key */
       if( !RatIsAbsInfinity(slackrow->rhs) && !RatIsAbsInfinity(slackrow->lhs) && !RatIsEqual(slackrow->lhs, slackrow->rhs) && aggrinfo->negslackweights[i] >= 0 )
          key += 1;

       SCIPcertificatePrintProofMessage(certificate, " %d ", key);
       RatMultReal(tmpval, tmpval, mirinfo->scale);
       SCIPcertificatePrintProofRational(certificate, tmpval, 10);
   }

   SCIPdebugMessage("Correcting for integer slacks  ( needed for v >= 0 part ) \n");
   for( i = 0; i < mirinfo->nslacks; i++ )
   {
       size_t key;
       SCIP_ROWEXACT* slackrow;
       slackrow = SCIProwGetRowExact(mirinfo->slackrows[i]);
       SCIP_Longint upar;

       assert(slackrow != NULL);
       assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));
       if( mirinfo->slackroundeddown[i] )
         RatSetReal(tmpval, 0);
       else
       {
         // value = weight * scale * sign -> compute (1 - fr)/(1-f0) where fr is the fractionality of value
         RatSetReal(tmpval, mirinfo->slackweight[i]);
         RatMultReal(tmpval, tmpval, mirinfo->slackscale[i]);
         RatMultReal(tmpval, tmpval, mirinfo->slacksign[i]);
         RatRoundInteger(&upar, tmpval, SCIP_R_ROUND_UPWARDS);
         RatDiffReal(tmpval, tmpval, upar);
         RatNegate(tmpval, tmpval);
         RatMultReal(tmpval, tmpval, mirinfo->slacksign[i]);
         RatDiv(tmpval, tmpval, oneminusf0);
         RatDebugMessage("tmpval 1 %g \n", RatApproxReal(tmpval));
         RatSetReal(tmpval, mirinfo->slackusedcoef[i]);;
         RatDiffReal(tmpval, tmpval, mirinfo->slackcoefficients[i]);
         RatDebugMessage("tmpval 2 %g (slacksign %d, splitcoef %g, cutval %g) \n", RatApproxReal(tmpval), mirinfo->slacksign[i], mirinfo->slackcoefficients[i], mirinfo->slackusedcoef[i]);
       }

       key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) slackrow);
       /* for ranged rows, the key always corresponds to the >= part of the row;
          therefore we need to increase it by one to get the correct key */
       if( !RatIsAbsInfinity(slackrow->rhs) && !RatIsAbsInfinity(slackrow->lhs) && !RatIsEqual(slackrow->lhs, slackrow->rhs) && mirinfo->slackweight[i] >= 0 )
          key += 1;

       SCIPcertificatePrintProofMessage(certificate, " %d ", key);
       RatMultReal(tmpval, tmpval, mirinfo->scale);
       SCIPcertificatePrintProofRational(certificate, tmpval, 10);
   }

   SCIPcertificatePrintProofMessage(certificate, " } -1 \n");

   SCIPdebugMessage("Verifying right part of split disjunction, multipliers -f/(1-f) and 1/1-f \n");

   /* print the mir cut with proof (-f/1-f) * (\xi \ge \lfloor \beta + 1 \rfloor) + (1/1-f)(\xi - \nu \le \beta) */
   SCIPcertificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs);

   certificatePrintWeakDerStart(certificate, prob, SCIProwIsLocal(row));
   SCIPcertificatePrintProofMessage(certificate, " %d ", 1 + aggrinfo->naggrrows + aggrinfo->nnegslackrows + mirinfo->nslacks);

   /* (-f/1-f) * (\xi \ge \lfloor \beta + 1 \rfloor) */
   SCIPcertificatePrintProofMessage(certificate, "%d ", rightdisjunctionindex);
   RatSet(tmpval, mirinfo->frac);
   RatNegate(tmpval, tmpval); /* -f */
   RatDiv(tmpval, tmpval, oneminusf0); /* -f/(1-f) */
   /* multiply with scaling factor that was used in cut derivation */
   RatMultReal(tmpval, tmpval, mirinfo->scale);
   SCIPcertificatePrintProofRational(certificate, tmpval, 10);

   RatDivReal(tmpval, tmpval, mirinfo->scale);

   SCIPdebugMessage("Correcting for negative continous slacks \n");
   /** @todo exip: really necessary? if so this is only for rounding error correction */
   /* we also have to add the correct multipliers for the negative continuous slacks that were used here */
   for( i = 0; i < aggrinfo->nnegslackrows; i++ )
   {
      size_t key;
      SCIP_ROWEXACT* slackrow;
      slackrow = SCIProwGetRowExact(aggrinfo->negslackrows[i]);

      RatSetReal(value, -aggrinfo->negslackweights[i]);
      RatDiv(value, value, oneminusf0);
      RatAddReal(value, value, aggrinfo->substfactor[i]);
      RatDebugMessage("adding %q times row (negative continous slacks) (%g aggweight %g substfactor): ", value, -aggrinfo->negslackweights[i] / RatApproxReal(oneminusf0), aggrinfo->substfactor[i]);
      SCIPdebug(SCIProwExactPrint(slackrow, set->scip->messagehdlr, NULL));

      assert(slackrow != NULL);
      assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));

      key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) slackrow);
      /* for ranged rows, the key always corresponds to the >= part of the row;
         therefore we need to increase it by one to get the correct key */
      if( !RatIsAbsInfinity(slackrow->rhs) && !RatIsAbsInfinity(slackrow->lhs) && !RatIsEqual(slackrow->lhs, slackrow->rhs) && aggrinfo->negslackweights[i] >= 0 )
         key += 1;

      SCIPcertificatePrintProofMessage(certificate, " %d ", key);
      RatMultReal(value, value, mirinfo->scale);
      SCIPcertificatePrintProofRational(certificate, value, 10);
   }

   SCIPdebugMessage("Correcting for %d integer slacks \n", mirinfo->nrounddownslacks);
   /* we also have to add the correct multipliers for the integer slacks that were used here */
   for( i = 0; i < mirinfo->nslacks; i++ )
   {
      size_t key;
      SCIP_ROWEXACT* slackrow;
      slackrow = SCIProwGetRowExact(mirinfo->slackrows[i]);

      if( mirinfo->slackroundeddown[i] )
      {
         SCIPdebugMessage("Rounded down intger slack on row %s \n", mirinfo->slackrows[i]->name);
         RatSetReal(value, mirinfo->slackweight[i]);
         RatMultReal(value, value, mirinfo->slackscale[i]);
         RatMultReal(value, value, mirinfo->slacksign[i]);
         RatDiffReal(value, value, mirinfo->slackcoefficients[i] * (-mirinfo->slacksign[i]));
         RatMultReal(value, value, mirinfo->slacksign[i]);
         RatDiv(value, value, oneminusf0);
      }
      else
      {
         SCIPdebugMessage("Rounded up intger slack on row %s \n", mirinfo->slackrows[i]->name);
         RatSetReal(value, mirinfo->slackweight[i]);
         RatMultReal(value, value, mirinfo->slackscale[i]);
         RatMultReal(value, value, mirinfo->slacksign[i]);
         RatDiffReal(value, value, (mirinfo->slackcoefficients[i] * -mirinfo->slacksign[i]) - 1); // fr exactly
         RatDiff(value, value, mirinfo->frac); // fr - f0
         RatDiv(value, value, oneminusf0); // (fr-f0) / (1-f0)
         RatAddReal(value, value, (mirinfo->slackcoefficients[i] * -mirinfo->slacksign[i]) - 1); // (down(ar) + (fr-f0) / (1-f0)
         RatMultReal(value, value, mirinfo->slacksign[i]);
         RatDebugMessage("Exact coefficient %q(%g), used coefficient %g\n", value, RatApproxReal(value), mirinfo->slackusedcoef[i]);

         RatAddReal(value, value, mirinfo->slackusedcoef[i]);
      }

      RatDebugMessage("adding %q(%g) times row: ", value, RatApproxReal(value));
      SCIPdebug(SCIProwExactPrint(slackrow, set->scip->messagehdlr, NULL));

      assert(slackrow != NULL);
      assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));

      key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) slackrow);
      /* for ranged rows, the key always corresponds to the >= part of the row;
         therefore we need to increase it by one to get the correct key */
      if( !RatIsAbsInfinity(slackrow->rhs) && !RatIsAbsInfinity(slackrow->lhs) && !RatIsEqual(slackrow->lhs, slackrow->rhs) && mirinfo->slackroundeddown[i] >= 0 )
         key += 1;

      SCIPcertificatePrintProofMessage(certificate, " %d ", key);
      RatMultReal(value, value, mirinfo->scale);
      SCIPcertificatePrintProofRational(certificate, value, 10);
   }

   SCIPdebugMessage("Adding %d aggregation rows \n", aggrinfo->naggrrows);
   /* we also have to add the correct multipliers for the aggregation rows that were used here */
   for( i = 0; i < aggrinfo->naggrrows; i++ )
   {
      size_t key;
      SCIP_ROWEXACT* aggrrow;
      aggrrow = SCIProwGetRowExact(aggrinfo->aggrrows[i]);

      RatSetReal(value, aggrinfo->weights[i]);
      RatDiv(value, value, oneminusf0);

      SCIPdebugMessage("adding (%g/%g) = %g times row: ", aggrinfo->weights[i], RatApproxReal(oneminusf0), RatApproxReal(value));
      SCIPdebug(SCIProwExactPrint(aggrrow, set->scip->messagehdlr, NULL));

      assert(aggrrow != NULL);
      assert(SCIPhashmapExists(certificate->rowdatahash, (void*) aggrrow));

      key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) aggrrow);
      /* for ranged rows, the key always corresponds to the >= part of the row;
         therefore we need to increase it by one to get the correct key */
      if( !RatIsAbsInfinity(aggrrow->rhs) && !RatIsAbsInfinity(aggrrow->lhs) && !RatIsEqual(aggrrow->lhs, aggrrow->rhs) && aggrinfo->weights[i] >= 0 )
         key += 1;

      SCIPcertificatePrintProofMessage(certificate, " %d ", key);
      RatMultReal(value, value, mirinfo->scale);
      SCIPcertificatePrintProofRational(certificate, value, 10);
   }

   SCIPcertificatePrintProofMessage(certificate, " } -1\n");

   /* print the unsplitting of the split disjunction */
   SCIPcertificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs);
   SCIPcertificatePrintProofMessage(certificate, " { uns %d %d  %d %d  } -1\n", certificate->indexcounter - 3, leftdisjunctionindex,
         certificate->indexcounter - 2, rightdisjunctionindex);

   /* if we rounded the rhs we need to still certify that part */
   if( !SCIPsetIsInfinity(set, mirinfo->unroundedrhs) )
   {
      /* print the row with the rounded rhs */
      SCIPcertificatePrintRow(set, certificate, rowexact, SCIPsetInfinity(set));
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }

   SCIP_CALL( SCIPhashmapInsert(certificate->rowdatahash, SCIProwGetRowExact(row), (void*)((size_t)certificate->indexcounter - 1)) );

   SCIP_CALL( SCIPcertificateFreeAggrInfo(set, certificate, lp, aggrinfo, row) );

   RatFreeBuffer(set->buffer, &value);
   RatFreeBuffer(set->buffer, &oneminusf0);
   RatFreeBuffer(set->buffer, &tmpval);

   /* remove the mirinfo, move last element to the now freed up one */
   arraypos = mirinfo->arpos;
   SCIP_CALL( SCIPhashmapRemove(certificate->mirinfohash, (void*) row) );
   RatFreeBlock(certificate->blkmem, &(mirinfo->rhs));
   RatFreeBlock(certificate->blkmem, &(mirinfo->frac));
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->splitcoefficients), mirinfo->nsplitvars);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->varinds), mirinfo->nsplitvars);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->upperused), mirinfo->nsplitvars);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->localbdused), mirinfo->nsplitvars);

   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackrows), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slacksign), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackcoefficients), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackweight), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackscale), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackusedcoef), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackroundeddown), mirinfo->nslacks);

   BMSfreeBlockMemory(certificate->blkmem, &mirinfo);
   if( arraypos != certificate->nmirinfos - 1 )
   {
      certificate->mirinfo[arraypos] = certificate->mirinfo[certificate->nmirinfos - 1];
      certificate->mirinfo[arraypos]->arpos = arraypos;
   }
   certificate->nmirinfos--;

   return SCIP_OKAY;
}

/** @todo exip: refactor this so duplicates and redundant bounds do not get printed */
/** prints a variable bound to the problem section of the certificate file and returns line index */
SCIP_RETCODE SCIPcertificatePrintBoundCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           boundname,          /**< name of the bound constraint */
   SCIP_VAR*             var,                /**< variable to print the bound cons for */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   )
{
   SCIP_Longint varindex;

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   varindex = SCIPvarGetCertificateIndex(var);

   if( !isorigfile )
   {
      certificate->indexcounter++;
      certificate->conscounter++;

#ifndef NDEBUG
      certificate->lastinfo->isbound = TRUE;
      certificate->lastinfo->boundtype = isupper ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER;
      certificate->lastinfo->varindex = SCIPvarGetCertificateIndex(var);
      RatSet(certificate->lastinfo->boundval, boundval);
#endif

      if( isupper )
         SCIPvarSetUbCertificateIndexGlobal(var, certificate->indexcounter - 1);
      else
         SCIPvarSetLbCertificateIndexGlobal(var, certificate->indexcounter - 1);

      if( boundname == NULL )
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "B%d %c ", certificate->indexcounter - 1, (isupper ? 'L' : 'G'));
      else
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", boundname, (isupper ? 'L' : 'G'));
   }
   else
   {
      certificate->indexcounter_ori++;
      if( boundname == NULL )
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "B%d %c ", certificate->indexcounter_ori - 1, (isupper ? 'L' : 'G'));
      else
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", boundname, (isupper ? 'L' : 'G'));
   }

   SCIPcertificatePrintProblemRational(certificate, isorigfile, boundval, 10);
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " 1 %d 1\n", varindex);

   return SCIP_OKAY;
}

/** installs updated node data in parent node */
SCIP_RETCODE SCIPcertificateUpdateParentData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_Rational*        newbound            /**< pointer to value of new bound, NULL if infeasible */
   )
{
   SCIP_CERTNODEDATA* nodedataparent;
   SCIP_CERTNODEDATA* nodedata;

   assert(node != NULL);
   assert(fileindex >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   /* Retrieve node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( RatIsLT(newbound, nodedata->derbound_inherit) )
      return SCIP_OKAY;

   /* If the node is the root node, then when only update the index and bound */
   if( SCIPnodeGetParent(node) == NULL )
   {
      certificate->derindex_root = fileindex;
      if( newbound == NULL )
         certificate->rootinfeas = TRUE;
      else
         RatSet(certificate->rootbound, newbound);
      return SCIP_OKAY;
   }

   /* Retrieve parent node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, SCIPnodeGetParent(node)));
   nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

   /* First ensure whether the node is left/right child of node parent, then left/rightfilled tells us if,
    * a bound has already been derived for this node, if yes then depending on the value of the bound we update it.
    * Else it is the first time this node is processed.Therefore, we also have to set left/rightfilled to TRUE
    */
   if( certificateIsLeftNode(certificate, node) )
   {
      nodedataparent->leftfilled = TRUE;
      if( RatIsGT(newbound, nodedataparent->derbound_left) )
      {
         nodedataparent->derindex_left = fileindex;
         RatSet(nodedataparent->derbound_left, newbound);
      }
      if( RatIsNegInfinity(newbound) )
      {
         nodedataparent->derindex_left = fileindex;
         nodedataparent->leftinfeas = TRUE;
      }
   }
   else
   {
      nodedataparent->rightfilled = TRUE;
      if( RatIsGT(newbound, nodedataparent->derbound_right) )
      {
         nodedataparent->derindex_right = fileindex;
         RatSet(nodedataparent->derbound_right, newbound);
      }
      if( RatIsNegInfinity(newbound) )
      {
         nodedataparent->rightinfeas = TRUE;
         nodedataparent->derindex_right = fileindex;
      }
   }

   return SCIP_OKAY;
}

/** Print a dual bound from an exact lp solution */
SCIP_RETCODE SCIPcertificatePrintDualboundExactLP(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             usefarkas           /**< should an infeasibility proof be printed? */
   )
{
   SCIP_Rational** vals;
   SCIP_Rational* val;
   SCIP_Rational* lowerbound;
   SCIP_Rational* farkasrhs;
   SCIP_Rational* tmp;
   SCIP_Longint* ind;
   int len;
   int i;
   unsigned long key;

   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
   SCIPlpExactGetObjval(lpexact, set, prob, tmp);

   /* only print line if bound improved */
   if( !usefarkas && RatIsLT(tmp, SCIPnodeGetLowerboundExact(node)) )
   {
      RatFreeBuffer(set->buffer, &tmp);
      return SCIP_OKAY;
   }

   /* if at root node set, objintegral flag */
   if( SCIPnodeGetParent(node) == NULL )
      certificate->objintegral = SCIPprobIsObjIntegral(prob);

   assert(lpexact!= NULL);
   assert(certificate->transfile != NULL);

   /* if needed extend vals array */
   if( lpexact->ncols + lpexact->nrows > certificate->valssize )
   {
      SCIP_CALL( RatReallocBlockArray(certificate->blkmem, &(certificate->vals), certificate->valssize, lpexact->nrows + lpexact->ncols + 50) );
      certificate->valssize =  lpexact->ncols + lpexact->nrows + 50;
   }

   vals = certificate->vals;

   SCIP_CALL( RatCreateBuffer(set->buffer, &farkasrhs) );

   SCIPsetAllocBufferArray(set, &ind, (lpexact->nrows + lpexact->ncols));

   len = 0;
   for( i = 0; i < lpexact->ncols; ++i )
   {
      SCIP_COLEXACT* col = lpexact->cols[i];
      SCIP_VAR* var = SCIPcolExactGetVar(col);

      if( usefarkas )
         val = col->farkascoef;
      else
         val = col->redcost;

      assert(!RatIsAbsInfinity(val));

      if( !RatIsZero(val) )
      {
         if( usefarkas )
            RatNegate(vals[len], val);
         else
            RatSet(vals[len], val);

         ind[len] = RatIsNegative(vals[len]) ? SCIPvarGetUbCertificateIndexLocal(var) : SCIPvarGetLbCertificateIndexLocal(var);

         RatDebugMessage("Column %d for var %s has index %l and farkas coef %q \n", col->index, col->var->name, ind[len], val);

         /* update farkasrhs */
         if( usefarkas )
         {
            val = RatIsNegative(vals[len]) ? SCIPvarGetUbLocalExact(var) : SCIPvarGetLbLocalExact(var);
            RatAddProd(farkasrhs, vals[len], val);
         }
         len++;
      }
   }

   for( i = 0; i < lpexact->nrows; ++i )
   {
      SCIP_ROWEXACT* row;
      row = lpexact->rows[i];
      val = usefarkas ? row->dualfarkas : row->dualsol;
      assert(!RatIsAbsInfinity(val));

      if( !RatIsZero(val) )
      {
         RatSet(vals[len], val);
         key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) row);

         if( key == 0 && SCIProwGetOrigintype(SCIProwExactGetRow(row)) == SCIP_ROWORIGINTYPE_SEPA )
         {
            SCIP_CALL( SCIPcertificatePrintMirCut(set, lpexact->fplp, certificate, prob, SCIProwExactGetRow(row), 'L') );
            key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) row);
         }
         else if( key == 0 && SCIProwExactGetNNonz(row) == 1 )
         {
            SCIP_VAR* var = SCIPcolExactGetVar(SCIProwExactGetCols(row)[0]);
            key = ind[len] = RatIsNegative(val) ? SCIPvarGetUbCertificateIndexLocal(var) : SCIPvarGetLbCertificateIndexLocal(var);
         }

         assert(key != 0);

         ind[len] = key;
         /* if we have a ranged row, and the dual corresponds to the upper bound,
          * the index for the rhs-constraint is one larger in the certificate */
         if( !RatIsEqual(row->lhs, row->rhs) && !RatIsAbsInfinity(row->lhs) && RatIsNegative(val) )
             ind[len] += 1;

         RatDebugMessage("Row (index %d, %s has index %l and farkas coef %q ", row->index, row->fprow->name, ind[len], val);
         SCIPdebug(SCIProwExactPrint(row, set->scip->messagehdlr, NULL) );

         /* update farkasrhs */
         if( usefarkas )
         {
            val = RatIsPositive(vals[len]) ? row->lhs : row->rhs;
            RatAddProd(farkasrhs, vals[len], val);
            RatDiffProd(farkasrhs, vals[len], row->constant);
         }

         len++;
      }
   }

   SCIP_CALL( RatCreateBuffer(set->buffer, &lowerbound) );
   if( usefarkas )
   {
      RatSetString(lowerbound, "inf");
      /* Scale proof to have RHS = 1 */
      // for( i = 0; i < len; ++i)
      //    RatDiv(vals[i], vals[i], farkasrhs);
   }
   else
   {
      /* vipr does not accept infinity, so in case of objlimit, get the objval from the lpi */
      if( !RatIsInfinity(lpexact->lpobjval) )
         RatSet(lowerbound, lpexact->lpobjval);
      else
      {
         SCIPlpiExactGetObjval(lpexact->lpiexact, lowerbound);
      }
   }

   SCIPcertificatePrintDualbound(certificate, NULL, lowerbound, len, ind, vals);
   SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound);
   SCIPcertificateSetInheritanceData(certificate, node, certificate->indexcounter - 1, lowerbound);

   RatFreeBuffer(set->buffer, &lowerbound);
   SCIPsetFreeBufferArray(set, &ind);
   RatFreeBuffer(set->buffer, &farkasrhs);
   RatFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** Print a dual bound from the pseudo solution
 *  in case of a bound change (branching), this happens before the bound change is processed;
 *  therefore we add the option to give on varindex, boundchgindex pair to pass directly to the method
 */
SCIP_RETCODE  SCIPcertificatePrintDualboundPseudo(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_NODE*            node,               /**< current node */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             lowerchanged,       /**< to the modified indices address a change in lb or ub? */
   int                   modifiedvarindex,   /**< index of modified variable, or -1 */
   SCIP_Longint          boundchangeindex,   /**< index of unprocessed bound change in the certificate, or -1 */
   SCIP_Real             psval               /**< the pseudo obj value (or inf to use exact lp value) */
   )
{
   SCIP_VAR** vars;
   SCIP_Rational* pseudoobjval;
   SCIP_Rational** bounds;
   SCIP_Longint* dualind;
   int nvars;
   int duallen;
   int i;
   int nnonzeros;

   assert(certificate != NULL);
   assert(prob != NULL);
   assert(set != NULL);
   assert((modifiedvarindex >= 0 && boundchangeindex >= 0) || (modifiedvarindex == -1 && boundchangeindex == -1) );

   /* only print if not -infinity and certificate is active */
   if( !set->exact_enabled || !SCIPcertificateIsActive(set, certificate) || SCIPsetIsInfinity(set, -psval) )
      return SCIP_OKAY;

   if( psval < SCIPnodeGetLowerbound(node) )
      return SCIP_OKAY;

   /* if at root node set, objintegral flag */
   if( SCIPnodeGetParent(node) == NULL )
      certificate->objintegral = SCIPprobIsObjIntegral(prob);

   vars = SCIPprobGetVars(prob);
   nvars = SCIPprobGetNVars(prob);
   SCIP_CALL( RatCreateBuffer(set->buffer, &pseudoobjval) );

   /* infinity means we use the exact lp value */
   if( SCIPsetIsInfinity(set, psval) )
      SCIPlpExactGetPseudoObjval(lpexact, set, prob, pseudoobjval);
   else
      RatSetReal(pseudoobjval, psval);
   duallen = SCIPprobGetNObjVars(prob, set);
   SCIP_CALL( RatCreateBufferArray(set->buffer, &bounds, duallen) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualind, duallen) );
   /* computes pseudo objective value (with bound change if necessary) */
   /*lint --e{838}*/

   nnonzeros = 0;
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Rational* obj = SCIPvarGetObjExact(vars[i]);
      if( !RatIsZero(obj) )
      {
         RatSet(bounds[nnonzeros], obj);

         assert(!RatIsAbsInfinity(bounds[nnonzeros]));

         /* retrieve the line in the certificate of the bound */
         if( RatIsPositive(obj) )
         {
            dualind[nnonzeros] = SCIPvarGetLbCertificateIndexLocal(vars[i]);
            if( lowerchanged && modifiedvarindex == SCIPvarGetCertificateIndex(vars[i]) )
               dualind[nnonzeros] = boundchangeindex;
         }
         else
         {
            dualind[nnonzeros] = SCIPvarGetUbCertificateIndexLocal(vars[i]);
            if( !lowerchanged && modifiedvarindex == SCIPvarGetCertificateIndex(vars[i]) )
               dualind[nnonzeros] = boundchangeindex;
         }

         nnonzeros++;
      }
   }
   assert(nnonzeros == duallen);
   /* print pseudo solution into certificate file */
   SCIPcertificatePrintDualbound(certificate, NULL, pseudoobjval, duallen,
      dualind, bounds);

   SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1,
      pseudoobjval) );

   SCIPcertificateSetInheritanceData(certificate, node, certificate->indexcounter - 1, pseudoobjval);

   SCIPsetFreeBufferArray(set, &dualind);
   RatFreeBufferArray(set->buffer, &bounds, nnonzeros);
   RatFreeBuffer(set->buffer, &pseudoobjval);

   return SCIP_OKAY;
}

/** prints the bound that a node inherits from its parent to the certificate */
SCIP_RETCODE SCIPcertificatePrintInheritedBound(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   )
{

   SCIP_CERTNODEDATA* nodedata;
   SCIP_Rational* lowerbound;

   assert(node != NULL);

   /* check whether certificate output should be created */
   if( certificate->transfile == NULL || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get the current node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( nodedata->inheritedbound && nodedata->assumptionindex_self != - 1 )
   {
      SCIP_Longint ind[1];
      SCIP_Rational* val;

      ind[0] = nodedata->derindex_inherit;

      RatCreateBuffer(set->buffer, &lowerbound);
      RatCreateBuffer(set->buffer, &val);

      RatSet(lowerbound, nodedata->derbound_inherit);
      RatSetInt(val, 1, 1);

      SCIPcertificatePrintDualbound(certificate, NULL, lowerbound, 1, ind, &val);
      SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound);

      RatFreeBuffer(set->buffer, &lowerbound);
      RatFreeBuffer(set->buffer, &val);
   }

   return SCIP_OKAY;
}

/** prints dual bound to proof section */
SCIP_Longint SCIPcertificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   SCIP_Longint*         ind,                /**< index array */
   SCIP_Rational**       val                 /**< array of dual multipliers */
   )
{
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return 0;

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   if( linename == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, "L%d ", certificate->indexcounter - 1);
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "%s ", linename);
   }

   if( RatIsInfinity(lowerbound) )
   {
      SCIPcertificatePrintProofMessage(certificate, "G 1 0");
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "G ");
      SCIPcertificatePrintProofRational(certificate, lowerbound, 10);
      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
   }

   if( val == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, " { lin ... } -1\n");
   }
   else
   {
      int i;

      SCIPcertificatePrintProofMessage(certificate, " { lin %d", len);

      for( i = 0; i < len; i++ )
      {
         /** @todo exip: perform line breaking before exceeding maximum line length */
         assert(!RatIsAbsInfinity(val[i]));
         SCIPcertificatePrintProofMessage(certificate, " %d ", ind[i]);
         SCIPcertificatePrintProofRational(certificate, val[i], 10);
      }
      SCIPcertificatePrintProofMessage(certificate, " } -1\n");
   }

   /* print rounding derivation */
   if( !RatIsNegInfinity(lowerbound) && certificate->objintegral && !RatIsIntegral(lowerbound) )
   {
      certificate->indexcounter++;
      certificate->lastinfo->isbound = FALSE;

      SCIPcertificatePrintProofMessage(certificate, "R%d G ", certificate->indexcounter - 1);
      updateFilesize(certificate, 4.0 + ceil(log10(certificate->indexcounter - 1 + 1)));
      RatRound(lowerbound, lowerbound, SCIP_R_ROUND_UPWARDS);

      SCIPcertificatePrintProofRational(certificate, lowerbound, 10);

      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }

   return (certificate->indexcounter - 1);
}

/** update the parent certificate node data when branching, print branching into certificate if not already present */
SCIP_RETCODE SCIPcertificatePrintBranching(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP informations */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR*             branchvar,          /**< the variable that gets branched on */
   SCIP_BOUNDTYPE        boundtype,          /**< the bounding type */
   SCIP_Real             newbound            /**< the new bound */
   )
{
   SCIP_CERTNODEDATA* nodedataparent;
   SCIP_CERTNODEDATA* nodedata;
   SCIP_Rational* branchbound;

   assert(node != NULL);
   assert(stat != NULL);

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE || lp->diving )
      return SCIP_OKAY;

   /* check whether output should be created */
   if ( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_CALL( RatCreateBuffer(set->buffer, &branchbound) );
   RatSetReal(branchbound, newbound);

   assert(RatIsIntegral(branchbound));

   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( branchvar != NULL )
   {
      nodedata->assumptionindex_self = printBoundAssumption(set, certificate, branchvar,
         branchbound, boundtype);
   }

   if( SCIPnodeGetParent(node) != NULL && certificate->transfile != NULL )
   {
      nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

      if( branchvar != NULL )
      {
         if( boundtype == SCIP_BOUNDTYPE_UPPER )
            nodedataparent->assumptionindex_right = nodedata->assumptionindex_self;
         if( boundtype == SCIP_BOUNDTYPE_LOWER )
            nodedataparent->assumptionindex_left = nodedata->assumptionindex_self;
      }
   }

   RatFreeBuffer(set->buffer, &branchbound);

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateNewNodeData(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   assert(stat != NULL );
   assert(node != NULL );

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* check whether output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemory(certificate->blkmem, &nodedata) );

   nodedata->derindex_left = -1;
   nodedata->derindex_right = -1;
   nodedata->assumptionindex_left = -1;
   nodedata->assumptionindex_right = -1;
   SCIP_CALL( RatCreateString(certificate->blkmem, &nodedata->derbound_left, "-inf") );
   SCIP_CALL( RatCreateString(certificate->blkmem, &nodedata->derbound_right, "-inf") );
   SCIP_CALL( RatCreateString(certificate->blkmem, &nodedata->derbound_inherit, "-inf") );
   nodedata->assumptionindex_self = -1;
   nodedata->leftinfeas = FALSE;
   nodedata->leftfilled = FALSE;
   nodedata->rightinfeas = FALSE;
   nodedata->rightfilled = FALSE;
   nodedata->inheritedbound = TRUE;
   nodedata->derindex_inherit = -1;
   if( SCIPnodeGetParent(node) != NULL )
   {
      SCIP_NODE* parent = SCIPnodeGetParent(node);
      SCIP_CERTNODEDATA* parentdata;
      assert(SCIPhashmapExists(certificate->nodedatahash, parent));
      parentdata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, (void*) parent);
      assert(parentdata != NULL);

      nodedata->derindex_inherit = parentdata->derindex_inherit;
      RatSet(nodedata->derbound_inherit, parentdata->derbound_inherit);
   }

   /* link the node to its nodedata in the corresponding hashmap */
   SCIP_CALL( SCIPhashmapSetImage(certificate->nodedatahash, node, (void*)nodedata) );

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificatePrintAggrrow(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_Bool             local               /**< true if local bound information can be used */
   )
{
   int i;
   SCIP_Rational* tmpval;
   SCIP_ROWEXACT* rowexact;
   SCIP_VAR** vars;



   SCIP_CALL( RatCreateBuffer(set->buffer, &tmpval) );
   vars = SCIPprobGetVars(prob);

   SCIPdebugMessage("printing certificate for aggrrow: ");
   SCIPdebug(SCIPaggrRowPrint(set->scip, aggrrow, NULL));

   SCIPcertificatePrintProofMessage(certificate, "L%d %c ", certificate->indexcounter, 'L');

   RatSetReal(tmpval, SCIPaggrRowGetRhs(aggrrow));

   SCIPcertificatePrintProofRational(certificate, tmpval, 10);

   SCIPcertificatePrintProofMessage(certificate, " %d", SCIPaggrRowGetNNz(aggrrow));

   for( i = 0; i < SCIPaggrRowGetNNz(aggrrow); i++ )
   {
      int varindex;
      /** @todo exip: perform line breaking before exceeding maximum line length */

      varindex = SCIPvarGetCertificateIndex(vars[SCIPaggrRowGetInds(aggrrow)[i]]);
      RatSetReal(tmpval, SCIPaggrRowGetValue(aggrrow, i));

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIPcertificatePrintProofRational(certificate, tmpval, 10);
   }

   certificatePrintWeakDerStart(certificate, prob, local);
   SCIPcertificatePrintProofMessage(certificate, " %d", naggrrows);
   for( i = 0; i < naggrrows; i++ )
   {
      size_t key;
      rowexact = SCIProwGetRowExact(aggrrows[i]);

      SCIPdebugMessage("adding %g times row: ", weights[i]);
      SCIPdebug(SCIProwExactPrint(rowexact, set->scip->messagehdlr, NULL));
      RatSetReal(tmpval, weights[i]);

      assert(rowexact != NULL);
      assert(SCIPhashmapExists(certificate->rowdatahash, (void*) rowexact));

      key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) rowexact);
      /* for ranged rows, the key always corresponds to the >= part of the row;
         therefore we need to increase it by one to get the correct key */
      if( !RatIsAbsInfinity(rowexact->rhs) && !RatIsAbsInfinity(rowexact->lhs) && !RatIsEqual(rowexact->lhs, rowexact->rhs) && weights[i] >= 0 )
         key += 1;

      SCIPcertificatePrintProofMessage(certificate, " %d ", key);
      SCIPcertificatePrintProofRational(certificate, tmpval, 10);
   }

   SCIPcertificatePrintProofMessage(certificate, " } -1\n");

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   RatFreeBuffer(set->buffer, &tmpval);

   return SCIP_OKAY;
}

/** free all aggregation information */
SCIP_RETCODE SCIPcertificateClearAggrinfo(
   SCIP*                 scip                /**< global SCIP data structure */
   )
{
   int i;
   SCIP_CERTIFICATE* certificate;

   certificate = SCIPgetCertificate(scip);

   if( certificate == NULL || certificate->aggrinfo == NULL )
      return SCIP_OKAY;

   assert(certificate != NULL);

   for( i = certificate->naggrinfos - 1; i >= 0; i-- )
   {
      SCIP_CALL( SCIPcertificateFreeAggrInfo(scip->set, certificate, scip->lp, certificate->aggrinfo[i], NULL) );
   }

   BMSfreeBlockMemoryArray(certificate->blkmem, &certificate->aggrinfo, certificate->aggrinfosize);
   SCIPhashmapRemoveAll(certificate->aggrinfohash);
   SCIPhashmapFree(&certificate->aggrinfohash);
   certificate->naggrinfos = 0;
   certificate->aggrinfohash = NULL;

   return SCIP_OKAY;
}

/** free aggregation information for row */
SCIP_RETCODE SCIPcertificateFreeAggrInfo(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate structure */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_AGGREGATIONINFO* aggrinfo,           /**< SCIP aggregation info */
   SCIP_ROW*             row                 /**< row that should be freed, or NULL if not needed */
   )
{
   int arraypos;
   int i;

   /* remove the (no longer needed aggrinfo), move last element to now freed spot */
   arraypos = aggrinfo->arpos;

   if( row != NULL )
   {
      SCIP_CALL( SCIPhashmapRemove(certificate->aggrinfohash, (void*) row) );
   }

   for( i = 0; i < aggrinfo->naggrrows; i++ )
   {
      SCIProwRelease(&(aggrinfo->aggrrows[i]), certificate->blkmem, set, lp);
   }
   for( i = 0; i < aggrinfo->nnegslackrows; i++ )
   {
      SCIProwRelease(&(aggrinfo->negslackrows[i]), certificate->blkmem, set, lp);
   }

   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->substfactor), aggrinfo->nnegslackrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->negslackweights), aggrinfo->nnegslackrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->negslackrows), aggrinfo->nnegslackrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->weights), aggrinfo->naggrrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->aggrrows), aggrinfo->naggrrows);
   SCIPaggrRowFree(set->scip, &(aggrinfo->aggrrow));
   BMSfreeBlockMemory(certificate->blkmem, &aggrinfo);
   if( arraypos != certificate->naggrinfos - 1 )
   {
      certificate->aggrinfo[arraypos] = certificate->aggrinfo[certificate->naggrinfos - 1];
      certificate->aggrinfo[arraypos]->arpos = arraypos;
   }
   certificate->naggrinfos--;

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateNewAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   )
{
   int i;
   SCIP_AGGREGATIONINFO* info;
   SCIP_CERTIFICATE* certificate;

   assert(scip != NULL );

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);

   /* check whether output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBlockMemory(scip, &info) );

   SCIPdebugMessage("adding aggrinfo, with %d rows to certficate \n", naggrrows);

   SCIP_CALL( SCIPaggrRowCopy(scip, &(info->aggrrow), aggrrow) );

   info->naggrrows = naggrrows;
   info->nnegslackrows = nnegslackrows;
   info->fileindex = certificate->indexcounter - 1;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->aggrrows), naggrrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->weights), naggrrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->negslackrows), nnegslackrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->negslackweights), nnegslackrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->substfactor), nnegslackrows) );
   for( i = 0; i < naggrrows; i++ )
   {
      SCIPdebugMessage("adding row %s with weight %g to aggrinfo \n", SCIProwGetName(aggrrows[i]), weights[i]);
      info->aggrrows[i] = aggrrows[i];
      info->weights[i] = weights[i];
      SCIProwCapture(aggrrows[i]);
   }
   for( i = 0; i < nnegslackrows; i++ )
   {
      SCIPdebugMessage("adding (implicitly) row %s with negative slack multiplier %g to aggrinfo \n", SCIProwGetName(negslackrows[i]), negslackweights[i]);
      info->negslackrows[i] = negslackrows[i];
      info->negslackweights[i] = negslackweights[i];
      SCIProwCapture(negslackrows[i]);
   }

   /* link the node to its nodedata in the corresponding hashmap */
   if( certificate->aggrinfosize == certificate->naggrinfos )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &certificate->aggrinfo, certificate->aggrinfosize, certificate->aggrinfosize + 100) );
      certificate->aggrinfosize += 100;
   }
   certificate->aggrinfo[certificate->naggrinfos] = info;
   info->arpos = certificate->naggrinfos;
   certificate->naggrinfos++;
   certificate->workingaggrinfo = TRUE;

   return SCIP_OKAY;
}

/** create a new split info structure for the current cut */
SCIP_RETCODE SCIPcertificateNewMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_MIRINFO* mirinfo;
   SCIP_CERTIFICATE* certificate;

   assert(scip != NULL );

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);

   /* check whether output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   if( certificate->mirinfosize == certificate->nmirinfos )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &certificate->mirinfo, certificate->mirinfosize, certificate->mirinfosize + 100) );
      certificate->mirinfosize += 100;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, &(certificate->mirinfo[certificate->nmirinfos])) );

   mirinfo =  certificate->mirinfo[certificate->nmirinfos];

   SCIPdebugMessage("adding mirinfo, with to certficate \n");

   SCIP_CALL( RatCreateBlock(SCIPblkmem(scip), &(mirinfo->rhs)) );
   SCIP_CALL( RatCreateBlock(SCIPblkmem(scip), &(mirinfo->frac)) );
   mirinfo->nlocalvars = 0;
   mirinfo->nslacks = 0;
   mirinfo->nrounddownslacks = 0;
   mirinfo->nsplitvars = SCIPgetNVars(scip);
   mirinfo->arpos = certificate->nmirinfos;
   mirinfo->scale = 1.0;
   mirinfo->unroundedrhs = SCIPinfinity(scip);

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(mirinfo->splitcoefficients), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(mirinfo->upperused), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(mirinfo->localbdused), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackrows), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackcoefficients), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slacksign), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackweight), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackscale), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackusedcoef), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackroundeddown), SCIPgetNVars(scip)) );
   mirinfo->varinds = NULL;

   certificate->nmirinfos++;
   certificate->workingmirinfo = TRUE;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPstoreCertificateActiveAggregationInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that aggregation-info is stored for */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_AGGREGATIONINFO* aggrinfo;

   assert(SCIPisExactSolve(scip));

   if( !SCIPisExactSolve(scip) || !SCIPisCertificateActive(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);
   assert(certificate->workingaggrinfo);

   aggrinfo = certificate->aggrinfo[certificate->naggrinfos - 1];
   certificate->workingaggrinfo = FALSE;

   assert(aggrinfo != NULL);

   SCIP_CALL( SCIPhashmapSetImage(certificate->aggrinfohash, (void*) row, (void*) aggrinfo) );

   return SCIP_OKAY;
}

/** stores the active mir information in the certificate data structures for a row */
SCIP_RETCODE SCIPstoreCertificateActiveMirInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that mirinfo is stored for */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_MIRINFO* mirinfo;
   int i;
   int csplit;

   assert(SCIPisExactSolve(scip));

   if( !SCIPisExactSolve(scip) || !SCIPisCertificateActive(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);
   assert(certificate->workingmirinfo);

   certificate->workingmirinfo = FALSE;

   mirinfo = certificate->mirinfo[certificate->nmirinfos - 1];

   assert(mirinfo != NULL);

   assert(mirinfo->nsplitvars == SCIPgetNVars(scip));

   csplit = 0;

   /* make the mirinfo sparse again */
   SCIPallocBlockMemoryArray(scip, &(mirinfo->varinds), mirinfo->nsplitvars);
   for( i = 0; i < mirinfo->nsplitvars; i++ )
   {
      if( mirinfo->splitcoefficients[i] != 0.0 || mirinfo->upperused[i] || mirinfo->localbdused[i] )
      {
         mirinfo->splitcoefficients[csplit] = mirinfo->splitcoefficients[i];
         mirinfo->upperused[csplit] = mirinfo->upperused[i];
         mirinfo->localbdused[csplit] = mirinfo->localbdused[i];
         mirinfo->varinds[csplit] = i;
         csplit++;
      }
   }

   mirinfo->nsplitvars = csplit;
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->splitcoefficients), SCIPgetNVars(scip), mirinfo->nsplitvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->varinds), SCIPgetNVars(scip), mirinfo->nsplitvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->upperused), SCIPgetNVars(scip), mirinfo->nsplitvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(mirinfo->localbdused), SCIPgetNVars(scip), mirinfo->nsplitvars) );

   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackrows, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slacksign, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackcoefficients, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackweight, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackscale, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackusedcoef, SCIPgetNVars(scip), mirinfo->nslacks) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &mirinfo->slackroundeddown, SCIPgetNVars(scip), mirinfo->nslacks) );

   SCIP_CALL( SCIPhashmapSetImage(certificate->mirinfohash, (void*) row, (void*) mirinfo) );

   return SCIP_OKAY;
}

/** frees the active mir information */
SCIP_RETCODE SCIPfreeCertificateActiveMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_MIRINFO* mirinfo;

   assert(SCIPisExactSolve(scip));

   if( !SCIPisExactSolve(scip) || !SCIPisCertificateActive(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);

   mirinfo = certificate->mirinfo[certificate->nmirinfos - 1];

   assert(mirinfo != NULL);

   /* if the mirinfo is used it gets tranformed into sparse format, don't free it in that case */
   if( !certificate->workingmirinfo )
      return SCIP_OKAY;

   assert(mirinfo->varinds == NULL);

   RatFreeBlock(SCIPblkmem(scip), &(mirinfo->rhs));
   RatFreeBlock(SCIPblkmem(scip), &(mirinfo->frac));

   SCIPfreeBlockMemoryArray(scip, &(mirinfo->splitcoefficients), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->upperused), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->localbdused), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackrows), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackcoefficients), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slacksign), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackweight), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackscale), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackusedcoef), SCIPgetNVars(scip));
   SCIPfreeBlockMemoryArray(scip, &(mirinfo->slackroundeddown), SCIPgetNVars(scip));
   SCIPfreeBlockMemory(scip, &mirinfo);
   certificate->nmirinfos--;
   certificate->workingmirinfo = FALSE;

   return SCIP_OKAY;
}

/** frees the active aggregation information */
SCIP_RETCODE SCIPfreeCertificateActiveAggregationInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_AGGREGATIONINFO* aggrinfo;

   assert(SCIPisExactSolve(scip));

   if( !SCIPisExactSolve(scip) || !SCIPisCertificateActive(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   assert(certificate != NULL);


   /* if the mirinfo is used it gets tranformed into sparse format, don't free it in that case */
   if( !certificate->workingaggrinfo )
      return SCIP_OKAY;

   aggrinfo = certificate->aggrinfo[certificate->naggrinfos - 1];

   assert(aggrinfo != NULL);

   SCIP_CALL( SCIPcertificateFreeAggrInfo(scip->set, certificate, scip->lp, aggrinfo, NULL) );

   certificate->workingaggrinfo = FALSE;

   return SCIP_OKAY;
}

/** prints unsplitting information to proof section */
int SCIPcertificatePrintUnsplitting(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   )
{
   SCIP_CERTNODEDATA* nodedata;
   SCIP_Rational* lowerbound;
   SCIP_Bool infeas;

   assert(node != NULL);

   /* check whether certificate output should be created */
   if( certificate->transfile == NULL || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get the current node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, node);
   SCIP_CALL( RatCreateBuffer(set->buffer, &lowerbound) );
   infeas = FALSE;

   assert(nodedata != NULL);

   if( nodedata->leftfilled && nodedata->rightfilled )
   {
      if( nodedata->leftinfeas && nodedata->rightinfeas )
      {
         infeas = TRUE;
         RatSetString(lowerbound, "inf");
      }
      else if( nodedata->leftinfeas )
         RatSet(lowerbound, nodedata->derbound_right);
      else if( nodedata->rightinfeas )
         RatSet(lowerbound, nodedata->derbound_left);
      else
         RatMIN(lowerbound, nodedata->derbound_left, nodedata->derbound_right);

      if( RatIsInfinity(nodedata->derbound_left) && RatIsInfinity(nodedata->derbound_right) )
         infeas = TRUE;

      certificate->indexcounter++;
      certificate->lastinfo->isbound = FALSE;

      SCIPcertificatePrintProofMessage(certificate, "U%d ", certificate->indexcounter - 1);

      if( infeas )
      {
         SCIPcertificatePrintProofMessage(certificate, "G 1 0");
      }
      else
      {
         SCIPcertificatePrintProofMessage(certificate, "G ");
         SCIPcertificatePrintProofRational(certificate, lowerbound, 10);
         SCIPcertificatePrintProofMessage(certificate, " ");
         SCIPcertificatePrintProofMessage(certificate, "OBJ");
      }

      assert(nodedata->derindex_right < certificate->indexcounter - 1);
      assert(nodedata->derindex_left < certificate->indexcounter - 1);

      SCIPcertificatePrintProofMessage(certificate, " { uns %d %d  %d %d  } -1\n", nodedata->derindex_left, nodedata->assumptionindex_left,
         nodedata->derindex_right, nodedata->assumptionindex_right);
      SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound) );
   }
   else
   {
      SCIPdebugMessage("Node %lld is a leaf! \n", SCIPnodeGetNumber(node));
      /* if a leaf has an inherited bound, we need to print a bound for it and update the parent data
         don't do it if we interrupted the solve, e.g. due to timeout */
      if( nodedata->inheritedbound && nodedata->assumptionindex_self != - 1 )
      {
         SCIP_Longint ind[1];
         SCIP_Rational* val;

         ind[0] = nodedata->derindex_inherit;

         RatCreateBuffer(set->buffer, &val);

         RatSet(lowerbound, nodedata->derbound_inherit);
         RatSetInt(val, 1, 1);

         SCIPcertificatePrintDualbound(certificate, NULL, lowerbound, 1, ind, &val);
         SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound);

         RatFreeBuffer(set->buffer, &val);
      }
   }

   certificateFreeNodeData(certificate, node);

   RatFreeBuffer(set->buffer, &lowerbound);
   return (certificate->indexcounter - 1);
}

/** prints RTP section with lowerbound and upperbound range */
void SCIPcertificatePrintRtpRange(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective */
   SCIP_Rational*        upperbound          /**< pointer to upper bound on the objective */
   )
{
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "RTP range ");
   if( RatIsNegInfinity(lowerbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "-inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, isorigfile, lowerbound, 10);
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " ");
   if( RatIsInfinity(upperbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, isorigfile, upperbound, 10);
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");
 }

/** prints RTP section for infeasibility */
void SCIPcertificatePrintRtpInfeas(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile          /**< should the original solution be printed or in transformed space */
   )
{
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "RTP infeas\n");
 }
