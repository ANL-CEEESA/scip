/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objconshdlr.cpp,v 1.4 2004/02/04 17:27:30 bzfpfend Exp $"

/**@file   objconshdlr.cpp
 * @brief  C++ wrapper for constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objconshdlr.h"




/*
 * Data structures
 */

/** constraint handler data */
struct ConshdlrData
{
   scip::ObjConshdlr* objconshdlr;      /**< constraint handler object */
   Bool             deleteobject;       /**< should the constraint handler object be deleted when conshdlr is freed? */
};




/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
static
DECL_CONSFREE(consFreeObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_free(scip, conshdlr) );

   /* free conshdlr object */
   if( conshdlrdata->deleteobject )
      delete conshdlrdata->objconshdlr;

   /* free conshdlr data */
   delete conshdlrdata;
   SCIPconshdlrSetData(conshdlr, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of constraint handler (called when problem solving starts) */
static
DECL_CONSINIT(consInitObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_init(scip, conshdlr) );

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called when problem solving exits) */
static
DECL_CONSEXIT(consExitObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_exit(scip, conshdlr) );

   return SCIP_OKAY;
}


/** solving start notification method of constraint handler (called when presolving was finished) */
static
DECL_CONSSOLSTART(consSolstartObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_solstart(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
DECL_CONSDELETE(consDeleteObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_delete(scip, conshdlr, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
DECL_CONSTRANS(consTransObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_trans(scip, conshdlr, sourcecons, targetcons) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
DECL_CONSINITLP(consInitlpObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_initlp(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** separation method of constraint handler */
static
DECL_CONSSEPA(consSepaObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_sepa(scip, conshdlr, conss, nconss, nusefulconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_enfolp(scip, conshdlr, conss, nconss, nusefulconss, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_enfops(scip, conshdlr, conss, nconss, nusefulconss, objinfeasible, result) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_check(scip, conshdlr, conss, nconss, sol, 
                  checkintegrality, checklprows, result) );

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
DECL_CONSPROP(consPropObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_prop(scip, conshdlr, conss, nconss, nusefulconss, result) );

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
DECL_CONSPRESOL(consPresolObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_presol(scip, conshdlr, conss, nconss, nrounds,
                  nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
                  nnewdelconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
                  nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
                  ndelconss, nupgdconss, nchgcoefs, nchgsides, result) );

   return SCIP_OKAY;
}


/** conflict variable resolving method of constraint handler */
static
DECL_CONSRESCVAR(consRescvarObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_rescvar(scip, conshdlr, cons, infervar) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_lock(scip, conshdlr, cons, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_unlock(scip, conshdlr, cons, nunlockspos, nunlocksneg) );

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
DECL_CONSACTIVE(consActiveObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_active(scip, conshdlr, cons) );

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
DECL_CONSDEACTIVE(consDeactiveObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_deactive(scip, conshdlr, cons) );

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
static
DECL_CONSENABLE(consEnableObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_enable(scip, conshdlr, cons) );

   return SCIP_OKAY;
}


/** constraint disabling notification method of constraint handler */
static
DECL_CONSDISABLE(consDisableObj)
{  /*lint --e{715}*/
   CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->objconshdlr != NULL);

   /* call virtual method of conshdlr object */
   CHECK_OKAY( conshdlrdata->objconshdlr->scip_disable(scip, conshdlr, cons) );

   return SCIP_OKAY;
}




/*
 * constraint handler specific interface methods
 */

/** creates the constraint handler for the given constraint handler object and includes it in SCIP */
RETCODE SCIPincludeObjConshdlr(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjConshdlr* objconshdlr,      /**< constraint handler object */
   Bool             deleteobject        /**< should the constraint handler object be deleted when conshdlr is freed? */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create obj constraint handler data */
   conshdlrdata = new CONSHDLRDATA;
   conshdlrdata->objconshdlr = objconshdlr;
   conshdlrdata->deleteobject = deleteobject;

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, objconshdlr->scip_name_, objconshdlr->scip_desc_, 
                  objconshdlr->scip_sepapriority_, objconshdlr->scip_enfopriority_, objconshdlr->scip_checkpriority_,
                  objconshdlr->scip_sepafreq_, objconshdlr->scip_propfreq_, objconshdlr->scip_needscons_,
                  consFreeObj, consInitObj, consExitObj, consSolstartObj,
                  consDeleteObj, consTransObj, consInitlpObj,
                  consSepaObj, consEnfolpObj, consEnfopsObj, consCheckObj, 
                  consPropObj, consPresolObj, consRescvarObj,
                  consLockObj, consUnlockObj,
                  consActiveObj, consDeactiveObj, 
                  consEnableObj, consDisableObj,
                  conshdlrdata) );

   return SCIP_OKAY;
}
