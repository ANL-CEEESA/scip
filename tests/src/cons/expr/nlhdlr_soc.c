/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_soc.c
 * @brief  tests quadratic nonlinear handler methods
 * @author Fabian Wegscheider
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

/* XXX: need the consdata struct because we don't have getNlhdlrs or findNlhdlrs; I don't add those function because I'm unsure
 * we actually need them
 */
#include "scip/cons_expr.c"
#include "scip/cons_expr_nlhdlr_soc.c"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;

static SCIP_CONSEXPR_NLHDLR* nlhdlr;
static SCIP_CONSHDLR* conshdlr;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   int h;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get quadratic handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   nlhdlr = SCIPfindConsExprNlhdlr(conshdlr, "soc");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 1.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -0.9, 0.7, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   //BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test suite */
TestSuite(nlhdlrsoc, .init = setup, .fini = teardown);

static
void checkData(
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata,
   SCIP_VAR**            vars,
   SCIP_Real*            coefs,
   SCIP_Real*            offsets,
   SCIP_Real*            transcoefs,
   int*                  transcoefidx,
   int*                  nnonzeroes,
   int                   nvars,
   int                   nterms,
   int                   ntranscoefs,
   SCIP_Real             constant
   )
{
   int i;

   cr_assert_not_null(nlhdlrexprdata->vars);
   cr_assert_not_null(nlhdlrexprdata->coefs);
   cr_assert_not_null(nlhdlrexprdata->offsets);
   cr_assert_not_null(nlhdlrexprdata->nnonzeroes);
   cr_assert_not_null(nlhdlrexprdata->transcoefs);
   cr_assert_not_null(nlhdlrexprdata->transcoefsidx);

   cr_expect_eq(nlhdlrexprdata->nvars, nvars);
   cr_expect_eq(nlhdlrexprdata->nterms, nterms);
   cr_expect_eq(nlhdlrexprdata->ntranscoefs, ntranscoefs);
   cr_expect_eq(nlhdlrexprdata->constant, constant);

   for( i = 0; i < nvars; ++i )
   {
      cr_assert_not_null(nlhdlrexprdata->vars[i]);
      cr_expect_eq(nlhdlrexprdata->vars[i], vars[i], );
   }

   for( i = 0; i < nterms; ++i )
      cr_expect_eq(nlhdlrexprdata->coefs[i], coefs[i]);

   for( i = 0; i < nterms; ++i )
      cr_expect_eq(nlhdlrexprdata->offsets[i], offsets[i]);

   for( i = 0; i < nterms; ++i )
      cr_expect_eq(nlhdlrexprdata->nnonzeroes[i], nnonzeroes[i]);

   for( i = 0; i < ntranscoefs; ++i )
      cr_expect_eq(nlhdlrexprdata->transcoefs[i], transcoefs[i]);

   for( i = 0; i < ntranscoefs; ++i )
      cr_expect_eq(nlhdlrexprdata->transcoefsidx[i], transcoefidx[i]);
}

/* detects ||x|| < 1 as soc expression */
Test(nlhdlrsoc, detectandfree1, .description = "detects simple norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   int i;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*) "(<x>^2 + <y>^2 + <z>^2)^0.5", NULL, &expr) );

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "norm", expr, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* detect */
   SCIP_CALL( SCIPinitlpCons(scip, cons, &infeasible) );
   cr_assert(!infeasible);

   SCIP_CALL( SCIPclearCuts(scip) ); /* we have to clear the separation store */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* setup expected data */
   SCIP_VAR* vars[4] = {x, y, z, SCIPgetConsExprExprAuxVar(expr)};
   SCIP_Real coefs[4] = {1.0, 1.0, 1.0, 1.0};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[4] = {1.0, 1.0, 1.0, 1.0};
   int transcoefsidx[4] = {0, 1, 2, 3};
   int nnonzeroes[4] = {1, 1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 4, 0.0);

   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* free expr and cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* detects ||x|| - y <= 0 as soc expression */
Test(nlhdlrsoc, detectandfree2, .description = "detects simple norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* normexpr;
   SCIP_Bool infeasible;
   int i;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*) "(<x>^2 + <y>^2 + <z>^2)^0.5 - <w>", NULL, &expr) );
   normexpr = SCIPgetConsExprExprChildren(expr)[0];

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "norm", expr, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* detect */
   SCIP_CALL( SCIPinitlpCons(scip, cons, &infeasible) );
   cr_assert(!infeasible);

   SCIP_CALL( SCIPclearCuts(scip) ); /* we have to clear the separation store */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* find the nlhdlr expr data */
   for( i = 0; i < normexpr->nenfos; ++i )
   {
      if( normexpr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = normexpr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* setup expected data */
   SCIP_VAR* vars[4] = {x, y, z, SCIPgetConsExprExprAuxVar(normexpr)};
   SCIP_Real coefs[4] = {1.0, 1.0, 1.0, 1.0};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[4] = {1.0, 1.0, 1.0, 1.0};
   int transcoefsidx[4] = {0, 1, 2, 3};
   int nnonzeroes[4] = {1, 1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 4, 0.0);

   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* free expr and cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* detects sqrt( 2*(x + 1)^2 + 3*(y + sin(x) + 2)^2 ) as soc expression */
Test(nlhdlrsoc, detectandfree3, .description = "detects more complex norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* normexpr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_Bool infeasible;
   SCIP_Bool changed;
   int i;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr,
         (char*) "(8 + 2*(<x> + 1)^2 + 3*(sin(<y>) - 2)^2)^0.5 + 2*(<w> - 1)", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed) );
   cr_expect(changed);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;
   normexpr = SCIPgetConsExprExprChildren(expr)[1];

   SCIPdismantleConsExprExpr(scip, normexpr);

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "soc", expr, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* detect */
   SCIP_CALL( SCIPinitlpCons(scip, cons, &infeasible) );
   cr_assert(!infeasible);

   SCIP_CALL( SCIPclearCuts(scip) ); /* we have to clear the separation store */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* find the nlhdlr expr data */
   for( i = 0; i < normexpr->nenfos; ++i )
   {
      if( normexpr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = normexpr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* setup expected data */
   SCIP_VAR* sinauxvar = SCIPgetConsExprExprAuxVar(normexpr->children[0]->children[1]->children[0]->children[2]);
   SCIP_VAR* vars[3] = {x, sinauxvar, SCIPgetConsExprExprAuxVar(normexpr)};
   SCIP_Real coefs[3] = {2.0, 3.0, 1.0};
   SCIP_Real offsets[3] = {1.0, -2.0, 0.0};
   SCIP_Real transcoefs[3] = {1.0, 1.0, 1.0};
   int transcoefsidx[3] = {0, 1, 2};
   int nnonzeroes[3] = {1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 3, 8.0);

   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* free expr and cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}