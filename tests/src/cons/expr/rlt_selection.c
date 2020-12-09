/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   rlt.c
 * @brief  tests rlt cut selection
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr.c"
#include "scip/sepastore.h"
#include "scip/lp.h"
#include "scip/scip.h"
#include "scip/var.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1o;
static SCIP_VAR* x2o;
static SCIP_VAR* x3o;
static SCIP_VAR* x4o;
static SCIP_VAR* y12o;
static SCIP_VAR* b1o;
static SCIP_VAR* b2o;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* y12;
static SCIP_VAR* b1;
static SCIP_VAR* b2;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* include rlt separator */
   SCIP_CALL( SCIPincludeSepaRlt(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* get rlt separator */
   sepa = SCIPfindSepa(scip, "rlt");
   assert(sepa != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1o, "x1", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2o, "x2", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3o, "x3", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4o, "x4", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y12o, "y12", 2.0, 4.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b1o, "b1", 0, 1, 1, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b2o, "b2", 0, 1, 1, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1o) );
   SCIP_CALL( SCIPaddVar(scip, x2o) );
   SCIP_CALL( SCIPaddVar(scip, x3o) );
   SCIP_CALL( SCIPaddVar(scip, x4o) );
   SCIP_CALL( SCIPaddVar(scip, y12o) );
   SCIP_CALL( SCIPaddVar(scip, b1o) );
   SCIP_CALL( SCIPaddVar(scip, b2o) );


   /* get SCIP into SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, x1o, &x1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x2o, &x2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x3o, &x3) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x4o, &x4) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y12o, &y12) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b1o, &b1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b2o, &b2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4o) );
   SCIP_CALL( SCIPreleaseVar(scip, &y12o) );
   SCIP_CALL( SCIPreleaseVar(scip, &b1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &b2o) );
   cr_assert(x1 != NULL);
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

static
void checkCut(SCIP_ROW* cut, SCIP_VAR** vars, SCIP_Real* vals, int nvars, SCIP_Real lhs, SCIP_Real rhs)
{
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Bool found;
   int i;
   int j;

   cr_assert(cut != NULL);
   cr_expect_eq(SCIProwGetNNonz(cut), nvars, "\nExpected %d nonz, got %d", nvars, SCIProwGetNNonz(cut));
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(cut), lhs));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(cut), rhs));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];
      found = FALSE;

      for( j = 0; j < nvars; ++j )
      {
         if( var == vars[j] )
         {
            cr_expect(SCIPisEQ(scip, coef, vals[j]));
            found = TRUE;
         }
      }

      if( !found )
         cr_expect(FALSE, "found an unknown variable");
   }
}

Test(rlt_selection, sepadata, .init = setup, .fini = teardown, .description = "test creation and freeing of separator data")
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeasible;
   SCIP_SEPADATA* sepadata;
//   BILINVARDATA** blvardatas;
   /* TODO update this to use the new structures */

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(sepadata->conshdlr != NULL);

   /* create a cons with some bilinear expressions */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <t_x1>*<t_x2> + <t_x1>*<t_x3> + <t_x4>*<t_x2> + <t_x4>^2 <= 1", TRUE, TRUE,
                 TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   SCIP_CALL( SCIPaddCons(scip, cons) ); /* adds locks */

   /* creates auxvars and creates disaggregation variables and row */
   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeasible) );

   SCIP_CALL( SCIPcollectConsExprBilinTerms(scip, conshdlr, &cons, 1) );

   SCIP_CALL( createSepaData(scip, sepadata) );

   cr_expect_eq(sepadata->nbilinvars, 4, "\nExpected 4 bilinear vars, got %d", sepadata->nbilinvars);

   cr_expect_eq(sepadata->varssorted[0], x4, "\nExpected varssorted[0] to be x4, got %s", SCIPvarGetName(sepadata->varssorted[0]));
   cr_expect_eq(sepadata->varssorted[1], x1, "\nExpected varssorted[1] to be x1, got %s", SCIPvarGetName(sepadata->varssorted[1]));
   cr_expect_eq(sepadata->varssorted[2], x2, "\nExpected varssorted[2] to be x2, got %s", SCIPvarGetName(sepadata->varssorted[2]));
   cr_expect_eq(sepadata->varssorted[3], x3, "\nExpected varssorted[3] to be x3, got %s", SCIPvarGetName(sepadata->varssorted[3]));

//   blvardatas = sepadata->bilinvardatas;
//
//   cr_expect_eq(blvardatas[2]->nvarbilinvars, 2, "\nExpected 2 bilinear vars for x2, got %d", blvardatas[2]->nvarbilinvars);
//   cr_expect_eq(blvardatas[3]->nvarbilinvars, 1, "\nExpected 1 bilinear vars for x3, got %d", blvardatas[3]->nvarbilinvars);
//   cr_expect_eq(blvardatas[1]->nvarbilinvars, 2, "\nExpected 2 bilinear vars for x1, got %d", blvardatas[1]->nvarbilinvars);
//   cr_expect_eq(blvardatas[0]->nvarbilinvars, 2, "\nExpected 2 bilinear vars for x4, got %d", blvardatas[0]->nvarbilinvars);
//
//   cr_expect_eq(blvardatas[0]->varbilinvars[0], x4,
//         "\nBilinear var 0 for x4 should be x4, got %s", SCIPvarGetName(blvardatas[0]->varbilinvars[0]));
//   cr_expect_eq(blvardatas[0]->varbilinvars[1], x2,
//         "\nBilinear var 1 for x4 should be x2, got %s", SCIPvarGetName(blvardatas[0]->varbilinvars[1]));
//   cr_expect_eq(blvardatas[1]->varbilinvars[0], x3,
//         "\nBilinear var 0 for x1 should be x3, got %s", SCIPvarGetName(blvardatas[1]->varbilinvars[0]));
//   cr_expect_eq(blvardatas[1]->varbilinvars[1], x2,
//         "\nBilinear var 1 for x1 should be x2, got %s", SCIPvarGetName(blvardatas[1]->varbilinvars[1]));
//   cr_expect_eq(blvardatas[2]->varbilinvars[0], x4,
//         "\nBilinear var 0 for x2 should be x4, got %s", SCIPvarGetName(blvardatas[2]->varbilinvars[0]));
//   cr_expect_eq(blvardatas[2]->varbilinvars[1], x1,
//         "\nBilinear var 1 for x2 should be x1, got %s", SCIPvarGetName(blvardatas[2]->varbilinvars[1]));
//   cr_expect_eq(blvardatas[3]->varbilinvars[0], x1,
//         "\nBilinear var 0 for x3 should be x1, got %s", SCIPvarGetName(blvardatas[3]->varbilinvars[0]));

   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIPfreeBuffer(scip, &sepadata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* frees the disaggregation row */
   SCIP_CALL( SCIPclearCuts(scip) );
}

Test(rlt_selection, projection, .init = setup, .fini = teardown, .description = "test projection of problem")
{
   SCIP_ROW** rows;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   RLT_SIMPLEROW* projrows;
   SCIP_Bool allcst;

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 3) );

   /* create test row1: -10 <= 4x1 - 7x2 + x3 <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x1, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x2, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x3, 1.0) );

   /* specify solution (only x3 is not at bound) */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   vars[0] = x1; vals[0] = 5.0;
   vars[1] = x2; vals[1] = -6.0;
   vars[2] = x3; vals[2] = 2.0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 3, vars, vals) );
   cr_assert(SCIProwGetNNonz(rows[0]) == 3);

   SCIP_CALL( createProjRows(scip, rows, 1, sol, &projrows, TRUE, &allcst) );

   /* check results */

   /* the projected cut should be: -72 <= x3 <= -57 */
   cr_assert_eq(projrows[0].nnonz, 1, "\nExpected 1 non-zero in the projected row, got %d", projrows[0].nnonz);
   cr_assert_eq(projrows[0].coefs[0], 1.0, "\nExpected coef 0 in projected row 0 to be 1.0, got %f", projrows[0].coefs[0]);
   cr_assert_eq(projrows[0].vars[0], x3, "\nExpected var 0 in projected row 0 to be x3, got %s", SCIPvarGetName(projrows[0].vars[0]));
   cr_assert_eq(projrows[0].cst, 0.0, "\nExpected the const in projected row to be 0.0, got %f", projrows[0].cst);
   cr_assert_eq(projrows[0].lhs, -72.0, "\nExpected the lhs in projected row to be -72.0, got %f", projrows[0].lhs);
   cr_assert_eq(projrows[0].rhs, -57.0, "\nExpected the rhs in projected row to be -57.0, got %f", projrows[0].rhs);

   /* free memory */
   freeProjRows(scip, &projrows, 1);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseRow(scip, &rows[0]) );
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}

Test(rlt_selection, compute_projcut, .init = setup, .fini = teardown, .description = "test projected cut computation")
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   RLT_SIMPLEROW projrow;
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* cut;
   SCIP_ROW* row;
   SCIP_Bool success;
   SCIP_Real cut_val;

   /* create row: -10 <= x1 + 2x2 - x3 <= 20 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, "test_row", -10.0, 20.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x1, 1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x2, 2.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x3, -1.0) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 3) );

   /* specify solution (none of the variables is at bound) */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   vars[0] = x1; vals[0] = 0.0;
   vars[1] = x2; vals[1] = -5.0;
   vars[2] = x3; vals[2] = 2.0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 3, vars, vals) );

   /* fill in sepadata */
   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(sepadata->conshdlr != NULL);
   sepadata->maxusedvars = 4;

   /* create projected LP with row -10 <= x1 + 2x2 - x3 <= 20 */
   SCIP_CALL( createProjRow(scip, &(projrow), row, sol, FALSE) );

   /* compute a cut with x1, lb and lhs */
   SCIP_CALL( computeRltCut(scip, sepa, sepadata, &cut, NULL, &projrow, sol, NULL, NULL, x1, &success, TRUE, TRUE,
         FALSE, FALSE, TRUE) );
   assert(success);

   /* the cut should be -8 <= 8x1 */
   cut_val = 8.0;
   checkCut(cut, &x1, &cut_val, 1, -8.0, SCIPinfinity(scip));

   /* free memory */
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   freeProjRow(scip, &projrow);
   SCIPfreeBuffer(scip, &sepadata);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

Test(rlt_selection, compute_clique_cuts, .init = setup, .fini = teardown, .description = "test cut computation when cliques are present")
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* cut;
   SCIP_ROW** rows;
   SCIP_Bool success;
   SCIP_Bool infeasible;
   SCIP_VAR* clique_vars[2];
   SCIP_Bool clique_vals[2];
   int nbdchgs;
   SCIP_VAR** cut_vars;
   SCIP_Real* cut_vals;

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

   /* create test row1: -10 <= b2 <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], b2, 1.0) );

   scip->stat->nnz = 1;

   /* specify solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   vars[0] = b1;
   vals[0] = 1;
   vars[1] = b2;
   vals[1] = 0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 2, vars, vals) );

   /* fill in sepadata */
   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(sepadata->conshdlr != NULL);

   /*add a clique (1-b1) + (1-b2) <= 1*/
   clique_vars[0] = b1;
   assert(clique_vars[0] != NULL);
   clique_vars[1] = b2;
   assert(clique_vars[1] != NULL);
   clique_vals[0] = 0;
   clique_vals[1] = 0;
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /* compute a cut with b1, lb and lhs */
   SCIP_CALL( computeRltCut(scip, sepa, sepadata, &cut, rows[0], NULL, sol, NULL, NULL, b1, &success, TRUE, TRUE,
         FALSE, FALSE, FALSE) );
   assert(success);

   /* the cut should be 1 <= 11b1 + b2 */
   cut_vars = (SCIP_VAR*[4]) {b2, b1};
   cut_vals = (SCIP_Real[4]) {1.0, 11.0};
   checkCut(cut, cut_vars, cut_vals, 2, 1.0, SCIPinfinity(scip));

   /* free memory */
   if( cut != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }
   SCIP_CALL( SCIPreleaseRow(scip, &rows[0]) );
   SCIPfreeBuffer(scip, &sepadata);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}

///* test row marking */
//Test(rlt_selection, mark, .init = setup, .fini = teardown, .description = "test row marking")
//{
//   SCIP_ROW** rows;
//   SCIP_SOL* sol;
//   SCIP_VAR** vars;
//   SCIP_Real* vals;
//   SCIP_SEPADATA* sepadata;
//   SCIP_Bool success, infeasible;
//   SCIP_CONSEXPR_EXPR* expr;
//   SCIP_CONS* cons;
//   int c;
////   SCIP_LP* lp;
//   int* row_marks;
//   int* row_idcs;
//   int nmarked;
//   SCIP_HASHMAP* row_to_pos;
//
//   SCIPallocBufferArray(scip, &rows, 1);
//   SCIPallocBufferArray(scip, &vars, 6);
//   SCIPallocBufferArray(scip, &vals, 6);
//   SCIP_CALL( SCIPallocCleanBufferArray(scip, &row_marks, 1) );
//   SCIP_CALL( SCIPallocBufferArray(scip, &row_idcs, 1) );
//
////   SCIPlpCreate(&lp, scip->set, scip->messagehdlr, scip->stat, "lp");
//
//   /* create a row: -10 <= 4x1 - 7x2 + x3 <= 5 */
//   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, TRUE, FALSE) );
//   SCIPaddRow(scip, rows[0], FALSE, &infeasible);
////   SCIPlpAddRow(lp, SCIPblkmem(scip), scip->set, scip->eventqueue, scip->eventfilter, rows[0], 0);
//   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x1, 4.0) );
//   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x2, -7.0) );
//   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x3, 1.0) );
//   SCIPaddRow(scip, rows[0], FALSE, &infeasible);
//
//   SCIP_CALL( SCIPhashmapCreate(&row_to_pos, SCIPblkmem(scip), 1) );
//   SCIP_CALL( SCIPhashmapSetImageInt(row_to_pos, (void*)(size_t)SCIProwGetIndex(rows[0]), 0) );
//
//   /* create a cons with some bilinear expressions */
//   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <t_x1>*<t_x2> + <t_x1>*<t_x3> + <t_x3>*<t_x2> <= 1", TRUE, TRUE,
//                            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
//   cr_assert(success);
//   expr = SCIPgetExprConsExpr(scip, cons);
//   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
//   {
//      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], NULL) );
//   }
//
//   /* specify solution */
//   SCIPcreateSol(scip, &sol, NULL);
//   vars[0] = x1; vals[0] = 5.0;   /* [-1,5] */
//   vars[1] = x2; vals[1] = -5.0;  /* [-6,-3] */
//   vars[2] = x3; vals[2] = 2.0;   /* [1,3] */
//
//   vars[3] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]); vals[3] = -25.0; /* y12 = x1*x2 */
//   vars[4] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[1]); vals[4] = 10.0; /* y13 = x1*x3 */
//   vars[5] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[2]); vals[5] = -9.0; /* y23 > x2*x3 */
//   SCIP_CALL( SCIPsetSolVals(scip, sol, 6, vars, vals) );
//   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
//   {
//      SCIPevalConsExprExpr(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], sol, 0);
//   }
//
//   /* fill in sepadata */
//   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
//   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
//   cr_assert(sepadata->conshdlr != NULL);
//   SCIP_CALL( createSepaData(scip, sepadata) );
//   sepadata->maxusedvars = 4;
//   sepadata->maxncuts = 10;
//   sepadata->maxunknownterms = 100;
//
//   SCIPinfoMessage(scip, NULL, "\nvarssorted: ");
//   for( int i = 0; i < sepadata->nbilinvars; ++i )
//   {
//      SCIPinfoMessage(scip, NULL, "%s; prior = %d", SCIPvarGetName(sepadata->varssorted[i]), sepadata->varpriorities[i]);
//   }
//
//   /* mark rows */
//
//   row_marks[0] = 0;
//
//   /* multiply by x1 */
//   markRowsXj(scip, sepadata, conshdlr, sol, 0, FALSE, row_to_pos, row_marks, row_idcs, &nmarked);
//
//   /* no products involving x1 are violated => no mark should have been added */
//   cr_assert_eq(row_marks[0], 0, "\nExpected row_marks[0] = 0 for x1, got %d", row_marks[0]);
//
//
//   /* multiply by x2 */
//   markRowsXj(scip, sepadata, conshdlr, sol, 1, FALSE, row_to_pos, row_marks, row_idcs, &nmarked);
//   /* TODO this currently fails because the row doesn't get properly added to LP */
//   cr_assert_eq(row_marks[0], 1, "\nExpected row_marks[0] = 1 for x2, got %d", row_marks[0]);
//
//   row_marks[0] = 0;
//
//   /* multiply by x3 */
//   markRowsXj(scip, sepadata, conshdlr, sol, 2, FALSE, row_to_pos, row_marks, row_idcs, &nmarked);
//   cr_assert_eq(row_marks[0], 2, "\nExpected row_marks[0] = 2 for x3, got %d", row_marks[0]);
//
//   /* free memory */
//   SCIPclearCuts(scip);
//   SCIP_CALL( freeSepaData(scip, sepadata) );
//   SCIPfreeBuffer(scip, &sepadata);
//   SCIPfreeSol(scip, &sol);
//   SCIPreleaseCons(scip, &cons);
//   SCIPreleaseRow(scip, &rows[0]);
////   SCIPlpFree(&lp, SCIPblkmem(scip), scip->set, scip->eventqueue, scip->eventfilter);
//   row_marks[0] = 0;
//   SCIPfreeBufferArray(scip, &row_idcs);
//   SCIPfreeCleanBufferArray(scip, &row_marks);
//   SCIPfreeBufferArray(scip, &vals);
//   SCIPfreeBufferArray(scip, &vars);
//   SCIPfreeBufferArray(scip, &rows);
//}
