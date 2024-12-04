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

/**@file   examples/GenDisj/src/cmain.c
 * @brief  main file for general disjunction example
 * @author Suresh Bolusani
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdbool.h>

#define MAX_LINE_LEN 1024
#define MAX_VARS_PER_LINE 100

/** Helper function to read an initial solution file */
static
SCIP_RETCODE readSolutionFile(SCIP* scip, const char* filename, SCIP_HASHMAP* var_values, SCIP_SOL* solution)
{
   SCIP_FILE* file;
   SCIP_Bool error;
   int lineno;
   SCIP_Bool unknownvariablemessage;
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   error = FALSE;
   lineno = 0;
   unknownvariablemessage = FALSE;

   /* read the file */
   while( !SCIPfeof(file) && !(error) )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      char format[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real value;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* there are some lines which may precede the solution information */
      if( SCIPstrncasecmp(buffer, "solution status:", 16) == 0 || SCIPstrncasecmp(buffer, "objective value:", 16) == 0 ||
         SCIPstrncasecmp(buffer, "Log started", 11) == 0 || SCIPstrncasecmp(buffer, "Variable Name", 13) == 0 ||
         SCIPstrncasecmp(buffer, "All other variables", 19) == 0 || strspn(buffer, " \n\r\t\f") == strlen(buffer) ||
         SCIPstrncasecmp(buffer, "NAME", 4) == 0 || SCIPstrncasecmp(buffer, "ENDATA", 6) == 0 ||    /* allow parsing of SOL-format on the MIPLIB 2003 pages */
         SCIPstrncasecmp(buffer, "=obj=", 5) == 0 )    /* avoid "unknown variable" warning when reading MIPLIB SOL files */
         continue;

      /* parse the line */
      (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%ds %%%ds %%%ds\n", SCIP_MAXSTRLEN, SCIP_MAXSTRLEN, SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, varname, valuestring, objstring);
      if( nread < 2 )
      {
         SCIPerrorMessage("Invalid input line %d in solution file <%s>: <%s>.\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> in line %d of solution file <%s>\n",
               varname, lineno, filename);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the value */
      if( SCIPstrncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( SCIPstrncasecmp(valuestring, "+inf", 4) == 0 || SCIPstrncasecmp(valuestring, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( SCIPstrncasecmp(valuestring, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else if( SCIPstrncasecmp(valuestring, "unknown", 7) == 0 )
      {
         value = SCIP_UNKNOWN;
      }
      else
      {
         /* coverity[secure_coding] */
         nread = sscanf(valuestring, "%lf", &value);
         if( nread != 1 )
         {
            SCIPerrorMessage("Invalid solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
               valuestring, varname, lineno, filename);
            error = TRUE;
            break;
         }
      }
      SCIP_CALL(SCIPcaptureVar(scip, var)); // Capture for hashmap storage
      SCIP_CALL(SCIPhashmapInsert(var_values, (void*)var, (void*)(size_t)value));

      SCIP_CALL(SCIPsetSolVal(scip, solution, var, value));
   }

   SCIPfclose(file);
   return SCIP_OKAY;
}

/** reads parameters */
static
SCIP_RETCODE readParams(
      SCIP*                 scip,               /**< SCIP data structure */
      const char*           filename            /**< parameter file name, or NULL */
      )
{
   if( SCIPfileExists(filename) )
   {
      SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
   }
   else
      SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);

   return SCIP_OKAY;
}

/** read problem */
static
SCIP_RETCODE readProblem(
      SCIP*                 scip,               /**< SCIP data structure */
      const char*           filename            /**< input file name */
      )
{
   /********************
    * Problem Creation *
    ********************/

   if ( filename != NULL )
   {
      if ( SCIPfileExists(filename) )
      {
         SCIPinfoMessage(scip, NULL, "read problem <%s> ...\n\n", filename);
         SCIP_CALL( SCIPreadProb(scip, filename, NULL) );

         return SCIP_OKAY;
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "problem file <%s> not found.\n", filename);

         return SCIP_NOFILE;
      }
   }

   return SCIP_NOFILE;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runSCIP(
      int                   argc,               /**< number of shell parameters */
      char**                argv                /**< array with shell parameters */
      )
{
   if (argc < 5)
   {
      printf("Usage: %s <inst_file.mps(.gz)> <param_file.set> <gen_disj_file.txt> <sol_file.sol(.gz)>\n", argv[0]);
      return SCIP_ERROR;
   }

   const char* inst_filename = argv[1];
   const char* param_filename = argv[2];
   const char* gendisj_filename = argv[3];
   const char* sol_filename = argv[4];

   SCIP* scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* version information */
   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* read parameter file */
   SCIP_CALL( readParams(scip, param_filename) );

   /* read instance file */
   SCIP_CALL( readProblem(scip, inst_filename) );

   /* create a hashmap for storing the initial solution */
   SCIP_HASHMAP* var_values;
   SCIP_CALL(SCIPhashmapCreate(&var_values, SCIPblkmem(scip), 256));

   /* create an empty solution */
   SCIP_SOL* solution;
   SCIP_CALL(SCIPcreateSol(scip, &solution, NULL));

   /* read the initial solution file */
   SCIP_CALL(readSolutionFile(scip, sol_filename, var_values, solution));
   SCIPinfoMessage(scip, NULL, "Solution file '%s' loaded successfully.\n", sol_filename);

   /* open the general disjunctions file */
   SCIP_FILE* gendisj_file = SCIPfopen(gendisj_filename, "r");
   if( gendisj_file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", gendisj_filename);
      SCIPprintSysError(gendisj_filename);
      SCIP_CALL(SCIPfree(&scip));
      return SCIP_NOFILE;
   }

   char* token;
   int constraint_idx = 0;

   SCIP_Bool error;
   int lineno;

   error = FALSE;
   lineno = 0;

   /* read the gendisj file */
   while( !SCIPfeof(gendisj_file) && !(error) )
   {
      char buffer[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      int num_existing_vars = 0;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), gendisj_file) == NULL )
         break;
      lineno++;

      token = strtok(buffer, " \t\n");

      /* parse the number of existing variables */
      if( token != NULL )
         num_existing_vars = atoi(token);

      /* parse the variables */
      SCIP_VAR* vars[num_existing_vars + 1];
      SCIP_Real vals[num_existing_vars + 1];
      SCIP_Real new_var_value = 0.0;
      int var_count = 0;

      while( (token = strtok(NULL, " \t\n")) != NULL && (var_count < num_existing_vars) )
      {
         SCIP_Real coeff = 1.0;

         /* check for a leading '-' sign */
         if (token[0] == '-')
         {
            coeff = -1.0;
            token++; /* strip the '-' to get the actual variable name */
         }

         char* varname = strdup(token);

         var = SCIPfindVar(scip, varname);
         if (var == NULL)
         {
            SCIPinfoMessage(scip, NULL, "Variable '%s' not found. Skipping line.\n", varname);
            free(varname);
            goto cleanup_and_continue;
         }
         free(varname);

         vars[var_count] = var;
         vals[var_count] = -1.0 * coeff;
         var_count++;

         /* compute the value of the new variable based on the initial solution */
         void* var_value = SCIPhashmapGetImage(var_values, (void*)var);
         new_var_value += coeff * (var_value ? (SCIP_Real)(size_t)var_value : 0.0);
      }

      /* create a new variable `zdisjN` */
      char new_var_name[SCIP_MAXSTRLEN];
      SCIPsnprintf(new_var_name, SCIP_MAXSTRLEN, "zdisj%d", constraint_idx + 1);
      SCIP_VAR* new_var;

      // TODO: SCIP_VARTYPE_IMPLINT??
      SCIP_CALL(SCIPcreateVarBasic(scip, &new_var, new_var_name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER));
      SCIP_CALL(SCIPaddVar(scip, new_var));

      vars[var_count] = new_var;
      vals[var_count] = 1.0;
      var_count++;

      /* create and add the constraint */
      char cons_name[SCIP_MAXSTRLEN];
      SCIPsnprintf(cons_name, SCIP_MAXSTRLEN, "c%d", SCIPgetNConss(scip) + 1);
      SCIP_CONS* cons;

      SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons, cons_name, var_count, vars, vals, 0.0, 0.0));
      SCIP_CALL(SCIPaddCons(scip, cons));
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );

      /* release the constraint */
      SCIP_CALL(SCIPreleaseCons(scip, &cons));

      /* update the new variable's value in the solution */
      SCIP_CALL(SCIPsetSolVal(scip, solution, new_var, new_var_value));
      SCIPinfoMessage(scip, NULL, "\nAdded variable: %s = %.2f\n", new_var_name, new_var_value);

      constraint_idx++;

cleanup_and_continue:
      continue;
   }

   SCIPfclose(gendisj_file);

   /* pass the initial solution to SCIP */
   SCIP_Bool stored;
   SCIP_CALL(SCIPaddSolFree(scip, &solution, &stored));

   /* solve the problem */
   SCIP_CALL(SCIPsolve(scip));

   SCIP_CALL(SCIPfree(&scip));
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
      int                   argc,               /**< number of arguments from the shell */
      char**                argv                /**< array of shell arguments */
      )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
