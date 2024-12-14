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
#include <ctype.h>
#include <zlib.h>
#include <stdbool.h>

/** Helper function to read an initial solution file */
static
SCIP_RETCODE readSolutionFile(SCIP* scip, const char* filename, SCIP_HASHMAP* solvals, SCIP_SOL* sol)
{
   SCIP_FILE* file;
   SCIP_Bool error;
   int lineno;
   SCIP_Bool unknownvar;

   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         SCIPinfoMessage(scip, NULL, "\nreading initial solution file <%s>\n", filename);
         file = SCIPfopen(filename, "r");

         if( file == NULL )
         {
            SCIPerrorMessage("unable to read initial solution file <%s>\n", filename);

            return SCIP_READERROR;
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "initial solution file <%s> not found\n", filename);

         return SCIP_NOFILE;
      }
   }
   else
      return SCIP_NOFILE;

   error = FALSE;
   lineno = 0;
   unknownvar = FALSE;

   /* read the file */
   while( !SCIPfeof(file) && !(error) )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char varvalstr[SCIP_MAXSTRLEN];
      char varobjstr[SCIP_MAXSTRLEN];
      char format[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real val;
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
      nread = sscanf(buffer, format, varname, varvalstr, varobjstr);
      if( nread < 2 )
      {
         SCIPerrorMessage("invalid input line %d in solution file <%s>: <%s>.\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvar )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> in line %d of solution file <%s>\n",
                  varname, lineno, filename);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvar = TRUE;
         }
         continue;
      }

      /* cast the value */
      if( SCIPstrncasecmp(varvalstr, "inv", 3) == 0 )
         continue;
      else if( SCIPstrncasecmp(varvalstr, "+inf", 4) == 0 || SCIPstrncasecmp(varvalstr, "inf", 3) == 0 )
         val = SCIPinfinity(scip);
      else if( SCIPstrncasecmp(varvalstr, "-inf", 4) == 0 )
         val = -SCIPinfinity(scip);
      else if( SCIPstrncasecmp(varvalstr, "unknown", 7) == 0 )
      {
         val = SCIP_UNKNOWN;
      }
      else
      {
         /* coverity[secure_coding] */
         nread = sscanf(varvalstr, "%lf", &val);
         if( nread != 1 )
         {
            SCIPerrorMessage("invalid solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
                  varvalstr, varname, lineno, filename);
            error = TRUE;
            break;
         }
      }
      SCIP_CALL( SCIPhashmapInsert(solvals, (void*)var, (void*)(size_t)val) );

      SCIP_CALL( SCIPsetSolVal(scip, sol, var, val) );
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
   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n", filename);
         SCIPinfoMessage(scip, NULL, "===========================\n\n");
         SCIP_CALL( SCIPreadParams(scip, filename) );
         SCIP_CALL( SCIPwriteParams(scip, NULL, FALSE, TRUE) );
         SCIPinfoMessage(scip, NULL, "\n");

         return SCIP_OKAY;
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);

         return SCIP_NOFILE;
      }
   }

   return SCIP_NOFILE;
}

/** read problem */
static
SCIP_RETCODE readProblem(
      SCIP*                 scip,               /**< SCIP data structure */
      const char*           filename            /**< input file name */
      )
{
   if ( filename != NULL )
   {
      if ( SCIPfileExists(filename) )
      {
         SCIPinfoMessage(scip, NULL, "\nreading problem <%s> ...\n\n", filename);
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
   char* probname = NULL;
   char* gendisjname = NULL;
   char* settingsname = NULL;
   char* setfilenametowrite = NULL;
   char* soluname = NULL;
   char* logname = NULL;
   int randomseed;
   int permutationseed;
   SCIP_Bool randomseedread;
   SCIP_Bool permutationseedread;
   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool onlyversion;
   SCIP_Real primalreference;
   SCIP_Real dualreference;
   SCIP_Real duallimit;
   const char* dualrefstring;
   const char* primalrefstring;
   const char* duallimitstring;
   int i;

   quiet = FALSE;
   paramerror = FALSE;
   onlyversion = FALSE;
   randomseedread = FALSE;
   randomseed = 0;
   permutationseedread = FALSE;
   permutationseed = 0;
   primalreference = SCIP_UNKNOWN;
   dualreference = SCIP_UNKNOWN;
   duallimit = SCIP_UNKNOWN;
   primalrefstring = NULL;
   dualrefstring = NULL;
   duallimitstring = NULL;

   for( i = 1; i < argc; ++i )
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            printf("missing log filename after parameter '-l'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = TRUE;
      else if( strcmp(argv[i], "-v") == 0 )
         onlyversion = TRUE;
      else if( strcmp(argv[i], "--version") == 0 )
         onlyversion = TRUE;
      else if( strcmp(argv[i], "-s") == 0 )
      {
         i++;
         if( i < argc )
            settingsname = argv[i];
         else
         {
            printf("missing settings filename after parameter '-s'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-f") == 0 )
      {
         i++;
         if( i < argc )
            probname = argv[i];
         else
         {
            printf("missing problem filename after parameter '-f'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-g") == 0 )
      {
         i++;
         if( i < argc )
            gendisjname = argv[i];
         else
         {
            printf("missing general disjunctions filename after parameter '-g'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-x") == 0 )
      {
         i++;
         if( i < argc )
            soluname = argv[i];
         else
         {
            printf("missing solution filename after parameter '-x'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-r") == 0 )
      {
         /*read a random seed from the command line */
         i++;
         if( i < argc && isdigit(argv[i][0]) )
         {
            randomseed = atoi(argv[i]);
            randomseedread = TRUE;
         }
         else
         {
            printf("Random seed parameter '-r' followed by something that is not an integer\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-p") == 0 )
      {
         /*read a permutation seed from the command line */
         i++;
         if( i < argc && isdigit(argv[i][0]) )
         {
            permutationseed = atoi(argv[i]);
            permutationseedread = TRUE;
         }
         else
         {
            printf("Permutation seed parameter '-p' followed by something that is not an integer\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-d") == 0 )
      {
         /*read the dual bound limit from the command line */
         i++;
         if( i < argc )
         {
            duallimitstring = argv[i];
         }
         else
         {
            printf("Dual limit parameter '-d' used incorrectly\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-o") == 0 )
      {
         if( i >= argc - 2 )
         {
            printf("wrong usage of reference objective parameter '-o': -o <primref> <dualref>\n");
            paramerror = TRUE;
         }
         else
         {
            /* do not parse the strings directly, the settings could still influence the value of +-infinity */
            primalrefstring = argv[i + 1];
            dualrefstring = argv[i + 2];
         }
         i += 2;
      }
      else if( strcmp(argv[i], "-w") == 0 )
      {
         i++;
         if( i < argc )
            setfilenametowrite = argv[i];
         else
         {
            printf("missing settings filename for writing after parameter '-w'\n");
            paramerror = TRUE;
         }
      }
      else
      {
         printf("invalid parameter <%s>\n", argv[i]);
         paramerror = TRUE;
      }
   }

   if( !paramerror )
   {
      SCIP* scip = NULL;
      SCIP_HASHMAP* initsolvals = NULL;
      SCIP_SOL* initsol = NULL;
      SCIP_FILE* gendisjfile = NULL;

      /* initialize SCIP */
      SCIP_CALL( SCIPcreate(&scip) );

      /* include default SCIP plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

      /***********************************
       * create log file message handler *
       ***********************************/

      if( quiet )
      {
         SCIPsetMessagehdlrQuiet(scip, quiet);
      }

      if( logname != NULL )
      {
         SCIPsetMessagehdlrLogfile(scip, logname);
      }

      /***********************************
       * Version and library information *
       ***********************************/

      SCIPprintVersion(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPprintExternalCodes(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      if( onlyversion )
      {
         SCIPprintBuildOptions(scip, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
         return SCIP_OKAY;
      }

      /*****************
       * Load settings *
       *****************/

      if( settingsname != NULL )
      {
         SCIP_CALL( readParams(scip, settingsname) );
      }

      /******************
       * Write settings *
       ******************/

      if( setfilenametowrite != NULL )
      {
         SCIP_RETCODE retcode;

         retcode =  SCIPwriteParams(scip, setfilenametowrite, TRUE, FALSE);

         if( retcode == SCIP_FILECREATEERROR )
         {
            SCIPdialogMessage(scip, NULL, "error creating file  <%s>\n", setfilenametowrite);
         }
         else
         {
            SCIP_CALL( retcode );
            SCIPdialogMessage(scip, NULL, "saved parameter file <%s>\n", setfilenametowrite);
         }
      }

      /************************************
       * Change random seed, if specified *
       ***********************************/
      if( randomseedread )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "randomization/randomseedshift", randomseed) );
      }

      /******************************************
       * Change permutation seed, if specified *
       *****************************************/
      if( permutationseedread )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "randomization/permutationseed", permutationseed) );
      }

      /*********************
       * Solve the problem *
       ********************/

      if( probname != NULL )
      {
         SCIP_Bool validatesolve = FALSE;

         if( primalrefstring != NULL && dualrefstring != NULL )
         {
            char *endptr;
            if( ! SCIPparseReal(scip, primalrefstring, &primalreference, &endptr) ||
                  ! SCIPparseReal(scip, dualrefstring, &dualreference, &endptr) )
            {
               printf("error parsing primal and dual reference values for validation: %s %s\n", primalrefstring, dualrefstring);
               return SCIP_ERROR;
            }
            else
               validatesolve = TRUE;
         }

         /* read instance file */
         SCIP_CALL( readProblem(scip, probname) );

         /* read the solution file and create a SCIP_SOL */
         if( soluname != NULL )
         {
            /* create a hashmap for storing the initial solution */
            SCIP_CALL( SCIPhashmapCreate(&initsolvals, SCIPblkmem(scip), SCIPgetNVars(scip)) );

            /* create an empty solution */
            SCIP_CALL( SCIPcreateSol(scip, &initsol, NULL) );

            /* read the initial solution file */
            SCIP_CALL( readSolutionFile(scip, soluname, initsolvals, initsol) );
         }

         /* open the general disjunctions file */
         if( gendisjname != NULL )
         {
            if ( SCIPfileExists(gendisjname) )
            {
               SCIPinfoMessage(scip, NULL, "\nopening general disjunctions file <%s> ...\n", gendisjname);
               gendisjfile = SCIPfopen(gendisjname, "r");

               if( gendisjfile == NULL )
               {
                  SCIPerrorMessage("unable to read general disjunctions file <%s>\n", gendisjname);
                  SCIP_CALL( SCIPfree(&scip) );
                  BMScheckEmptyMemory();

                  return SCIP_READERROR;
               }
            }
            else
            {
               SCIPinfoMessage(scip, NULL, "general disjunctions file <%s> not found.\n", gendisjname);

               return SCIP_NOFILE;
            }
         }

         /* read the general disjunctions file, create and add general disjunction constraints and variables, and update
          * initsol and initsolvls
          */
         if( gendisjfile != NULL )
         {
            SCIP_Bool error = FALSE;
            int consid = 0;

            /* read the gendisj file */
            SCIPinfoMessage(scip, NULL, "reading general disjunctions...\n");
            while( !SCIPfeof(gendisjfile) && !(error) )
            {
               char buffer[SCIP_MAXSTRLEN];

               /* get next line */
               if( SCIPfgets(buffer, (int) sizeof(buffer), gendisjfile) == NULL )
                  break;

               char* token;
               int noldvarsindisj = 0;

               token = strtok(buffer, " \t\n");

               /* parse the number of existing variables in the disjunction */
               if( token != NULL )
                  noldvarsindisj = atoi(token);

               /* parse the variables */
               SCIP_Real disjvarval = 0.0;
               int nvarsindisj = 0;

               /* create an empty constraint */
               char disjconsname[SCIP_MAXSTRLEN];
               SCIPsnprintf(disjconsname, SCIP_MAXSTRLEN, "c%d", SCIPgetNConss(scip) + 1);
               SCIP_CONS* cons;
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, disjconsname, 0, NULL, NULL, 0.0, 0.0) );

               while( ((token = strtok(NULL, " \t\n")) != NULL) )
               {
                  SCIP_Real coeff = 1.0;

                  /* check for a leading '-' sign */
                  if (token[0] == '-')
                  {
                     coeff = -1.0;
                     token++; /* strip the '-' to get the actual variable name */
                  }

                  char* varname = strdup(token);
                  SCIP_VAR* var;

                  var = SCIPfindVar(scip, varname);
                  if (var == NULL)
                  {
                     SCIPinfoMessage(scip, NULL, "variable '%s' not found. Skipping line.\n", varname);
                     free(varname);

                     goto cleanup_and_continue;
                  }

                  /* add the existing variable and coefficient to the new constraint */
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, -1.0 * coeff) );

                  if( initsolvals != NULL )
                  {
                     /* compute the value of the new disjunction variable based on the initial solution */
                     void* varval = SCIPhashmapGetImage(initsolvals, (void*)var);
                     disjvarval += coeff * (varval ? (SCIP_Real)(size_t)varval : 0.0);
                  }

                  nvarsindisj++;
                  free(varname);
               }
               assert(nvarsindisj == noldvarsindisj);

               /* create a new variable `zdisjN` */
               char disjvarname[SCIP_MAXSTRLEN];
               SCIPsnprintf(disjvarname, SCIP_MAXSTRLEN, "zdisj%d", consid + 1);
               SCIP_VAR* disjvar;

               /* create and add a new auxiliary variable representing the general disjunction */
               // TODO: SCIP_VARTYPE_IMPLINT??
               SCIP_CALL( SCIPcreateVarBasic(scip, &disjvar, disjvarname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
               SCIP_CALL( SCIPaddVar(scip, disjvar) );

               /* add the new variable and its coefficient to the new constraint */
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, disjvar, 1.0) );

               /* add the constraint */
               SCIP_CALL( SCIPaddCons(scip, cons) );

               /* print the general disjunction */
               /*
                  SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                  SCIPinfoMessage(scip, NULL, "\n");
                  */

               if( initsol != NULL )
               {
                  /* update the disjunction variable's value in the solution */
                  SCIP_CALL( SCIPsetSolVal(scip, initsol, disjvar, disjvarval) );
               }

               /* release the variable and constraint */
               SCIP_CALL( SCIPreleaseVar(scip, &disjvar) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );

               consid++;

cleanup_and_continue:
               continue;
            }
            SCIPinfoMessage(scip, NULL, "added <%d> general disjunctions.\n", consid);
            SCIPfclose(gendisjfile);
         }

         if( initsol != NULL )
         {
            /* pass the initial solution to SCIP */
            SCIP_Bool stored;
            SCIP_CALL( SCIPaddSolFree(scip, &initsol, &stored) );
         }
         assert(initsol == NULL);

         /* set the dual bound limit, if exists */
         if( duallimitstring != NULL )
         {
            char *endptr;
            if( ! SCIPparseReal(scip, duallimitstring, &duallimit, &endptr) )
            {
               printf("error parsing dual limit value: %s\n", duallimitstring);
               return SCIP_ERROR;
            }
            else
               SCIP_CALL( SCIPsetRealParam(scip, "limits/dual", duallimit) );
         }

         /* solve the problem */
         SCIP_CALL( SCIPsolve(scip) );

         /* display statistics */
         SCIPinfoMessage(scip, NULL, "\nStatistics\n");
         SCIPinfoMessage(scip, NULL, "==========\n\n");
         SCIP_CALL( SCIPprintStatistics(scip, NULL) );

         /* validate the solve */
         if( validatesolve )
         {
            SCIP_CALL( SCIPvalidateSolve(scip, primalreference, dualreference, SCIPfeastol(scip), FALSE, NULL, NULL, NULL) );
         }
      }
      else
      {
         printf("missing problem filename; provide it!\n");
         SCIPinfoMessage(scip, NULL, "\n");
         return SCIP_ERROR;
      }

      if( initsolvals != NULL )
      {
         /* free hash map */
         SCIPhashmapFree(&initsolvals);
      }

      if( scip != NULL )
      {
         /* free SCIP */
         SCIP_CALL( SCIPfree(&scip) );
      }
      BMScheckEmptyMemory();
   }
   else
   {
      printf("\nsyntax: %s [-l <logfile>] [-q] [-s <settings>] [-r <randseed>] [-f <problem>] [-g <gendisjunctions>]\n"
            "  -v, --version : print version and build options\n"
            "  -l <logfile>  : copy output into log file\n"
            "  -q            : suppress screen messages\n"
            "  -s <settings> : load parameter settings (.set) file\n"
            "  -f <problem>  : load and solve problem file\n"
            "  -o <primref> <dualref> : pass primal and dual objective reference values for validation at the end of the solve\n"
            "  -g <gendisjunctions>: load general disjunctions file, create and add general disjunction constraints\n"
            "  -w <setfilenametowrite> : write settings to this file\n"
            "  -d <duallimit> : dual bound limit (limits/dual parameter value)\n"
            "  -r <randseed> : nonnegative integer to be used as random seed. "
            "Has priority over random seed specified through parameter settings (.set) file\n",
         argv[0]);
      printf("\n");
   }

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

   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
