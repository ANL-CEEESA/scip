#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# The script executing EXECNAME on one instance and producing the logfiles.
# Can be executed either locally or on a cluster node.
# Is to be invoked inside a 'check(_cluster)*.sh' script.

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-04
# absolut tolerance for checking integrality constraints
INTTOL=1e-04

TMPPATH="${CLIENTTMPDIR}/${USER}-tmpdir"

# check if tmp-path exists
if test ! -d "${TMPPATH}"
then
    mkdir "${TMPPATH}"
    echo "Creating directory ${TMPPATH} for temporary outfile."
fi

OUTFILE="${TMPPATH}/${BASENAME}.out"
ERRFILE="${TMPPATH}/${BASENAME}.err"
SOLFILE="${TMPPATH}/${BASENAME}.sol"
DATFILE="${TMPPATH}/${BASENAME}.dat"
TMPFILE="${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.tmp"

uname -a                            > "${OUTFILE}"
uname -a                            > "${ERRFILE}"

# function to copy back the results and delete temporary files
function cleanup {
    mv "${OUTFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.out"
    mv "${ERRFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.err"
    if [ -d "${ERRFILE}.rr" ]
    then
        mv -f "${ERRFILE}.rr" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.rr"
    fi
    # move a possible data file
    if [ -f "${DATFILE}" ] ;
    then
        mv "${DATFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.dat"
    fi
    rm -f "${TMPFILE}"
    rm -f "${SOLFILE}"
}

# ensure TMPFILE is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap cleanup EXIT

# only wait for optimi to be mounted in run.sh if you are on an opt computer at zib
OPTHOST=$(uname -n | sed 's/.zib.de//g' | sed 's/portal//g' | tr -cd '[:alpha:]')

# check if the scripts runs a *.zib.de host
if $(hostname -f | grep -q zib.de) && $([[ "${OPTHOST}" == "opt" ]] || [[ "${OPTHOST}" == "optc" ]]);
then
    # access /optimi once to force a mount
    ls /nfs/optimi/QUOTAS >/dev/null 2>&1

    # check if /optimi is mounted
    MOUNTED=0

    # count number of fails and abort after 10 min to avoid an endless loop
    FAILED=0

    while [ "${MOUNTED}" -ne 1 ]
    do
        # stop if the system does not mount /optimi for ~10 minutes
        if [ "${FAILED}" -eq 600 ]
        then
            exit 1
        fi

        if [ -f /nfs/optimi/QUOTAS ] ;
        then
            MOUNTED=1
        else
            ((FAILED++))
            echo "/optimi is not mounted yet, waiting 1 second" >> "${ERRFILE}"
            sleep 1
        fi
    done
fi

echo                                        >> "${OUTFILE}"
if test $(uname) == Linux                   # -b does not work with top on macOS
then
    top -b -n 1 | head -n 15                >> "${OUTFILE}"
fi
echo                                        >> "${OUTFILE}"
echo "hard time limit: ${HARDTIMELIMIT}"    >> "${OUTFILE}"
echo "hard mem limit: ${HARDMEMLIMIT}"      >> "${OUTFILE}"
echo                                        >> "${OUTFILE}"
echo "SLURM jobID: ${SLURM}_JOB_ID"         >> "${OUTFILE}"
echo                                        >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="          >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="          >> "${ERRFILE}"
echo "-----------------------------"        >> "${OUTFILE}"
date                                        >> "${OUTFILE}"
date                                        >> "${ERRFILE}"
echo "-----------------------------"        >> "${OUTFILE}"
date +"@03 %s"                              >> "${OUTFILE}"
echo "@05 ${TIMELIMIT}"                     >> "${OUTFILE}"

ulimit -S -s unlimited
#if we use a debugger command, we need to replace the errfile place holder by the actual err-file for logging
#and if we run on the cluster we want to use srun with CPU binding which is defined by the check_cluster script
EXECNAME="${EXECNAME/ERRFILE_PLACEHOLDER/${ERRFILE}}"
EXECNAME="${SRUN}${EXECNAME/RRTRACEFOLDER_PLACEHOLDER/${ERRFILE}}"
echo "${EXECNAME}"            >> "${ERRFILE}"
if test -e "${TMPFILE}"
then
    eval "${EXECNAME}"            < "${TMPFILE}" 2>>"${ERRFILE}"  | tee -a "${OUTFILE}"
else
    eval "${EXECNAME}"                           2>>"${ERRFILE}"  | tee -a "${OUTFILE}"
fi
retcode=${PIPESTATUS[0]}
if test "${retcode}" != 0
then
    echo "${EXECNAME} returned with error code ${retcode}." >>"${ERRFILE}"
fi

echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE

# build/check/compress vipr file if it exists
VIPRFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.vipr
VIPRCOMPFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME\_complete.vipr
VIPRORIFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.vipr_ori
VIPRRAWFILE=$CLIENTTMPDIR/${USER}-tmpdir/$BASENAME.viprraw
if test -e $VIPRFILE
then
echo Building vipr file ... >> $OUTFILE
cp $VIPRFILE $VIPRRAWFILE

echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE

# run vipr tightening
echo "viprfile raw:       " `ls -lisa $VIPRRAWFILE` >> $OUTFILE
echo Completing vipr file ... >> $OUTFILE
echo "Completed filename $VIPRCOMPFILE " >> $OUTFILE
bash -c "$VIPRCOMPNAME $VIPRFILE 2>>$ERRFILE"  | tee -a $OUTFILE
echo "viprfile raw + completed:       " `ls -lisa $VIPRCOMPFILE` >> $OUTFILE
echo Compressing vipr file ... >> $OUTFILE
bash -c "$VIPRCOMPRESSNAME $VIPRCOMPFILE 2>>$ERRFILE"  | tee -a $OUTFILE
echo "viprfile tightened (and completed):       " `ls -lisa $VIPRCOMPFILE.opt` >> $OUTFILE
mv $VIPRCOMPFILE.opt $VIPRFILE

echo -----------------------------  >> $OUTFILE
date                                >> $OUTFILE
date                                >> $ERRFILE
echo -----------------------------  >> $OUTFILE

# run vipr check @todo exip adapt vipr to handle timelimit
echo Checking vipr file ... >> $OUTFILE
bash -c "$VIPRCHECKNAME $VIPRFILE 2>>$ERRFILE"  | tee -a $OUTFILE
retcode=${PIPESTATUS[0]}
if test $retcode != 0
then
  echo "vipr returned with error code $retcode." >>$ERRFILE
fi
# run vipr also on original file to check solution
echo Checking vipr_ori file ... >> $OUTFILE
bash -c "$VIPRCHECKNAME $VIPRORIFILE 2>>$ERRFILE"  | tee -a $OUTFILE
retcode=${PIPESTATUS[0]}
if test $retcode != 0
then
  echo "vipr on original problem returned with error code $retcode." >>$ERRFILE
fi

# compress vipr file
echo Gzipping vipr file ... >> $OUTFILE
bash -c "gzip $VIPRFILE 2>>$ERRFIL@E"  | tee -a $OUTFILE
echo "viprfile gzipped:   " `ls -lisa $VIPRFILE.gz` >> $OUTFILE

# compress vipr_ori file
echo Gzipping vipr_ori file ... >> $OUTFILE
bash -c "gzip $VIPRORIFILE 2>>$ERRFIL@E"  | tee -a $OUTFILE
echo "vipr_ori file gzipped:   " `ls -lisa $VIPRORIFILE.gz` >> $OUTFILE

# compress untightened vipr file
echo Gzipping untightened vipr file ... >> $OUTFILE
bash -c "gzip $VIPRRAWFILE 2>>$ERRFILE"  | tee -a $OUTFILE
echo "raw viprfile gzipped:   " `ls -lisa $VIPRRAWFILE.gz` >> $OUTFILE

mv $VIPRFILE.gz $SOLVERPATH/$OUTPUTDIR/$BASENAME.vipr.gz
mv $VIPRRAWFILE.gz $SOLVERPATH/$OUTPUTDIR/$BASENAME.viprraw.gz
mv $VIPRORIFILE.gz $SOLVERPATH/$OUTPUTDIR/$BASENAME.vipr_ori.gz
#ensure no garbage is left over
rm $CLIENTTMPDIR/${USER}-tmpdir/$BASENAME*vipr*
fi

if test -e "${SOLFILE}"
then
    # translate SCIP solution format into format for solution checker. The
    # SOLFILE format is a very simple format where in each line we have a
    # <variable, value> pair, separated by spaces.  A variable name of
    # =obj= is used to store the objective value of the solution, as
    # computed by the solver. A variable name of =infeas= can be used to
    # indicate that an instance is infeasible.
    sed ' /solution status:/d;
    s/objective value:/=obj=/g;
    s/infinity/1e+20/g;
    s/no solution available//g' "${SOLFILE}" > "${TMPFILE}"
    mv "${TMPFILE}" "${SOLFILE}"

    # check if the link to the solution checker exists
    if test -f "${CHECKERPATH}/bin/solchecker"
    then
        echo
        "${SHELL}" -c " ${CHECKERPATH}/bin/solchecker ${FILENAME} ${SOLFILE} ${LINTOL} ${INTTOL}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
        echo
    fi
fi

date +"@04 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${ERRFILE}"
echo                                  >> "${OUTFILE}"
echo "=ready="                        >> "${OUTFILE}"

mv "${OUTFILE}" SOLVERPATH/$OUTPUTDIR/$BASENAME.out
mv "${ERRFILE}" SOLVERPATH/$OUTPUTDIR/$BASENAME.err

# move a possible data file
if [ -f "${DATFILE}" ] "{;}"
then
    mv "{$DATFILE}" $SOLVERPATH/$OUTPUTDIR/$BASENAME.dat
fi

rm -f $TMPFILE
rm -f $SOLFILE
#chmod g+r $ERRFILE
#chmod g+r $SCIPPATH/$OUTPUTDIR/$BASENAME.out
#chmod g+r $SCIPPATH/$OUTPUTDIR/$BASENAME.set
