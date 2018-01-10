#!/bin/bash

# run a NikhefORCA job on the stoomboot cluster

callcommand="$0 $*"

#PROG=$1
#echo "prog = ${PROG}"

#TARGETDIR=$2
#echo "output directory = '${TARGETDIR}'"

#ARGS=$3
#echo "args = '${ARGS}'"

echo ""

CALLDIR=`pwd` # directory from which this script was called
QUEUE=generic # can also be e.g. express, generic, long

# create script
SCRIPTNAME=tmp_script_$$.sh
cat >> $SCRIPTNAME << EOF
#!/bin/bash -l

export PATH=/project/datagrid/anaconda/bin:$PATH
##source activate py27root5 #doesnt seem to work
source activate py27root5
export PATH=/project/datagrid/anaconda/bin:$PATH

# print working directory
cd \$TMPDIR
echo "Working directory = " 
pwd
echo ""

echo "Directory contents: "
ls
echo ""

# copy necessary input files
echo "Copying input files"
cp $CALLDIR/src/NeutrinoFlux.py . # the program itself
cp $CALLDIR/create_fluxes.py . # the program itself
cp $CALLDIR/src/Statistics.py . # the program itself
cp $CALLDIR/config.m4 .  # config files

echo "Directory contents: "
ls
echo ""

# run the program
echo "Running"

m4 -D ER=0 -D N=6 -D Z=20 -D TF=False config.m4 > ERES30_F_OPT.cfg
python create_fluxes.py ERES30_F_OPT.cfg ERES30_F_OPT
python Statistics.py ERES30_F_OPT_H0H1_m0.1_zmax20_n6_alpha2_ZdecayTrue.txt

echo "Directory contents: "
ls
echo ""

# copy output files to output directory
echo "Copying output files."
cp *ERES30_F_OPT* ${CALLDIR}/.

EOF

echo "created ${SCRIPTNAME}, now submitting to queue."

# change access rights
chmod u=rwx ${SCRIPTNAME}

# submit script
qsub -N CnuB_fluxes_ERES30_F_OPT -q ${QUEUE} -V -j oe  -m be -M rasam@nikhef.nl ./$SCRIPTNAME

# remove original script
rm ${SCRIPTNAME}


