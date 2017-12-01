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
QUEUE=generic # can also be e.g. express

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
for i in 0 #20 30 100
do
    m4 -D ER=$i -D N=2 -D Z=20 -D TF=False config.m4 > topdown_Eres_$i_F.cfg
    python create_fluxes.py topdown_Eres_$i.cfg topdown_Eres_$i_F
    python Statistics.py topdown_Eres_$i_F

    m4 -D ER=$i -D N=2 -D Z=20 -D TF=True config.m4 > topdown_Eres_$i_T.cfg
    python create_fluxes.py topdown_Eres_$i.cfg topdown_Eres_$i_T
    python Statistics.py topdown_Eres_$i_T
done

echo "Directory contents: "
ls
echo ""

# copy output files to output directory
echo "Copying output files."
cp scenario1* ${CALLDIR}/.

EOF
  
echo "created ${SCRIPTNAME}, now submitting to queue."

# change access rights
chmod u=rwx ${SCRIPTNAME}

# submit script
qsub -N CnuB_fluxes_topdown -q ${QUEUE} -V -j oe  -m be -M rasam@nikhef.nl ./$SCRIPTNAME

# remove original script
rm ${SCRIPTNAME}


