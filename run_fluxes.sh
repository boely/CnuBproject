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
QUEUE=express # can also be e.g. generic

# create script
SCRIPTNAME=tmp_script_$$.sh
cat >> $SCRIPTNAME << EOF
#!/bin/bash -l

export PATH=/project/datagrid/anaconda/bin:$PATH
##source activate py27root5 #doesnt seem to work
source /project/anatares/root_v5.34.03/bin/thisroot.sh

# register starting time of script
date1=\$(date +"%s")

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
#cp $CALLDIR/${PROG} . # the program itself
cp $CALLDIR/src/NeutrinoFlux.py . # the program itself
cp $CALLDIR/src/Statistics.py . # the program itself
cp $CALLDIR/create_fluxes.py . # the program itself
cp $CALLDIR/*.cfg .  # config files
cp $CALLDIR/config.m4 .  # config files

echo "Directory contents: "
ls
echo ""

# run the program
echo "Running"
for i in 1 5 10 20
do
    m4 -D ZI=$i config.m4 > scenario1_Eres_$i.cfg
    python create_fluxes.py scenario1_Eres_$i.cfg scenario1_Eres_$i
    python Statistics.py scenario1_Eres_$i
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
qsub -N CnuB_fluxes -q ${QUEUE} -V -j oe  ./$SCRIPTNAME

# remove original script
rm ${SCRIPTNAME}


