#! /bin/bash

#module load charmm
#export LC_CTYPE=en_US.UTF-8
#export LC_ALL=en_US.UTF-8

ids=`cat $1`
CWD=`pwd`
for id in $ids; do
    echo $id
    idl=`echo $id | tr '[:upper:]' '[:lower:]'|rev|cut -d '/' -f 1 |rev`
    seg=`echo $idl | sed  's/\./_seg\./'`
    outa=`echo $idl | sed 's/\./_proa\./'`
    outx=`echo $idl | sed 's/\./_prox\./'`
    outdir=`echo $idl | cut -d '.' -f 1`
    outdir=gbsw_dynamics/$outdir/
    mkdir $outdir
    mkdir ${outdir}/proa
    mkdir ${outdir}/prox
    mkdir ${outdir}/dimer

    convpdb.pl -segnames  $id > ${outdir}$seg
    sed -i 's/HSE/HSD/' ${outdir}$seg
    convpdb.pl -readseg  -chain A ${outdir}$seg > ${outdir}/proa/$outa
    convpdb.pl -readseg  -chain X ${outdir}$seg > ${outdir}/prox/$outx
    cd $outdir/dimer

    cp ../../../input.charmm .
    sed -i "s|XXA|\.\./proa/$outa|" input.charmm
    sed -i "s|XXX|\.\./prox/$outx|" input.charmm

    charmm < input.charmm > out.log
    cd $CWD
done
exit
cd ../proa/
cp ../../../input_PROA.charmm .
sed -i "s/XXA/$outa/" input_PROA.charmm

charmm < input_PROA.charmm > out.log

cd ../prob/
cp ../../../input_PROB.charmm .
sed -i "s/XXB/$outb/" input_PROB.charmm

charmm < input_PROB.charmm > out.log
