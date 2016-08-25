#! /bin/bash
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /netapp/home/mmravic/CHAMP/a5/CHARMM_GBSW/logs                       #-- output directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=4G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=4G,scratch=4G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=18:00:00                #-- runtime limit (see above; this requests 24 hours)
#$ -t 1-163                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)

#module load charmm
#export LC_CTYPE=en_US.UTF-8
#export LC_ALL=en_US.UTF-8
### ./estimate_GBSW.sh ../model_list.txt 
inx=1

#taskID=3
taskID=$SGE_TASK_ID


ids=`cat $1`
CWD=`pwd`
for id in $ids; do
    #echo $id $inx

    if [ $inx -eq $taskID ]
        then

    echo $id $inx "!!!!!"
   # break

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
    cp $id $outdir
    perl /netapp/home/mmravic/CHAMP/bin/convpdb.pl -segnames  $id > ${outdir}$seg
    sed -i 's/HSE/HSD/' ${outdir}$seg
    perl /netapp/home/mmravic/CHAMP/bin/convpdb.pl -readseg  -chain A ${outdir}$seg > ${outdir}/proa/$outa
    perl /netapp/home/mmravic/CHAMP/bin/convpdb.pl -readseg  -chain X ${outdir}$seg > ${outdir}/prox/$outx
    cd $outdir/dimer

    cp /netapp/home/mmravic/CHAMP/a5/CHARMM_GBSW/input.charmm .
    sed -i "s|XXA|\.\./proa/$outa|" input.charmm
    sed -i "s|XXX|\.\./prox/$outx|" input.charmm

    charmm < input.charmm > out.log
    cd $CWD
    break
    fi

    inx=$(( inx + 1))

done
exit
#cd ../proa/
#cp ../../../input_PROA.charmm .
#sed -i "s/XXA/$outa/" input_PROA.charmm

#charmm < input_PROA.charmm > out.log

#cd ../prob/
#cp ../../../input_PROB.charmm .
#sed -i "s/XXB/$outb/" input_PROB.charmm

#charmm < input_PROB.charmm > out.log
