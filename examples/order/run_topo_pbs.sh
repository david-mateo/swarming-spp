mail='itsame@mario.com'
# Location of the source code
srcdir=$(pwd)/
#Location where to put the executables
exes=$(pwd)/exes/
mkdir -p $exes
#Location where to put the results
outdir=$(pwd)/logs/
mkdir -p $outdir

nois=`seq 0.05 0.05 1.0`

for noi in $nois ; do
        exe=topo_n$noi
        if [ $noi == 0.55 ]; then
            sendmail='-m ea -M'$mail
        else
            sendmail=''
        fi
        prog_name=topo_$noi
        out_name=topo_n$noi
        make vicsek_topo eta=$noi
        mv vicsek_topo $exes/$exe
        qsub $sendmail -N $prog_name -o $outdir$out_name.res -e $outdir$out_name.err -l nodes=1:ppn=1 -l walltime=4:00:00 << EOF
cd $exes
./$exe
EOF
done
