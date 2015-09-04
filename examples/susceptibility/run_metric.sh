# Location of the source code
srcdir=$(pwd)/
mail='david_mateo@sutd.edu.sg'
exes=$(pwd)/exes/
mkdir -p $exes
out_base=$(pwd)/logs/
mkdir -p $out_base

rs=`seq 0.2 0.2 3.0`
dens=1

#Location where to put the results
for den in $dens ; do
    outdir=${out_base}metric_den${den}/
    mkdir -p $outdir
    for r in $rs ; do
        exe=vicsek_metric_r$r
        if [ $r == 1.0 ]; then
            sendmail='-m ea -M'$mail
        else
            sendmail=''
        fi
        prog_name=celfcorr_r${r}
        out_name=correlation_r${r}
        make $exe
        mv $exe $exes
        qsub $sendmail -N $prog_name -o $outdir$out_name.res -e $outdir$out_name.err -l nodes=1:ppn=1 -l walltime=30:00:00 << EOF
cd $exes
./$exe
EOF
    done
done
