# Location of the source code
srcdir=$(pwd)/
mail='david_mateo@sutd.edu.sg'
exes=$(pwd)/exes/
mkdir -p $exes
out_base=$(pwd)/logs/
mkdir -p $out_base

ks=`seq 0 2 60`
dens=1

#Location where to put the results
for den in $dens ; do
    outdir=${out_base}topo_den${den}/
    mkdir -p $outdir
    for k in $ks ; do
        if [ $k -lt 10 ] ; then kk=0$k ; else kk=$k ; fi
        exe=vicsek_topo_k$k
        if [ $k == 16 ]; then
            sendmail='-m ea -M'$mail
        else
            sendmail=''
        fi
        prog_name=celfcorr_k${kk}
        out_name=correlation_k${kk}
        make $exe
        mv $exe $exes
        qsub $sendmail -N $prog_name -o $outdir$out_name.res -e $outdir$out_name.err -l nodes=1:ppn=1 -l walltime=30:00:00 << EOF
cd $exes
./$exe
EOF
    done
done
