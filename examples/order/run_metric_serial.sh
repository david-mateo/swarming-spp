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
        exe=${exes}/metric_n$noi
        out_name=${outdir}/metric_n$noi
        make vicsek_metric eta=$noi
        mv vicsek_metric $exe
        echo 'Computing eta='$noi
        $exe > ${out_name}.res 2> ${out_name}.err
done

