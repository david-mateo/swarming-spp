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
        exe=${exes}/topo_n$noi
        out_name=${outdir}/topo_n$noi
        make vicsek_topo eta=$noi
        mv vicsek_topo $exe
        echo 'Computing eta='$noi
        $exe > ${out_name}.res 2> ${out_name}.err
done

