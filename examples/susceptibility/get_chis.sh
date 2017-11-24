# Process correlations to obtain the susceptibility at each sampling time.
# This scripts takes the result from running vicsek_*.cpp as stdin. To get
# the susceptibilites e.g. in Vicsek topologic, run
#       ./vicsek_topo_k10 | ./get_chis.sh > susceptibilites.dat
# (or, to keep the correlations, do
#       ./vicsek_topo_k10 > correla_histo.dat
#       cat correla_histo.dat | ./get_chis.sh > susceptibilites.dat
# )
#
# Note: The factor of 2 is there because libspp computes the sum over i>j instead of i!=j.
awk '
/# Number of agents/ {nag=$5}
/#Iteration/ {q=0. ; chi=0. ; iter=$2 ; output=1 }
!/#/ {
    if($1==""){
        if(output) {output=0 ; print iter, 2*chi/nag}
    }else{
        q+= $2 ; if(q>chi) chi=q ;
    }
}' 
