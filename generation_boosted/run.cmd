# LO
import model loop_sm-no_b_mass
define p = g u c d s u~ c~ d~ s~ b b~
define j = g u c d s u~ c~ d~ s~ b b~
generate p p > t t~
output -f
launch
1
2
4
done
set mxx_min_pdg {6:800}
set nevents 10000
set iseed <RND>
./madspin.top.card
./delphes.card
done
output -f
launch
1
2
4
done
set mxx_min_pdg {6:800}
set nevents 10000
set iseed <RND>
./madspin.antitop.card
./delphes.card
done
exit