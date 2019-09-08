#/usr/bin/env sh

rm -f output_dummy.xyz #code will otherwise append
#runs mcvnt with a set of supplied settings
#order is: seed temp dataxyz nstep stridetrj stridelog mcstep outputf 
./mcnvt \
    1238 \
    0.23 \
    input_dummy.xyz \
    10000000 \
    5000 \
    100 \
    0.01 \
    output_dummy.xyz
