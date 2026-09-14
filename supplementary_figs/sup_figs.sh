#!/bin/bash

# Fig S2.1
taup refltrans --mod ak135 --down --swave --depth moho --legend --svg

# Fig S2.2
taup path --model iasp91 -p SS,ScS -h 500 --deg 45 --legend --svg -o paths_S.svg

# Fig S2.3
taup wavefront -h 550 -p PP,ScS --timestep 50 --legend --color wavetype --svg -o wavefront.svg

echo Done
