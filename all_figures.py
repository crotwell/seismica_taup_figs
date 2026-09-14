#!/usr/bin/env python
import subprocess
import sys
from contextlib import chdir
import taup

sys.path.append("fig1")
from fig1 import fig1
sys.path.append("fig2")
from fig2 import fig2
sys.path.append("fig3")
from fig3 import fig3
sys.path.append("fig4")
from fig4 import fig4
sys.path.append("fig6")
from fig6 import fig6



def main():
    taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
    taup_path="~/Code/seis/TauP/build/install/TauP/bin/taup"

    with taup.TauPServer(taup_path=taup_path) as taupserver:
        print("#### Figure 1 ####")
        with chdir("fig1"):
            fig1(taupserver)
        print("#### Figure 2 ####")
        with chdir("fig2"):
            fig2(taupserver)
        print("#### Figure 3 ####")
        with chdir("fig3"):
            fig3(taupserver)
        print("#### Figure 4 ####")
        with chdir("fig4"):
            fig4(taupserver)
        print("#### Figure 6 ####")
        with chdir("fig6"):
            fig6(taupserver)

        # supplementary figures:
        with chdir("supplementary_figs"):
            subprocess.run("./sup_figs.sh", shell=True)
        print("#### Done ####")

if __name__ == '__main__':
    main()
