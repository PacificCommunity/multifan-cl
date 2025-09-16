# Compilation scripts
by Yukio Takeuchi

As of January 2019, these scripts allow to compile OpenBlas, QD, ADMB and MFCL.

**Requirements**

OpenBLAS-0.2.20.tar.gz
qh.h
qd-2.3.22.tar.gz

An internet connection and `git` are also required.

**Build**

```./build_openblas4mfcl.sh
./build_qd4mfcl.sh
git clone https://github.com/PacificCommunity/ofp-sam-admb.git
cd ofp-sam-admb
git checkout e06e6fe54718bf4520b8d10de2b375da8efd2329
cd ..
git clone https://github.com/PacificCommunity/ofp-sam-mfcl.git
./buildMFCL.sh
./buildMFCLd.sh```

**Output**
Produced executables are located in:

`./devvsn11/<os>/<date><qd_opt>_QD_<openblas_opt>_OpenBLAS/mfclsdbg64`

Where:
`<os>` is the name of the current Operating System (ie: Ubuntu16.04)
`<date>` is the date of compilation (ie: 20190115)
`<qd_opt>` are QD flags (ie: O0)
`<openblas_opt>` are OpenBLAS flags (ie: default_dynamic)

For example:
`./devvsn11/Ubuntu16.04/20190115_O0_QD_default_dynamic_OpenBLAS/mfclsdbg64`
