### PrivHZChain: Privacy-Preserving Monitoring of Controlled Substance Supply Chains Using Hierarchical Blockchain and ZKP

Monitoring the Supply Chain of Controlled Substances with Privacy-Preserving Hierarchical Blockchain

#### Building the Program

This code uses MITSUNARI Shigeo's MCL for elliptic curve gruop operations. You can clone it by following commands:


###### install [herumi/mcl](https://github.com/herumi/mcl)
```
git clone https://github.com/herumi/mcl;
cd mcl;
make -j4;

mkdir build
cd build
cmake ..
make
```

Make sure that the `MCL` library path is well recognized by the program. You can add the path by setting the `LD_LIBRARY_PATH` as the following command line:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mcl/lib
```

After cloning the `MCL` library, you can compile the code by the following command line. Note that `test.cpp` contains all the test codes.

```
g++ -o test.exe MonitoringCSA/test.cpp -O2 -Iinclude -I../mcl/include -lmcl -L../mcl/lib
```

#### Reproducing Experimental Data

The `test.cpp` contains the following test codes:

- Core Functions
    - ZKMP_test: Membership proof generation and verification.
    - ZKNMP_test: Non-membership proof generation and verification.
    - ZKSP_test: Set split proof generation and verification.
    - ZKIPP_test: Zero-knowledge inner pairing proof generation and verification.

- Applications 
    - OwnPf_test: Ownership proof generation and verification.
    - NonOwnPf_test: Non-ownership proof generation and verification.
    - SetupTest: Measuing setup time for preparing the parameter of accumulator (by the global manager)
    - TxGenTest: (Entry/Transfer/Exit) Transaction proof generation.
    - TxVrfyTest: (Entry/Transfer/Exit) Transaction verification.
    - TxAggTest: Summary Transaction proof generation and verification.
    
#### Enjoy!