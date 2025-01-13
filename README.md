# tinyfactor
A collection and benchmark of integer factorization methods for at most 128 bit numbers

To make ggnfs-mpqs:  
cd ggnfs_mpqs  
make liblasieve.a  
cd..  
make mpqstest  

To make ggnfs-mpqs3:  
cd ggnfs_mpqs  
make liblasieve.a  
cd..  
make mpqs3test  

To make the other tests:  
make all <USE_AVX2=1> <USE_BMI2=1> <SKYLAKEX=1> <ICELAKE=1> <COMPILER=icc>  


