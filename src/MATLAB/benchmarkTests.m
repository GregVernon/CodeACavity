clear
close all
clc

cd Native
FCN = @() lidCavity(20,10,1e-5,false);
t(1) = timeit(FCN);
cd ..

cd Built-in
FCN = @() lidCavity(20,10,1e-5,false);
t(2) = timeit(FCN);
cd ..

cd SimpleVectorization
FCN = @() lidCavity(20,10,1e-5,false);
t(3) = timeit(FCN);
cd ..

% cd SimpleVectorization_GPU
% FCN = @() lidCavity(20,10,1e-5,false);
% t(4) = timeit(FCN);
% cd ..

cd AdvancedVectorization_GPU
FCN = @() lidCavity(20,10,1e-5,false);
t(5) = timeit(FCN);
cd ..

cd LinearAlgebra
FCN = @() lidCavity(20,10,1e-5,false);
t(6) = timeit(FCN);
cd ..

cd LinearAlgebra_GPU
FCN = @() lidCavity(20,10,1e-5,false);
t(7) = timeit(FCN);
cd ..

cd LinearAlgebra_CPU_GPU
FCN = @() lidCavity(20,10,1e-5,false);
t(8) = timeit(FCN);
cd ..
