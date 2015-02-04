% SolvePoiseuille
clear all
close all
clc

myFlow1 = FVM_Flow();
myFlow1.InitializeProblem;

for i = 1:2e3
    i
    myFlow1.CopyNewToOld();
    myFlow1.ClearB();
    if (mod(i,10) == 1)
        myFlow1.ClearA_and_Inverse();
    	myFlow1.CalcJacobian_and_Inverse();
        myFlow1.WriteToFile();
    end
    myFlow1.CalcRHS();
    myFlow1.SolveAndUpdate();
    myFlow1.Visualize();
    cenerVorticity = myFlow1.CenterOmega()
end