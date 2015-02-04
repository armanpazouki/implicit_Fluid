clear all
clc
format compact
gamma = 0.5;
dT = .001;
bodyForce = [.1; 0; 0];
%dx = .001;
global h;
global m;
global cMax;
global cMin;
h = .0002;
% dx = 1.1 * h;           % dx < 1.1 works good --> results in density = initial density
% dy = 1.15 * h;
% dz = 1.05 * h;
dx = 1.0 * h;
dy = 1.0 * h;
dz = 1.0 * h;
global etha;
etha = .1 * h;
eps = .0001;
rho=1180;
m = rho * dx * dy * dz;

v = [.14; .14; .14];
dV = 0*v;%[.1; .5; .4];

rPM = [rho; 0; .05]; %rho pres mu
dRPM = [0;0;0]; %[0; 15; 0];

nX = 7;
nY = 7;
nZ = 7;
numMesh = [nX,nY,nZ];

% middlepoint = floor([1147 / numMesh / numMesh, mod(1147, numMesh * numMesh) / numMesh, mod(mod(1147, numMesh * numMesh), numMesh)])
posVec = zeros(nX, nY, nZ, 3);
velVec = zeros(nX, nY, nZ, 3);
rPMVec = zeros(nX, nY, nZ, 3);

cMin = [0,0,0];
cMax = cMin + [dx, dy, dz].*numMesh;
hashVal = [];

%% initialize here
% counter = 0;
% for i = 1:nX
%     for j = 1:nY
%         for k = 1:nZ
%             posVec(i,j,k,:) = Perturb([(i-1)*dx; (j-1)*dy; (k-1)*dz], .01 * [dx; dy; dz]);
%             velVec(i,j,k,:) = Perturb(v, dV);
%             rPMVec(i,j,k,:) = Perturb(rPM, dRPM);
%             counter = counter + 1;  
%            % hashVal[counter] = 
%         end
%     end
% end
%% initialize from file
ifd = fopen('../C/Executable/povFiles/fluid0.csv', 'r');
sphData = fscanf(ifd, '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f', [11 inf]);
fclose(ifd);
size(sphData)
sphData = sphData';
counter = 0;
for i = 1:nX
    for j = 1:nY
        for k = 1:nZ
            counter = counter + 1; 
            posVec(i,j,k,:) = sphData(counter, 1:3);
            velVec(i,j,k,:) = sphData(counter, 4:6);
            rPMVec(i,j,k,:) = sphData(counter, 8:10);
           % hashVal[counter] = 
        end
    end
end
%%  Calc Jacobian Matlab
sumA_r = zeros(3);
sumA_v = zeros(3);
sumA_p = zeros(3,1);

sumB_r = zeros(3);
sumB_v = zeros(3);
sumB_p = zeros(3,1);

sumC_r = zeros(1,3);
sumC_v = zeros(1,3);
sumC_p = zeros(1);

sumU = zeros(3,1);
sumU = sumU - bodyForce;
sumL = zeros(1);


JacobianTotal = zeros(4*counter);
ResidualTotal = zeros(4*counter, 1);

countA = 0;
for i = 1:nX
    for j = 1:nY
        for k = 1:nZ
            countA = countA + 1;
            countA4 = 4 * countA - 3;
            r_a = [posVec(i,j,k,1); posVec(i,j,k,2); posVec(i,j,k,3)];
            v_a = [velVec(i,j,k,1); velVec(i,j,k,2); velVec(i,j,k,3)];
            rpm_a = [rPMVec(i,j,k,1); rPMVec(i,j,k,2); rPMVec(i,j,k,3)];
            rho_a = rpm_a(1);
            p_a = rpm_a(2);
            mu_a = rpm_a(3);
                      
            sumA_r = zeros(3);
            sumA_v = zeros(3);
            sumA_p = zeros(3,1);
            sumB_r = zeros(3);
            sumB_v = zeros(3);
            sumB_p = zeros(3,1);
            sumC_r = zeros(1,3);
            sumC_v = zeros(1,3);
            sumC_p = zeros(1);
            
            countB = 0;
            for i2 = 1:nX
                for j2 = 1:nY
                    for k2 = 1:nZ
                        countB = countB + 1;
                        countB4 = 4 * countB - 3;
                        
                        if (countB == countA)
                            continue
                        end
                        
                        r_b = [posVec(i2,j2,k2,1); posVec(i2,j2,k2,2); posVec(i2,j2,k2,3)];
                        v_b = [velVec(i2,j2,k2,1); velVec(i2,j2,k2,2); velVec(i2,j2,k2,3)];
                        rpm_b = [rPMVec(i2,j2,k2,1); rPMVec(i2,j2,k2,2); rPMVec(i2,j2,k2,3)];
                        rho_b = rpm_b(1);
                        p_b = rpm_b(2);
                        mu_b = rpm_b(3);
                        
                        r_ab = Distance(r_a, r_b);
                        v_ab = v_a - v_b;
                        
                        dw = dW(r_ab, h);
                        d2w = d2W(r_ab, h);
                        qr = q_r(r_ab, h);
                        qrr = q_rr(r_ab, h);
                        w = W(r_ab, h);
                        gradW = dw * qr;    
                        Jw = dw*qrr + d2w * qr' * qr;            
                        denomSquare = (norm(r_ab))^2 + etha^2;  
                 
                        %************************ Calc Summations
                        dumSumAr = m * ( p_a / (rho_a^2) + p_b / (rho_b^2) ) * Jw;
                        dumSumBr = m * (mu_a + mu_b) / (0.5*(rho_a + rho_b))^2 ...
                            *1.0/denomSquare * (v_ab * r_ab' * (Jw - 2 * (r_ab' * gradW') / denomSquare * eye(3))...
                            + v_ab * gradW);
                        dumSumCr = m / rho_b * v_ab' * Jw;
                        
                        dumSumAv = zeros(3);
                        dumSumBv = m * (mu_a + mu_b) / (0.5*(rho_a + rho_b))^2 ...
                            *1.0/denomSquare * (r_ab' * gradW') * eye(3);
                        dumSumCv = m / rho_b * gradW;
                        
                        dumSumAp = m / rho_a^2 * gradW';
                        dumSumBp = zeros(3,1);
                        dumSumCp = zeros(1);
                        
                        %************************ calc off-diagonals
                        JacobianTotal(countA4:countA4+3, countB4:countB4+3) = ...
                            [-gamma * dT * (dumSumAv - dumSumBv) - (gamma * dT)^2*(dumSumAr - dumSumBr), (dumSumAp - dumSumBp);
                            -(dumSumCv) - (gamma * dT) * (dumSumCr), 0];
                        %*******************************************                        
                        if (countA == 6 & countB == 13)
                            r_a'
                            r_b'
                            r_ab'
%                             q = norm(r_ab,2) / h
%                             Jw
%                             dw
%                             d2w
%                             qr
%                             'multU2'
%                             m * (mu_a + mu_b) / (0.5*(rho_a + rho_b))^2 *1.0/denomSquare
%                             m
%                             denomSquare
%                             1/denomSquare
%                             'rpm_a'
%                             rpm_a(1)
%                             rpm_a(2)
%                             rpm_a(3)
%                             dumSumAr(1)
%                             'm * ( p_a / (rho_a^2) + p_b / (rho_b^2) )'
%                             m * ( p_a / (rho_a^2) + p_b / (rho_b^2) )
%                             -((dumSumCv) + (gamma * dT) * (dumSumCr))

                        end
                        %*******************************************
                                           
                                                
                        sumA_r = sumA_r + dumSumAr;
                        sumB_r = sumB_r + dumSumBr;
                        sumC_r = sumC_r + dumSumCr;
                        
                        sumA_v = sumA_v + dumSumAv;
                        sumB_v = sumB_v + dumSumBv;
                        sumC_v = sumC_v + dumSumCv;

                        sumA_p = sumA_p + dumSumAp;
                        sumB_p = sumB_p + dumSumBp;
                        sumC_p = sumC_p + dumSumCp;
                        
                        sumU = sumU + U(r_ab, v_ab, rpm_a, rpm_b);
                        sumL = sumL + L(r_ab, v_ab, rpm_a, rpm_b);
                    end
                end
            end
            JacobianTotal(countA4:countA4+3, countA4:countA4+3) = ...
                [eye(3) + gamma * dT * (sumA_v - sumB_v) + (gamma * dT)^2*(sumA_r - sumB_r), (sumA_p - sumB_p);
                (sumC_v) + (gamma * dT) * (sumC_r), 0];
            ResidualTotal(countA4:countA4+3) = [sumU; sumL];
        end
    end
end
%% compare with C code
% note: the data to be loaded should be in sorted array format, not
% original indices
path = '../C/Executable/coo_data/';
coo_row = importdata(strcat(path, 'coo_row.txt'));
coo_col = importdata(strcat(path, 'coo_col.txt'));
coo_val = importdata(strcat(path, 'coo_val.txt'));
% coo_val(79392)
numComp = max(coo_row) + 1;
mySparseMat = sparse(coo_row + 1, coo_col + 1, coo_val,numComp,numComp);
myFullMat = full(mySparseMat);
% spy(JacobianTotal,'o');
% figure
% spy(myFullMat,'o');

% svdDecomposition = [(svd(JacobianTotal))';(svd(myFullMat))']
%% find the max off components
difMat = (myFullMat - JacobianTotal);
% [indMax_I,indMax_J] = find(difMat==max(difMat(:)))
% [indMin_I,indMin_J] = find(difMat==min(difMat(:)))
absDifMat = abs(difMat);
absJacobianTotal = abs(JacobianTotal);

figure
bb = absDifMat > (.00001 * absJacobianTotal);
spy(bb, 'o');

'maxAbs'
maxAbs = max(absDifMat(:))
% difMatRefineryIndex = (absDifMat > .001 * maxAbs);
% difMatRefined = zeros(size(difMat));
% difMatRefined(difMatRefineryIndex) = difMat(difMatRefineryIndex);
% spy(difMatRefined, 'o');
% refinedCompareWithJacobian = difMatRefined > (.00001 * absJacobianTotal);
% figure
% spy(refinedCompareWithJacobian,'*')
% 
% myDifff = zeros(size(difMat));
% myDifff(difMatRefineryIndex) = 1;
% figure
% spy(myDifff,'*')
%% solve using bicgstab
clear delX1 delX2 flag1 flag2
maxit=50;
err=1e-5;
% delXAnalitical = JacobianTotal\ResidualTotal;
[delX1, flag1] = bicgstab(JacobianTotal, ResidualTotal, err, maxit);
flag1
[delX2, flag2] = bicgstab(myFullMat, ResidualTotal, err, maxit);
flag2
'max(JacobianTotal*delX1 - ResidualTotal)'
max(JacobianTotal*delX1 - ResidualTotal)
'max(myFullMat*delX2 - ResidualTotal)'
max(myFullMat*delX2 - ResidualTotal)
max(ResidualTotal - delX2)

% save('jacobianMatlab.txt', 'JacobianTotal');
% save('jacobianC.txt', 'myFullMat');
% save('ResidualTotal.txt', 'ResidualTotal');

%% find max difference
myTolerance = 1e-3;
JacobianTotalRefined = round(JacobianTotal / myTolerance) * myTolerance;
myFullMatRefined = round(myFullMat / myTolerance) * myTolerance;
diffMatRounded = JacobianTotalRefined;
nonZeroIndex = (JacobianTotalRefined ~= 0);
diffMatRounded(nonZeroIndex) = (JacobianTotalRefined(nonZeroIndex) - myFullMatRefined(nonZeroIndex))./JacobianTotalRefined(nonZeroIndex);
myTolerance = 1e-2;
diffMatRounded = round(diffMatRounded / myTolerance) * myTolerance;
bb = max(diffMatRounded);
%% rank of constraint matrix
% myConstraintMat = JacobianTotal(4:4:864,:);
% 'size myConstraintMat'
% size(myConstraintMat)
% 'rank myConstraintMat'
% rank(myConstraintMat)
% spy(myConstraintMat, 'o');
% axis equal
% myTolerance = 1e-1;
% refineMat = round(myConstraintMat / myTolerance) * myTolerance;
% myConstB = unique(refineMat, 'rows');
% figure 
% spy(refineMat, 'o')
%% corner cells
% % (i-1)*(nY*nZ) + (j-1)*(nZ) + (k) 
% myDependencyMatInd = [
%     (1-1)*(nY*nZ) + (1-1)*(nZ) + (1),
%     (1-1)*(nY*nZ) + (1-1)*(nZ) + (nZ),
%     (1-1)*(nY*nZ) + (nY-1)*(nZ) + (1),
%     (1-1)*(nY*nZ) + (nY-1)*(nZ) + (nZ),
%     (nX-1)*(nY*nZ) + (1-1)*(nZ) + (1),
%     (nX-1)*(nY*nZ) + (1-1)*(nZ) + (nZ),
%     (nX-1)*(nY*nZ) + (nY-1)*(nZ) + (1),
%     (nX-1)*(nY*nZ) + (nY-1)*(nZ) + (nZ)];
% myDependencyMat = JacobianTotal(myDependencyMatInd,:);
%% load A, b, and res
path = '../C/Executable/coo_data/';
A = importdata(strcat(path, 'A.dat'));
b = importdata(strcat(path, 'b.dat'));
x = importdata(strcat(path, 'xUnknown.dat'));
x2 = A\b;
max(abs(x-x2))
                


         