% Flow solver based on index 2 DAE and FDM
% Unfortunately it is rank deficient. Therefore, instead of A\b or
% inverse(A)*b we used linsolve

function PoiseuilleSolver()
    clear all
    close all
    clc
    format longG
    global numGrid Re bF L h k numFluidNode
    numGrid = 19;
    numFluidNode = numGrid - 2;
    Re = 1;
    bF = 0;%0.5;
    L = 1;
    h = L / (numGrid - 1);
    k = .01; %delT

    fig1 = figure;
    set(fig1, 'position', [50, 50, 500, 400]);
    fig2 = figure;
    set(fig2, 'position', [600, 80, 500, 400]);
    fig3 = figure;           
    set(fig3, 'position', [1200, 100, 500, 400]);
    A = zeros(2*numGrid^2+2*numGrid, 2*numGrid^2);  % number of nodes: numGrid^2, 
                                                    % 2 unknown per node: pressure and u
                                                    % the remaining components are constraints (boundary)                                                             
    b = zeros(2*numGrid^2 + 2*numGrid,1);
    x = zeros(2*numGrid^2,1);
    UP = zeros(2*numGrid^2);
    U = zeros(numGrid);
    P = zeros(numGrid);
    for tCount = 1:200
        for i=2:numGrid-1
            for j=2:numGrid-1
                % N-S
                 A(MapU(i,j), MapU(i,j)) = (1/k + 2/h^2/Re);
                 A(MapU(i,j), MapU(i+1,j)) = (-1/h^2/Re);
                 A(MapU(i,j), MapU(i-1,j)) = (-1/h^2/Re);
                 A(MapU(i,j), MapP(i,j+1)) = 1/(h);
%                     A(MapU(i,j), MapP(i,j+1)) = 1/(2*h);
                 A(MapU(i,j), MapP(i,j)) = -1/(h);  
%                     A(MapU(i,j), MapP(i,j-1)) = -1/(2*h);
                % Continuity
                 A(MapP(i,j), MapU(i,j+1)) = 1;
%                     A(MapP(i,j), MapU(i,j+1)) = 1;
                 A(MapP(i,j), MapU(i,j)) = -1;
%                     A(MapP(i,j), MapU(i,j-1)) = -1;
                % forcing term
                 b(MapU(i,j)) = bF + 1/k * UP(MapU(i,j));
                % b(MapP(i,j)) is zero and we know it
            end
        end
        %% boudndary
        % periodic x-direction
        for i=2:numGrid-1
             % velocity
             A(MapU(i,1), MapU(i,1)) = 1;
             A(MapU(i,1), MapU(i,numGrid - 1)) = -1;
             b(MapU(i,1)) = 0;
             A(MapU(i,numGrid), MapU(i,numGrid)) = 1;
             A(MapU(i,numGrid), MapU(i,2)) = -1;
             b(MapU(i,numGrid)) = 0;
             
             % pressure
             A(MapP(i,1), MapP(i,2)) = 1;
             A(MapP(i,1), MapP(i,numGrid - 1)) = 0;%-1;
             b(MapP(i,1)) = (numGrid - 2); %bF * (numGrid - 2);
             A(MapP(i,numGrid), MapP(i,numGrid-1)) = 1;
             A(MapP(i,numGrid), MapP(i,2)) = 0;% -1;
             b(MapP(i,numGrid)) = 0;%-bF * (numGrid - 2);
        end
        for j = 2:numGrid-1
             % velocity
             A(MapU(1,j), MapU(1,j)) = 1;
             A(MapU(1,j), MapU(3,j)) = 1;
             b(MapU(1,j)) = 0;
             A(MapU(numGrid,j), MapU(1,j)) = 1;
             A(MapU(numGrid,j), MapU(3,j)) = 1;
             b(MapU(numGrid,j)) = 0;
             
             % pressure
             A(MapP(1,j), MapP(1,j)) = 1;
             A(MapP(1,j), MapP(3,j)) = -1;
             b(MapP(1,j)) = 0;
             A(MapP(numGrid,j), MapP(1,j)) = 1;
             A(MapP(numGrid,j), MapP(3,j)) = -1;
             b(MapP(numGrid,j)) = 0;
        end
        % corners
        % velocity
        A(MapU(1,1), MapU(1,1)) = 1;
        b(MapU(1,1)) = 0;
        A(MapU(1,numGrid), MapU(1,numGrid)) = 1;
        b(MapU(1,numGrid)) = 0;
        A(MapU(numGrid,1), MapU(numGrid,1)) = 1;
        b(MapU(numGrid,1)) = 0;
        A(MapU(numGrid,numGrid), MapU(numGrid,numGrid)) = 1;
        b(MapU(numGrid,numGrid)) = 0;
        % pressure
        A(MapP(1,1), MapP(1,1)) = 1;
        b(MapP(1,1)) = 0;
        A(MapP(1,numGrid), MapP(1,numGrid)) = 1;
        b(MapP(1,numGrid)) = 0;
        A(MapP(numGrid,1), MapP(numGrid,1)) = 1;
        b(MapP(numGrid,1)) = 0;
        A(MapP(numGrid,numGrid), MapP(numGrid,numGrid)) = 1;
        b(MapP(numGrid,numGrid)) = 0;    
        %%  
        fprintf('step %d rankA %d\n', tCount, rank(A));
%         x = A\b;
size(A)
size(b)
% x = pinv(A)*b;
x = linsolve(A,b);
% svd(A(1:numGrid,:));
%         disp(A);
'max diff'
max(max(A*x-b))
        UP = x;
        [U, P] = DeMapUP(UP);
        figure(fig1);
%         for i=1:numGrid
%             U2(i,:) = (i-numGrid/2).^2 * ones(numGrid,1);
%         end
%         pcolor(U2);
        pcolor(U);
        colorbar
        figure(fig2);
        pcolor(P);
        colorbar
        figure(fig3);
        spy(A, 'o');
    end
end
    



function ind = MapU(i,j)
global numGrid
    ind = (j-1) * numGrid + i;
end

function ind = MapP(i,j)
global numGrid
    ind = numGrid^2 + (j-1) * numGrid + i;
end

function out = GhostBoundary(i,j,UP)
global numGrid Re bF L h k 
    out = 0;
    if (i > 1 & i < numGrid & j >= 1 & j < numGrid)
        out = 0;
    end
%     if (j == 1) % minus boundary, pressure
%         p_pm = UP(MapP(i, numGrid)) + bF * (L+h) % should be (L + gridSize)
%         out = out + 1 / (2 * h) * p_pm;
%     end
    if (j == numGrid) % plus boundary, pressure
        p_pp = UP(MapP(i, 1)) - bF * (L+h); 
        out = out - 1 / (2 * h) * p_pp;
    end
    if (i == 1)        % wall
        uGhost = 0;% -UP(MapU(2, j));
        out = out + uGhost / (h^2 * Re);
        % you can also add pressure term herein (i.e. dp/dy = 0 term)
    end
    if (i == numGrid)
        uGhost = 0;%-UP(MapU(numGrid-1, j));
        out = out + uGhost / (h^2 * Re);
    end
end

function [U, P] = DeMapUP(UP)
global numGrid Re bF L h k 
    U = zeros(numGrid);
    P = zeros(numGrid);
    for s = 1:numGrid^2
        j = floor((s-1)/numGrid)+1;
        i = mod(s-1, numGrid)+1;
        U(i, j) = UP(s);
    end
    for l = numGrid^2+1:2*numGrid^2
        s = l - numGrid^2;
        j = floor((s-1)/numGrid)+1;
        i = mod(s-1, numGrid)+1;
        P(i, j) = UP(l);
    end
end