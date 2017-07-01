% Flow solver based on index 2 DAE, FVM, and staggered grid

%the code has issue somewhere. increase the domain size and you will see
%that. the boundaries are not working well

classdef FVM_Flow < handle
    properties
        ni;
        nj;
        dx;
        dy;
        dt;
        lx;
        ly;
        u;
        v;
        p;
        psi;
        u_old;
        v_old;
       
        
        A;
        invA;
        b;
        fp;
        gu;
        
        rho;
        nue;
        fx;
        fy;
        
        % boundaries
        uT; % top
        uB; % bottom
        uL; % left
        uR; % right
        vT; % top
        vB; % bottom
        vL; % left
        vR; % right
        
        fgP;
        fgU;
        
        fileCounter;
        
        
    end %properties
    methods
        function obj = FVM_FLow()
        end
        
        function InitializeProblem(obj)
            obj.ni = 41; % number of nodes along x-direction
            obj.nj = 41; % number of nodes along y-direction
            obj.lx = 1;
            obj.ly = 1;
            
            obj.dx = obj.lx / (obj.ni - 1);
            obj.dy = obj.ly / (obj.nj - 1);
            obj.dt = .05;
            obj.u = zeros(obj.ni, obj.nj - 1);
            obj.v = zeros(obj.ni - 1, obj.nj);
            obj.p = zeros(obj.ni-1, obj.nj-1);
            obj.psi = zeros(obj.ni-1, obj.nj-1);
            obj.u_old = obj.u;
            obj.v_old = obj.v;

            numUVPcomps = obj.ni * (obj.nj - 1) + (obj.ni - 1) * obj.nj + (obj.ni - 1) * (obj.nj - 1);
            obj.A = zeros(numUVPcomps + 1, numUVPcomps);
            obj.b = zeros(numUVPcomps + 1, 1);
            obj.fp = zeros(obj.ni * (obj.nj - 1) + (obj.ni - 1) * obj.nj, (obj.ni - 1) * (obj.nj - 1));
            obj.gu = zeros((obj.ni - 1) * (obj.nj - 1), obj.ni * (obj.nj - 1) + (obj.ni - 1) * obj.nj);

            obj.rho = 1;
            obj.nue = .01;
            obj.fx = 0;
            obj.fy = 0;%-9.81;

            % boundaries
            obj.uT = 1; % top
            obj.uB = 0; % bottom
            obj.uL = 0; % left
            obj.uR = 0; % right
            obj.vT = 0; % top
            obj.vB = 0; % bottom
            obj.vL = 0; % left
            obj.vR = 0; % right
            
            obj.fgP = figure;
            obj.fgU = figure;
            
            obj.fileCounter = 0;
        end
        
        %% Clear A
        function ClearA_and_Inverse(obj)
            obj.A = zeros(size(obj.A));
            obj.invA = zeros(size(obj.A));
            obj.fp = zeros(size(obj.fp));
            obj.gu = zeros(size(obj.gu));
        end
        
        %% Clear B
        function ClearB(obj)
            obj.b = zeros(size(obj.A));
        end
            
        %% RHS for : psi_q * (-del_q) = psi
        function rhsU(obj, i, j)
            row = obj.mapU(i,j);
            obj.b(row) = obj.b(row) + ...
                (obj.U(i,j) - obj.u_old(i,j)) / obj.dt + ... 
                .25 /obj.dx * ((obj.U(i+1,j) + obj.U(i,j))^2 - (obj.U(i,j) + obj.U(i-1,j))^2) + ...
                .25 /obj.dy * ((obj.U(i,j+1) + obj.U(i,j)) * (obj.V(i-1,j+1) + obj.V(i,j+1)) - ...
                               (obj.U(i,j-1) + obj.U(i,j)) * (obj.V(i-1,j) + obj.V(i,j))) + ...
                (1.0 / obj.rho) * 1.0 / obj.dx *(obj.P(i,j) - obj.P(i-1,j)) + ...
                (-obj.nue) * (...
                    1.0 / (obj.dx)^2 *(obj.U(i+1,j) - 2 * obj.U(i,j) + obj.U(i-1,j)) + ...
                    1.0 / (obj.dy)^2 *(obj.U(i,j+1) - 2 * obj.U(i,j) + obj.U(i,j-1)) ...
                ) + ...
                (-obj.fx);
        end
        function rhsV(obj, i, j)
            row = obj.mapV(i,j); 
            obj.b(row) = obj.b(row) + ...  
                (obj.V(i,j) - obj.v_old(i,j)) / obj.dt + ... 
                .25 /obj.dx * ((obj.U(i+1,j) + obj.U(i+1,j-1)) * (obj.V(i+1,j) + obj.V(i,j)) - ...
                               (obj.U(i,j) + obj.U(i,j-1)) * (obj.V(i,j) + obj.V(i-1,j))) + ...
                .25 /obj.dy * ((obj.V(i,j+1) + obj.V(i,j))^2 - (obj.V(i,j-1) + obj.V(i,j))^2) + ...
                (1.0 / obj.rho) * 1.0 / obj.dy *(obj.P(i,j) - obj.P(i,j-1)) + ...
                (-obj.nue) * (...
                    1.0 / (obj.dx)^2 *(obj.V(i+1,j) - 2 * obj.V(i,j) + obj.V(i-1,j)) + ...
                    1.0 / (obj.dy)^2 *(obj.V(i,j+1) - 2 * obj.V(i,j) + obj.V(i,j-1)) ...
                ) + ...
                (-obj.fy);
        end
        function rhsContinuity(obj, i, j)
            row = obj.mapP(i,j); 
            obj.b(row) = obj.b(row) + ... 
                1.0 / obj.dx *(obj.U(i+1,j) - obj.U(i,j)) + ...
                1.0 / obj.dy *(obj.V(i,j+1) - obj.V(i,j));
        end        
        %% Derivatives
        function uu_x(obj, i, j) % appears in u equation, i.e. obj.mapU(i,j) row of A
%             inout = .25 /obj.dx * (...
%                 (obj.U(i+1,j) + obj.U(i,j))^2 - (obj.U(i,j) + obj.U(i-1,j))^2 ...
%                 );
              row = obj.mapU(i,j); 
              obj.A(row, obj.mapU(i+1,j)) = obj.A(row, obj.mapU(i+1,j)) + obj.multIJ_uEq(i+1,j) * 2.0 * .25 / obj.dx * (obj.U(i+1,j) + obj.U(i,j));
              obj.A(row, obj.mapU(i-1,j)) = obj.A(row, obj.mapU(i-1,j)) - obj.multIJ_uEq(i-1,j) * 2.0 * .25 / obj.dx * (obj.U(i,j) + obj.U(i-1,j));
        end
        function uv_y(obj, i, j) % appears in u equation, i.e. obj.mapU(i,j) row of A
%             inout = .25 /obj.dy * (...
%                 (obj.U(i,j+1) + obj.U(i,j)) * (obj.V(i-1,j+1) + obj.V(i,j+1))   ...
%                 -(obj.U(i,j-1) + obj.U(i,j)) * (obj.V(i-1,j) + obj.V(i,j)) ...
%                 );
              row = obj.mapU(i,j); 
              obj.A(row, obj.mapU(i,j+1)) = obj.A(row, obj.mapU(i,j+1)) + obj.multIJ_uEq(i,j+1) * .25 / obj.dy * (obj.V(i-1,j+1) + obj.V(i,j+1));              
              obj.A(row, obj.mapU(i,j)) = obj.A(row, obj.mapU(i,j)) + obj.multIJ_uEq(i,j) * .25 / obj.dy * ( (obj.V(i-1,j+1) + obj.V(i,j+1)) - (obj.V(i-1,j) + obj.V(i,j)) );
              obj.A(row, obj.mapU(i,j-1)) = obj.A(row, obj.mapU(i,j-1)) - obj.multIJ_uEq(i,j-1) * .25 / obj.dy * (obj.V(i-1,j) + obj.V(i,j)); 
              
              obj.A(row, obj.mapV(i-1,j+1)) = obj.A(row, obj.mapV(i-1,j+1)) + obj.multIJ_uEq(i-1,j+1) * .25 / obj.dy * (obj.U(i,j+1) + obj.U(i,j));
              obj.A(row, obj.mapV(i,j+1)) = obj.A(row, obj.mapV(i,j+1)) + obj.multIJ_uEq(i,j+1) * .25 / obj.dy * (obj.U(i,j+1) + obj.U(i,j));
              obj.A(row, obj.mapV(i-1,j)) = obj.A(row, obj.mapV(i-1,j)) - obj.multIJ_uEq(i-1,j) * .25 / obj.dy * (obj.U(i,j-1) + obj.U(i,j));
              obj.A(row, obj.mapV(i,j)) = obj.A(row, obj.mapV(i,j)) - obj.multIJ_uEq(i,j) * .25 / obj.dy * (obj.U(i,j-1) + obj.U(i,j));        
        end
        function vu_x(obj, i, j) % appears in v equation, i.e. obj.mapV(i,j) row of A
%             inout = .25 /obj.dx * (...
%                 (obj.U(i+1,j) + obj.U(i+1,j-1)) * (obj.V(i+1,j) + obj.V(i,j)) ...
%                 -(obj.U(i,j) + obj.U(i,j-1)) * (obj.V(i,j) + obj.V(i-1,j)) ...
%                 );
              row = obj.mapV(i,j); 
              obj.A(row, obj.mapU(i+1,j)) = obj.A(row, obj.mapU(i+1,j)) + obj.multIJ_vEq(i+1,j) * .25 / obj.dx * (obj.V(i+1,j) + obj.V(i,j));
              obj.A(row, obj.mapU(i+1,j-1)) = obj.A(row, obj.mapU(i+1,j-1)) + obj.multIJ_vEq(i+1,j-1) * .25 / obj.dx * (obj.V(i+1,j) + obj.V(i,j));
              obj.A(row, obj.mapU(i,j)) = obj.A(row, obj.mapU(i,j)) - obj.multIJ_vEq(i,j) * .25 / obj.dx * (obj.V(i,j) + obj.V(i-1,j));
              obj.A(row, obj.mapU(i,j-1)) = obj.A(row, obj.mapU(i,j-1)) - obj.multIJ_vEq(i,j-1) * .25 / obj.dx * (obj.V(i,j) + obj.V(i-1,j));
              
              obj.A(row, obj.mapV(i+1,j)) = obj.A(row, obj.mapV(i+1,j)) + obj.multIJ_vEq(i+1,j) * .25 / obj.dx * (obj.U(i+1,j) + obj.U(i+1,j-1));
              obj.A(row, obj.mapV(i,j)) = obj.A(row, obj.mapV(i,j)) + obj.multIJ_vEq(i,j) * .25 / obj.dx * ( (obj.U(i+1,j) + obj.U(i+1,j-1)) - (obj.U(i,j) + obj.U(i,j-1)) );
              obj.A(row, obj.mapV(i-1,j)) = obj.A(row, obj.mapV(i-1,j)) - obj.multIJ_vEq(i-1,j) * .25 / obj.dx * (obj.U(i,j) + obj.U(i,j-1));
        end
        function vv_y(obj, i, j) % appears in v equation, i.e. obj.mapV(i,j) row of A
%             inout = .25 /obj.dy * (...
%                 (obj.V(i,j+1) + obj.V(i,j))^2 - (obj.V(i,j_1) + obj.V(i,j))^2 ...
%                 );
              row = obj.mapV(i,j); 
              obj.A(row, obj.mapV(i,j+1)) =  obj.A(row, obj.mapV(i,j+1)) + obj.multIJ_vEq(i,j+1) * 2.0 * .25 / obj.dy * (obj.V(i,j+1) + obj.V(i,j));
              obj.A(row, obj.mapV(i,j-1)) =  obj.A(row, obj.mapV(i,j-1)) - obj.multIJ_vEq(i,j-1) * 2.0 * .25 / obj.dy * (obj.V(i,j-1) + obj.V(i,j));
        end
        
        function u_xx(obj, i, j) % appears in u equation, i.e. obj.mapU(i,j) row of A
%             inout = (-obj.nue) * 1.0 / (obj.dx)^2 *(obj.U(i+1,j) - 2 * obj.U(i,j) + obj.U(i-1,j));
              row = obj.mapU(i,j); 
              obj.A(row, obj.mapU(i+1,j)) = obj.A(row, obj.mapU(i+1,j)) + obj.multIJ_uEq(i+1,j) * (-obj.nue) * 1.0 / (obj.dx)^2;  
              obj.A(row, obj.mapU(i,j)) = obj.A(row, obj.mapU(i,j)) - obj.multIJ_uEq(i,j) * (-obj.nue) * 2.0 / (obj.dx)^2; 
              obj.A(row, obj.mapU(i-1,j)) = obj.A(row, obj.mapU(i-1,j)) + obj.multIJ_uEq(i-1,j) * (-obj.nue) * 1.0 / (obj.dx)^2;
        end
        function u_yy(obj, i, j) % appears in u equation, i.e. obj.mapU(i,j) row of A
%             inout = (-obj.nue) * 1.0 / (obj.dy)^2 *(obj.U(i,j+1) - 2 * obj.U(i,j) + obj.U(i,j-1));
              row = obj.mapU(i,j); 
              obj.A(row, obj.mapU(i,j+1)) = obj.A(row, obj.mapU(i,j+1)) + obj.multIJ_uEq(i,j+1) * (-obj.nue) * 1.0 / (obj.dy)^2; 
              obj.A(row, obj.mapU(i,j)) = obj.A(row, obj.mapU(i,j)) - obj.multIJ_uEq(i,j) * (-obj.nue) * 2.0 / (obj.dy)^2; 
              obj.A(row, obj.mapU(i,j-1)) = obj.A(row, obj.mapU(i,j-1)) + obj.multIJ_uEq(i,j-1) * (-obj.nue) * 1.0 / (obj.dy)^2; 
        end
        function v_xx(obj, i, j) % appears in v equation, i.e. obj.mapV(i,j) row of A
%             inout = (-obj.nue) * 1.0 / (obj.dx)^2 *(obj.V(i+1,j) - 2 * obj.V(i,j) + obj.V(i-1,j));
              row = obj.mapV(i,j); 
              obj.A(row, obj.mapV(i+1,j)) = obj.A(row, obj.mapV(i+1,j)) + obj.multIJ_vEq(i+1,j) * (-obj.nue) * 1.0 / (obj.dx)^2;  
              obj.A(row, obj.mapV(i,j)) = obj.A(row, obj.mapV(i,j)) - obj.multIJ_vEq(i,j) * (-obj.nue) * 2.0 / (obj.dx)^2; 
              obj.A(row, obj.mapV(i-1,j)) = obj.A(row, obj.mapV(i-1,j)) + obj.multIJ_vEq(i-1,j) * (-obj.nue) * 1.0 / (obj.dx)^2;
        end
        function v_yy(obj, i, j) % appears in v equation, i.e. obj.mapV(i,j) row of A
%             inout = (-obj.nue) * 1.0 / (obj.dy)^2 *(obj.V(i,j+1) - 2 * obj.V(i,j) + obj.V(i,j-1));
              row = obj.mapV(i,j); 
              obj.A(row, obj.mapV(i,j+1)) = obj.A(row, obj.mapV(i,j+1)) + obj.multIJ_vEq(i,j+1) * (-obj.nue) * 1.0 / (obj.dy)^2; 
              obj.A(row, obj.mapV(i,j)) = obj.A(row, obj.mapV(i,j)) - obj.multIJ_vEq(i,j) * (-obj.nue) * 2.0 / (obj.dy)^2; 
              obj.A(row, obj.mapV(i,j-1)) = obj.A(row, obj.mapV(i,j-1)) + obj.multIJ_vEq(i,j-1) * (-obj.nue) * 1.0 / (obj.dy)^2; 
        end
        
        function p_x(obj, i, j)
%             inout = ( 1.0 / obj.rho) * 1.0 / obj.dx *(obj.P(i,j) - obj.P(i-1,j));
              row = obj.mapU(i,j); 
              obj.A(row, obj.mapP(i,j)) = obj.A(row, obj.mapP(i,j)) + ( 1.0 / obj.rho) * 1.0 / obj.dx; 
              obj.A(row, obj.mapP(i-1,j)) = obj.A(row, obj.mapP(i-1,j)) - ( 1.0 / obj.rho) * 1.0 / obj.dx;
              
              obj.fp(row, obj.mapP_Forfp(i,j)) = obj.fp(row, obj.mapP_Forfp(i,j)) - ( 1.0 / obj.rho) * 1.0 / obj.dx; 
              obj.fp(row, obj.mapP_Forfp(i-1,j)) = obj.fp(row, obj.mapP_Forfp(i-1,j)) + ( 1.0 / obj.rho) * 1.0 / obj.dx; 
        end
        function p_y(obj, i, j)
%             inout = ( 1.0 / obj.rho) * 1.0 / obj.dy *(obj.P(i,j) - obj.P(i,j-1));
              row = obj.mapV(i,j); 
              obj.A(row, obj.mapP(i,j)) = obj.A(row, obj.mapP(i,j)) + ( 1.0 / obj.rho) * 1.0 / obj.dy; 
              obj.A(row, obj.mapP(i,j-1)) = obj.A(row, obj.mapP(i,j-1)) - ( 1.0 / obj.rho) * 1.0 / obj.dy;
              
              obj.fp(row, obj.mapP_Forfp(i,j)) = obj.fp(row, obj.mapP_Forfp(i,j)) + ( 1.0 / obj.rho) * 1.0 / obj.dy; 
              obj.fp(row, obj.mapP_Forfp(i,j-1)) = obj.fp(row, obj.mapP_Forfp(i,j-1)) - ( 1.0 / obj.rho) * 1.0 / obj.dy;               
        end
        function u_t(obj, i, j)
              row = obj.mapU(i,j); 
              obj.A(row, obj.mapU(i,j)) = obj.A(row, obj.mapU(i,j)) + 1.0 / obj.dt;
        end
        function v_t(obj, i, j)
              row = obj.mapV(i,j); 
              obj.A(row, obj.mapV(i,j)) = obj.A(row, obj.mapV(i,j)) + 1.0 / obj.dt;
        end
        function u_x(obj, i, j)
%             inout = 1.0 / obj.dx *(obj.U(i+1,j) - obj.U(i,j));
              row = obj.mapP(i,j); 
              obj.A(row, obj.mapU(i,j)) = obj.A(row, obj.mapU(i,j)) - 1.0 / obj.dx;  
              obj.A(row, obj.mapU(i+1,j)) = obj.A(row, obj.mapU(i+1,j)) + 1.0 / obj.dx;
              
              row_gu = obj.mapP_Forfp(i,j);
              obj.gu(row_gu, obj.mapU(i,j)) = obj.gu(row_gu, obj.mapU(i,j)) - 1.0 / obj.dx;
              obj.gu(row_gu, obj.mapU(i+1,j)) = obj.gu(row_gu, obj.mapU(i+1,j)) + 1.0 / obj.dx;
        end
        function v_y(obj, i, j)
%             inout = 1.0 / obj.dy *(obj.V(i,j+1) - obj.V(i,j));
              row = obj.mapP(i,j); 
              obj.A(row, obj.mapV(i,j)) = obj.A(row, obj.mapV(i,j)) - 1.0 / obj.dy;  
              obj.A(row, obj.mapV(i,j+1)) = obj.A(row, obj.mapV(i,j+1)) + 1.0 / obj.dy;    
              
              row_gu = obj.mapP_Forfp(i,j);
              obj.gu(row_gu, obj.mapV(i,j)) = obj.gu(row_gu, obj.mapV(i,j)) - 1.0 / obj.dy;
              obj.gu(row_gu, obj.mapV(i,j+1)) = obj.gu(row_gu, obj.mapV(i,j+1)) + 1.0 / obj.dy;
        end

        %%
        function SetOnePressureA(obj)
            i = 1;%obj.ni - 1;
            j = 1;%obj.nj - 1;
            row = obj.mapP(obj.ni - 1, obj.nj - 1) + 1;
            obj.A(row, obj.mapP(i,j)) = 1;
        end
        %%
        function SetOnePressureB(obj)
            row = obj.mapP(obj.ni - 1, obj.nj - 1) + 1;
            obj.b(row) = 0;
        end
        %%
        function CopyNewToOld(obj)
            obj.u_old = obj.u;
            obj.v_old = obj.v;
        end
        
        %%
        function CalcJacobian_and_Inverse(obj)
            for i=1:obj.ni
                for j = 1:obj.nj-1
                    % x-momentum
                    obj.u_t(i,j);
                    obj.uu_x(i,j);
                    obj.uv_y(i,j);
                    obj.p_x(i,j);
                    obj.u_xx(i,j);
                    obj.u_yy(i,j);
                end
            end
            for i=1:obj.ni-1
                for j = 1:obj.nj
                    % y-momentum
                    obj.v_t(i,j);
                    obj.vu_x(i,j);
                    obj.vv_y(i,j);
                    obj.p_y(i,j);
                    obj.v_xx(i,j);
                    obj.v_yy(i,j);
                end
            end
            for i=1:obj.ni-1
                for j = 1:obj.nj-1
                    % continuity
                    obj.u_x(i, j);
                    obj.v_y(i, j);
                end 
            end
            %set one pressure
            obj.SetOnePressureA()
            % calce InvA
            obj.invA = pinv(obj.A);         
        end      
        %%
        function CalcRHS(obj)
            for i=1:obj.ni
                for j = 1:obj.nj-1
                    obj.rhsU(i, j);
                end
            end
            for i=1:obj.ni-1
                for j = 1:obj.nj
                    obj.rhsV(i, j);
                end
            end
            for i=1:obj.ni-1
                for j = 1:obj.nj-1
                    obj.rhsContinuity(i, j);
                end 
            end
            %set one pressure
            obj.SetOnePressureB()
        end
        
        %% Retrieve U, V, and P values
        function out = U(obj,i,j)
            if (i < 1) out = 2 * obj.uL - obj.u(i + 2, j); 
            elseif (i > obj.ni) out = 2 * obj.uR - obj.u(i - 2, j); 
            elseif (j < 1) out = 2 * obj.uB - obj.u(i, j + 1);
            elseif (j > obj.nj - 1) out = 2 * obj.uT - obj.u(i, j - 1);
            else out = obj.u(i,j);
            end
        end
        function out = V(obj,i,j)
            if (i < 1) out = 2 * obj.vL - obj.v(i + 1, j); 
            elseif (i > obj.ni - 1) out = 2 * obj.vR - obj.v(i - 1, j); 
            elseif (j < 1) out = 2 * obj.vB - obj.v(i, j + 2);
            elseif (j > obj.nj) out = 2 * obj.vT - obj.v(i, j - 2);
            else out = obj.v(i,j);
            end
        end
        function out = P(obj,i,j)
            if (i < 1) out = obj.p(i + 1, j); 
            elseif (i > obj.ni - 1) out = obj.p(i - 1, j); 
            elseif (j < 1) out = obj.p(i, j + 1);
            elseif (j > obj.nj - 1) out = obj.p(i, j - 1);
            else out = obj.p(i,j);
            end
        end
        %% solve and update            
        function SolveAndUpdate(obj)
            obj.CheckRank();
            delX = obj.invA * obj.b;
            for i=1:obj.ni
                for j = 1:obj.nj-1
                    obj.u(i, j) = obj.u(i, j) - delX(obj.mapU(i,j));
                end
            end
            for i=1:obj.ni-1
                for j = 1:obj.nj
                    obj.v(i, j) = obj.v(i, j) - delX(obj.mapV(i,j));
                end
            end
            for i=1:obj.ni-1
                for j = 1:obj.nj-1
                    obj.p(i, j) = obj.p(i, j) - delX(obj.mapP(i,j));
                end
            end
        end
        %% Check Rank
        function CheckRank(obj)
            rankA = rank(obj.A);
            sizeA = size(obj.A, 2);
            GuFp = obj.gu*obj.fp;
            sizeGuFp = size(GuFp);
            rankGuFp = rank(GuFp);
            fprintf('sizeA %d rankA %d sizeGuFp %d rankGuFp %d \n', rankA, sizeA, sizeGuFp, rankGuFp);
        end
        %% Streamline
        function Visualize(obj)
            for i=1:obj.ni-1
                for j=1:obj.nj-1
                    x(i,j) = (i-0.5)*obj.dx;
                    y(i,j) = (j-0.5)*obj.dy;
                end
            end
            figure(obj.fgU);

            quiver(x, y, obj.u(1:obj.ni - 1,:), obj.v(:,1:obj.nj-1),2,'k');
%             norm([obj.u(1:obj.ni - 1,:), obj.v(:,1:obj.nj-1)],'inf')

            
%             figure(obj.fgP);
%             surf(obj.p);
%             x = [obj.dx/2 : obj.dx : obj.lx - obj.dx/2];
%             y = [obj.dy/2 : obj.dy : obj.ly - obj.dy/2];

            mU = obj.u(1:obj.ni - 1,:);
            mV = obj.v(:,1:obj.nj-1);
            
            [U2,V2,x2,y2]=obj.reformUV(mU, mV);

%             figure
            
            y00=[.1*obj.nj*obj.dy:2*obj.dy:obj.nj*obj.dy];
            x00=.5*obj.ni*obj.dx*ones(size(y00));
            x01 = [.001*obj.ni*obj.dx:1*obj.dx:.9*obj.ni*obj.dx];
            y01=.7*obj.nj*obj.dy*ones(size(x01));
%             h = streamline(x2,y2,U2,V2,x01,y01,[.05,20000]);
%             set(h,'Color','k');
            
%             streamslice(x2,y2,U2,V2);
%             obj.CalcStream();
%             [nx,ny]=size(obj.psi);
%             psi2 = obj.psi(3:nx-2,3:ny-2);
%             contour(obj.dx*[3:nx-2], obj.dy*[3:ny-2], psi2', 50);
%             xlim([0,obj.dx*obj.ni]);
%             ylim([0,obj.dy*obj.nj]);
            
            
%             contourf(sqrt((obj.u(1:obj.ni-1, 1:obj.nj-1)).^2 + (obj.v(1:obj.ni-1, 1:obj.nj-1)).^2), 20);
        end
        
        %% mapping functions
        function out = mapU(obj, iOrig, jOrig)
            i = iOrig;
            j = jOrig;
            if (iOrig < 1) i = iOrig + 2;
            elseif (iOrig > obj.ni) i = iOrig - 2;
            end
            if (jOrig < 1) j = jOrig + 1;
            elseif (jOrig > obj.nj - 1) j = jOrig - 1;
            end
            out = (i - 1) * (obj.nj-1) + j;
        end
        function out = mapV(obj, iOrig, jOrig)
            i = iOrig;
            j = jOrig;
            if (iOrig < 1) i = iOrig + 1;
            elseif (iOrig > obj.ni - 1) i = iOrig - 1;
            end
            if (jOrig < 1) j = jOrig + 2;
            elseif (jOrig > obj.nj) j = jOrig - 2;
            end
            start = obj.ni * (obj.nj - 1);
            out = start + (i - 1) * obj.nj + j;
        end
        function out = mapP(obj, iOrig, jOrig)
            i = iOrig;
            j = jOrig;
            if (iOrig < 1) i = iOrig + 1;
            elseif (iOrig > obj.ni-1) i = iOrig - 1;
            end
            if (jOrig < 1) j = jOrig + 1;
            elseif (jOrig > obj.nj-1) j = jOrig - 1;
            end
            start = obj.ni * (obj.nj - 1) + (obj.ni - 1) * obj.nj; % velocities
            out = start + (i - 1) * (obj.nj-1) + j;
        end
        function out = mapP_Forfp(obj, iOrig, jOrig)
            i = iOrig;
            j = jOrig;
            if (iOrig < 1) i = iOrig + 1;
            elseif (iOrig > obj.ni-1) i = iOrig - 1;
            end
            if (jOrig < 1) j = jOrig + 1;
            elseif (jOrig > obj.nj-1) j = jOrig - 1;
            end
            out = (i - 1) * (obj.nj-1) + j;
        end
        function out = multIJ_uEq(obj, i, j)
            out = 1;
            if ((i < 1) | (i > obj.ni) | (j < 1) | (j > obj.nj - 1))  
                out = -1;
            end
        end 
        function out = multIJ_vEq(obj, i, j)
            out = 1;
            if ((i < 1) | (i > obj.ni - 1) | (j < 1) | (j > obj.nj))  
                out = -1;
            end
        end
        %% Check Rank
        function CalcStream(obj)
            obj.psi=zeros(size(obj.psi));
            for i=1:obj.ni-1
                for j=2:obj.nj-1
                    obj.psi(i,j) = obj.psi(i,j-1) + obj.dy * obj.u(i,j-1);
                end
            end
        end
        %% write to file
        function WriteToFile(obj)
            obj.fileCounter = obj.fileCounter + 1;
            fileName1 = strcat('csvFiles/gridData', num2str(obj.fileCounter));
            fileName1 = strcat(fileName1,'.csv');
            fileID1 = fopen(fileName1,'w');
            for i=1:obj.ni-1
                for j=1:obj.nj-1
                    fprintf(fileID1, '%d, %d, %d, %f, %f, %f, %f, %f, %f, %f,\n', i, j, 0, i * obj.dx, ...
                        j * obj.dy, 0, obj.u(i, j), obj.v(i, j), 0, obj.p(i,j));
                end
            end
            fclose(fileID1);
        end
        
        %%
        function [U2,V2,x2,y2]=reformUV(obj, U, V)
            nx = size(U,1);
            ny = size(U,2);
            U2 = ones(ny,nx);
            V2 = ones(ny,nx);
            for i=1:ny
                for j=1:nx
                    U2(i,j) = U(j, ny-i+1);
                    V2(i,j) = V(j, ny-i+1);
                    x2(i,j) = j*obj.dx;
                    y2(i,j) = (ny-i+1)*obj.dy;
                end
            end
        end
        %% calc center vorticity
        function omega = CenterOmega(obj)
            i2 = 0;
            j2 = 0;
            dx = 0;
            dy = 0;
            if (mod(obj.ni,2) == 0)
                i2 = obj.ni/2;
                if (mod(obj.nj,2) == 0)
                    j2 = obj.nj/2;
                    dudy = .5 / obj.dy * (.5 * (obj.u(i2,j2+1) + obj.u(i2+1,j2+1)) ...
                        - .5 * (obj.u(i2,j2-1) + obj.u(i2+1,j2-1)));
                    dvdx = .5 / obj.dx * (.5 * (obj.v(i2+1,j2) + obj.v(i2+1,j2+1)) ...
                        - .5 * (obj.v(i2-1,j2) + obj.v(i2-1,j2+1)));
                else
                    j2 = (obj.nj+1)/2;
                    dudy = 1.0 / obj.dy * (.5 * (obj.u(i2,j2) + obj.u(i2+1,j2)) ...
                        - .5 * (obj.u(i2,j2-1) + obj.u(i2+1,j2-1)));
                    dvdx = .5 / obj.dx * (obj.v(i2+1,j2) - obj.v(i2-1,j2));
                end
            else
                i2 = (obj.ni+1)/2;
                if (mod(obj.nj,2) == 0)
                    j2 = obj.nj/2;
                    dudy = .5 / obj.dy * (obj.u(i2,j2+1) - obj.u(i2,j2-1));
                    dvdx = 1.0 / obj.dx * (.5 * (obj.v(i2,j2) + obj.v(i2,j2+1)) ...
                        - .5 * (obj.v(i2-1,j2) + obj.v(i2-1,j2+1)));
                else
                    j2 = (obj.nj+1)/2;
                    dudy = 1.0 / obj.dy * (obj.u(i2,j2) - obj.u(i2,j2-1));
                    dvdx = 1.0 / obj.dx * (obj.v(i2,j2) - obj.v(i2-1,j2));
                end
            end
            omega = dvdx - dudy;
        end
       
    end %methods
end % classdef
    