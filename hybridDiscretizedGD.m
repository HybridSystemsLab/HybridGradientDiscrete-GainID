classdef hybridDiscretizedGD < HybridSystem

    methods 
        function this = hybridDiscretizedGD()
            state_dim = 6; % (z,tau,thetahat,tauS)
            this = this@HybridSystem(state_dim);
        end

        function xdot = flowMap(this, x, t, j)
            global A B Ts Fw Va Vb K
            
            z = x(1:2);
            tau = x(3);
            thetahat = x(4:5);
            tauS = x(6);
            
            if tau < Ts
               v = Fw*Vb;
            else
               v = Fw*Va;
            end
            u = -K*z + v;

            zdot = A*z + B*u;
            tdot = 1;
            thetahatdot = 0*thetahat;
            tauSdot = 1;

            xdot = [zdot; tdot; thetahatdot; tauSdot];
        end

        function xplus = jumpMap(this, x, t, j)
            global lambda Ts Fw Va Vb K gammac gammad s
            
            z = x(1:2);
            tau = x(3);
            thetahat = x(4:5);
            tauS = x(6);
            
            if tauS >= s        % take a sample during flows
                zplus = z;
                tauplus = tau;
                tauSplus = 0;
                
                if tau < Ts
                   v = Fw*Vb;
                else
                   v = Fw*Va;
                end
                u = -K*z + v;
                y = u - v;
                psi = -z;
                thetahatplus = thetahat + s*gammac*psi*(y - psi'*thetahat);
                
            elseif tau >= Ts    % jump due to the jump set
                zplus = [z(1); -lambda*z(2)];
                tauplus = 0;
                tauSplus = tauS;

                if tauplus < Ts
                   vplus = Fw*Vb;
                else
                   vplus = Fw*Va;
                end
                uplus = -K*zplus + vplus;
                yplus = uplus - vplus;
                psiplus = -zplus;
                thetahatplus = thetahat + gammad*psiplus*(yplus - psiplus'*thetahat)/(1 + gammad*norm(psiplus)^2);
                
            else
                zplus = z;
                tauplus = tau;
                thetahatplus = thetahat;
                tauSplus = tauS;
            end

            
            xplus = [zplus; tauplus; thetahatplus; tauSplus];
        end
        
        function inC = flowSetIndicator(this, x)
            inC = 1;
        end

        function inD = jumpSetIndicator(this, x)
            global zMax s
            
            z = x(1:2);
            tau = x(3);
            thetahat = x(4:5);
            tauS = x(6);
            
            inD1 = tauS >= s;                   % take a sample during flows
            inD2 = z(1) >= zMax && z(2)>= 0;    % jump due to jump set
            inD = inD1 || inD2;
        end
    end
end
