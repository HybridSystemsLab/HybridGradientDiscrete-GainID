classdef hybridDiscretizedGD < HybridSystem
    
    properties(SetAccess = immutable)
        K 
        Fw 
        zMax 
        lambda 
        Va 
        Vb 
        Ts 
        A 
        B 
        theta 
        gammac 
        gammad 
        s
    end

    methods 
        function this = hybridDiscretizedGD(parameters)
            state_dim = 6; % (z,tau,thetahat,tauS)
            this = this@HybridSystem(state_dim);
            
            this.K = parameters.K;
            this.Fw = parameters.Fw;
            this.zMax = parameters.zMax;
            this.lambda = parameters.lambda;
            this.Va = parameters.Va;
            this.Vb = parameters.Vb;
            this.Ts = parameters.Ts;
            this.A = parameters.A;
            this.B = parameters.B;
            this.theta = parameters.theta;
            this.gammac = parameters.gammac;
            this.gammad = parameters.gammad;
            this.s = parameters.s;
        end

        function xdot = flowMap(this, x, t, j)
            z = x(1:2);
            tau = x(3);
            thetahat = x(4:5);
            tauS = x(6);
            
            if tau < this.Ts
               v = this.Fw*this.Vb;
            else
               v = this.Fw*this.Va;
            end
            u = -this.K*z + v;

            zdot = this.A*z + this.B*u;
            tdot = 1;
            thetahatdot = 0*thetahat;
            tauSdot = 1;

            xdot = [zdot; tdot; thetahatdot; tauSdot];
        end

        function xplus = jumpMap(this, x, t, j)
            z = x(1:2);
            tau = x(3);
            thetahat = x(4:5);
            tauS = x(6);
            
            if tauS >= this.s        % take a sample during flows
                zplus = z;
                tauplus = tau;
                tauSplus = 0;
                
                if tau < this.Ts
                   v = this.Fw*this.Vb;
                else
                   v = this.Fw*this.Va;
                end
                u = -this.K*z + v;
                y = u - v;
                psi = -z;
                thetahatplus = thetahat + this.s*this.gammac*psi*(y - psi'*thetahat);
                
            elseif tau >= this.Ts    % jump due to the jump set
                zplus = [z(1); -this.lambda*z(2)];
                tauplus = 0;
                tauSplus = tauS;

                if tauplus < this.Ts
                   vplus = this.Fw*this.Vb;
                else
                   vplus = this.Fw*this.Va;
                end
                uplus = -this.K*zplus + vplus;
                yplus = uplus - vplus;
                psiplus = -zplus;
                thetahatplus = thetahat + this.gammad*psiplus*(yplus - psiplus'*thetahat)/(1 + this.gammad*norm(psiplus)^2);
                
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
            z = x(1:2);
            tau = x(3);
            thetahat = x(4:5);
            tauS = x(6);
            
            inD1 = tauS >= this.s;                   % take a sample during flows
            inD2 = z(1) >= this.zMax && z(2)>= 0;    % jump due to jump set
            inD = inD1 || inD2;
        end
    end
end
