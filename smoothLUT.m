classdef smoothLUT
%SMOOTHLUT generates uniform knots for a fast and smooth bicubic B-spline 
% interpolation of regularly spaced two dimensional data using the 
% algorithm described in [1].
%
% Syntax:
%   myLUT = smoothLUT(xk,yk,zk)
%   myLUT = smoothLUT(xk,yk,zk,r)
%   myLUT = smoothLUT(xk,yk,zk,r,G)
%   myLUT = smoothLUT(xk,yk,zk,r,G,lb,ub)
%
%
% Description:
%   myLUT = smoothLUT(xk,yk,zk) generates knot values c for a smooth 
%   bicubic B-spline interpolation of the data xk, yk and zk. The knot
%   vectors equal the given xk- and yk-vectors, just expanded by one
%   knot to each side. The resulting knots are stored in the property c and
%   can therefore be accessed by myLUT.c .
%
%   myLUT = smoothLUT(xk,yk,zk,r) generates knot values c such that the 
%   weight of the smoothness over the accuracy is r. The smoothness is 
%   hereby represented by the integral of the squared second derivative, 
%   the accuracy is respresented by the sum of the squared deviation
%   between the given and the interpolated data values.
%
%   myLUT = smoothLUT(xk,yk,zk,r,G) introduces the individual accuarcy
%   weight matrix G which allows to lower the accuracy weight of every data
%   point.
%
%   myLUT = smoothLUT(xk,yk,zk,r,G,lb,ub) makes the z-values of the 
%   interpolated table stay within the range of the lower bound lb and the 
%   upper bound ub.
%
%
% Examples:
%   Creation of knots 
%     [~,~,zk]=peaks(21);
%     test=smoothLUT([0:1:20]',[0:1:20],zk)
%     test.showInterpolation(30,30)
%
%   Smoothing of noisy data
%     [~,~,zk]=peaks(21);
%     zk=randn(21)+zk;
%     test=smoothLUT([0:1:20]',[0:1:20],zk)
%     test.showInterpolation(30,30,'derivatives',true)
%     test=smoothLUT([0:1:20]',[0:1:20],zk,2)
%     test.showInterpolation(30,30,'derivatives',true)
%
%   Local smoothing of noisy data
%     [~,~,zk]=peaks(21);
%     G=[0.3*ones(10,21); ones(11,21)];
%     test=smoothLUT([0:1:20]',[0:1:20],randn(21)+zk,[],G)
%     test.showInterpolation(30,30,'derivatives',true)
%
%   Addition of bounds in z-direction
%     [~,~,zk]=peaks(21);
%     test=smoothLUT([0:1:20]',[0:1:20],zk,[],[],-2,2)
%     test.showInterpolation(30,30)
%
%
% [1] R. Mitze, D. Dillkötter, S. Gros, A. Schild and M. Mönnigmann. Fast and
%  smooth surface B-spline interpolation for regularly spaced data used in
%  system modeling to make MPC real-time feasible. Proceedings of the
%  European Control Conference 2018 (ECC18), Limassol: 667-672, 2018
%
% 
% See also Lut2Vhdl

%  AUTHORS
%    
%   2018    Ruth Mitze and Martin Mönnigmann:
%           Ruhr-Universität Bochum
%           Systems Theory and Automatic Control
%   mailto: ruth.mitze@rub.de 
%   mailto: martin.moennigmann@rub.de

%  LICENSE
%    
%    This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation; either version 3 of the License, or (at
%  your option) any later version.
%    This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
%  General Public License for more details.
%    You should have received a copy of the GNU Lesser General Public
%  License along with this library; if not, write to the  Free Software
%  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%  MA 02110-1335 USA
    
    properties (SetAccess = private)
        xk         % x-data xi, nx1 equidistant numerical vector
        yk         % y-data yj, 1xm equidistant numerical vector
        zk         % z-data zij, nxm numerical matrix
        r  = .1    % overall smoothness weight, positive scalar, default=0.1
        G  = 1     % individual accuracy weights, nxm numerical matrix [0,...,1], default=1
        lb = []    % lower z-bound for the interpolated table, scalar, default=[] (deactivated)
        ub = []    % upper z-bound for the interpolated table, scalar, default=[] (deactivated)
        c          % knot values, (n+2)x(m+2) numerical matrix
    end
    
    properties (SetAccess = private, GetAccess = private)
        hx % step size xk-data
        hy % step size yk-data
        x0 % first entry of xk-data
        y0 % first entry of yk-data
        n  % number of elements in xk-data
        m  % number of elements in yk-data
    end
    
    methods
        function obj = set.xk(obj,new_xk)
            if isnumeric(new_xk) && isreal(new_xk) && length(new_xk(1,:))==1 && length(new_xk(:,1))>1 && max(new_xk(2:end)-new_xk(1:end-1))==min(new_xk(2:end)-new_xk(1:end-1))
                obj.xk = new_xk;
                obj.n  = length(new_xk);
                obj.hx = new_xk(2)-new_xk(1);
                obj.x0 = new_xk(1);
            else
                error('xk must be an equidistant numerical coloumn vector');
            end
        end
        
        function obj = set.yk(obj,new_yk)
            if isnumeric(new_yk) && isreal(new_yk) && length(new_yk(1,:))>1 && length(new_yk(:,1))==1 && max(new_yk(2:end)-new_yk(1:end-1))==min(new_yk(2:end)-new_yk(1:end-1))
                obj.yk = new_yk;
                obj.m  = length(new_yk);
                obj.hy = new_yk(2)-new_yk(1);
                obj.y0 = new_yk(1);
            else
                error('yk must be an equidistant numerical line vector');
            end
        end
        
        function obj = set.zk(obj,new_zk)
            if isnumeric(new_zk) && isreal(new_zk) && length(new_zk(1,:))>1 && length(new_zk(:,1))>1
                obj.zk = new_zk;
            else
                error('zk must be a numerical matrix');
            end
        end
        
        function obj = set.r(obj,new_r)
            if isempty(new_r)
            elseif isnumeric(new_r) && isscalar(new_r) && isreal(new_r) && new_r>=0
                obj.r = new_r;
            else
                error('r must be a positive scalar, or empty');
            end
        end
        
        function obj = set.G(obj,new_G)
            if isempty(new_G)
            elseif isnumeric(new_G) && isreal(new_G) && min(min(new_G))>=0 && max(max(new_G))<=1
                obj.G = new_G;
            else
                error('G must be a numerical matrix with entries between zero and one, a scalar, or empty');
            end
        end
        
        function obj = set.lb(obj,new_lb)
            if isempty(new_lb)
            elseif isnumeric(new_lb) && isscalar(new_lb) && isreal(new_lb)
                obj.lb = new_lb;
            else
                error('lb must be a scalar, or empty');
            end
        end
        
        function obj = set.ub(obj,new_ub)
            if isempty(new_ub)
            elseif isnumeric(new_ub) && isscalar(new_ub) && isreal(new_ub)
                obj.ub = new_ub;
            else
                error('ub must be a scalar, or empty');
            end
        end
        
        
        function obj = smoothLUT(new_xk,new_yk,new_zk,new_r,new_G,new_lb,new_ub)
            %SMOOTHLUT constructs an instance of class smoothLUT
            %
            % 
            % See also smoothLUT
            if nargin<3
                error('not enough input arguments (at least 3: xk,yk,zk)')
            elseif nargin>7
                error('too many input arguments (not more than 7: xk,yk,zk,r,G,lb,ub)')
            else
                if min([length(new_xk) length(new_yk)]~=size(new_zk))
                    error('dimensions of xk (nx1), yk (1xm) and zk (nxm) do not match')
                else
                    obj.xk = new_xk;
                    obj.yk = new_yk;
                    obj.zk = new_zk;
                    if nargin>3
                        obj.r = new_r;
                        if nargin>4
                            if max(size(new_G)~=[1 1]) && max(size(new_G)~=size(new_zk)) && ~isempty(new_G)
                                error('individual weight matrix G must be a matrix of size(zk), a scalar, or empty')
                            else
                                obj.G = new_G;
                                if nargin>5
                                    obj.lb = new_lb;
                                    if nargin>6
                                        if ~isempty(new_lb) && ~isempty(new_ub) && new_lb>=new_ub
                                            error('lower bound lb must be smaller than upper bound ub')
                                        else
                                            obj.ub = new_ub;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % compute matrices A, Axx, Axy and Ayy (see Sec. IIIA in [1])
            A             = getMatrixA(obj);
            [Axx,Axy,Ayy] = getMatrixA2(obj);
            
            zk_transp = (obj.zk)';
            zk_row    = zk_transp(:);
            
            % set individual weighting matrix G
            if ~isscalar(obj.G)
                indivWeights = obj.G';
                indivWeights = indivWeights(:);
                indivWeights = diag(indivWeights);
            else
                indivWeights = obj.G;
            end
            
            % set lower an upper boundaries
            if isempty(obj.lb)
                LB = [];
            else
                LB = obj.lb*ones((obj.n+2)*(obj.m+2),1);
            end
            if isempty(obj.ub)
                UB = [];
            else
                UB = obj.ub*ones((obj.n+2)*(obj.m+2),1);
            end
            
            % set cost function matrices
            P      = A'*(indivWeights)'*indivWeights*A+obj.r*(Axx+2*Axy+Ayy);
            Q      = 2*zk_row'*(indivWeights)'*indivWeights*A;
            
            % compute knots with quadratic program (see Eq. 18 in [1])
            c_row  = quadprog(2*P,(-Q)',[],[],[],[],LB,UB);
            c_line = c_row';
            obj.c  = (reshape(c_line,[obj.m+2,obj.n+2]))';
        end
        
        
        function z = determineValue(obj,x,y,varargin)
            %DETERMINEVALUE determines the interpolated value at a given position.
            % 
            % Syntax:
            %   myLUT.DETERMINEVALUE(x,y)
            %
            %
            % Examples:
            %   Determine interpolated value
            %     [~,~,zk]=peaks(21);
            %     test=smoothLUT([0:1:20]',[0:1:20],zk)
            %     test.determineValue(15,11)
            %
            %   Determine interpolated value and its derivatives up to 2nd order
            %     [~,~,zk]=peaks(21);
            %     test=smoothLUT([0:1:20]',[0:1:20],zk)
            %     test.determineValue(15,11,'derivatives',true)
            % 
            %
            % Input Arguments:
            %   x must be a scalar in the range of xk
            %   y must be a scalar in the range of yk
            %
            %
            % Name-Value Pair Arguments:
            %   'derivatives' (default = false):
            %       false: determines interpolated value 
            %              ans = s(x,y)
            %       true:  determines interpolated value and derivatives up to 2nd order 
            %              ans = [s(x,y), ds(x,y)/dx, ds(x,y)/dy, d2s(x,y)/x2, d2s(x,y)/dy2, d2s(x,y)/dxdy]
            
            if x<obj.x0 || x>=obj.x0+(obj.n-1)*obj.hx
                error('x-position must lie in the span of xk')
            elseif y<obj.y0 || y>=obj.y0+(obj.m-1)*obj.hy
                error('y-position must lie in the span of yk')
            else

                p                  = inputParser;
                defaultDerivatives = false;
                addParameter(p,'derivatives',defaultDerivatives,@islogical);
                parse(p,varargin{:});
                
                % find intervals (see Eq. 20 in [1])
                k = floor((x-obj.x0)/obj.hx)+2;
                l = floor((y-obj.y0)/obj.hy)+2;

                % compute basis functions
                [Bx,Bx_prime,Bx_primeprime] = getBasisfct(obj,x,k,obj.hx,obj.x0,p.Results.derivatives);
                [By,By_prime,By_primeprime] = getBasisfct(obj,y,l,obj.hy,obj.y0,p.Results.derivatives);

                % relevant knots
                c_relevant = obj.c(k+(-1:2),l+(-1:2));

                % compute interpolated value (see Eq. 19 in [1])
                z(1) = Bx*c_relevant*By';
                if p.Results.derivatives
                    z(2) = Bx_prime*c_relevant*By';
                    z(3) = Bx*c_relevant*By_prime';
                    z(4) = Bx_primeprime*c_relevant*By';
                    z(5) = Bx*c_relevant*By_primeprime';
                    z(6) = Bx_prime*c_relevant*By_prime';
                end
            end
        end
        
        function showInterpolation(obj,amountX,amountY,varargin)
            %SHOWINTERPOLATION plots interpolated look-up table.
            %
            % Syntax:
            %   myLUT.SHOWINTERPOLATION(nx,ny)
            % 
            %
            % Examples:
            %   Plot interpolated table
            %     [~,~,zk]=peaks(21);
            %     test=smoothLUT([0:1:20]',[0:1:20],zk)
            %     test.showInterpolation(30,30)
            %
            %   Plot interpolated table and derivatives up to 2nd order
            %     [~,~,zk]=peaks(21);
            %     test=smoothLUT([0:1:20]',[0:1:20],zk)
            %     test.showInterpolation(30,30,'derivatives',true)
            %
            %
            % Input Arguments:
            %   nx: number of plotted points in x-direction, positive integer
            %   ny: number of plotted points in y-direction, positive integer
            %
            %
            % Name-Value Pair Arguments:
            %   'derivatives' (default = false):
            %       false: plots manifold
            %       true:  plots manifold and derivatives up to 2nd order
            
            if ~isnumeric(amountX) || ~isscalar(amountX) || ~isreal(amountX) || amountX<1 || mod(amountX,1)
                error('number of plotted points in x-direction must be a positive integer')
            elseif ~isnumeric(amountY) || ~isscalar(amountY) || ~isreal(amountY) || amountY<1 || mod(amountY,1)
                error('number of plotted points in y-direction must be a positive integer')
                
            else
                p=inputParser;
                defaultDerivatives = false;
                addParameter(p,'derivatives',defaultDerivatives,@islogical);
                parse(p,varargin{:});
                
                % grid over (x,y)-space, determine z-values at grid points
                x = obj.x0:obj.hx*(obj.n-1)/amountX:obj.x0+obj.hx*(obj.n-1);
                y = obj.y0:obj.hy*(obj.m-1)/amountY:obj.y0+obj.hy*(obj.m-1);
                z = zeros(amountX,amountY,6);
                for i=1:amountX
                    for j=1:amountY
                        z(i,j,:) = determineValue(obj,x(i),y(j),'derivatives',p.Results.derivatives);
                    end
                end                
                [x_map,y_map]   = meshgrid(x(1:end-1),y(1:end-1));
                x_map           = x_map';
                y_map           = y_map';
                [xk_map,yk_map] = meshgrid(obj.xk,obj.yk);
                xk_map          = xk_map';
                yk_map          = yk_map';

                % initialize figure
                figure()
                if p.Results.derivatives
                    set(gcf,'Position',[0 0 2000 1000])
                end
                
                % plot
                if p.Results.derivatives
                    subplot(2,3,1)
                end
                hold on
                mesh(xk_map,yk_map,obj.zk,'EdgeColor','r')
                hidden off
                surf(x_map,y_map,z(:,:,1))
                xlabel('x')
                ylabel('y')
                zlabel('z')
                view(3)

                if p.Results.derivatives
                    subplot(2,3,2)
                    surf(x_map,y_map,z(:,:,2))
                    xlabel('x')
                    ylabel('y')
                    zlabel('dz/dx')

                    subplot(2,3,3)
                    surf(x_map,y_map,z(:,:,3))
                    xlabel('x')
                    ylabel('y')
                    zlabel('dz/dy')

                    subplot(2,3,4)
                    surf(x_map,y_map,z(:,:,4))
                    xlabel('x')
                    ylabel('y')
                    zlabel('d^2z/dx^2')

                    subplot(2,3,5)
                    surf(x_map,y_map,z(:,:,5))
                    xlabel('x')
                    ylabel('y')
                    zlabel('d^2z/dy^2')

                    subplot(2,3,6)
                    surf(x_map,y_map,z(:,:,6))
                    xlabel('x')
                    ylabel('y')
                    zlabel('d^2z/dxdy')
                end
            end
        end
        
        function write2file(obj,fileName)
            %WRITE2FILE saves knot data in .mat-file
            % 
            % Syntax:
            %   myLUT.WRITE2FILE('fileName')
            %
            %
            % Desctiption:
            %   knotvector in x-direction: x_start:x_stepsize:x_end
            %   knotvector in y-direction: y_start:y_stepsize:y_end
            %   knot values (z-direction): c
            %
            %
            % Examples:
            %   Save knot data
            %     [~,~,zk]=peaks(21);
            %     test=smoothLUT([0:1:20]',[0:1:20],zk)
            %     test.write2file('testFile')
            %
            %   Access knot data
            %     load('testFile.mat')
            %     x_knotvector=x_start:x_stepsize:x_end;
            %     y_knotvector=y_start:y_stepsize:y_end;
            %     c_knotvalues=c;
            %
            % 
            % Input Arguments:
            %   fileName must be a char or a string
            
            if ischar(fileName) || isstring(fileName)
                x_start    = obj.x0-obj.hx;
                x_stepsize = obj.hx;
                x_end      = obj.x0+obj.n*obj.hx;
                y_start    = obj.y0-obj.hy;
                y_stepsize = obj.hy;
                y_end      = obj.y0+obj.m*obj.hy;
                c          = obj.c;

                save(fileName,'x_start','x_stepsize','x_end','y_start','y_stepsize','y_end','c')
            else
                error('filename must be a char or a string')
            end
        end
    end
    
    methods (Access = private)
        function [A] = getMatrixA(obj) 
            % s(x_St,y_St)=A*c with A=B^T.*B (see Sec. IIIA in [1])
            A     = zeros(obj.n*obj.m,(obj.n+2)*(obj.m+2));
            Zeile = 1;
            for j = 0:obj.n-1
                for i = 1:obj.m
                    A(Zeile,j*(obj.m+2)+i:j*(obj.m+2)+i+2)                         = [1 4 1]/36;
                    A(Zeile,j*(obj.m+2)+i+(obj.m+2):j*(obj.m+2)+i+(obj.m+2)+2)     = [4 16 4]/36;
                    A(Zeile,j*(obj.m+2)+i+2*(obj.m+2):j*(obj.m+2)+i+2*(obj.m+2)+2) = [1 4 1]/36;
                    Zeile = Zeile+1;
                end
            end
        end
        
        function [Axx,Axy,Ayy] = getMatrixA2(obj) 
            % (see Appendix in [1])
            
            % integral (2.derivative)^2 dxdy = c'*A_uv*c with Auv=int_0_1 (B_3(h),k_Abl)^2 dh
            Aquadr   = [1/7 129/140 3/7 1/140; 129/140 297/35 933/140 3/7; 3/7 933/140 297/35 129/140; 1/140 3/7 129/140 1/7]/36;
            Aiquadr  = [1/5 7/30 -2/5 -1/30; 7/30 17/15 -29/30 -2/5; -2/5 -29/30 17/15 7/30; -1/30 -2/5 7/30 1/5]/4;
            Aiiquadr = [1/3 -1/2 0 1/6; -1/2 1 -1/2 0; 0 -1/2 1 -1/2; 1/6 0 -1/2 1/3];
            
            Axx = zeros((obj.n+2)*(obj.m+2),(obj.n+2)*(obj.m+2));
            for i=1:obj.n-1 % sum over x
                for j=1:obj.m-1 % sum over y
                    for u=1:4 % outer Kronecker
                        for v=1:4 % inner Kronecker
                            Axx((i-1+u-1)*(obj.m+2)+j:(i-1+u-1)*(obj.m+2)+j+3,(i-1+v-1)*(obj.m+2)+j:(i-1+v-1)*(obj.m+2)+j+3) = Axx((i-1+u-1)*(obj.m+2)+j:(i-1+u-1)*(obj.m+2)+j+3,(i-1+v-1)*(obj.m+2)+j:(i-1+v-1)*(obj.m+2)+j+3)+Aiiquadr(u,v)*Aquadr;
                        end
                    end
                end
            end
            
            Ayy = zeros((obj.n+2)*(obj.m+2),(obj.n+2)*(obj.m+2));
            for i=1:obj.n-1
                for j=1:obj.m-1
                    for u=1:4
                        for v=1:4
                            Ayy((i-1+u-1)*(obj.m+2)+j:(i-1+u-1)*(obj.m+2)+j+3,(i-1+v-1)*(obj.m+2)+j:(i-1+v-1)*(obj.m+2)+j+3) = Ayy((i-1+u-1)*(obj.m+2)+j:(i-1+u-1)*(obj.m+2)+j+3,(i-1+v-1)*(obj.m+2)+j:(i-1+v-1)*(obj.m+2)+j+3)+Aquadr(u,v)*Aiiquadr;
                        end
                    end
                end
            end
            
            Axy = zeros((obj.n+2)*(obj.m+2),(obj.n+2)*(obj.m+2));
            for i=1:obj.n-1
                for j=1:obj.m-1
                    for u=1:4
                        for v=1:4
                            Axy((i-1+u-1)*(obj.m+2)+j:(i-1+u-1)*(obj.m+2)+j+3,(i-1+v-1)*(obj.m+2)+j:(i-1+v-1)*(obj.m+2)+j+3) = Axy((i-1+u-1)*(obj.m+2)+j:(i-1+u-1)*(obj.m+2)+j+3,(i-1+v-1)*(obj.m+2)+j:(i-1+v-1)*(obj.m+2)+j+3)+Aiquadr(u,v)*Aiquadr;
                        end
                    end
                end
            end
        end
        
        function [B,B_prime,B_primeprime] = getBasisfct(~,x,k,h,x0,derivatives) 
            % values of basisfunctions at x (shifted by previous knot xk,
            % stretched by step size h, cubical degree), (see Sec. IIIB in [1])
            % 'derivatives' = true - additionally provides derivated basisfunctions
                
            % subexpressions
            p0  = (x-(x0+h*(k-2)))/h;
            p1  = 1-p0;
            p2  = 2-p0;
            p3  = 1+p0;
            p4  = p0*p0;
            p5  = p1*p1;
            p6  = p2*p2;
            p7  = p3*p3;
            p8  = p4*p0;
            p9  = p5*p1;
            p10 = p6*p2;
            p11 = p7*p3;
            
            % basis functions
            B = [p9 p10-4*p9 p11-4*p8 p8]/6;
            if derivatives
                B_prime      = -[p5 p6-4*p5 4*p4-p7 -p4]/2/h;
                B_primeprime = [p1 p2-4*p1 p3-4*p0 p0]/h/h;
            else
                B_prime      = [];
                B_primeprime = [];
            end
        end
    end
    
end

