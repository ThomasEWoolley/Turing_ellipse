classdef EllipsePDE< pde
    % diskpdeo2: alternative disk pdeo, generates radially symmetric mesh
    %  pdeo=diskpdeo2(R,nr,nphi)
    % R=radius, nr,nphi=#of discret.points in r and phi,
    % based on pts2dom meshing (points2domain, for octave)
    methods(Access = public)
        function o=EllipsePDE(A,B,nphi,hmax) % constructor
            a=sqrt(A^2-B^2);
            uc=atanh(B/A);

            o.grid=grid2D;  v=linspace(0,2*pi,nphi); v(end)=[];  x=[]; y=[];

                x=[x a*cosh(uc)*cos(v)];
                y=[y a*sinh(uc)*sin(v)];



            % plot(x,y)

            o.grid.freeGeometry([x;y]); % use point2domain meshing
            h=1;
            while h > hmax;
                o.grid.refineMesh;
                h = h/2;
            end
            o.fem=lagrange12D;

            % o.grid.pts2dom([x;y]); % use point2domain meshing
            % o.fem=lagrange12D;
        end
    end
    methods(Access = protected) % only here since required by pde-class
        function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
    end
end
