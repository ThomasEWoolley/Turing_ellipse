classdef Polys4< pde
    % diskpdeo2: alternative disk pdeo, generates radially symmetric mesh
    %  pdeo=diskpdeo2(R,nr,nphi)
    % R=radius, nr,nphi=#of discret.points in r and phi,
    % based on pts2dom meshing (points2domain, for octave)
    methods(Access = public)
        function o=Polys4(R,nr,nphi,hmax) % constructor


            o.grid=grid2D; r=linspace(0,R,nr); x=[]; y=[];
            for i=nr;
                for j=[0:nphi-1]
                    z=r(i)*exp(1i*j*2*pi/nphi);
                    x=[x real(z)];
                    y=[y imag(z)];

                    if rem(nphi,2)==1
%                         if j==(nphi+1)/2
                            x=[x real(r(i)*exp(1i*j*2*pi/nphi)+r(i)*exp(1i*(j+1)*2*pi/nphi))/2];
                            y=[y imag(r(i)*exp(1i*j*2*pi/nphi)+r(i)*exp(1i*(j+1)*2*pi/nphi))/2];
%                         end
                    end


                end
            end

% plot(x,y)

            o.grid.freeGeometry([x;y]); % use point2domain meshing
            h=1;
            while h > hmax;
                o.grid.refineMesh;
                h = h/2;
            end
            o.fem=lagrange12D;
        end
    end
    methods(Access = protected) % only here since required by pde-class
        function r=df(~,~,~); r=0; end % rather use p.fuha.sG in pderesi
    end
end
