function [gama1,gama2]=VanLaar(A12,A21,xA)
gama1 = exp(A12.*(A21.*(1-xA)./(A12.*xA + A21.*(1 - xA))).^2) ;
gama2 = exp(A21.*(A12.*xA./(A12.*xA + A21.*(1-xA))).^2) ;
end