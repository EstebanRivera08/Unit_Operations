function f = Pvap(A,B,C,T)
    f = 10.^(A - B./(T + C)) ;
end