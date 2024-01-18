function E = elas_matrix(Elasmm)
E = zeros(2,2,12);
for i = 1:12
    E(1,:,i) =  Elasmm{1,i};
    E(2,:,i) =  Elasmm{2,i};
end
end