function [ErruL2,ErrpL2,uL2,pL2] = error_V2_Q2_k0(node,elem,uh,ph,info,pde,M)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s 

%K = pde.K; % coefficient matrix

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s


%% Get auxiliary data
% exact solution
ue = pde.zeta_exact;  pe = pde.phi_exact;
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);
div = info.div;

%% Compute L2 error
ErruL2 = 0; ErrpL2 = 0;
uL2 = 0; pL2 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    %xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel); 
    
    m11 = @(x,y) 1+0*x;
    m21 = @(x,y) 0+0*x;
    m12 = @(x,y) 0+0*x;
    m22 = @(x,y) 1+0*x;
    m1 = {m11, m21}; 
    m2 = {m12, m22};
    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration     
    am1 = @(x,y) 0+0*x;  
    am2 = @(x,y) 0+0*x; 
    for i = 1:2
        am1 = @(x,y) am1(x,y) + a(i)*m1{i}(x,y);
        am2 = @(x,y) am2(x,y) + a(i)*m2{i}(x,y);
    end
    am = @(x,y) [am1(x,y), am2(x,y)];
    damk = @(x,y) sum(div{iel}.*chi{iel});

    M1 = pde.M_1;

    errdiv = @(x,y) (pde.div_z([x,y])-damk(x,y)).^2;  
    errf = @(x,y) (M1([x,y])*(ue([x,y])-am(x,y))')'*(ue([x,y])-am(x,y))';
    errp = @(x,y) (pe([x,y])-ph(iel)).^2;
    u = @(x,y) (M1([x,y])*am(x,y)')'*am(x,y)';  
    p = @(x,y) ph(iel).^2;
    divu = @(x,y) damk(x,y).^2;
    % elementwise error
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    ErruL2 = ErruL2 + sum(integralTri(errf,4,nodeT,elemT)) + M*integralTri(errdiv,4,nodeT,elemT);%
    ErrpL2 = ErrpL2 + (1/M)*integralTri(errp,4,nodeT,elemT) + pde.theta*integralTri(errp,4,nodeT,elemT);
    uL2 = uL2 + sum(integralTri(u,4,nodeT,elemT)) + M*integralTri(divu,4,nodeT,elemT);
    pL2 = pL2 + (1/M)*integralTri(p,4,nodeT,elemT) + pde.theta*integralTri(p,4,nodeT,elemT);
end
h = mean(aux.diameter);
ErruL2 = sqrt(ErruL2);  
ErrpL2 = sqrt(ErrpL2);
uL2 = sqrt(uL2);
pL2 = sqrt(pL2);
