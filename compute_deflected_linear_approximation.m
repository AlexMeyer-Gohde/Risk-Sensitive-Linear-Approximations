function [deflect_,oo_, M_]=compute_deflected_linear_approximation(M_,options_,oo_,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute_deflected_linear_approximation.m
%
% This file produces risk-sensitive linear approximations. To calculate the
% approximation that is accurate to second order in the perturbation
% parameter, make sure you set order=3 when calling stoch_simul in your
% .mod file to have Dynare calculate the terms from a third order perturbation
%
% Call this file by placing the line:
% [deflect_,oo_, M_]=compute_deflected_linear_approximation(M_,options_,oo_,type);
% after the stoch_simul command in your Dynare .mod file
%
% The input type is a string equal to:
%
%  'stochastic_steady_state'  for the approximation around the approximated
%  stochastic steady state
%
% 'ergodic' for the approximation around the approximated ergodic mean
%
%
%THIS VERSION: Aug. 8, 2016
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,'ergodic')==0 && strcmp(type,'stochastic_steady_state')==0
    disp('Incompatible Type')
    deflect_=[];
    return
end
if options_.order>=3
 %% Recover the certainty equivalent decision rule   
    options_.order=1;
    if isempty(options_.qz_criterium)
        options_.qz_criterium = 1+1e-6;%not always set by dynare for higher order perturbations
    end
    exo_simul=oo_.exo_simul;oo_.exo_simul=zeros(M_.maximum_lag+2,M_.exo_nbr);% these two lines are needed to prevent dynare from using the simulation to evaluate the derivatives
    exo_det_simul=oo_.exo_det_simul;exo_det_simul=[];
    [numeric_version]=return_dynare_version(dynare_version);
    if numeric_version<4.3 %Starting in Dynare 4.3.0, dr1 is removed. 
        [oo.dr,waste1,waste2,waste3,waste4] = dr1(oo_.dr,0,M_,options_,oo_);
    elseif numeric_version>=4.3 %replace with stochastic_solvers in 4.3
        [oo.dr,waste1] =stochastic_solvers(oo_.dr,0,M_,options_,oo_);   
    else
        disp('Error, no certainty evuivalent dr calculated');
    end;
    oo_.dr.ghx=oo.dr.ghx;
    oo_.dr.ghu=oo.dr.ghu;
    oo_.exo_simul=exo_simul;oo_.exo_det_simul=exo_det_simul;
    options_.order=3;
end


if options_.order>=2
    [numeric_version]=return_dynare_version(dynare_version);
    if numeric_version >= 4.4 
        oo_.dr.nstatic = M_.nstatic;
        oo_.dr.npred = M_.nspred; % note M_.nspred = M_.npred+M_.nboth;
        oo_.dr.nboth = M_.nboth;
        oo_.dr.nfwrd = M_.nfwrd;
    end
    if isempty(options_.qz_criterium)==1
        options_.qz_criterium=1.000001;
    end
    [state_var, waste] =  lyapunov_symm(oo_.dr.ghx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:),oo_.dr.ghu(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:)*M_.Sigma_e*oo_.dr.ghu(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:)',2-options_.qz_criterium,options_.lyapunov_complex_threshold);
end

if options_.order>=1
    deflect_.y=oo_.dr.ys(oo_.dr.order_var);
    deflect_.y_x=oo_.dr.ghx;
    deflect_.y_u=oo_.dr.ghu;
end
if options_.order>=2
    if strcmp(type,'ergodic')==1
        x_2=(eye(oo_.dr.npred)-oo_.dr.ghx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:))\(oo_.dr.ghuu(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:)*vec(M_.Sigma_e)+oo_.dr.ghxx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:)*vec(state_var)+oo_.dr.ghs2(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:));
        y_2=[oo_.dr.ghx(1:oo_.dr.nstatic,:)*x_2+oo_.dr.ghuu(1:oo_.dr.nstatic,:)*vec(M_.Sigma_e)+oo_.dr.ghxx(1:oo_.dr.nstatic,:)*vec(state_var)+oo_.dr.ghs2(1:oo_.dr.nstatic,:);x_2;oo_.dr.ghx(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*x_2+oo_.dr.ghuu(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*vec(M_.Sigma_e)+oo_.dr.ghxx(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*vec(state_var)+oo_.dr.ghs2(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)];
        deflect_.y=deflect_.y+0.5*y_2;
    elseif strcmp(type,'stochastic_steady_state')==1
        x_2=(eye(oo_.dr.npred)-oo_.dr.ghx(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:))\(oo_.dr.ghs2(oo_.dr.nstatic+1:oo_.dr.nstatic+oo_.dr.npred,:));
        y_2=[oo_.dr.ghx(1:oo_.dr.nstatic,:)*x_2+oo_.dr.ghs2(1:oo_.dr.nstatic,:);x_2;oo_.dr.ghx(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)*x_2+oo_.dr.ghs2(oo_.dr.nstatic+oo_.dr.npred+1:M_.endo_nbr,:)];
        deflect_.x_2=x_2;
        deflect_.y=deflect_.y+0.5*y_2;
    end
end
if options_.order==3
        deflect_.y_x=deflect_.y_x+0.5*(oo_.dr.ghxx*kron(x_2,eye(oo_.dr.npred))+2*(oo_.dr.g_1(:,1:oo_.dr.npred)-oo_.dr.ghx));
        deflect_.y_u=deflect_.y_u+0.5*(oo_.dr.ghxu*kron(x_2,eye(M_.exo_nbr))+2*(oo_.dr.g_1(:,1+oo_.dr.npred:end)-oo_.dr.ghu));
end
deflect_.y=deflect_.y(oo_.dr.inv_order_var);
deflect_.y_y=[zeros(M_.endo_nbr,oo_.dr.nstatic),deflect_.y_x,zeros(M_.endo_nbr,oo_.dr.nfwrd)];
deflect_.y_y=deflect_.y_y(oo_.dr.inv_order_var,oo_.dr.inv_order_var);
deflect_.y_e=deflect_.y_u(oo_.dr.inv_order_var,:);
end
function [numeric_version]=return_dynare_version(dynare_version)
dynare_ver=dynare_version;
[dynare_ver_1,dynare_ver_remain] = strtok(dynare_ver, '.');
[dynare_ver_2,dynare_ver_remain2] = strtok(dynare_ver_remain, '.');
[dynare_ver_3,dynare_ver_remain3] = strtok(dynare_ver_remain2, '.');
numeric_version=str2num([dynare_ver_1 '.' dynare_ver_2 dynare_ver_3]);
end
