function [deflect_]=compute_deflected_irf(M_,options_,deflect_,var_list_,plots)
if options_.irf==0
    disp('The option irf has been set to zero. To see IRFS, please set irf to a positive integer')
    return
end
if isempty(var_list_)==1
    var_list_=M_.endo_names;
end
[waste, variable_select] = ismember(cellstr(var_list_), cellstr(M_.endo_names));
deflect_irf{1}=deflect_.y_e;
%;%linear irf at impact
for i=2:options_.irf%iterating through the recursion for the linear irf
   deflect_irf{i}=deflect_.y_y*deflect_irf{i-1};
end
temp = cat(3,deflect_irf{:});
for jj=1:M_.exo_nbr
deflect_.irf((jj-1)*M_.endo_nbr+[1:M_.endo_nbr],:)=squeeze(temp(:,jj,:));
end

if plots==1
TIME=(-1:options_.irf-1);
for j=1:M_.exo_nbr
    figure;
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(variable_select)-1);0.25, 0.25,0.25])
    plot(TIME,[zeros(length(variable_select),1) deflect_.irf((j-1)*(M_.endo_nbr)+variable_select,:)*(M_.Sigma_e(j,j))^(1/2)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''Impulse Responses to a Shock in %s'')',M_.exo_names(j,:)))
    legend(M_.endo_names(variable_select,:),'Location','Best');
    ylabel('% Deviations from Steady-State');
    xlabel('Years since Shock Realization');
    hold on
    plot(TIME, 0*TIME,'k')
    hold off
%     figure;
%     clf('reset')
%     for i=1:length(variable_select)
%         subplot(ceil(length(variable_select)^(1/2)),round(length(variable_select)^(1/2)),i); plot(TIME, 0*TIME,'k', TIME, [zeros(1,1) deflect_.irf((j-1)*(M_.endo_nbr)+variable_select(i),:)*(M_.Sigma_e(j,j))^(1/2)],'k:.','MarkerEdgeColor','auto','MarkerSize',8);
%         legend(M_.endo_names(variable_select(i),:),'Location','Best')
%         ylabel('% Deviations from Steady-State');
%         xlabel('Years since Shock Realization');
%         eval(sprintf('title({''Impulse Response to a Shock in %s''})',M_.exo_names(j,:)))
%     end
end
end
end