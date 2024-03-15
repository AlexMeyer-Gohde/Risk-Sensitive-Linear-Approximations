function [deflect_]=compute_deflected_moments(M_,deflect_,var_list_)

if cond(eye(M_.endo_nbr)-deflect_.y_y)<1e14
    deflect_.info_cond=0;
Gamma_1_spec=[];
Gamma_1_spec_hp=[];
Xcorr_spec=[];
Xcorr_spec_hp=[];
[waste, variable_select] = ismember(cellstr(var_list_), cellstr(M_.endo_names));

CORRELATION_SELECT=variable_select;
plots=0;
grid_size=2^10;
PERIOD=4;
HP_LAMBDA=1600*(PERIOD/4)^4; %Ravn, Morten and Harald Uhlig (2002) "On Adjusting the HP-Filter for the Frequency of Observations"
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;


s_Y=cell(1,grid_size);
s_Y_hp=cell(1,grid_size);
s_yy=zeros(M_.endo_nbr^2,grid_size);
s_yy_hp=zeros(M_.endo_nbr^2,grid_size);
Gamma_1_spec=cell(1,length(s_Y)+1);
Gamma_1_spec_hp=cell(1,length(s_Y)+1);
Xcorr_spec=cell(1,M_.endo_nbr);
Xcorr_spec_hp=cell(1,M_.endo_nbr);

hp = 4*HP_LAMBDA*(1 - cos(grid_points)).^2 ./ (1 + 4*HP_LAMBDA*(1 - cos(grid_points)).^2);



for n = 1 : grid_size
e_plus=exp( ((-1)^(1/2))*grid_points(n));
e_minus=exp( (-(-1)^(1/2))*grid_points(n));
s_Y{n}=((2*pi)^(-1))*inv(eye(M_.endo_nbr)-deflect_.y_y*e_minus)*deflect_.y_e*M_.Sigma_e*deflect_.y_e'*inv(eye(M_.endo_nbr)-deflect_.y_y'*e_plus);
s_Y_hp{n}=hp(n)^2*s_Y{n};
end



for n = 1 : grid_size
s_yy(:,n)=s_Y{n}(:);
s_yy_hp(:,n)=s_Y_hp{n}(:);
end


s_yy=real(ifft(conj(s_yy')))*2*pi;
s_Y=cell(1,grid_size);


s_yy_hp=real(ifft(conj(s_yy_hp')))*2*pi;
s_Y_hp=cell(1,grid_size);



for n = 1 : grid_size
s_temp = reshape(s_yy(n,:),M_.endo_nbr,M_.endo_nbr);
s_Y{n}= s_temp;
s_hp_temp = reshape(s_yy_hp(n,:),M_.endo_nbr,M_.endo_nbr);
s_Y_hp{n}=s_hp_temp;
end


for j=1:length(s_Y)+1
if j-length(s_Y)/2<=0
Gamma_1_spec{j}=s_Y{j+length(s_Y)/2};
Gamma_1_spec_hp{j}=s_Y_hp{j+length(s_Y)/2};
else
Gamma_1_spec{j}=s_Y{j-length(s_Y)/2};
Gamma_1_spec_hp{j}=s_Y_hp{j-length(s_Y)/2};
end
end

deflect_.standard_deviations=diag(Gamma_1_spec{length(s_Y)/2+1}).^(1/2);
deflect_.standard_deviations_hp=diag(Gamma_1_spec_hp{length(s_Y)/2+1}).^(1/2);
% 
% Standard_Deviations=[char(M_.endo_names(variable_select,:)),repmat(char(32),[length(variable_select),3]), num2str(deflect_.standard_deviations(variable_select),'% 0.5f')]
% disp(' ')
% 
% Standard_Deviations_HP=[char(M_.endo_names(variable_select,:)),repmat(char(32),[length(variable_select),3]), num2str(deflect_.standard_deviations_hp(variable_select),'% 0.5f')]
% disp(' ')
% 
% Means=[char(M_.endo_names(variable_select,:)),repmat(char(32),[length(variable_select),3]), num2str(deflect_.y(variable_select),'% 0.5f')]

[mm,nn]=size(deblank(char(M_.endo_names(variable_select,:))));
disp([repmat(char(32),[1,nn]),repmat(char(32),[1,10]),'Mean',repmat(char(32),[1,5]), 'Std Dev',repmat(char(32),[1,5]),'Std Dev HP'])
disp([char(deblank(M_.endo_names(variable_select,:))), repmat(char(32),[length(variable_select),7]), num2str(deflect_.y(variable_select),'% 0.10f'),...
    repmat(char(32),[length(variable_select),5]),num2str(deflect_.standard_deviations(variable_select),'% 0.5f'),repmat(char(32),...
    [length(variable_select),5]),num2str(deflect_.standard_deviations_hp(variable_select),'% 0.5f'),repmat(char(32),[length(variable_select),3])])




for j=1:length(Gamma_1_spec)
for i=1:M_.endo_nbr
Xcorr_spec{i}(:,j)=Gamma_1_spec{j}(:,i)./((Gamma_1_spec{length(s_Y)/2+1}(i,i).^(1/2))*(diag(Gamma_1_spec{length(s_Y)/2+1}).^(1/2)));
Xcorr_spec_hp{i}(:,j)=Gamma_1_spec_hp{j}(:,i)./((Gamma_1_spec_hp{length(s_Y)/2+1}(i,i).^(1/2))*(diag(Gamma_1_spec_hp{length(s_Y)/2+1}).^(1/2)));
end
end
% [mm,nn]=size(char(M_.endo_names(variable_select,:)));
% disp(['Correlations with ', char(M_.endo_names(2,:))])
% [repmat(char(32),[1,1+nn]),repmat(char(32),[1,2]),'j-3',repmat(char(32),[1,3]),'j-2',repmat(char(32),[1,3]),'j-1',repmat(char(32),[1,4]),'j',repmat(char(32),[1,4]),'j+1',repmat(char(32),[1,3]),'j+2',repmat(char(32),[1,3]),'j+3',repmat(char(32),[1,1]);
%    repmat(char(45),[1,1+nn+6*7]);
% char(M_.endo_names(variable_select,:)),repmat(char(32),[length(variable_select),1]),num2str(Xcorr_spec{}(variable_select),'% 0.3f')   ]
% [] 
% char(M_.endo_names(variable_select,:)),repmat(char(32),[length(variable_select),1]), num2str(deflect_.standard_deviations(variable_select),'% 0.3f')]
% disp(' ')
% Standard_Deviations_HP=[char(M_.endo_names(variable_select,:)),repmat(char(32),[length(variable_select),3]), num2str(deflect_.standard_deviations_hp(variable_select),'% 0.5f')]
% disp(' ')
%%%%%%Begin correlation plots%%%%%%%%%%%


if exist('CORRELATION_SELECT','var')==0
CORRELATION_SELECT=1:(M_.endo_nbr);
end
if exist('CORRELATION_HORIZON','var')==0
CORRELATION_HORIZON=20;
end

if plots==1
for j=1:length(CORRELATION_SELECT)
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(CORRELATION_SELECT)-1);0.25, 0.25,0.25])
    CROSS=plot(-CORRELATION_HORIZON:CORRELATION_HORIZON,Xcorr_spec{CORRELATION_SELECT(j)}(CORRELATION_SELECT,grid_size/2+1-CORRELATION_HORIZON:grid_size/2+1+CORRELATION_HORIZON)',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''Frequency-Domain Cross-Correlations of Selected Variables at t+j with %s at t'')',M_.endo_names(CORRELATION_SELECT(j),:)))
    legend(M_.endo_names(CORRELATION_SELECT,:),'Location','Best');
    ylabel('Correlation Coefficient');
    xlabel('j');
    axis([-Inf,Inf,-1,1])
    pause;
end
for j=1:length(CORRELATION_SELECT)
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(CORRELATION_SELECT)-1);0.25, 0.25,0.25])
    CROSS=plot(-CORRELATION_HORIZON:CORRELATION_HORIZON,Xcorr_spec_hp{CORRELATION_SELECT(j)}(CORRELATION_SELECT,grid_size/2+1-CORRELATION_HORIZON:grid_size/2+1+CORRELATION_HORIZON)',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''HP-Filtered Frequency-Domain Cross-Correlations of Selected Variables at t+j with %s at t'')',M_.endo_names(CORRELATION_SELECT(j),:)))
    legend(M_.endo_names(CORRELATION_SELECT,:),'Location','Best');
    ylabel('Correlation Coefficient');
    xlabel('j');
    axis([-Inf,Inf,-1,1])
    pause;
end
end
deflect_.Gamma_1_spec=Gamma_1_spec;
deflect_.Gamma_1_spec_hp=Gamma_1_spec_hp;
deflect_.Xcorr_spec=Xcorr_spec;
deflect_.Xcorr_spec_hp=Xcorr_spec_hp;
else
    deflect_.info_cond=1;
end