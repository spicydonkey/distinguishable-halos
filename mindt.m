% halo counts min dt analyser - diagnostics for detector ringing

function h=mindt(zxy)
vz=9.81*0.41;       % z-vel of atom at detection

% get minimum diff of detection time in all shots
min_dt=cellfun(@(x) min(diff(sortrows(x(:,1))))/vz,zxy,'UniformOutput',false);

% simple histogram
deadtime=100e-9;    % deadtime of TDC

bedge_dt=[0,linspace(deadtime,20e-6,100)];  % bin edge composed of deadtime zone and onwards
% bedge_dt=linspace(0,2000e-9,50);

hhist=figure();
hold on;
% for ii=1:size(min_dt,2)
%     [N{ii},edges{ii}]=histcounts(vertcat(min_dt{:,ii}),bedge_dt);
%     cents{ii}=(edges{ii}(1:end-1)+edges{ii}(2:end))/2;
% 
%     bar(cents{ii},N{ii},'LineStyle','none');
% end
for ii=1:size(min_dt,2)
    h{ii}=histogram(vertcat(min_dt{:,ii}),bedge_dt,'LineStyle','none');
end

box on;
axis tight;
xlabel('min $\Delta t$ in shot');
ylabel('no. of shots');
legend({'$m_F=0$','$m_F=1$'});

end
