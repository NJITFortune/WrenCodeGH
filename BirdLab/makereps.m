for i = 1:length(fem14jan13.syls);
    
fem14jan13.e1dat{i} = fem14jan13.e1spikes(find(fem14jan13.e1spikes > fem14jan13.syls(i,1) & fem14jan13.e1spikes < fem14jan13.syls(i,2))) - fem14jan13.syls(i,1);
fem14jan13.e2dat{i} = fem14jan13.e2spikes(find(fem14jan13.e2spikes > fem14jan13.syls(i,1) & fem14jan13.e2spikes < fem14jan13.syls(i,2))) - fem14jan13.syls(i,1);
fem14jan13.e3dat{i} = fem14jan13.e3spikes(find(fem14jan13.e3spikes > fem14jan13.syls(i,1) & fem14jan13.e3spikes < fem14jan13.syls(i,2))) - fem14jan13.syls(i,1);
fem14jan13.e4dat{i} = fem14jan13.e4spikes(find(fem14jan13.e4spikes > fem14jan13.syls(i,1) & fem14jan13.e4spikes < fem14jan13.syls(i,2))) - fem14jan13.syls(i,1);
    
end;   
