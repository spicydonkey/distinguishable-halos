function collated_data = collate_shots(cell_data)
% Collate all shots

collated_data=cell(1,size(cell_data,2));    % initilialise

for ii=1:length(collated_data)
    % vertical concatenation (each data point must be a 1xN vector) in
    % appropriate cell
    collated_data{ii}=vertcat(cell_data{:,ii});
end