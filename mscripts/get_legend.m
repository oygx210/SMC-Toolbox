function [leg] = get_legend(data_str)
%GET_LEGEND Summary of this function goes here
%   Construct the plot legend for an data array by adding the
%   incremental indices

data=evalin('base',data_str);
if strcmp(class(data),'timeseries')
    sz=size(data.Data,2);
else
    sz=size(data,2);
end
if sz>1
    leg=strcat(data_str,string(num2cell(1:sz)));
else
    leg={data_str};
end

end

