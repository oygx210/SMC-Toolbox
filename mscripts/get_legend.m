function [leg] = get_legend(data_str,indices)
%GET_LEGEND Summary of this function goes here
%   Construct the plot legend for an data array by adding the
%   incremental indices

data=evalin('base',data_str);

if strcmp(class(data),'timeseries')
    sz=size(data.Data,2);
else
    sz=size(data,2);
end

if nargin>1
    if max(indices)>sz
        error(['max index ',num2str(max(indices)),' is greater than size ',num2str(sz)]);
    else
        idx=indices;
    end
else
    idx=1:sz;
end

if sz>1
    leg=strcat(data_str,string(num2cell(idx)));
else
    leg={data_str};
end

end

