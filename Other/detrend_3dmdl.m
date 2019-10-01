function out_mdl= detrend_3dmdl(mdl, options)
asGroup= false;
if nargin >1
    if ~isfield(options, 'as_group')
        asGroup = true;
    end
end

data= mdl.elem_data;
out_mdl= mdl;
sz_data= size(out_mdl.elem_data);
rows= sz_data(1);
ts= sz_data(2);
db= struct;
db.idx= [];
for i= 1:rows
    if ~isnan(data(i, :))
        db.idx= [db.idx; i];
    end
end

valid= length(db.idx);
db.means= zeros(valid, 1);
db.trends= zeros(valid, ts);
db.dts= zeros(valid, ts);
for i =1:valid
    pixel_idx= db.idx(i);
    temp= data(pixel_idx, :);
    if ~asGroup
        out_mdl.elem_data(pixel_idx, :)= (detrend(temp) + mean(temp))';
    else
        db.dts(i, :)= temp;
        db.means(i)= mean(temp);
        trend= temp- detrend(temp);
        db.trends(i, :)= trend;
    end
end

if asGroup % detrend based on mean of signal trends
    mean_trend= mean(db.trends);
    for i =1:valid
        pixel_idx= db.idx(i);
        temp= data(pixel_idx, :);
        out_mdl.elem_data(pixel_idx, :)= (temp(:)'- mean_trend) + db.means(i);
    end
end
    
end