function b_resampled=time_sunc(a_time,b,b_time,order_b)
%assumes data a is taken exactly at the time indicated, integrated over the
%period before. 

if (nargin<3)
    order_b =0;
end

b = b-mean(b(:));
b_shape = size(b);
colons = repmat({':'},1,length(b_shape)-1);%seriously, matlab? This is the best you can do?
br_shape = b_shape;
br_shape(1)=length(a_time);
br_shape(end+1)=order_b+1;
b_resampled = zeros(br_shape);

if a_time(1)>0
    a_time = [0;a_time];
end

for q=1:(length(a_time)-1)
    this_slice = (b_time>a_time(q))&(b_time<a_time(q+1));
    datas = b(this_slice,colons{:});
    l_data = length(datas);
    for w=0:(order_b)
        basis = repmat(((1:l_data)').^w,[1,b_shape(2:end)]);
        basis = basis/sum(basis);
        weight=sum(datas.*basis);
        b_resampled(q,colons{:},w+1)=weight/sum(this_slice);
        datas = datas-weight*basis;
    end
end

if order_b == 0
    b_resampled=squeeze(b_resampled);
end
    