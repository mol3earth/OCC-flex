function broken_axis = axisbreaker(axis)

Ay = axis(:,2);
SegRid.Idx = find(sign(diff(Ay)) > -1);
temp2 = [];
if isempty(SegRid.Idx)
    broken_axis = axis;
else
    if isempty(SegRid.Idx)
        SegRid.Idx = length(Ay);
    else
        SegRid.CrossOvers = Ay(SegRid.Idx-1);
        SegRid.Idx = [1 SegRid.Idx' length(Ay)];
        jk=0;
        for n = 1:(length(SegRid.CrossOvers)+1)
            temp1 = axis*NaN;
            temp1(SegRid.Idx(n)+jk:SegRid.Idx(n+1),:) = axis(SegRid.Idx(n)+jk:SegRid.Idx(n+1),:);        
            temp2 = [temp2  temp1];
            jk = 1;
        end
        broken_axis = temp2;
    end
end