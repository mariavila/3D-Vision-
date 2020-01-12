function [ fD1 , fS1 , fD2 , fS2] = Classify_Dyn_Stat_Ref( f1, f2)
    dif = CalcDistance(f1, f2)
    idx_dyn = find(dif>=130&dif<200)
    fD1 = f1(:, idx_dyn);
    fD2 = f2(:, idx_dyn);
    idx_stat = find(dif<130);
    fS1 = f1(:, idx_stat);
    fS2 = f2(:, idx_stat);
end

function euclideanDistance = CalcDistance(p1, p2)
s_p1 = size(p1);
euclideanDistance = zeros(s_p1(2),1);
for i = 1:s_p1(2)
    x1 = p1(1,i)/p1(3,i);
    x2 = p2(1,i)/p2(3,i);
    y1 = p1(2,i)/p1(3,i);
    y2 = p2(2,i)/p2(3,i);
    euclideanDistance(i) = sqrt((x2-x1)^2+(y2-y1)^2);
end
end