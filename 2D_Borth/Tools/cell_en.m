
function U_en=cell_en(Un,en)
U_en=cell(1,length(Un));
for i=1:length(Un)
    U_en{1,i}=Un{1,i}(:,en);
end
end

