%Функция расчёта режима
for i=1:size(model.BusTable,1)
    if model.BusTable(i).Vinst~=0
        model.BusTable(i).V=model.BusTable(i).Vinst;
    else
    model.BusTable(i).V=model.BusTable(i).Unom
    end
end
Y = zeros(size(model.BusTable,1));
Y=CreateYTable(model,Y);
model=Iteration(model,Y); 