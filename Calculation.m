%Функция расчёта режима. Не реализована заданная точность расчёта уст.
%режима, исходная привязка напряжений V к Unom, СХН
Y = zeros(size(model.BusTable,1));
Y=CreateYTable(model,Y);
model=Iteration(model,Y); 