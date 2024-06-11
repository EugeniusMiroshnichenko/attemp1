%Функция создания матрицы проводимостей
function Ytable = CreateYTable(model,Y)
  for i=1:size(model.LineTable,1) 
        st=model.LineTable(i).Nstartin;
        en=model.LineTable(i).Nendin;
        if model.LineTable(i).Type==1
            Y(st,st)=Y(st,st)-1/(model.LineTable(i).R+1i*model.LineTable(i).X)-model.LineTable(i).G-1i*model.LineTable(i).B; %Возможно не тот знак перед B
            Y(en,en)=Y(en,en)-1/((model.LineTable(i).R+1i*model.LineTable(i).X)*((model.LineTable(i).Ktm)^2));     
        else
        Y(st,st)=Y(st,st)-1/(model.LineTable(i).R+1i*model.LineTable(i).X)-model.LineTable(i).G/2-1i*model.LineTable(i).B/2;
        Y(en,en)=Y(en,en)-1/((model.LineTable(i).R+1i*model.LineTable(i).X)*((model.LineTable(i).Ktm)^2))-model.LineTable(i).G/2-1i*model.LineTable(i).B/2; %Возможно не тот знак перед B
        end
        Y(st,en)= 1/((model.LineTable(i).R+1i*model.LineTable(i).X)*model.LineTable(i).Ktm);
        Y(en,st)=Y(st,en);
  end
  Ytable=Y;
end