
%Одна итерация методом Ньютона
function model1 = Iteration(model,Y)
   discrepancy=CalculateDisperancy(model,Y);
   jacobian=CalculateJacobian(model,Y); 
        %Вектор изменения
    dx(:,1)=inv(jacobian)*(-discrepancy(:,1));
    dx(:,2)=discrepancy(:,2);
    dx(:,3)=discrepancy(:,3);
   %Присвоение новых значений переменных
   for i=1:size(dx,1)
        if dx(i,3)==1 % соответсвует невязке по реактивной мощности 
            model.BusTable(dx(i,2)).V=model.BusTable(dx(i,2)).V+dx(i,1); 
        end
        if dx(i,3)==0 % соответсвует невяке по активной мощности
            model.BusTable(dx(i,2)).D=model.BusTable(dx(i,2)).D+dx(i,1);
        end        
   end
  model1=model;
  maxDisperancy= 0;
  for i=1:size(discrepancy,1)
      if discrepancy(i,1) > maxDisperancy && discrepancy(i,3)==0;
          maxDisperancy=discrepancy(i,1);
      end
  end
  if maxDisperancy>0.1
      model1=Iteration(model1,Y);
  end

end


function disperancy1 = CalculateDisperancy (model,Y)
 %Вектор невязок
    for i=1:size(model.BusTable,1)
       discrepancy(i,1)=0;
       discrepancy(i,2)=i; %Пометка номера узла
       discrepancy(i,3)=1; %Пометка реактивной мощности
       discrepancy(i+size(model.BusTable,1),1)=0;
       discrepancy(i+size(model.BusTable,1),2)=i;%Пометка номера узла
       discrepancy(i+size(model.BusTable,1),3)=0;%Пометка активной мощности
           for j=1:size(model.BusTable,1)
           discrepancy(i+size(model.BusTable,1),1)=discrepancy(i+size(model.BusTable,1),1)+real(Y(i,j))*model.BusTable(i).V*model.BusTable(j).V*cos(model.BusTable(i).D-model.BusTable(j).D)+imag(Y(i,j))*model.BusTable(i).V*model.BusTable(j).V*sin(model.BusTable(i).D-model.BusTable(j).D);
           if model.BusTable(i).Type==3 
           discrepancy(i,1)=discrepancy(i,1)+imag(Y(i,j))*model.BusTable(i).V*model.BusTable(j).V*cos(model.BusTable(i).D-model.BusTable(j).D)-real(Y(i,j))*model.BusTable(i).V*model.BusTable(j).V*sin(model.BusTable(i).D-model.BusTable(j).D)
           end              
           end       
           discrepancy(i+size(model.BusTable,1),1)=discrepancy(i+size(model.BusTable,1),1)-model.BusTable(i).Pload*CalculatePSLC(model,i)+model.BusTable(i).Pgen;
           if model.BusTable(i).Type==3 
           discrepancy(i,1)=discrepancy(i,1)+model.BusTable(i).Qload*CalculateQSLC(model,i);
           end 
    end
    %Исключение базовых и генераторных узлов в векторе невязок
    x=1; deleted=0; i=1;  
    while (x<=size(model.BusTable,1)-deleted)
        if model.BusTable(i).Type==4
            discrepancy(x+size(model.BusTable,1),:)=[];
            discrepancy(x,:)=[];         
            deleted=deleted+1;
            x=x-1;
        end
        if model.BusTable(i).Type==0
            discrepancy(x,:)=[];           
            deleted=deleted+1;
            x=x-1;
        end 
        x=x+1;
        i=i+1;
    end
    disperancy1=discrepancy;
end

function jacobian1 = CalculateJacobian (model,Y)
    % Диагональные элементы матрицы Якоби
    for i=1:size(model.BusTable,1)
        jacobian(i,i+size(model.BusTable,1))=0;
         jacobian(i+size(model.BusTable,1),i+size(model.BusTable,1))=0;
        if model.BusTable(i).Type==3
            jacobian(i,i)=2*model.BusTable(i).V*imag(Y(i,i));  %dQ/dV
            jacobian(i+size(model.BusTable,1),i)=2*model.BusTable(i).V*real(Y(i,i)); %dP/dV
        end
        

        for j=1:size(model.BusTable,1)    
            if i~=j
                if model.BusTable(i).Type==3
                jacobian(i,i)=jacobian(i,i)+model.BusTable(j).V*(imag(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)-real(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)); %dQ/dV
                jacobian(i+size(model.BusTable,1),i)=jacobian(i+size(model.BusTable,1),i)+model.BusTable(j).V*(real(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)+imag(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)); %dP/dV
                jacobian(i,i+size(model.BusTable,1))=jacobian(i,i+size(model.BusTable,1))-model.BusTable(i).V*model.BusTable(j).V*(imag(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)+real(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)); %dQ/ddelta
                end   
                jacobian(i+size(model.BusTable,1),i+size(model.BusTable,1))=jacobian(i+size(model.BusTable,1),i+size(model.BusTable,1))-model.BusTable(i).V*model.BusTable(j).V*(real(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)-imag(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)); %dP/ddelta
            end
        end          
    end
    % Недиагональные элементы матрицы Якоби
    for i=1:size(model.BusTable,1)
        for j=1:size(model.BusTable,1)
                if i~=j
                jacobian(i,j)= model.BusTable(i).V*(imag(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)-real(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)); %dQ/dV nondiagonal
                jacobian(i+size(model.BusTable,1),j)= model.BusTable(i).V*(real(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)+imag(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)); %dP/dV nondiagonal
                jacobian(i,j+size(model.BusTable,1))= model.BusTable(i).V*model.BusTable(j).V*(imag(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)+real(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)); %dQ/ddelta nondiagonal
                jacobian(i+size(model.BusTable,1),j+size(model.BusTable,1))= model.BusTable(i).V*model.BusTable(j).V*(real(Y(i,j))*sin(model.BusTable(i).D-model.BusTable(j).D)-imag(Y(i,j))*cos(model.BusTable(i).D-model.BusTable(j).D)); %dP/ddelta nondiagonal
                
                end
        end
        
    end
    %Исключение базового и генерирующего узла в матрице Якоби
    x=1; deleted=0; i=1;
    while (x<=size(model.BusTable,1)-deleted)
        if model.BusTable(i).Type==4
            jacobian(x+size(model.BusTable,1),:)=[];
            jacobian(x,:)=[]; 
            jacobian(:,x+size(model.BusTable,1))=[];
            jacobian(:,x)=[]; 
            deleted=deleted+1;
            x=x-1;
        end
        if model.BusTable(i).Type==0
            jacobian(x,:)=[];          
            jacobian(:,x)=[]; 
            deleted=deleted+1;
            x=x-1;
        end 
        x=x+1;
        i=i+1;
    end
    jacobian1=jacobian;

end
%Расчёт коэф. для СХН по активной мощности
function kpslc = CalculatePSLC(model,i)
    if  model.BusTable(i).Plslc~=0;
        k=model.BusTable(i).V/model.BusTable(i).Unom;       
        kpslc=0.83-0.3*k+0.47*k^2;       
    else
        kpslc=1;
    end
end
%Расчёт коэф. для СХН по реактивной мощности
function kqslc = CalculateQSLC(model,i)
    if model.BusTable(i).Type==3 && model.BusTable(i).Plslc~=0;
        k=model.BusTable(i).V/model.BusTable(i).Unom;
        if k < 0.815 
            kqslc=0.721+0.158*k;
        elseif k >1.2 
            kqslc=1.49;
        else
            kqslc=3.7-7*k+4.3*k^2;
        end
    else
        kqslc=1;
    end
end
