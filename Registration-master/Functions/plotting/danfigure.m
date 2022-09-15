 
function num = danfigure(num)

try
    set(0,'currentfigure',num);
catch
    try
        figure(num,'position',[391 182 1485 796]);
    catch
        num = figure('position',[391 182 1485 796]);
    end
end