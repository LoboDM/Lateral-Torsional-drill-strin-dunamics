function [ Map_Regimens,colour ]= Regimen(Whirl,StickSlip,Impact)
% Regimen  Evaluates combinations of regimes and separates a specific
%          color for each one.
% 
%   [ Map_Regimens,colour ]= Regimen(Whirl,StickSlip,Impact)
%   with Mapa_Regimes as the map with all regimes and cor as the color of 
%   each regime on the map.
% 
%  Inputs:
%  Precessao    -> Precession orientation flag                                   
%  StickSlip    -> Stick-slip occurrence flag
%  Impacto      -> Impact/contact occurrence flag
%
%  Outputs:
%  Mapa_Regimes -> Map with all regimes
%  cor          -> Color of each regime on the map
% 
%  LAST MODIFIED: 12/05/2020 BY LUCAS VOLPI, JORDAN BARBOZA AND DANIEL LOBO
%  CREATED BY LAVI (COPPE-UFRJ) FOR PETROBRAS


[m,n] = size(Whirl);
Map_Regimens = zeros(m,n);
 for j=1:n
    for i=1:m
        if Whirl(i,j) == 0 && StickSlip(i,j) == 0 %&& ...
                %Impact(i,j) == 1
        % Forward, no Stick-slip and no impact
        Map_Regimens(i,j) = 0;   
        elseif Whirl(i,j) == 0 && StickSlip(i,j) == 1 %&& ...
              %  Impact(i,j) == 1
        % Forward, Stick-slip and no impact
        Map_Regimens(i,j) = 1;
        elseif Whirl(i,j) == 1 && StickSlip(i,j) == 0 %&& ...
              %  Impact(i,j) == 0
        % Backward, no Stick-slip and permanent contact
        Map_Regimens(i,j) = 2;
        elseif Whirl(i,j) == 1 && StickSlip(i,j) == 1 %&& ...
             %   Impact(i,j) == 2
        % Backward, Stick-slip and impact
        Map_Regimens(i,j) = 3;
%         elseif Whirl(i,j) == 1 && StickSlip(i,j) == 0 && ...
%                 Impact(i,j) == 1
%         % Backward, no Stick-slip and no contact
%         Map_Regimens(i,j) = 4;
%         elseif Whirl(i,j) == 1 && StickSlip(i,j) == 1 && ...
%                 Impact(i,j) == 1
%         % Backward, Stick-slip and no contact
%         Map_Regimens(i,j) = 5;
%         elseif Whirl(i,j) == 1 && StickSlip(i,j) == 1 && ...
%                 Impact(i,j) == 0
%         % Backward, no Stick-slip and no contact
%         Map_Regimens(i,j) = 6;
        end
    end
end


Ind = max(max(Map_Regimens))+1;
colour = zeros(Ind,3);

for i = 1:(Ind)
    colour(i,:) = [(1 - i/Ind) (1-i/Ind) 1];
end

end