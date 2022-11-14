function H_s = Hs_extract(H_grid,z_grid,theta_grid,theta,z)
    
%     theta_aux = rem(rem(theta,2*pi) + 2*pi,2*pi);
%     
%     theta_grid(:,end+1) = theta_grid(:,end) + theta_grid(:,2);
%     z_grid(:,end+1) = z_grid(:,end);
%     H_grid(:,end+1) = H_grid(:,1);
%     
%     H_s = interp2(theta_grid,z_grid,H_grid,theta_aux,z,'nearest',1e10);

H_map = H_grid(:);
z_map = z_grid(:);
theta_map = theta_grid(:);

H_s = zeros(size(theta));
for ii = 1:length(theta)
    theta_aux = rem(rem(theta(ii),2*pi) + 2*pi,2*pi);

    aux0 = find(theta_map<=theta_aux, 1, 'last' );

    aux1 = find(z_map<=z(ii));
    H_s(ii) = H_map(max(aux1(aux1<=aux0)));
    
end