function out = getcumdrive2(r1,r2,mode)

if mode == 1
    out = (max(0,r1)).^2 + r2.^2; % half-wave rectified r1 + full-wave rectified r2
elseif mode == 2
    out = r1.^2 + (max(0,r2)).^2; % full-wave rectified r1 + half-wave rectified r2 
elseif mode == 3
    out = r1.^2 + r2.^2; % full-wave rectified r1 + full wave rectified r2 
elseif mode == 4
    out = (max(0,r1)).^2 + (max(0,r2)).^2; % half-wave rectified r1 + half-wave rectified r2 
end

end


