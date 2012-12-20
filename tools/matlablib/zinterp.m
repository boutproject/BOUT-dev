function result = zinterp(v, zind, nz)

v = squeeze(v);
if ( size(v, 1) ~= 1)
    v = v';
end
if ( size(zind, 1) ~= 1)
    zind = zind';
end

z0 = round(zind);
p = zind - z0;

z0 = z0 - (p<0);
p = p + (p<0);

z0 = mod(z0, nz-1);
z0 = z0 + (z0 == 0)*(nz-1);

zp = mod(z0+1, nz-1);
zp = zp + (zp==0)*(nz-1);

zm = mod(z0-1+nz, nz-1);
zm = zm + (zm==0)*(nz-1);

result = 0.5*p.*(p-1).*v(zm) + (1-p.*p).*v(z0) + 0.5*p.*(p+1).*v(zp);

end



