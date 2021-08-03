function b = interpnan(b)
nanx = isnan(b);t = 1:numel(b);
b(nanx) = interp1(t(~nanx),b(~nanx),t(nanx));
end