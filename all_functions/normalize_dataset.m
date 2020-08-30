function datanew = normalize_dataset (data, mn, st)

datanew = (data - repmat(mn, size(data, 1), 1)) ./ repmat(st, size(data, 1), 1);

return