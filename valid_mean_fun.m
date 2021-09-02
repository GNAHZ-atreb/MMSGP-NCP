function [vm_iter, vm_cput, vm_err, vm_fe] = valid_mean_fun(iter, cput, err, fe, suc)

ind = find(suc); num_suc = sum(suc);

vm_iter = sum(iter(ind))/num_suc;

vm_cput = sum(cput(ind))/num_suc;

vm_err = sum(err(ind))/num_suc;

vm_fe = sum(fe(ind))/num_suc;

end