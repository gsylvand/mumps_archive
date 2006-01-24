function printmumpsstat(id)
%
% printmumpsstat(id)
% print mumps info
%

disp(['After analysis : Estimated operations                ' num2str(id.RINFO(1))]);
disp(['After analysis : Estimated space for factors         ' int2str(id.INFO(3))]);
disp(['After analysis : Estimated integer space             ' int2str(id.INFO(4))]);
disp(['After analysis : Estimated max front size            ' int2str(id.INFO(5))]);
disp(['After analysis : Number of node in the tree          ' int2str(id.INFO(6))]);
disp(['After analysis : Minimum value of MAXIS              ' int2str(id.INFO(7))]);
disp(['After analysis : Minimum value of MAXS               ' int2str(id.INFO(8))]);
disp(['After analysis : Estimated total size (Mbytes)       ' int2str(id.INFO(15))]);

disp(['After factorization : Assembly operations            ' num2str(id.RINFO(2))]);
disp(['After factorization : Elimination operations         ' num2str(id.RINFO(3))]);
disp(['After factorization : Real space to store LU         ' int2str(id.INFO(9))]);
disp(['After factorization : Integer space                  ' int2str(id.INFO(10))]);
disp(['After factorization : Largest front size             ' int2str(id.INFO(11))]);
disp(['After factorization : Number of off-diagonal pivots  ' int2str(id.INFO(12))]);
disp(['After factorization : Number of delayed pivots       ' int2str(id.INFO(13))]);
disp(['After factorization : Number of memory compresses    ' int2str(id.INFO(14))]);
disp(['After factorization : Total size needed (Mbytes)     ' int2str(id.INFO(16))]);
