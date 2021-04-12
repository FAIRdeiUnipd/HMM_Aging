function NII = create_4D_nii(REFnii,matrix)

REFnii = load_untouch_nii(REFnii);
REFnii.img = matrix;
REFnii.hdr.dime.dim(1) = 4;
REFnii.hdr.dime.dim(5) = size(matrix,4);
REFnii.hdr.dime.pixdim(5) = 0;
REFnii.hdr.dime.datatype = 16;
REFnii.hdr.dime.glmax = max(matrix(:));
REFnii.hdr.dime.glmin = min(matrix(:));
REFnii.hdr.dime.cal_max = max(matrix(:));
REFnii.hdr.dime.cal_min = min(matrix(:));

NII = REFnii;