Nell=readFileNifti('t1_1mm.nii.gz')
mrVista 3

NellNumRightCenter=Nell;
NellNumRightCenter.data=zeros(size(Nell.data));
NellNumRightCenter.fname='NellNumRightCenter.nii';
m=2; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); NellNumRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(NellNumRightCenter)

NellNumLeftCenter=Nell;
NellNumLeftCenter.data=zeros(size(Nell.data));
NellNumLeftCenter.fname='NellNumLeftCenter.nii';
m=1; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); NellNumLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(NellNumLeftCenter)

NellSizeRightCenter=Nell;
NellSizeRightCenter.data=zeros(size(Nell.data));
NellSizeRightCenter.fname='NellSizeRightCenter.nii';
m=4; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); NellSizeRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(NellSizeRightCenter)

NellSizeLeftCenter=Nell;
NellSizeLeftCenter.data=zeros(size(Nell.data));
NellSizeLeftCenter.fname='NellSizeLeftCenter.nii';
m=3; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); NellSizeLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(NellSizeLeftCenter)


%do some things in MNI tools
!nii2mnc NellNumLeftCenter.nii
!nii2mnc NellNumRightCenter.nii
!nii2mnc NellSizeLeftCenter.nii
!nii2mnc NellSizeRightCenter.nii

!mincresample -like model.mnc -transformation nell3.xfm NellNumLeftCenter.mnc NellNumLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation nell3.xfm NellNumRightCenter.mnc NellNumRightCenterResampled.mnc
!mincresample -like model.mnc -transformation nell3.xfm NellSizeLeftCenter.mnc NellSizeLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation nell3.xfm NellSizeRightCenter.mnc NellSizeRightCenterResampled.mnc

!mnc2nii -nii NellNumLeftCenterResampled.mnc
!mnc2nii -nii NellNumRightCenterResampled.mnc
!mnc2nii -nii NellSizeLeftCenterResampled.mnc
!mnc2nii -nii NellSizeRightCenterResampled.mnc


%Back to matlab
NellNumLeftCenterResampled=readFileNifti('NellNumLeftCenterResampled.nii')
NellNumRightCenterResampled=readFileNifti('NellNumRightCenterResampled.nii')
NellSizeLeftCenterResampled=readFileNifti('NellSizeLeftCenterResampled.nii')
NellSizeRightCenterResampled=readFileNifti('NellSizeRightCenterResampled.nii')

[tmp, ii]=max(NellNumLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(NellNumLeftCenterResampled.data), (ii));
NellNumLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(NellNumRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(NellNumRightCenterResampled.data), (ii));
NellNumRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(NellSizeLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(NellSizeLeftCenterResampled.data), (ii));
NellSizeLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(NellSizeRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(NellSizeRightCenterResampled.data), (ii));
NellSizeRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]


Dumo=readFileNifti('t1_1mm.nii.gz')
mrVista 3

DumoNumRightCenter=Dumo;
DumoNumRightCenter.data=zeros(size(Dumo.data));
DumoNumRightCenter.fname='DumoNumRightCenter.nii';
m=2; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); DumoNumRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(DumoNumRightCenter)

DumoNumLeftCenter=Dumo;
DumoNumLeftCenter.data=zeros(size(Dumo.data));
DumoNumLeftCenter.fname='DumoNumLeftCenter.nii';
m=1; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); DumoNumLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(DumoNumLeftCenter)

DumoSizeRightCenter=Dumo;
DumoSizeRightCenter.data=zeros(size(Dumo.data));
DumoSizeRightCenter.fname='DumoSizeRightCenter.nii';
m=4; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); DumoSizeRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(DumoSizeRightCenter)

DumoSizeLeftCenter=Dumo;
DumoSizeLeftCenter.data=zeros(size(Dumo.data));
DumoSizeLeftCenter.fname='DumoSizeLeftCenter.nii';
m=3; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); DumoSizeLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(DumoSizeLeftCenter)


%do some things in MNI tools
!nii2mnc DumoNumLeftCenter.nii
!nii2mnc DumoNumRightCenter.nii
!nii2mnc DumoSizeLeftCenter.nii
!nii2mnc DumoSizeRightCenter.nii

!gunzip t1_1mm.nii.gz
!nii2mnc t1_1mm.nii
!gzip t1_1mm.nii
!mritotal t1_1mm.mnc Dumo.xfm

!mincresample -like model.mnc -transformation Dumo.xfm DumoNumLeftCenter.mnc DumoNumLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Dumo.xfm DumoNumRightCenter.mnc DumoNumRightCenterResampled.mnc
!mincresample -like model.mnc -transformation Dumo.xfm DumoSizeLeftCenter.mnc DumoSizeLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Dumo.xfm DumoSizeRightCenter.mnc DumoSizeRightCenterResampled.mnc

!mnc2nii -nii DumoNumLeftCenterResampled.mnc
!mnc2nii -nii DumoNumRightCenterResampled.mnc
!mnc2nii -nii DumoSizeLeftCenterResampled.mnc
!mnc2nii -nii DumoSizeRightCenterResampled.mnc


%Back to matlab
DumoNumLeftCenterResampled=readFileNifti('DumoNumLeftCenterResampled.nii');
DumoNumRightCenterResampled=readFileNifti('DumoNumRightCenterResampled.nii');
DumoSizeLeftCenterResampled=readFileNifti('DumoSizeLeftCenterResampled.nii');
DumoSizeRightCenterResampled=readFileNifti('DumoSizeRightCenterResampled.nii');

[tmp, ii]=max(DumoNumLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(DumoNumLeftCenterResampled.data), (ii));
DumoNumLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(DumoNumRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(DumoNumRightCenterResampled.data), (ii));
DumoNumRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(DumoSizeLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(DumoSizeLeftCenterResampled.data), (ii));
DumoSizeLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(DumoSizeRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(DumoSizeRightCenterResampled.data), (ii));
DumoSizeRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]



Harv=readFileNifti('t1_1mm.nii.gz')
mrVista 3

HarvNumRightCenter=Harv;
HarvNumRightCenter.data=zeros(size(Harv.data));
HarvNumRightCenter.fname='HarvNumRightCenter.nii';
m=2; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); HarvNumRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(HarvNumRightCenter)

HarvNumLeftCenter=Harv;
HarvNumLeftCenter.data=zeros(size(Harv.data));
HarvNumLeftCenter.fname='HarvNumLeftCenter.nii';
m=1; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); HarvNumLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(HarvNumLeftCenter)

HarvSizeRightCenter=Harv;
HarvSizeRightCenter.data=zeros(size(Harv.data));
HarvSizeRightCenter.fname='HarvSizeRightCenter.nii';
m=4; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); HarvSizeRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(HarvSizeRightCenter)

HarvSizeLeftCenter=Harv;
HarvSizeLeftCenter.data=zeros(size(Harv.data));
HarvSizeLeftCenter.fname='HarvSizeLeftCenter.nii';
m=3; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); HarvSizeLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(HarvSizeLeftCenter)


%do some things in MNI tools
!nii2mnc HarvNumLeftCenter.nii
!nii2mnc HarvNumRightCenter.nii
!nii2mnc HarvSizeLeftCenter.nii
!nii2mnc HarvSizeRightCenter.nii

!gunzip t1_1mm.nii.gz
!nii2mnc t1_1mm.nii
!gzip t1_1mm.nii
!mritotal t1_1mm.mnc Harv.xfm

!mincresample -like model.mnc -transformation Harv.xfm HarvNumLeftCenter.mnc HarvNumLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Harv.xfm HarvNumRightCenter.mnc HarvNumRightCenterResampled.mnc
!mincresample -like model.mnc -transformation Harv.xfm HarvSizeLeftCenter.mnc HarvSizeLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Harv.xfm HarvSizeRightCenter.mnc HarvSizeRightCenterResampled.mnc

!mnc2nii -nii HarvNumLeftCenterResampled.mnc
!mnc2nii -nii HarvNumRightCenterResampled.mnc
!mnc2nii -nii HarvSizeLeftCenterResampled.mnc
!mnc2nii -nii HarvSizeRightCenterResampled.mnc


%Back to matlab
HarvNumLeftCenterResampled=readFileNifti('HarvNumLeftCenterResampled.nii');
HarvNumRightCenterResampled=readFileNifti('HarvNumRightCenterResampled.nii');
HarvSizeLeftCenterResampled=readFileNifti('HarvSizeLeftCenterResampled.nii');
HarvSizeRightCenterResampled=readFileNifti('HarvSizeRightCenterResampled.nii');

[tmp, ii]=max(HarvNumLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(HarvNumLeftCenterResampled.data), (ii));
HarvNumLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(HarvNumRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(HarvNumRightCenterResampled.data), (ii));
HarvNumRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(HarvSizeLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(HarvSizeLeftCenterResampled.data), (ii));
HarvSizeLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(HarvSizeRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(HarvSizeRightCenterResampled.data), (ii));
HarvSizeRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]


Klei=readFileNifti('t1_1mm.nii.gz')
mrVista 3

KleiNumRightCenter=Klei;
KleiNumRightCenter.data=zeros(size(Klei.data));
KleiNumRightCenter.fname='KleiNumRightCenter.nii';
m=2; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); KleiNumRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(KleiNumRightCenter)

KleiNumLeftCenter=Klei;
KleiNumLeftCenter.data=zeros(size(Klei.data));
KleiNumLeftCenter.fname='KleiNumLeftCenter.nii';
m=1; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); KleiNumLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(KleiNumLeftCenter)

KleiSizeRightCenter=Klei;
KleiSizeRightCenter.data=zeros(size(Klei.data));
KleiSizeRightCenter.fname='KleiSizeRightCenter.nii';
m=4; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); KleiSizeRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(KleiSizeRightCenter)

KleiSizeLeftCenter=Klei;
KleiSizeLeftCenter.data=zeros(size(Klei.data));
KleiSizeLeftCenter.fname='KleiSizeLeftCenter.nii';
m=3; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); KleiSizeLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(KleiSizeLeftCenter)


%do some things in MNI tools
!nii2mnc KleiNumLeftCenter.nii
!nii2mnc KleiNumRightCenter.nii
!nii2mnc KleiSizeLeftCenter.nii
!nii2mnc KleiSizeRightCenter.nii

!gunzip t1_1mm.nii.gz
!nii2mnc t1_1mm.nii
!gzip t1_1mm.nii
!mritotal t1_1mm.mnc Klei.xfm
!mincresample -like model.mnc -transformation Klei.xfm t1_1mm.mnc t1_1mm_Resampled.mnc

!mincresample -like model.mnc -transformation Klei.xfm KleiNumLeftCenter.mnc KleiNumLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Klei.xfm KleiNumRightCenter.mnc KleiNumRightCenterResampled.mnc
!mincresample -like model.mnc -transformation Klei.xfm KleiSizeLeftCenter.mnc KleiSizeLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Klei.xfm KleiSizeRightCenter.mnc KleiSizeRightCenterResampled.mnc

!mnc2nii -nii KleiNumLeftCenterResampled.mnc
!mnc2nii -nii KleiNumRightCenterResampled.mnc
!mnc2nii -nii KleiSizeLeftCenterResampled.mnc
!mnc2nii -nii KleiSizeRightCenterResampled.mnc


%Back to matlab
KleiNumLeftCenterResampled=readFileNifti('KleiNumLeftCenterResampled.nii');
KleiNumRightCenterResampled=readFileNifti('KleiNumRightCenterResampled.nii');
KleiSizeLeftCenterResampled=readFileNifti('KleiSizeLeftCenterResampled.nii');
KleiSizeRightCenterResampled=readFileNifti('KleiSizeRightCenterResampled.nii');

[tmp, ii]=max(KleiNumLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(KleiNumLeftCenterResampled.data), (ii));
KleiNumLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(KleiNumRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(KleiNumRightCenterResampled.data), (ii));
KleiNumRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(KleiSizeLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(KleiSizeLeftCenterResampled.data), (ii));
KleiSizeLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(KleiSizeRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(KleiSizeRightCenterResampled.data), (ii));
KleiSizeRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]


Frac=readFileNifti('t1_1mm.nii.gz')
mrVista 3

FracNumRightCenter=Frac;
FracNumRightCenter.data=zeros(size(Frac.data));
FracNumRightCenter.fname='FracNumRightCenter.nii';
m=2; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); FracNumRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(FracNumRightCenter)

FracNumLeftCenter=Frac;
FracNumLeftCenter.data=zeros(size(Frac.data));
FracNumLeftCenter.fname='FracNumLeftCenter.nii';
m=1; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); FracNumLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(FracNumLeftCenter)

FracSizeRightCenter=Frac;
FracSizeRightCenter.data=zeros(size(Frac.data));
FracSizeRightCenter.fname='FracSizeRightCenter.nii';
m=4; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); FracSizeRightCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(FracSizeRightCenter)

FracSizeLeftCenter=Frac;
FracSizeLeftCenter.data=zeros(size(Frac.data));
FracSizeLeftCenter.fname='FracSizeLeftCenter.nii';
m=3; for n=1:size(VOLUME{1}.ROIs(m).coords, 2); FracSizeLeftCenter.data(VOLUME{1}.ROIs(m).coords(3,n),217-VOLUME{1}.ROIs(m).coords(2,n),181-VOLUME{1}.ROIs(m).coords(1,n))=1; end
writeFileNifti(FracSizeLeftCenter)


%do some things in MNI tools
!nii2mnc FracNumLeftCenter.nii
!nii2mnc FracNumRightCenter.nii
!nii2mnc FracSizeLeftCenter.nii
!nii2mnc FracSizeRightCenter.nii

!gunzip t1_1mm.nii.gz
!nii2mnc t1_1mm.nii
!gzip t1_1mm.nii
!mritotal t1_1mm.mnc Frac.xfm
!mincresample -like model.mnc -transformation Frac2.xfm t1_1mm.mnc t1_1mm_Resampled.mnc

!mincresample -like model.mnc -transformation Frac2.xfm FracNumLeftCenter.mnc FracNumLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Frac2.xfm FracNumRightCenter.mnc FracNumRightCenterResampled.mnc
!mincresample -like model.mnc -transformation Frac2.xfm FracSizeLeftCenter.mnc FracSizeLeftCenterResampled.mnc
!mincresample -like model.mnc -transformation Frac2.xfm FracSizeRightCenter.mnc FracSizeRightCenterResampled.mnc

!mnc2nii -nii FracNumLeftCenterResampled.mnc
!mnc2nii -nii FracNumRightCenterResampled.mnc
!mnc2nii -nii FracSizeLeftCenterResampled.mnc
!mnc2nii -nii FracSizeRightCenterResampled.mnc


%Back to matlab
FracNumLeftCenterResampled=readFileNifti('FracNumLeftCenterResampled.nii');
FracNumRightCenterResampled=readFileNifti('FracNumRightCenterResampled.nii');
FracSizeLeftCenterResampled=readFileNifti('FracSizeLeftCenterResampled.nii');
FracSizeRightCenterResampled=readFileNifti('FracSizeRightCenterResampled.nii');

[tmp, ii]=max(FracNumLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(FracNumLeftCenterResampled.data), (ii));
FracNumLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(FracNumRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(FracNumRightCenterResampled.data), (ii));
FracNumRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(FracSizeLeftCenterResampled.data(:));;
[x,y,z]=ind2sub(size(FracSizeLeftCenterResampled.data), (ii));
FracSizeLeftCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

[tmp, ii]=max(FracSizeRightCenterResampled.data(:));;
[x,y,z]=ind2sub(size(FracSizeRightCenterResampled.data), (ii));
FracSizeRightCenterResampled.sto_xyz(1:3,4)'+[x,y,z]

