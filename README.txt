Date released 25/08/2023.
Matlab functions for the efficient privacy-preserving camera model identification.

-------------------------------------------------------------------
 Copyright
-------------------------------------------------------------------

Copyright (c) University of Technology of Troyes and Hangzhou Dianzi University.
All rights reserved.
This work should only be used for nonprofit purposes.

By downloading and/or using any of these files, you implicitly agree to all the
terms of the license, as specified in the document LICENSE.txt

-------------------------------------------------------------------
 Contents
-------------------------------------------------------------------

(1) Camera fingerprint estimation.
   In this package, we provide the code of camera fingerprint estimation. 
   See "cal_parameters.m" for more details.
   The detailed techniques have been described in [1]. 
   

(2) Camera model identification. 
   In the package, we provide the code of camera model identification, which is committed to verify whether an inquiry image is from the known camera.
   See "GLRT.m" for more details.
   The detailed techniques have been described in [2]. 
     

[1] Chen Y, Retraint F, Qiao T,
''Image splicing forgery detection using simplified generalized noise model''
Signal Processing: Image Communication, 2022.

[2] Chen Y, Qiao T, Retraint F, and Hu G, 
''Efficient Privacy-Preserving Forensic Method for Camera Model Identification''
IEEE Transactions on Information Forensics and Security, 2022.

Note: To benefit the users understanding , we also propode some sample images from the benchmark dataset Dresden.  
-------------------------------------------------------------------
 Requirements
-------------------------------------------------------------------

All the functions and scripts were tested on MATLAB R2016b,
the operation is not guaranteed with older version of MATLAB.

-------------------------------------------------------------------
 Feedback
-------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact Tong Qiao at tong.qiao@hdu.edu.cn
For other information you can see https://q-bigbig.github.io/pub.html
