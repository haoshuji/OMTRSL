clear
load news20_length1_OMTFP_res.mat
AP_te_omtRSL_FP = APTeOMTFP;   AP_tr_omtRSL_FP = APTrOMTFP;
mAP_te_omtRSL_FP = mAPTeOMTFP; mAP_tr_omtRSL_FP = mAPTrOMTFP;
time_omtRSL_FP = timOMTFP;
save('../25-Jan-2016/news20_length1_omtRSL_FP_results.mat', 'AP_te_omtRSL_FP','mAP_te_omtRSL_FP',...
    'AP_tr_omtRSL_FP','mAP_tr_omtRSL_FP','time_omtRSL_FP');