
evaluation_metrics = function(est_hubs, true_hubs, node_names){
  
  if(sum(is.na(est_hubs)) == 1){
    est_hubs = as.character(est_hubs)
  }
  
  true_hubs = as.character(true_hubs)
  
  binary_hub_est = node_names %in% est_hubs
  
  binary_hub_true = node_names %in% true_hubs
  
  TP = sum(binary_hub_est*binary_hub_true)
  
  FP = sum(binary_hub_est*!binary_hub_true)
  
  FN = sum((!binary_hub_est)*binary_hub_true)
  
  TN = sum((!binary_hub_est)*(!binary_hub_true))
  
  TP = as.numeric(TP)
  FP = as.numeric(FP)
  FN = as.numeric(FN)
  TN = as.numeric(TN)
  
  FDR = FP/(FP + TP)
  
  if(FP + TP == 0) FDR = 0
  
  TPR = TP/(TP + FN) 
  
  # TPR = power = sensitivity = prop. of correctly est. hub nodes
  
  Precision = TP/(TP + FP) 
  
  # under null (no true positives, no false positives), 
  # precision might not be meaningful metric for a good
  # screening method
  
  FPR = FP/(FP + TN) # aka fallout
  
  B = sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  
  MCC = (TP*TN - FP*FN)/B
  
  if(B == 0) MCC = 0
  
  results = list(TP = TP, FP = FP, FN = FN, TN = TN,
                 FDR = FDR,
                 TPR = TPR,
                 Precision = Precision,
                 FPR = FPR,
                 MCC = MCC)
  
  return(results)
  
}
