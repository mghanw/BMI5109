"""Discovery-fold response evidence for intervention-assisted orientation."""

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_ind


def response_table(values,meta,anchors,qc,fold):
    """Continuous mean-response tests; holdout rows never select seed directions."""
    use=meta.dictionary_fold.eq(fold)
    rows=[]
    for condition in sorted(meta.condition.unique()):
        control=values.loc[use & meta.condition.eq(condition) & meta.target.eq("NTC")]
        for target in anchors:
            treated=values.loc[use & meta.condition.eq(condition) & meta.target.eq(target)]
            _,pvalues=ttest_ind(treated,control,equal_var=False,axis=0)
            estimate=treated.mean()-control.mean()
            scale=control.std().clip(lower=1e-8)
            quality=qc[(qc.condition==condition)&(qc.target==target)].iloc[0]
            for j,destination in enumerate(values):
                if destination==target:
                    continue
                rows.append(dict(fold=fold,condition=condition,source=target,destination=destination,
                                 n_target=len(treated),n_control=len(control),
                                 estimate=float(estimate[destination]),
                                 standardized_estimate=float(estimate[destination]/scale[destination]),
                                 p_value=float(pvalues[j]),test_valid=bool(np.isfinite(pvalues[j])),
                                 target_log_expression_change=float(estimate[target]),
                                 knockdown_remaining_ratio=float(quality.knockdown_remaining_ratio),
                                 knockdown_pass=bool(quality.eligible_both_conditions and estimate[target]<0)))
    frame=pd.DataFrame(rows)
    frame["q_value"]=false_discovery_control(frame.p_value.fillna(1).to_numpy(),method="bh")
    frame.loc[~frame.test_valid,"q_value"]=np.nan
    frame["eligible_seed"]=frame.test_valid & frame.knockdown_pass & frame.q_value.lt(.05)
    return frame


def evidence_for_condition(evidence,condition):
    rows=evidence[evidence.condition.eq(condition)&evidence.eligible_seed]
    return rows[["source","destination","estimate","standardized_estimate","q_value"]].to_dict("records")
