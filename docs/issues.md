# Known issues in ARETE

## PopPUNK

We have experienced issues with PopPUNK in ARETE runs, primarily related to the distances and clusters generated and how these affect both the **subsampling** and **recombination** subworkflows.

Sometimes the distances generated are too small or the number of clusters changes between executions of the same dataset.
While we can't solve the latter since it is a result of how PopPUNK itself is defined, the former can be mitigated by adjusting the subsampling thresholds with `--core_similarity` (99.9 by default) and `--accessory_similarity` (99 by default), or disabling subsampling altogether (`--enable_subsetting false`).

A useful course of action is to first run your dataset through the [PopPUNK entry](https://beiko-lab.github.io/arete/usage/#poppunk-entry), and then choose the appropriate parameters for your final pipeline run.
The PopPUNK entry takes at most 2 hours to run, even on [large datasets](https://beiko-lab.github.io/arete/resource_profiles/).
Disabling PopPUNK in your execution is also simple to do with `--skip_poppunk`.

## rSPR

rSPR, which you can enable with `--run_rspr` or with the [rSPR entry](https://beiko-lab.github.io/arete/usage/#rspr-entry), is known to be **very slow**, especially with larger datasets.

The default of 3 days for rSPR runtimes should be enough for some runs, but for most larger datasets it won't be sufficient. In this case, if you do want to run rSPR, we suggest two possible routes:

- Increasing the default time allocation for the `RSPR_EXACT` processes. [Check out how](https://beiko-lab.github.io/arete/usage/#custom-resource-requests).
- Ignoring timeout errors altogether and finish the pipeline execution with whatever finished running in these 3 days. **This is the default course of action for ARETE**.

By choosing the second course of action, we ignore timeout errors generated with `RSPR_EXACT` and finish the execution of downstream processes, i.e. `RSPR_HEATMAP`, with whatever results we already have.
This process generates a heatmap of Tree size and Exact rSPR distance.

While `RSPR_HEATMAP` should execute with the results that were generated up to the timeout, we have heard from users that this process can still not run, even when results from `RSPR_EXACT` were generated.
This is unfortunate but it shouldn't be a big problem: The only output given by `RSPR_HEATMAP` is the aforementioned heatmap, which can also be generated externally by using our [`rspr_heatmap.py`](https://github.com/beiko-lab/arete/blob/master/bin/rspr_heatmap.py) script or your own downstream analysis.
