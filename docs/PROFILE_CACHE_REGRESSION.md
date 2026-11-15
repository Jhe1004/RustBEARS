# Profile cache regression and controlled discovery

The Psychotria BAYAREALIKE+J regression checks what happens **after a better constrained point is found**. It does not guarantee that an optimizer will discover every likelihood basin on every platform.

Windows CI run 43676080953 reproduced the native score of the recorded inferior basin but did not discover the better basin during the requested profile grid: `rebased_mle` remained false. That original failure is retained in the project CI evidence. The test revision does not claim to repair that optimizer discovery limitation, and changes no production optimizer or probability code.

The fixture uses the real 19-tip Psychotria data. Its inferior reference has `d=e=1e-12`, `j=0.10811626692419736`, and log likelihood about `-24.260680130280434`. The test first runs a genuine constrained optimization at the same fixed `d=1e-12`, initialized at the recorded better basin `e=5`, `j=0.11646904731717003`. It verifies convergence and independently recomputes the native likelihood of the returned parameters. The resulting score is about `-23.355689965854083`.

Only this already computed constrained result is replayed into the initial grid under the inferior reference. No score is fabricated. The production unrestricted rebase, subsequent constrained fits, endpoint searches, and independent endpoint re-evaluation all execute normally. The test preserves the `1e-4` endpoint residual requirements and verifies that the parameter box does not change. It additionally requires the shared upper-box ordinate to be evaluated under both the old and the new reference basin, so merely returning cached ordinates cannot satisfy the regression.

An independent in-memory mutation removed the reference parameters and score from the cache-generation identity, without editing production files. The revised test then reproduced the original failure exactly: endpoint deviance residual `0.0557264` exceeded `0.0001`. This confirms that controlled discovery still detects the historical stale-cache defect.

Local validation of the complete profile test file passed all 12 tests in 26.39 seconds on macOS. Windows validation of this revision remains a separate CI check; local success is not reported as Windows success. The analytical extinction/rebase tests continue to exercise the unmodified discovery path.
