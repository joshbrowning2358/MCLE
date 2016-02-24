Robust-Skew-t
=============

Repository to store research code for the robust skew-t fitting project.


While in development mode, I'll keep track of open tasks below:

- Compare speed of your new approach vs. the approach used with Mandy.
- Read up on the literature, particularly the papers Mandy sent.
- Define the scope of the project: what questions do we want to answer, what do we want to ignore?
- Determine when the non-robust model outperforms.  Is this a function of the "extremeness" of the outliers as well as the sample size?
- Write code (or find a package) for comparing two multivariate distributions.  This will be important when evaluating our estimates.  We could do this using an integral of the absolute value of the difference between the two densities.  Is there a way we could convert this into a one-dimensional integral?
- What effect is there when we use the wrong distribution (skew normal when the data is multivariate t?)
- Which constraint functions work best?  Can we determine some alternatives to the one currently in use?
- How can we optimally select k?  We currently have three approaches:
  * Flatness of the cdf
  * Minimum of the density
  * MAD (requires a bit of work to fully develop the idea).  However, maybe we can determine a theoretical relationship.
  * If a random number generator from the distribution is provided, as well as an expected number/percent of penalized observations, then we could simulate at the end to check if our parameters/calculated likelihoods lead to the expected percent of penalized observations.  If they do not, we could adjust k and refit.  However, ensuring convergence may be more challenging.
- Additional distributions: what else would be interesting?  Check with Mandy.
