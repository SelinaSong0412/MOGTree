#' MOGTree: Multi-Objective Gain Trees for Personalized Treatment Decisions
#'
#' `MOGTree` implements the core estimation steps used in the manuscript
#' "Multi-Objective Gain Trees for Making Personalized Treatment Decision Based
#' on Treatment-Outcome Contrast Profile". The package fits tree-based patient
#' partitions that target heterogeneity in treatment-outcome contrast profiles
#' and provides helper functions for propensity score estimation, regression
#' outcome modeling, and AIPW-based summaries across a grid of outcome weights.
#'
#' @docType package
#' @name MOGTree
#' @import data.tree
#' @import dplyr
#' @import ggforce
#' @import ggpattern
#' @import ggplot2
#' @importFrom grDevices chull
#' @importFrom nnet multinom
#' @importFrom randomForest randomForest
#' @importFrom stats ave cov lm predict qchisq quantile rnorm setNames var
#' @importFrom stringr str_c
#' @importFrom utils capture.output combn globalVariables
"_PACKAGE"

utils::globalVariables(
  c(
    "Custom.color",
    "Group",
    "Treatment",
    "i",
    "j",
    "pattern_type",
    "textsize",
    "value",
    "weight",
    "x",
    "y"
  )
)

#' Additional documentation for `DTRtree()`
#'
#' @name DTRtree
#' @rdname DTRtree
#' @return `DTRtree()` returns a data frame describing the fitted binary tree.
#'   The columns are `node`, `X`, `cutoff`, `mEy`, and `group`, where `X` and
#'   `cutoff` record the selected splitting variable and cutoff for each
#'   internal node and `group` indexes terminal nodes.
#' @examples
#' set.seed(1)
#' n <- 80
#' H <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n)
#' )
#' A <- sample(1:3, n, replace = TRUE)
#' Y1 <- 1 + 0.6 * (H$X1 > 0) + c(0.8, 0.2, -0.3)[A] + rnorm(n, sd = 0.3)
#' Y2 <- 0.5 - 0.4 * (H$X1 > 0) + c(-0.1, 0.5, 0.7)[A] + rnorm(n, sd = 0.3)
#' w <- seq(0, 1, by = 0.25)
#'
#' fit <- DTRtree(
#'   Ys = cbind(Y1, Y2),
#'   A = A,
#'   H = H,
#'   depth = 2,
#'   minsplit = 10,
#'   w_vec = cbind(w, 1 - w),
#'   weight_combine = "mean",
#'   lambda = 0.01
#' )
#'
#' fit
NULL

#' Additional documentation for `predict_leaf.DTR()`
#'
#' @name predict_leaf.DTR
#' @rdname predict_leaf.DTR
#' @return `predict_leaf.DTR()` returns an integer vector of terminal node ids,
#'   one per row of `newdata`.
#' @examples
#' set.seed(1)
#' n <- 60
#' H <- data.frame(
#'   X1 = rnorm(n),
#'   X2 = rnorm(n)
#' )
#' A <- sample(1:3, n, replace = TRUE)
#' Y1 <- 1 + 0.5 * (H$X1 > 0) + c(0.9, 0.3, -0.2)[A] + rnorm(n, sd = 0.3)
#' Y2 <- 0.5 - 0.3 * (H$X1 > 0) + c(-0.2, 0.4, 0.8)[A] + rnorm(n, sd = 0.3)
#' w <- seq(0, 1, by = 0.5)
#'
#' fit <- DTRtree(
#'   Ys = cbind(Y1, Y2),
#'   A = A,
#'   H = H,
#'   depth = 2,
#'   minsplit = 10,
#'   w_vec = cbind(w, 1 - w),
#'   weight_combine = "mean",
#'   lambda = 0.01
#' )
#'
#' predict_leaf.DTR(fit, H[1:6, , drop = FALSE])
NULL

#' Additional documentation for `M.propen()`
#'
#' @name M.propen
#' @rdname M.propen
#' @return `M.propen()` returns a propensity score matrix with one column per
#'   treatment level and one row per subject.
#' @examples
#' set.seed(2)
#' X <- data.frame(X1 = rnorm(30), X2 = rnorm(30))
#' A <- sample(1:3, 30, replace = TRUE)
#' M.propen(A, X)
NULL

#' Additional documentation for `Reg.mu()`
#'
#' @name Reg.mu
#' @rdname Reg.mu
#' @return `Reg.mu()` returns a list with components `mus.reg`, the estimated
#'   counterfactual mean matrix for the last treatment stage, and `RegModel`,
#'   the fitted linear model.
#' @examples
#' set.seed(3)
#' X <- data.frame(X1 = rnorm(40), X2 = rnorm(40))
#' A <- sample(1:3, 40, replace = TRUE)
#' Y <- 0.5 + 0.8 * X$X1 + c(0.7, 0.2, -0.1)[A] + rnorm(40, sd = 0.4)
#' Reg.mu(Y, A, X)
NULL

#' Additional documentation for `mus.AIPW()`
#'
#' @name mus.AIPW
#' @rdname mus.AIPW
#' @return `mus.AIPW()` returns an AIPW-estimated counterfactual mean matrix
#'   with one column per treatment option.
#' @examples
#' set.seed(4)
#' X <- data.frame(X1 = rnorm(50), X2 = rnorm(50))
#' A <- sample(1:3, 50, replace = TRUE)
#' Y <- 1 + 0.5 * X$X1 + c(0.6, 0.1, -0.4)[A] + rnorm(50, sd = 0.4)
#' pi_hat <- M.propen(A, X)
#' mu_reg <- Reg.mu(Y, A, X)$mus.reg
#' mus.AIPW(Y, A, pi_hat, mu_reg)
NULL

#' Internal tree-construction helpers
#'
#' These functions support split evaluation and tree construction inside
#' `DTRtree()`. They are documented so the source package has complete manual
#' pages, but they are not part of the stable user-facing API.
#'
#' @name mogtree-tree-helpers
#' @aliases Benefit_gain combine.mat
#' @usage Benefit_gain(w_vec, mu.hat.selected, left)
#' @usage combine.mat(a, b)
#' @keywords internal
NULL

#' Internal plotting and tree-display helpers
#'
#' These functions are retained to support the analysis scripts in the
#' repository. They are not part of the stable user-facing API.
#'
#' @name mogtree-plot-helpers
#' @aliases bin2ten ten2bin vis_tree my_vis my_vis_pair my_vis_elipse my_vis_mean my_vis_mean_pair my_vis_all_pair
#' @usage ten2bin(x, max_h)
#' @usage bin2ten(vec)
#' @usage vis_tree(max_h, treeout, number)
#' @usage my_vis(mu1.hat, mu2.hat, chosen, k, label, lwd = 3)
#' @usage my_vis_pair(mu1.hat, mu2.hat, chosen1, chosen2, k, label, lwd = 3)
#' @usage my_vis_elipse(mu1.hat, mu2.hat, chosen, k, ref_k, label, lwd = 3, textsize)
#' @usage my_vis_mean(mu1.hat, mu2.hat, chosen, k, label, lwd = 3)
#' @usage my_vis_mean_pair(mu1.hat, mu2.hat, chosen1, chosen2, i, j, label, k,
#'   point_size = 3, my_magick = c("bricks", "horizontalsaw", "verticalsaw",
#'   "hs_horizontal", "hs_vertical"))
#' @usage my_vis_all_pair(mu1.hat, mu2.hat, leaf_pred, label, k, alpha = 1,
#'   my_magick = c("bricks", "horizontalsaw", "verticalsaw", "hs_horizontal",
#'   "hs_vertical"))
#' @details Several plotting helpers expect the calling script to define
#'   objects such as `Custom.color` and `textsize`.
#' @keywords internal
NULL
