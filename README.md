# SAA 2025 Settlement Size

<!-- badges: start -->
[![Static Badge](https://img.shields.io/badge/Quarto-Presentation-74AADB?style=social&logo=Quarto&logoColor=74AADB&logoSize=auto&labelColor=f8f9fa)](https://quarto.org/docs/presentations/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
<!-- badges: end -->

Repository for my presentation at the Society for American Archaeology 90th
Annual Meeting in Denver. Slides made using [quarto](https://quarto.org/)
and R.

**Slides:** <https://kbvernon.github.io/saa_2025-bnn_dating/>

You can press "F" to enter Full Screen mode and "S" to enter Speaker view. The
hamburger button ☰ at the bottom left of the screen opens a menu that allows you
to easily navigate the presentation.

**Session:** "Multi-Scale Data and the History of Human Development in the US
Southwest"

**Title:** "A Bayesian neural network for relative dating in archaeology"

**Authors:** K. Blake Vernon, Scott Ortman

**Abstract**

The two most powerful forces driving long-term development in human societies
are climate change and demography, so it should come as no surprise that
researchers working at the bleeding edge of archaeological methods are laser
focused on finding new techniques to estimate key climate and demographic
quantities in the past. To aid in these efforts, specifically demographic
reconstruction, this project proposes a new spatio-temporal deep learning (DL)
framework for the relative-dating of human settlements. The core of this new
framework is a Bayesian Neural Network (BNN) model trained on multi-site
archaeological tree-ring date and ceramic assemblage data. In this context, a
BNN has two important advantages: (i) unlike other DL models, it is less
susceptible to overfitting, and (ii) unlike other archaeological cross-dating
methods, it provides a straightforward measure of uncertainty. We also show how
to incorporate a Moran Eigenvector Map into a BNN to account for spatial
autocorrelation. Training and testing of a BNN for relative dating are done
using the cyberSW database maintained by Archaeology Southwest, which provides
an extremely large ceramic assemblage dataset covering much of the US Southwest.
Importantly, the approach should be generalizable to any temporally diagnostic
artifact type for which count data exist.
