#####################################################
#       Packages that need to be installed          #
#####################################################
# From github - Note that these packages should be installed if you installed dynverse already. In this case, you only need to install dynplot, because the version in development have a fix for a bug
Sys.setenv(GITHUB_PAT = "b48dcc345ed2755770b0137d3ff94de1523a42c6")
install.packages("devtools")
library('devtools')
install_github("dynverse/dynfeature")
install_github("dynverse/dyno")
install_github("dynverse/dynplot", ref = "devel", force = T)
